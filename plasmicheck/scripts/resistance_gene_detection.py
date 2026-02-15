"""Resistance and reporter gene detection from GenBank annotations.

This module identifies resistance genes, selection markers, and reporter genes
from plasmid GenBank CDS features using curated regex patterns.
"""

from __future__ import annotations

import glob
import logging
import os
import re
from dataclasses import dataclass
from typing import Any

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

# ============================================================================
# Curated Pattern Lists
# ============================================================================

# Bacterial resistance genes (antibiotics)
BACTERIAL_RESISTANCE_PATTERNS: list[re.Pattern[str]] = [
    re.compile(r"\bAmpR\b", re.IGNORECASE),
    re.compile(r"\bampicillin\b", re.IGNORECASE),
    re.compile(r"\bbla\b", re.IGNORECASE),
    re.compile(r"\bKanR\b", re.IGNORECASE),
    re.compile(r"\bkanamycin\b", re.IGNORECASE),
    re.compile(r"\bChlR\b", re.IGNORECASE),
    re.compile(r"\bchloramphenicol\b", re.IGNORECASE),
    re.compile(r"\bTetR\b", re.IGNORECASE),
    re.compile(r"\btetracycline\b", re.IGNORECASE),
    re.compile(r"\bSpec\b", re.IGNORECASE),
    re.compile(r"\bStrep\b", re.IGNORECASE),
    re.compile(r"\bCamR\b", re.IGNORECASE),
    re.compile(r"\bGenR\b", re.IGNORECASE),
    re.compile(r"\bGentamicin\b", re.IGNORECASE),
]

# Mammalian selection markers
MAMMALIAN_RESISTANCE_PATTERNS: list[re.Pattern[str]] = [
    re.compile(r"\bPuroR\b", re.IGNORECASE),
    re.compile(r"\bPuro\b", re.IGNORECASE),
    re.compile(r"\bHygroR\b", re.IGNORECASE),
    re.compile(r"\bHygro\b", re.IGNORECASE),
    re.compile(r"\bHygromycin\b", re.IGNORECASE),
    re.compile(r"\bNeoR\b", re.IGNORECASE),
    re.compile(r"\bNeo\b", re.IGNORECASE),
    re.compile(r"\bNeomycin\b", re.IGNORECASE),
    re.compile(r"\bBlasticidinR\b", re.IGNORECASE),
    re.compile(r"\bBsd\b", re.IGNORECASE),
    re.compile(r"\bBsdR\b", re.IGNORECASE),
    re.compile(r"\bZeocinR\b", re.IGNORECASE),
    re.compile(r"\bZeocin\b", re.IGNORECASE),
    re.compile(r"\bBleoR\b", re.IGNORECASE),
]

# Reporter genes
REPORTER_PATTERNS: list[re.Pattern[str]] = [
    re.compile(r"\b[em]?GFP\b", re.IGNORECASE),  # GFP, eGFP, mGFP
    re.compile(r"\bEGFP\b", re.IGNORECASE),  # Explicit EGFP
    re.compile(r"\bmCherry\b", re.IGNORECASE),
    re.compile(r"\b[tdm]?RFP\b", re.IGNORECASE),  # RFP, tRFP, dRFP, mRFP
    re.compile(r"\bdsRed\b", re.IGNORECASE),
    re.compile(r"\bYFP\b", re.IGNORECASE),
    re.compile(r"\bCFP\b", re.IGNORECASE),
    re.compile(r"\bluciferase\b", re.IGNORECASE),
    re.compile(r"\blacZ\b", re.IGNORECASE),
]

# Affinity tags to exclude (not resistance/reporter genes)
EXCLUDE_PATTERNS: list[re.Pattern[str]] = [
    re.compile(r"\bHis[-_]?tag\b", re.IGNORECASE),
    re.compile(r"\b6His\b", re.IGNORECASE),
    re.compile(r"\bFLAG\b", re.IGNORECASE),
    re.compile(r"\bMyc\b", re.IGNORECASE),
    re.compile(r"\bV5\b", re.IGNORECASE),
    re.compile(r"\bHA\b", re.IGNORECASE),
]

# Combined default patterns
DEFAULT_PATTERNS: list[re.Pattern[str]] = (
    BACTERIAL_RESISTANCE_PATTERNS + MAMMALIAN_RESISTANCE_PATTERNS + REPORTER_PATTERNS
)


# ============================================================================
# Data Structure
# ============================================================================


@dataclass
class DetectedGene:
    """A resistance/reporter gene detected from GenBank annotations."""

    name: str  # Verbatim GenBank label (display name)
    sanitized_name: str  # For TSV row keys (spaces/special -> underscores)
    start: int  # 0-based inclusive
    end: int  # 0-based inclusive
    region: str  # "backbone" or "insert" (or "unknown" if no insert region)
    qualifier_source: str  # Which qualifier matched ("label", "product", "gene")


# ============================================================================
# Helper Functions
# ============================================================================


def sanitize_gene_name(name: str) -> str:
    """Replace non-alphanumeric chars (except _ and -) with underscores.

    Args:
        name: Gene name from GenBank annotation

    Returns:
        Sanitized name suitable for TSV column headers
    """
    return re.sub(r"[^a-zA-Z0-9_-]", "_", name)


def _matches_resistance_gene(
    text: str, patterns: list[re.Pattern[str]], exclude_patterns: list[re.Pattern[str]]
) -> bool:
    """Check if text matches any pattern and is not excluded.

    Args:
        text: Text to match (qualifier value)
        patterns: List of compiled regex patterns to match
        exclude_patterns: List of compiled regex patterns to exclude

    Returns:
        True if text matches a pattern and is not excluded
    """
    if not text:
        return False

    # Check exclusions first
    for exclude_pattern in exclude_patterns:
        if exclude_pattern.search(text):
            return False

    # Check if matches any pattern
    return any(pattern.search(text) for pattern in patterns)


def _identify_gene_name(
    qualifiers: dict[str, list[str]],
    patterns: list[re.Pattern[str]],
    exclude_patterns: list[re.Pattern[str]],
) -> tuple[str, str] | None:
    """Return (gene_name, qualifier_source) or None.

    Checks qualifiers in priority order: /label, /product, /gene.

    Args:
        qualifiers: Feature qualifiers dict
        patterns: List of compiled regex patterns to match
        exclude_patterns: List of compiled regex patterns to exclude

    Returns:
        (gene_name, qualifier_source) if a match is found, else None
    """
    # Priority order for qualifier sources
    qualifier_priority = [
        ("label", qualifiers.get("label", [])),
        ("product", qualifiers.get("product", [])),
        ("gene", qualifiers.get("gene", [])),
    ]

    for qualifier_source, values in qualifier_priority:
        for value in values:
            if _matches_resistance_gene(value, patterns, exclude_patterns):
                return (value, qualifier_source)

    return None


def _gene_overlaps_insert(
    gene_start: int, gene_end: int, insert_region: tuple[int, int] | None
) -> bool:
    """Check if gene overlaps insert region.

    Args:
        gene_start: Gene start position (0-based inclusive)
        gene_end: Gene end position (0-based inclusive)
        insert_region: (start, end) inclusive boundaries from cDNA_positions.txt, or None

    Returns:
        True if gene overlaps insert region, False if no overlap or insert_region is None
    """
    if insert_region is None:
        return False

    insert_start, insert_end = insert_region

    # Check for any overlap (conservative: partial overlap = excluded)
    return not (gene_end < insert_start or gene_start > insert_end)


# ============================================================================
# Main Functions
# ============================================================================


def detect_resistance_genes(
    genbank_file: str,
    insert_region: tuple[int, int] | None = None,
    custom_patterns: list[str] | None = None,
) -> list[DetectedGene]:
    """Parse GenBank CDS features, match against patterns, exclude insert-overlapping genes.

    Args:
        genbank_file: Path to .gb/.gbk file
        insert_region: (start, end) inclusive boundaries from cDNA_positions.txt, or None
        custom_patterns: Additional regex patterns from config.json (appended to defaults)

    Returns:
        List of DetectedGene objects for matched genes outside the insert region.
        Empty list if no genes match or file cannot be parsed.
    """
    if not os.path.exists(genbank_file):
        logging.warning(f"GenBank file not found: {genbank_file}")
        return []

    # Combine default patterns with custom patterns
    patterns: list[re.Pattern[str]] = list(DEFAULT_PATTERNS)
    if custom_patterns:
        patterns.extend([re.compile(pattern, re.IGNORECASE) for pattern in custom_patterns])

    detected_genes: list[DetectedGene] = []

    try:
        # Parse GenBank file
        for record in SeqIO.parse(genbank_file, "genbank"):  # type: ignore[no-untyped-call]
            # Process each CDS feature
            for feature in record.features:
                if feature.type != "CDS":
                    continue

                # Get qualifiers
                qualifiers: dict[str, list[Any]] = feature.qualifiers

                # Try to identify a resistance/reporter gene
                result = _identify_gene_name(qualifiers, patterns, EXCLUDE_PATTERNS)
                if result is None:
                    continue

                gene_name, qualifier_source = result

                # Extract coordinates
                # Handle CompoundLocation (origin-crossing features)
                if isinstance(feature.location, CompoundLocation):
                    # Take min start and max end across all parts
                    parts = feature.location.parts
                    gene_start = min(int(part.start) for part in parts)
                    gene_end = max(int(part.end) - 1 for part in parts)  # Convert to inclusive
                else:
                    gene_start = int(feature.location.start)  # 0-based
                    gene_end = int(feature.location.end) - 1  # Convert to 0-based inclusive

                # Check if gene overlaps insert region
                if _gene_overlaps_insert(gene_start, gene_end, insert_region):
                    logging.debug(
                        f"Excluding gene '{gene_name}' (overlaps insert region): "
                        f"{gene_start}-{gene_end}"
                    )
                    continue

                # Classify region
                region = "unknown" if insert_region is None else "backbone"

                # Create DetectedGene
                detected_gene = DetectedGene(
                    name=gene_name,
                    sanitized_name=sanitize_gene_name(gene_name),
                    start=gene_start,
                    end=gene_end,
                    region=region,
                    qualifier_source=qualifier_source,
                )
                detected_genes.append(detected_gene)

        # Sort by start position
        detected_genes.sort(key=lambda g: g.start)

        logging.info(f"Detected {len(detected_genes)} resistance/reporter genes in {genbank_file}")

    except Exception as e:
        logging.warning(f"Error parsing GenBank file {genbank_file}: {e}")
        return []

    return detected_genes


def locate_genbank_file(output_dir: str, plasmid_file: str | None = None) -> str | None:
    """Locate the GenBank file for gene detection.

    Checks (in order):
    1. .gb file in output_dir/references/ (new subfolder structure)
    2. .gb file in output_dir (flat structure, backward compatibility)
    3. Original plasmid_file if it's .gb/.gbk

    Args:
        output_dir: Directory to search for .gb files
        plasmid_file: Original plasmid file path (may be .gb, .gbk, .xdna, etc.)

    Returns:
        Path to GenBank file, or None if not found
    """
    # Check for .gb file in output_dir and references/ subdirectory
    if os.path.exists(output_dir):
        gb_files = glob.glob(os.path.join(output_dir, "*.gb"))
        gbk_files = glob.glob(os.path.join(output_dir, "*.gbk"))

        # Also check references/ subdirectory (new structure)
        refs_dir = os.path.join(output_dir, "references")
        if os.path.isdir(refs_dir):
            gb_files += glob.glob(os.path.join(refs_dir, "*.gb"))
            gbk_files += glob.glob(os.path.join(refs_dir, "*.gbk"))

        all_gb_files = gb_files + gbk_files

        if all_gb_files:
            # Return the first .gb file found
            return all_gb_files[0]

    # Check if plasmid_file is a GenBank file
    if (
        plasmid_file
        and plasmid_file.endswith((".gb", ".gbk", ".genbank"))
        and os.path.exists(plasmid_file)
    ):
        return plasmid_file

    return None
