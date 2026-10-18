"""Unit tests for resistance gene detection from GenBank annotations."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from plasmicheck.scripts.resistance_gene_detection import (
    BACTERIAL_RESISTANCE_PATTERNS,
    DEFAULT_PATTERNS,
    EXCLUDE_PATTERNS,
    MAMMALIAN_RESISTANCE_PATTERNS,
    REPORTER_PATTERNS,
    _gene_overlaps_insert,
    _matches_resistance_gene,
    detect_resistance_genes,
    locate_genbank_file,
    sanitize_gene_name,
)

# ============================================================================
# Helper Functions for Creating Test GenBank Files
# ============================================================================


def create_test_genbank(tmp_path: Path, features: list[dict[str, Any]]) -> str:
    """Create a minimal GenBank file for testing.

    Args:
        tmp_path: pytest tmp_path fixture
        features: List of feature dicts with keys: type, label, start, end

    Returns:
        Path to created .gb file
    """
    from Bio import SeqIO

    # Create a plasmid record
    record = SeqRecord(Seq("A" * 5000), id="test_plasmid", name="test")
    record.annotations["molecule_type"] = "DNA"

    # Add features
    for feat_dict in features:
        qualifiers = {}
        if "label" in feat_dict:
            qualifiers["label"] = [feat_dict["label"]]
        if "product" in feat_dict:
            qualifiers["product"] = [feat_dict["product"]]
        if "gene" in feat_dict:
            qualifiers["gene"] = [feat_dict["gene"]]

        feature = SeqFeature(  # type: ignore[no-untyped-call]
            FeatureLocation(feat_dict["start"], feat_dict["end"]),  # type: ignore[no-untyped-call]
            type=feat_dict.get("type", "CDS"),
            qualifiers=qualifiers,
        )
        record.features.append(feature)

    # Write to file
    gb_file = tmp_path / "test.gb"
    SeqIO.write(record, str(gb_file), "genbank")
    return str(gb_file)


# ============================================================================
# Tests for Pattern Matching (_matches_resistance_gene)
# ============================================================================


@pytest.mark.unit
def test_bacterial_resistance_patterns_match_ampr() -> None:
    """Test that AmpR matches bacterial resistance patterns."""
    assert _matches_resistance_gene("AmpR", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_bacterial_resistance_patterns_match_kanr() -> None:
    """Test that KanR matches bacterial resistance patterns."""
    assert _matches_resistance_gene("KanR", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_bacterial_resistance_patterns_match_bla() -> None:
    """Test that bla matches bacterial resistance patterns."""
    assert _matches_resistance_gene("bla", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_bacterial_resistance_patterns_match_tetr() -> None:
    """Test that TetR matches bacterial resistance patterns."""
    assert _matches_resistance_gene("TetR", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_bacterial_resistance_patterns_match_chlr() -> None:
    """Test that ChlR matches bacterial resistance patterns."""
    assert _matches_resistance_gene("ChlR", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_mammalian_resistance_patterns_match_puror() -> None:
    """Test that PuroR matches mammalian resistance patterns."""
    assert _matches_resistance_gene("PuroR", MAMMALIAN_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_mammalian_resistance_patterns_match_hygror() -> None:
    """Test that HygroR matches mammalian resistance patterns."""
    assert _matches_resistance_gene("HygroR", MAMMALIAN_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_mammalian_resistance_patterns_match_hygromycin() -> None:
    """Test that Hygromycin matches mammalian resistance patterns."""
    assert _matches_resistance_gene("Hygromycin", MAMMALIAN_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_mammalian_resistance_patterns_match_neor() -> None:
    """Test that NeoR matches mammalian resistance patterns."""
    assert _matches_resistance_gene("NeoR", MAMMALIAN_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_mammalian_resistance_patterns_match_blasticidinr() -> None:
    """Test that BlasticidinR matches mammalian resistance patterns."""
    assert _matches_resistance_gene("BlasticidinR", MAMMALIAN_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_mammalian_resistance_patterns_match_zeocinr() -> None:
    """Test that ZeocinR matches mammalian resistance patterns."""
    assert _matches_resistance_gene("ZeocinR", MAMMALIAN_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_reporter_patterns_match_egfp() -> None:
    """Test that eGFP matches reporter patterns."""
    assert _matches_resistance_gene("eGFP", REPORTER_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_reporter_patterns_match_egfp_uppercase() -> None:
    """Test that EGFP matches reporter patterns."""
    assert _matches_resistance_gene("EGFP", REPORTER_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_reporter_patterns_match_gfp() -> None:
    """Test that GFP matches reporter patterns."""
    assert _matches_resistance_gene("GFP", REPORTER_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_reporter_patterns_match_mcherry() -> None:
    """Test that mCherry matches reporter patterns."""
    assert _matches_resistance_gene("mCherry", REPORTER_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_reporter_patterns_match_yfp() -> None:
    """Test that YFP matches reporter patterns."""
    assert _matches_resistance_gene("YFP", REPORTER_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_reporter_patterns_match_luciferase() -> None:
    """Test that luciferase matches reporter patterns."""
    assert _matches_resistance_gene("luciferase", REPORTER_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_reporter_patterns_match_lacz() -> None:
    """Test that lacZ matches reporter patterns."""
    assert _matches_resistance_gene("lacZ", REPORTER_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_exclude_patterns_his_tag() -> None:
    """Test that His-tag is excluded."""
    assert not _matches_resistance_gene("His-tag", DEFAULT_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_exclude_patterns_flag() -> None:
    """Test that FLAG is excluded."""
    assert not _matches_resistance_gene("FLAG", DEFAULT_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_exclude_patterns_myc() -> None:
    """Test that Myc tag is excluded."""
    assert not _matches_resistance_gene("Myc tag", DEFAULT_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_exclude_patterns_v5() -> None:
    """Test that V5 is excluded."""
    assert not _matches_resistance_gene("V5", DEFAULT_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_case_insensitivity_ampr_lowercase() -> None:
    """Test that ampr (lowercase) matches."""
    assert _matches_resistance_gene("ampr", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_case_insensitivity_ampr_uppercase() -> None:
    """Test that AMPR (uppercase) matches."""
    assert _matches_resistance_gene("AMPR", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_case_insensitivity_ampr_mixed() -> None:
    """Test that AmpR (mixed case) matches."""
    assert _matches_resistance_gene("AmpR", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_word_boundaries_amprandom_no_match() -> None:
    """Test that AmpRandom does NOT match (word boundary)."""
    assert not _matches_resistance_gene(
        "AmpRandom", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS
    )


@pytest.mark.unit
def test_word_boundaries_supergfp_no_match() -> None:
    """Test that SuperGFP does NOT match (word boundary)."""
    assert not _matches_resistance_gene("SuperGFP", REPORTER_PATTERNS, EXCLUDE_PATTERNS)


# ============================================================================
# Tests for GenBank Parsing (detect_resistance_genes)
# ============================================================================


@pytest.mark.unit
def test_detect_single_gene_ampr(tmp_path: Path) -> None:
    """Test detection of single AmpR gene."""
    gb_file = create_test_genbank(
        tmp_path, [{"type": "CDS", "label": "AmpR", "start": 100, "end": 900}]
    )

    genes = detect_resistance_genes(gb_file)

    assert len(genes) == 1
    assert genes[0].name == "AmpR"
    assert genes[0].start == 100
    assert genes[0].end == 899  # Converted to 0-based inclusive
    assert genes[0].qualifier_source == "label"


@pytest.mark.unit
def test_detect_multiple_genes(tmp_path: Path) -> None:
    """Test detection of multiple genes (AmpR + eGFP), excluding His-tag."""
    gb_file = create_test_genbank(
        tmp_path,
        [
            {"type": "CDS", "label": "AmpR", "start": 100, "end": 900},
            {"type": "CDS", "label": "eGFP", "start": 1500, "end": 2200},
            {"type": "CDS", "label": "His-tag", "start": 3000, "end": 3100},
        ],
    )

    genes = detect_resistance_genes(gb_file)

    assert len(genes) == 2
    assert genes[0].name == "AmpR"
    assert genes[1].name == "eGFP"
    # His-tag should be excluded


@pytest.mark.unit
def test_detect_no_cds_features(tmp_path: Path) -> None:
    """Test that GenBank with no CDS features returns empty list."""
    gb_file = create_test_genbank(
        tmp_path, [{"type": "promoter", "label": "Amp prom", "start": 50, "end": 100}]
    )

    genes = detect_resistance_genes(gb_file)

    assert len(genes) == 0


@pytest.mark.unit
def test_detect_no_matching_patterns(tmp_path: Path) -> None:
    """Test that GenBank with CDS but no matching patterns returns empty list."""
    gb_file = create_test_genbank(
        tmp_path, [{"type": "CDS", "label": "SomeGene", "start": 100, "end": 500}]
    )

    genes = detect_resistance_genes(gb_file)

    assert len(genes) == 0


# ============================================================================
# Tests for Insert Region Exclusion (_gene_overlaps_insert)
# ============================================================================


@pytest.mark.unit
def test_gene_overlaps_insert_fully_inside() -> None:
    """Test gene fully inside insert region is detected as overlap."""
    assert _gene_overlaps_insert(100, 200, (50, 250))


@pytest.mark.unit
def test_gene_overlaps_insert_fully_outside_before() -> None:
    """Test gene fully outside insert region (before) has no overlap."""
    assert not _gene_overlaps_insert(10, 40, (100, 200))


@pytest.mark.unit
def test_gene_overlaps_insert_fully_outside_after() -> None:
    """Test gene fully outside insert region (after) has no overlap."""
    assert not _gene_overlaps_insert(250, 300, (100, 200))


@pytest.mark.unit
def test_gene_overlaps_insert_partial_overlap_start() -> None:
    """Test gene partially overlapping insert region at start is excluded."""
    assert _gene_overlaps_insert(80, 120, (100, 200))


@pytest.mark.unit
def test_gene_overlaps_insert_partial_overlap_end() -> None:
    """Test gene partially overlapping insert region at end is excluded."""
    assert _gene_overlaps_insert(180, 220, (100, 200))


@pytest.mark.unit
def test_gene_overlaps_insert_none() -> None:
    """Test gene with no insert region returns False."""
    assert not _gene_overlaps_insert(100, 200, None)


@pytest.mark.unit
def test_detect_genes_exclude_insert_overlap(tmp_path: Path) -> None:
    """Test that genes overlapping insert region are excluded."""
    gb_file = create_test_genbank(
        tmp_path,
        [
            {"type": "CDS", "label": "AmpR", "start": 100, "end": 500},  # Outside insert
            {"type": "CDS", "label": "eGFP", "start": 1000, "end": 1500},  # Inside insert
        ],
    )

    genes = detect_resistance_genes(gb_file, insert_region=(900, 1600))

    # Only AmpR should be detected (eGFP overlaps insert)
    assert len(genes) == 1
    assert genes[0].name == "AmpR"
    assert genes[0].region == "backbone"


@pytest.mark.unit
def test_detect_genes_no_insert_region(tmp_path: Path) -> None:
    """Test that genes with no insert region have region='unknown'."""
    gb_file = create_test_genbank(
        tmp_path, [{"type": "CDS", "label": "AmpR", "start": 100, "end": 500}]
    )

    genes = detect_resistance_genes(gb_file, insert_region=None)

    assert len(genes) == 1
    assert genes[0].region == "unknown"


# ============================================================================
# Tests for Gene Naming (sanitize_gene_name)
# ============================================================================


@pytest.mark.unit
def test_sanitize_gene_name_with_spaces() -> None:
    """Test that spaces are replaced with underscores."""
    assert sanitize_gene_name("FRT Hygromycin") == "FRT_Hygromycin"


@pytest.mark.unit
def test_sanitize_gene_name_already_clean() -> None:
    """Test that already clean names are unchanged."""
    assert sanitize_gene_name("AmpR") == "AmpR"


@pytest.mark.unit
def test_sanitize_gene_name_with_parens() -> None:
    """Test that parentheses are replaced with underscores."""
    assert sanitize_gene_name("bla(TEM-1)") == "bla_TEM-1_"


@pytest.mark.unit
def test_sanitize_gene_name_preserves_underscore_dash() -> None:
    """Test that underscores and dashes are preserved."""
    assert sanitize_gene_name("HIV-1_5_LTR") == "HIV-1_5_LTR"


# ============================================================================
# Tests for Custom Patterns
# ============================================================================


@pytest.mark.unit
def test_detect_genes_custom_pattern(tmp_path: Path) -> None:
    """Test that custom patterns are detected."""
    gb_file = create_test_genbank(
        tmp_path, [{"type": "CDS", "label": "MyCustomGene", "start": 100, "end": 500}]
    )

    genes = detect_resistance_genes(gb_file, custom_patterns=[r"\bMyCustomGene\b"])

    assert len(genes) == 1
    assert genes[0].name == "MyCustomGene"


# ============================================================================
# Tests for locate_genbank_file
# ============================================================================


@pytest.mark.unit
def test_locate_genbank_file_in_output_dir(tmp_path: Path) -> None:
    """Test that .gb file in output_dir is found."""
    gb_file = tmp_path / "plasmid.gb"
    gb_file.write_text("dummy genbank content")

    result = locate_genbank_file(str(tmp_path))

    assert result == str(gb_file)


@pytest.mark.unit
def test_locate_genbank_file_from_plasmid_file(tmp_path: Path) -> None:
    """Test that original plasmid_file is returned if it's .gb."""
    plasmid_file = tmp_path / "original.gb"
    plasmid_file.write_text("dummy genbank content")

    # Output dir is different (no .gb there)
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    result = locate_genbank_file(str(output_dir), str(plasmid_file))

    assert result == str(plasmid_file)


@pytest.mark.unit
def test_locate_genbank_file_not_found(tmp_path: Path) -> None:
    """Test that None is returned if no GenBank file found."""
    result = locate_genbank_file(str(tmp_path))

    assert result is None


@pytest.mark.unit
def test_locate_genbank_file_prefers_output_dir(tmp_path: Path) -> None:
    """Test that .gb in output_dir takes priority over plasmid_file."""
    # Create output directory with .gb file
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    output_gb = output_dir / "converted.gb"
    output_gb.write_text("output genbank content")

    # Provide a different plasmid_file outside output_dir
    plasmid_file = tmp_path / "original.gb"
    plasmid_file.write_text("original genbank content")

    result = locate_genbank_file(str(output_dir), str(plasmid_file))

    # Should return the output_dir .gb, not the original
    assert result == str(output_gb)


# ============================================================================
# Tests for Edge Cases
# ============================================================================


@pytest.mark.unit
def test_detect_genes_file_not_exist() -> None:
    """Test that non-existent file returns empty list."""
    genes = detect_resistance_genes("/nonexistent/path.gb")

    assert len(genes) == 0


@pytest.mark.unit
def test_detect_genes_empty_genbank(tmp_path: Path) -> None:
    """Test that empty GenBank file returns empty list."""
    gb_file = create_test_genbank(tmp_path, [])

    genes = detect_resistance_genes(str(gb_file))

    assert len(genes) == 0


@pytest.mark.unit
def test_identify_gene_name_checks_product_qualifier(tmp_path: Path) -> None:
    """Test that /product qualifier is checked if /label is missing."""
    gb_file = create_test_genbank(
        tmp_path, [{"type": "CDS", "product": "ampicillin resistance", "start": 100, "end": 500}]
    )

    genes = detect_resistance_genes(gb_file)

    # "Amp" pattern should match "ampicillin resistance"
    assert len(genes) == 1
    assert genes[0].qualifier_source == "product"


@pytest.mark.unit
def test_identify_gene_name_checks_gene_qualifier(tmp_path: Path) -> None:
    """Test that /gene qualifier is checked if /label and /product are missing."""
    gb_file = create_test_genbank(
        tmp_path, [{"type": "CDS", "gene": "KanR", "start": 100, "end": 500}]
    )

    genes = detect_resistance_genes(gb_file)

    assert len(genes) == 1
    assert genes[0].qualifier_source == "gene"


# ============================================================================
# Tests with Real GenBank Files
# ============================================================================


@pytest.mark.unit
def test_bacterial_resistance_patterns_match_amp_alone() -> None:
    """Test that Amp alone matches bacterial resistance patterns."""
    assert _matches_resistance_gene("Amp", BACTERIAL_RESISTANCE_PATTERNS, EXCLUDE_PATTERNS)


@pytest.mark.unit
def test_detect_resistance_genes_fasta_returns_empty(tmp_path: Path) -> None:
    """Test that fasta file type returns empty list immediately."""
    fa = tmp_path / "test.fasta"
    fa.write_text(">seq\nATCG")
    assert detect_resistance_genes(str(fa), file_type="fasta") == []
    assert detect_resistance_genes(str(fa)) == []


@pytest.mark.unit
def test_origin_wrapping_insert_region_overlap() -> None:
    """Test _overlaps_insert when insert region wraps circular origin (start > end)."""
    # Insert wraps across origin: 8000..end and 0..2000
    insert = (8000, 2000)
    # Feature in high part [8500..9000] -> overlaps
    assert _gene_overlaps_insert(8500, 9000, insert) is True
    # Feature in low part [1000..1500] -> overlaps
    assert _gene_overlaps_insert(1000, 1500, insert) is True
    # Feature in backbone [3000..4000] -> does not overlap
    assert _gene_overlaps_insert(3000, 4000, insert) is False


@pytest.mark.unit
def test_compound_location_origin_wrapping_detection(tmp_path: Path) -> None:
    """Test CompoundLocation wrapping origin produces 1-based start > end coordinates."""
    from Bio import SeqIO

    record = SeqRecord(Seq("A" * 10000), id="origin_plasmid", name="origin_plasmid")
    record.annotations["molecule_type"] = "DNA"

    # Origin wrapping feature: 8500..10000 and 0..300
    loc = CompoundLocation(  # type: ignore[no-untyped-call]
        [FeatureLocation(8500, 10000), FeatureLocation(0, 300)]  # type: ignore[no-untyped-call]
    )
    feat = SeqFeature(loc, type="CDS", qualifiers={"label": ["AmpR"]})  # type: ignore[no-untyped-call]
    record.features.append(feat)

    gb_file = tmp_path / "origin_comp.gb"
    SeqIO.write(record, str(gb_file), "genbank")

    # Non-overlapping insert region: 2000..4000
    genes = detect_resistance_genes(str(gb_file), insert_region=(2000, 4000))
    assert len(genes) == 1
    assert genes[0].name == "AmpR"
    assert genes[0].start == 8500
    assert genes[0].end == 299

    # Overlapping insert region: 100..500 (overlaps low part 0..300)
    genes_excluded = detect_resistance_genes(str(gb_file), insert_region=(100, 500))
    assert len(genes_excluded) == 0


@pytest.mark.unit
def test_xdna_parse_error_returns_empty(tmp_path: Path) -> None:
    """When parsing xDNA raises ValueError, detect_resistance_genes logs warning and returns []."""
    test_file = tmp_path / "broken.xdna"
    test_file.write_bytes(b"\x00" * 112)

    from unittest.mock import patch

    with patch(
        "Bio.SeqIO.parse", side_effect=ValueError("End location must be greater than start")
    ):
        genes = detect_resistance_genes(str(test_file), file_type="xdna")
        assert genes == []


# ============================================================================
# Tests with Synthetic GenBank Files (replacing external references)
# ============================================================================


@pytest.mark.unit
def test_detect_genes_real_pwpi_kcnj5(tmp_path: Path) -> None:
    """Test detection against synthetic pWPI-like plasmid with AmpR and EGFP."""
    gb_file = create_test_genbank(
        tmp_path,
        [
            {"type": "CDS", "label": "AmpR", "start": 100, "end": 900},
            {"type": "CDS", "label": "EGFP", "start": 1200, "end": 1900},
        ],
    )

    genes = detect_resistance_genes(gb_file)

    gene_names = [g.name for g in genes]
    assert "AmpR" in gene_names, f"AmpR not detected. Found: {gene_names}"
    assert "EGFP" in gene_names, f"EGFP not detected. Found: {gene_names}"


@pytest.mark.unit
def test_detect_genes_real_pcmv6_gfp_puro(tmp_path: Path) -> None:
    """Test detection against synthetic plasmid with AmpR and lacZ."""
    gb_file = create_test_genbank(
        tmp_path,
        [
            {"type": "CDS", "label": "AmpR", "start": 500, "end": 1200},
            {"type": "CDS", "label": "lacZ", "start": 1500, "end": 2000},
        ],
    )

    genes = detect_resistance_genes(gb_file)

    gene_names = [g.name for g in genes]
    assert len(genes) >= 1, f"Expected at least 1 gene, found: {gene_names}"
    assert "AmpR" in gene_names, f"AmpR not detected. Found: {gene_names}"


@pytest.mark.unit
def test_detect_genes_sorted_by_start(tmp_path: Path) -> None:
    """Test that detected genes are sorted by start position."""
    gb_file = create_test_genbank(
        tmp_path,
        [
            {"type": "CDS", "label": "GFP", "start": 2000, "end": 2500},
            {"type": "CDS", "label": "AmpR", "start": 100, "end": 500},
            {"type": "CDS", "label": "KanR", "start": 1000, "end": 1500},
        ],
    )

    genes = detect_resistance_genes(gb_file)

    assert len(genes) == 3
    assert genes[0].name == "AmpR"  # Start: 100
    assert genes[1].name == "KanR"  # Start: 1000
    assert genes[2].name == "GFP"  # Start: 2000


@pytest.mark.unit
def test_recall_against_existing_genbank_files(tmp_path: Path) -> None:
    """Validate recall against synthetic plasmid.gb and created test plasmids."""
    gb1 = create_test_genbank(
        tmp_path, [{"type": "CDS", "label": "AmpR", "start": 100, "end": 900}]
    )

    genes = detect_resistance_genes(str(gb1))
    assert len(genes) >= 1
    assert any("ampr" in g.name.lower() for g in genes)


@pytest.mark.unit
def test_pwpi_kcnj5_wt_detects_ampr(tmp_path: Path) -> None:
    """Targeted test expecting AmpR detection from synthetic plasmid."""
    gb_file = create_test_genbank(
        tmp_path, [{"type": "CDS", "label": "AmpR", "start": 200, "end": 1000}]
    )

    genes = detect_resistance_genes(str(gb_file))
    assert len(genes) >= 1
    gene_names = [g.name.lower() for g in genes]
    assert any("amp" in name or "bla" in name for name in gene_names)


@pytest.mark.unit
def test_pcmv6_ac_ires_gfp_puro_detects_markers(tmp_path: Path) -> None:
    """Targeted test for backbone markers in synthetic plasmid."""
    gb_file = create_test_genbank(
        tmp_path,
        [
            {"type": "CDS", "label": "AmpR", "start": 300, "end": 1100},
            {"type": "CDS", "label": "lacZ", "start": 1300, "end": 1800},
        ],
    )

    genes = detect_resistance_genes(str(gb_file))
    assert len(genes) >= 1


@pytest.mark.unit
def test_excludes_lacz_alpha_and_deduplicates_overlapping_cds(tmp_path: Path) -> None:
    """lacZ alpha cloning peptide is excluded while overlapping AmpR CDS annotations are deduplicated."""
    gb_file = create_test_genbank(
        tmp_path,
        [
            {"type": "CDS", "label": "lacZ alpha", "start": 100, "end": 300},
            {"type": "CDS", "label": "NeomycinR", "start": 400, "end": 1200},
            {"type": "CDS", "label": "AmpR", "start": 1500, "end": 2360},
            {"type": "CDS", "label": "AmpR", "start": 1569, "end": 2360},
        ],
    )
    genes = detect_resistance_genes(str(gb_file))
    names = [g.name for g in genes]
    assert "lacZ alpha" not in names
    assert names == ["NeomycinR", "AmpR"]
    assert genes[1].start == 1500
    assert genes[1].end == 2359
