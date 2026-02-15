# Phase 10: Resistance Gene Detection - Research

**Researched:** 2026-02-15
**Domain:** GenBank annotation parsing, resistance gene pattern matching, per-gene coverage metrics
**Confidence:** HIGH

## Summary

Phase 10 adds automatic resistance gene detection from plasmid GenBank annotations, computing per-gene coverage metrics (mean depth, breadth, coordinates) as an independent contamination signal. The standard approach uses Biopython SeqIO to parse CDS features from GenBank files, regex pattern matching against qualifiers (/label, /product, /gene) to identify resistance/reporter genes, and pysam.count_coverage() to compute per-gene metrics using the existing Phase 9 infrastructure.

The key technical challenge is that xDNA files are currently converted to GenBank then immediately deleted (line 68 in convert_plasmidfile_to_fasta.py). The solution is to either preserve the intermediate .gb file or parse xDNA directly with Biopython (SeqIO supports both formats). GenBank files (.gb/.gbk) are used directly without conversion.

**Primary recommendation:** Preserve the intermediate .gb file during xDNA conversion by modifying convert_plasmidfile_to_fasta.py to optionally skip deletion, then parse GenBank annotations in compare_alignments.py after coverage metrics computation.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Biopython | 1.75+ | GenBank/xDNA parsing | Industry standard for biological sequence formats, built-in SeqIO module handles GenBank CDS features with qualifiers |
| pysam | 0.17+ | Per-gene coverage computation | Already in use (Phase 9), count_coverage() provides per-base depth arrays for any genomic region |
| Python re | stdlib | Regex pattern matching | Standard library, case-insensitive word-boundary matching for gene names |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| numpy | 1.20+ | Coverage array aggregation | Already in use (Phase 9), compute mean/breadth from pysam depth arrays |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Regex patterns | BLAST against CARD/ResFinder | Higher recall but requires external databases, network access, slow; deferred to future phase |
| GenBank parsing | Custom text parsing | Fragile, won't handle join() locations or edge cases; Biopython is battle-tested |

**Installation:**
Already present in project dependencies (pyproject.toml). No new libraries required.

## Architecture Patterns

### Recommended Project Structure
```
plasmicheck/scripts/
├── coverage_metrics.py          # Phase 9 infrastructure (reuse)
├── resistance_gene_detection.py # NEW: gene detection module
└── compare_alignments.py        # Integration point (add gene detection call)
```

### Pattern 1: GenBank CDS Feature Parsing
**What:** Parse CDS features from GenBank records, extract qualifiers, match against curated regex patterns
**When to use:** All GenBank/xDNA plasmid files
**Example:**
```python
# Source: https://widdowquinn.github.io/2018-03-06-ibioic/01-introduction/03-parsing.html
from Bio import SeqIO

for record in SeqIO.parse("plasmid.gb", "genbank"):
    for feature in record.features:
        if feature.type == "CDS":
            # Qualifiers are a dict, values are lists
            label = feature.qualifiers.get("label", [""])[0]
            product = feature.qualifiers.get("product", [""])[0]
            gene = feature.qualifiers.get("gene", [""])[0]

            # Location: feature.location.start (0-based), feature.location.end (exclusive)
            # For join() locations: isinstance(feature.location, CompoundLocation)
            start = int(feature.location.start)
            end = int(feature.location.end) - 1  # Convert to inclusive
```

**Key insight:** Qualifiers are always lists (even single values), access with [0]. Locations use 0-based half-open intervals (start inclusive, end exclusive). Join locations (e.g., CDS crossing plasmid origin) are CompoundLocation objects.

### Pattern 2: Regex Pattern Matching for Gene Names
**What:** Case-insensitive, word-boundary regex matching against gene name patterns
**When to use:** Identifying resistance/reporter genes from CDS qualifiers
**Example:**
```python
# Source: https://docs.python.org/3/library/re.html
import re

# Curated patterns (bacterial resistance genes)
PATTERNS = [
    r'\bAmpR\b',           # Ampicillin resistance
    r'\bKanR\b',           # Kanamycin resistance
    r'\bbla\b',            # Beta-lactamase family (bla, blaTEM, blaCTX-M, etc.)
    r'\bPuroR\b',          # Puromycin resistance (mammalian)
    r'\bHygroR?\b',        # Hygromycin resistance
    r'\b[em]?GFP\b',       # GFP variants (GFP, eGFP, mGFP)
    r'\bmCherry\b',        # mCherry reporter
]

def matches_pattern(text: str, patterns: list[str]) -> str | None:
    """Check if text matches any pattern, return pattern name."""
    for pattern in patterns:
        if re.search(pattern, text, re.IGNORECASE):
            return pattern
    return None

# Usage: check all qualifiers
for qualifier_value in [label, product, gene]:
    if matches_pattern(qualifier_value, PATTERNS):
        # Found a resistance gene
        break
```

**Key insight:** Use `\b` word boundaries to avoid false positives (e.g., "AmpR" matches but "AmpRandom" does not). `re.IGNORECASE` handles case variations (AmpR, ampr, AMPR). Escape special characters with `re.escape()` if pattern contains regex metacharacters.

### Pattern 3: Per-Gene Coverage Computation
**What:** Compute mean depth and breadth for specific gene coordinates using pysam
**When to use:** After detecting resistance genes, to quantify gene-level coverage
**Example:**
```python
# Source: Phase 9 coverage_metrics.py (existing code)
from plasmicheck.scripts.coverage_metrics import get_depth_array, compute_coverage_stats

# Gene coordinates (0-based inclusive start, inclusive end from GenBank)
gene_start = 12714  # From GenBank feature
gene_end = 13373

# Convert to half-open interval for pysam
depth_array = get_depth_array(
    bam_path="plasmid_alignment.bam",
    contig="plasmid_contig",
    start=gene_start,
    end=gene_end + 1  # pysam uses half-open [start, end)
)

# Compute statistics
stats = compute_coverage_stats(depth_array, breadth_thresholds=[1, 5])
# Returns: {'mean_depth': 42.5, 'median_depth': 40.0, 'breadth_1x': 0.98, ...}
```

**Key insight:** Reuse existing Phase 9 infrastructure. GenBank coordinates are 1-based inclusive on both ends, but Biopython converts to 0-based half-open internally. For pysam, pass (start, end+1) to match half-open convention.

### Pattern 4: Insert Region Overlap Detection
**What:** Exclude genes overlapping the insert region (experimental target, not contamination signal)
**When to use:** After detecting genes, before computing coverage
**Example:**
```python
def gene_overlaps_insert(gene_start: int, gene_end: int, insert_region: tuple[int, int]) -> bool:
    """Check if gene overlaps insert region (inclusive boundaries)."""
    insert_start, insert_end = insert_region
    # No overlap: gene ends before insert starts, or gene starts after insert ends
    return not (gene_end < insert_start or gene_start > insert_end)

# Usage: skip genes in insert region
if not gene_overlaps_insert(gene_start, gene_end, (insert_start, insert_end)):
    # Compute coverage for this gene
    pass
```

### Pattern 5: Summary.tsv Output Format
**What:** Write per-gene metrics to summary.tsv with sanitized row keys
**When to use:** After computing gene coverage metrics
**Example:**
```python
# Sanitize gene name for TSV row key (replace spaces/special chars with underscores)
def sanitize_gene_name(name: str) -> str:
    """Sanitize gene name for TSV row keys."""
    return re.sub(r'[^a-zA-Z0-9_-]', '_', name)

# Write to summary.tsv
with open("summary.tsv", "a") as f:
    # Summary row: comma-separated list of detected genes
    gene_names = ["AmpR", "FRT Hygromycin", "EGFP"]
    f.write(f"ResistanceGenes\t{','.join(gene_names)}\n")

    # Per-gene metrics
    for gene_name, metrics in gene_metrics.items():
        sanitized = sanitize_gene_name(gene_name)
        f.write(f"ResGene_{sanitized}_MeanDepth\t{metrics['mean_depth']:.2f}\n")
        f.write(f"ResGene_{sanitized}_Breadth\t{metrics['breadth_1x']:.2f}\n")
        f.write(f"ResGene_{sanitized}_Start\t{metrics['start']}\n")
        f.write(f"ResGene_{sanitized}_End\t{metrics['end']}\n")
        f.write(f"ResGene_{sanitized}_Region\t{metrics['region']}\n")  # "backbone" or "insert"
```

**Key insight:** TSV row keys must be valid identifiers (no spaces). Display names in ResistanceGenes summary row keep original GenBank labels.

### Anti-Patterns to Avoid
- **Don't parse GenBank as text:** Use Biopython SeqIO, not regex on raw GenBank format. GenBank format is complex (join locations, multiline fields, escaping).
- **Don't compute coverage for insert-overlapping genes:** They're the experimental target, not contamination markers. Filter them out before coverage computation.
- **Don't use greedy regex:** Use `\b` word boundaries to avoid false positives (e.g., `\bGFP\b` not `.*GFP.*`).
- **Don't assume single qualifier:** Check /label, /product, /gene (in that order of priority). Some GenBank files use different qualifiers.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| GenBank parsing | Custom text parser | Biopython SeqIO | Handles edge cases: join() locations, circular plasmids, escaped quotes in qualifiers, multiline fields |
| Per-gene coverage | Manual BAM iteration | pysam.count_coverage() | Already in use (Phase 9), handles paired-end reads, quality thresholds, edge cases |
| Coverage statistics | Custom aggregation | coverage_metrics.py functions | Phase 9 infrastructure, tested, handles empty arrays, zero-division |
| Gene name sanitization | Manual string replacement | re.sub(r'[^a-zA-Z0-9_-]', '_', name) | Handles Unicode, special chars, consistent with existing sanitize_filename() pattern |

**Key insight:** xDNA → GenBank conversion already exists (convert_xdna_to_gb()), don't reimplement. Just preserve the intermediate .gb file instead of deleting it.

## Common Pitfalls

### Pitfall 1: xDNA Intermediate .gb File Deleted
**What goes wrong:** xDNA files are converted to .gb then immediately deleted (line 68, convert_plasmidfile_to_fasta.py). Gene detection needs this .gb file.
**Why it happens:** Original design only needed FASTA output, .gb was temporary.
**How to avoid:**
- Option A: Modify convert_plasmidfile_to_fasta() to optionally preserve .gb file (add `preserve_genbank=True` parameter)
- Option B: Parse xDNA directly with Biopython (SeqIO.parse(file, "xdna")) instead of converting
- Option C: Preserve .gb file in output_subfolder, delete after gene detection completes
**Warning signs:** FileNotFoundError when trying to open .gb file for gene detection
**Recommendation:** Option A (preserve_genbank parameter). Cleanest, backward-compatible, allows gene detection to run after FASTA conversion.

### Pitfall 2: Qualifier Values Are Lists
**What goes wrong:** Accessing `feature.qualifiers["label"]` returns a list `["AmpR"]`, not a string. Code expects string, crashes.
**Why it happens:** GenBank format allows multiple values per qualifier (rare but legal). Biopython always returns lists for consistency.
**How to avoid:** Always access with `[0]` and handle missing keys with `.get()`:
```python
label = feature.qualifiers.get("label", [""])[0]
```
**Warning signs:** TypeError: 'list' object is not callable, or unexpected list in string operations
**Recommendation:** Use `.get(key, [""])[0]` pattern consistently. Safer than try/except.

### Pitfall 3: GenBank Coordinates vs. Pysam Half-Open Intervals
**What goes wrong:** Biopython feature.location uses 0-based half-open [start, end). GenBank text format uses 1-based inclusive [start, end]. Pysam uses 0-based half-open. Mixing these causes off-by-one errors.
**Why it happens:** Different conventions across tools. Biopython converts GenBank → 0-based internally.
**How to avoid:**
- Biopython locations are already 0-based half-open: use `feature.location.start` and `feature.location.end` directly for pysam
- For inclusive end coordinate (summary.tsv output), use `feature.location.end - 1`
**Warning signs:** Coverage metrics off by 1 base, gene coordinates don't match IGV visualization
**Recommendation:** Document coordinate system in code comments. Test with known gene (e.g., AmpR in pWPI-KCNJ5-WT.gb at complement(4650..5309)).

### Pitfall 4: Circular Plasmid Features Crossing Origin
**What goes wrong:** CDS feature crosses plasmid origin (e.g., position 12000..500 on 12332bp plasmid). Simple start/end logic fails.
**Why it happens:** Circular plasmids can have features wrapping around position 0. GenBank represents these as join() locations.
**How to avoid:**
- Check if `isinstance(feature.location, CompoundLocation)`
- For CompoundLocation, use `feature.location.parts` to get sub-regions
- Most common case: two parts (e.g., [12000..12332] + [0..500])
- Concatenate coverage arrays from both parts before computing statistics
**Warning signs:** Negative start coordinates, start > end, missing coverage for genes known to be present
**Recommendation:** Test with existing plasmids. None of the 6 GenBank files in reference/plasmid/ appear to have origin-crossing CDS features, but handle defensively.

### Pitfall 5: Low Confidence Gene Name Matches
**What goes wrong:** Regex pattern too broad, matches affinity tags (His-tag, FLAG, Myc) or ambiguous names (e.g., "ampicillin" in a comment field).
**Why it happens:** GenBank qualifiers contain free text, inconsistent naming conventions.
**How to avoid:**
- Use word boundaries: `\bAmpR\b` not `Amp`
- Check multiple qualifiers in priority order: /label first (most curated), then /product, then /gene
- Curated pattern list with high-confidence patterns only
- Skip matches that also contain tag keywords ("His", "FLAG", "Myc", "V5")
**Warning signs:** Reporter genes detected in plasmids without fluorescent proteins, false positives on tag sequences
**Recommendation:** Start with conservative pattern list, log all matches during development, validate against existing 6 GenBank files (90% recall requirement).

### Pitfall 6: FASTA-Only Plasmids
**What goes wrong:** User provides FASTA plasmid file (no annotations). Gene detection fails or crashes.
**Why it happens:** FASTA format has no gene annotations, only sequence.
**How to avoid:**
- Check file extension before gene detection (.gb, .gbk, .xdna → has annotations; .fasta, .fa → no annotations)
- If FASTA-only, set `ResistanceGenesAvailable: false` in summary.tsv, log warning, skip gene detection
- If GenBank available but no genes detected, set `ResistanceGenesAvailable: true, ResistanceGenes: None`
**Warning signs:** FileNotFoundError, BioPython parse errors on FASTA files
**Recommendation:** Add file type check at start of gene detection function. Graceful degradation, not hard failure.

## Code Examples

Verified patterns from official sources and existing codebase:

### GenBank CDS Feature Parsing
```python
# Source: https://widdowquinn.github.io/2018-03-06-ibioic/01-introduction/03-parsing.html
# Verified: Biopython 1.75+ documentation
from Bio import SeqIO

def parse_genbank_cds_features(genbank_file: str) -> list[dict]:
    """Parse CDS features from GenBank file, extract qualifiers and coordinates.

    Returns:
        List of dicts with keys: 'label', 'product', 'gene', 'start', 'end'
    """
    features = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                # Qualifiers are dicts with list values
                label = feature.qualifiers.get("label", [""])[0]
                product = feature.qualifiers.get("product", [""])[0]
                gene = feature.qualifiers.get("gene", [""])[0]

                # Coordinates: 0-based half-open
                start = int(feature.location.start)
                end = int(feature.location.end) - 1  # Convert to inclusive

                features.append({
                    "label": label,
                    "product": product,
                    "gene": gene,
                    "start": start,
                    "end": end,
                })
    return features
```

### Resistance Gene Pattern Matching
```python
# Source: Curated from Addgene blog and CARD database conventions
# https://blog.addgene.org/plasmids-101-everything-you-need-to-know-about-antibiotic-resistance-genes
# https://en.vectorbuilder.com/resources/faq/optimize-drug-selection-markers.html
import re

# Tier 1: Drug resistance genes
BACTERIAL_RESISTANCE_PATTERNS = [
    r'\bAmpR\b',                  # Ampicillin resistance
    r'\bAmp\b',                   # Short form
    r'\bbla\b',                   # Beta-lactamase family
    r'\bKanR?\b',                 # Kanamycin resistance (KanR, Kan)
    r'\bChlR?\b',                 # Chloramphenicol resistance
    r'\bTetR?\b',                 # Tetracycline resistance
    r'\bSpec\b',                  # Spectinomycin resistance
    r'\bStrep\b',                 # Streptomycin resistance
]

MAMMALIAN_RESISTANCE_PATTERNS = [
    r'\bPuroR?\b',                # Puromycin resistance (PuroR, Puro)
    r'\bHygroR?\b',               # Hygromycin resistance (HygroR, Hygro)
    r'\bHygromycin\b',            # Full name
    r'\bNeoR?\b',                 # Neomycin/G418 resistance
    r'\bBlasticidinR?\b',         # Blasticidin resistance
    r'\bBsdR\b',                  # Blasticidin short form
    r'\bZeocinR?\b',              # Zeocin resistance
    r'\bBleoR\b',                 # Bleomycin/Zeocin resistance
]

# Tier 2: Reporter/fluorescent genes
REPORTER_PATTERNS = [
    r'\b[em]?GFP\b',              # GFP, eGFP, mGFP
    r'\bmCherry\b',               # mCherry
    r'\b[tdm]?RFP\b',             # RFP, tRFP, dRFP, mRFP
    r'\bdsRed\b',                 # dsRed
    r'\bYFP\b',                   # YFP
    r'\bCFP\b',                   # CFP
    r'\bluciferase\b',            # Luciferase
    r'\blacZ\b',                  # Beta-galactosidase
    r'\bEGFP\b',                  # Explicit EGFP
]

# Exclude patterns (affinity tags - not useful for contamination detection)
EXCLUDE_PATTERNS = [
    r'\bHis[-_]?tag\b',
    r'\b\d+His\b',                # 6His, etc.
    r'\bFLAG\b',
    r'\bMyc\b',
    r'\bV5\b',
    r'\bHA\b',
]

ALL_PATTERNS = BACTERIAL_RESISTANCE_PATTERNS + MAMMALIAN_RESISTANCE_PATTERNS + REPORTER_PATTERNS

def matches_resistance_gene(text: str) -> bool:
    """Check if text matches any resistance/reporter gene pattern."""
    if not text:
        return False

    # Check exclusions first
    for exclude_pattern in EXCLUDE_PATTERNS:
        if re.search(exclude_pattern, text, re.IGNORECASE):
            return False

    # Check resistance/reporter patterns
    for pattern in ALL_PATTERNS:
        if re.search(pattern, text, re.IGNORECASE):
            return True

    return False

def identify_gene_name(label: str, product: str, gene: str) -> str | None:
    """Identify gene name from qualifiers, return GenBank label verbatim."""
    # Priority: label > product > gene
    for qualifier_value in [label, product, gene]:
        if matches_resistance_gene(qualifier_value):
            return qualifier_value  # Return verbatim (original GenBank name)
    return None
```

### Per-Gene Coverage Computation (Reusing Phase 9)
```python
# Source: Phase 9 coverage_metrics.py (existing infrastructure)
from plasmicheck.scripts.coverage_metrics import get_depth_array, compute_coverage_stats

def compute_gene_coverage(
    bam_path: str,
    contig: str,
    gene_start: int,
    gene_end: int,
    breadth_thresholds: list[int]
) -> dict[str, float]:
    """Compute coverage metrics for a single gene region.

    Args:
        bam_path: Path to BAM alignment file
        contig: Reference contig name
        gene_start: 0-based inclusive start (from Biopython feature.location.start)
        gene_end: 0-based exclusive end (from Biopython feature.location.end)
        breadth_thresholds: Depth thresholds for breadth calculation

    Returns:
        Dict with keys: mean_depth, median_depth, breadth_1x, cv
    """
    # get_depth_array expects half-open interval [start, end)
    depth_array = get_depth_array(bam_path, contig, gene_start, gene_end)

    # Compute statistics using Phase 9 infrastructure
    stats = compute_coverage_stats(depth_array, breadth_thresholds)

    return stats
```

### Integration Point (compare_alignments.py)
```python
# Source: Phase 10 integration pattern
# Location: plasmicheck/scripts/compare_alignments.py, after coverage metrics computation

def compare_alignments(
    plasmid_bam: str,
    human_bam: str,
    output_basename: str,
    threshold: float = DEFAULT_THRESHOLD,
) -> None:
    # ... existing code ...

    # Phase 9: Compute region coverage metrics
    coverage_metrics, coverage_fallback = compute_region_coverage_metrics(
        plasmid_bam, insert_region, BREADTH_THRESHOLDS
    )

    # Phase 10: Detect resistance genes and compute per-gene coverage
    # Locate GenBank file (same directory as output_basename)
    genbank_file = locate_genbank_file(output_basename)

    if genbank_file:
        gene_metrics = detect_resistance_genes(
            genbank_file, plasmid_bam, insert_region, BREADTH_THRESHOLDS
        )
        resistance_genes_available = True
    else:
        gene_metrics = {}
        resistance_genes_available = False
        logging.warning("No GenBank file found, skipping resistance gene detection")

    # ... write summary.tsv ...
    with open(f"{output_basename}.summary.tsv", "w") as summary_file:
        # ... existing metrics ...

        # Resistance gene summary
        summary_file.write(f"ResistanceGenesAvailable\t{resistance_genes_available}\n")
        if resistance_genes_available:
            gene_names = list(gene_metrics.keys()) if gene_metrics else None
            gene_list = ",".join(gene_names) if gene_names else "None"
            summary_file.write(f"ResistanceGenes\t{gene_list}\n")

            # Per-gene metrics
            for gene_name, metrics in gene_metrics.items():
                sanitized = sanitize_gene_name(gene_name)
                summary_file.write(f"ResGene_{sanitized}_MeanDepth\t{metrics['mean_depth']:.2f}\n")
                summary_file.write(f"ResGene_{sanitized}_Breadth\t{metrics['breadth_1x']:.2f}\n")
                summary_file.write(f"ResGene_{sanitized}_Start\t{metrics['start']}\n")
                summary_file.write(f"ResGene_{sanitized}_End\t{metrics['end']}\n")
                summary_file.write(f"ResGene_{sanitized}_Region\t{metrics['region']}\n")
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual BLAST against CARD/ResFinder | GenBank annotation parsing | Ongoing (2024-2025) | Faster, no external dependencies, works offline |
| Single resistance marker (NeoR standard) | Multiple markers (BleoR, PuroR, HygroR) | 2021 research | BleoR shows ~10x higher expression than NeoR (relevant for detection sensitivity) |
| Beta-lactamase "bla" generic | Variant-specific (blaTEM, blaCTX-M) | Ongoing | Better specificity, but generic "bla" still common in plasmid annotations |
| LacZ for cloning selection | Fluorescent proteins (GFP, mCherry) | 2010s+ | Reporter genes now more common, both should be detected |

**Deprecated/outdated:**
- Text-based GenBank parsing: Use Biopython SeqIO (handles edge cases, maintained)
- Ampicillin-only resistance: Modern plasmids use diverse markers (Hygromycin, Puromycin, Blasticidin)

## Open Questions

Things that couldn't be fully resolved:

1. **Pattern completeness (90% recall target)**
   - What we know: Curated patterns cover common resistance genes (Amp, Kan, Puro, Hygro) and reporters (GFP, mCherry). 6 existing GenBank files in reference/plasmid/ contain AmpR, Hygromycin, EGFP.
   - What's unclear: Whether patterns achieve >=90% recall on existing files without manual validation. Edge cases like "FRT Hygromycin" (fusion construct) vs. "HygroR" (simple marker).
   - Recommendation: Test against all 6 GenBank files during implementation, log all detected genes, manually verify against GenBank file annotations. Adjust patterns to achieve 90% recall threshold.

2. **Config.json extension mechanism**
   - What we know: User should be able to add custom patterns for niche plasmids. Config.json structure needs "resistance_gene_patterns" section.
   - What's unclear: Merge strategy (append to built-in patterns vs. replace). Priority order if custom pattern conflicts with built-in.
   - Recommendation: Append custom patterns to built-in list (additive). Document in config.json with example pattern. Built-in patterns take precedence (check first).

3. **GenBank file location strategy**
   - What we know: xDNA files create .gb in output_subfolder (same directory as .fasta). Original .gb files are in reference/plasmid/.
   - What's unclear: Should gene detection use the intermediate .gb (guaranteed same coordinates as FASTA) or the original input file (may have different sanitization)?
   - Recommendation: Use intermediate .gb file in output_subfolder. Preserves alignment between gene coordinates and FASTA sequence after sanitization. Requires modifying convert_plasmidfile_to_fasta() to preserve .gb file.

4. **Insert region edge case handling**
   - What we know: Genes overlapping insert region should be excluded (experimental target). cDNA_positions.txt provides insert boundaries.
   - What's unclear: If a gene partially overlaps insert (e.g., 10% of bases), should it be excluded entirely or reported with caveat?
   - Recommendation: Exclude any gene with overlap (conservative). Partial overlaps likely indicate insert-proximal regulatory sequences (promoters), not independent contamination markers. Log excluded genes for debugging.

## Sources

### Primary (HIGH confidence)
- Biopython SeqIO documentation: https://widdowquinn.github.io/2018-03-06-ibioic/01-introduction/03-parsing.html
- Python re module documentation: https://docs.python.org/3/library/re.html
- Pysam count_coverage API: https://pysam.readthedocs.io/en/latest/api.html
- Existing codebase (Phase 9 coverage_metrics.py, compare_alignments.py)

### Secondary (MEDIUM confidence)
- Addgene resistance gene overview: https://blog.addgene.org/plasmids-101-everything-you-need-to-know-about-antibiotic-resistance-genes
- VectorBuilder selection marker guide: https://en.vectorbuilder.com/resources/faq/optimize-drug-selection-markers.html
- Reporter gene overview: https://en.vectorbuilder.com/resources/vector-academy/troubleshooting/reporter-gene.html
- Beta-lactamase variants (Nature Communications): https://www.nature.com/articles/s41467-024-49621-2
- Mammalian selection marker comparison (2021): https://pmc.ncbi.nlm.nih.gov/articles/PMC8258971/

### Tertiary (LOW confidence)
- Resistance gene naming conventions: https://academic.oup.com/jac/article/71/3/569/2363654 (notes inconsistency, no single standard)
- Circular plasmid join locations: https://github.com/biopython/biopython/issues/895 (Biopython GitHub issue)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Biopython and pysam are established, well-documented, already in use
- Architecture: HIGH - Phase 9 infrastructure proven, GenBank parsing is standard Biopython pattern
- Pitfalls: HIGH - Identified from codebase inspection (xDNA .gb deletion), Biopython documentation (qualifier lists, coordinates), existing GenBank files
- Pattern recall: MEDIUM - Curated patterns based on literature, but 90% recall needs validation against actual files

**Research date:** 2026-02-15
**Valid until:** 90 days (stable domain - Biopython API stable, resistance gene names change slowly)
