# Phase 15: Splice Junction Metadata - Research

**Researched:** 2026-02-16
**Domain:** BAM CIGAR parsing, Python dataclasses, file format extension
**Confidence:** HIGH

## Summary

This phase extends existing cDNA position extraction to include splice junction metadata from spliced alignment CIGAR strings. The technical implementation is straightforward: parse the primary alignment's CIGAR tuples (already available via pysam), count N operations (representing introns), and store additional fields in the existing human-readable text format.

**Architecture context:** PlasmiCheck already uses pysam extensively for BAM parsing, has an established text file format (cDNA_positions.txt), and a parser function (parse_insert_region) that downstream scripts depend on. The task is to extend this existing infrastructure without breaking backward compatibility.

**Primary recommendation:** Extend parse_insert_region() to return a dataclass (InsertRegionData) with optional fields for new metadata. This provides type safety, backward compatibility (old files work with None defaults), and clear API for downstream consumers.

## Standard Stack

The codebase already uses all required libraries. No new dependencies needed.

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pysam | >=0.22.1 | BAM file parsing, CIGAR access | Industry standard for SAM/BAM manipulation in Python |
| Python dataclasses | stdlib (3.10+) | Structured data for parser return | Type-safe alternative to tuples, mypy-compatible |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pandas | >=2.2.2 | summary.tsv and report table rendering | Already used for report generation |
| Jinja2 | >=3.0.0 | HTML template rendering | Already used for report templates |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| dataclass | TypedDict | TypedDict doesn't provide constructor validation or defaults |
| dataclass | NamedTuple | NamedTuple is immutable, harder to extend with optional fields |
| pysam | custom BAM parser | Reinventing the wheel, pysam handles edge cases |

**Installation:**
No new packages required. All dependencies already in pyproject.toml.

## Architecture Patterns

### Current Implementation Structure
```
plasmicheck/scripts/
├── spliced_alignment.py         # Performs spliced alignment, extracts cDNA positions
├── compare_alignments.py        # Reads cDNA_positions.txt via parse_insert_region()
├── generate_report.py           # Per-sample HTML reports (uses summary.tsv)
├── generate_summary_reports.py  # Multi-sample aggregation reports
└── run_pipeline.py              # Orchestrates full workflow
```

### Pattern 1: File Format Extension (Backward Compatible)
**What:** Add new lines to existing text file while preserving parse-ability of old files
**When to use:** Extending established file formats that existing tools depend on
**Example:**
```python
# Old format (3 lines):
cDNA start position in plasmid: 988
cDNA end position in plasmid: 2249
INSERT_REGION: (988, 2249)

# New format (11 lines):
cDNA start position in plasmid: 988
cDNA end position in plasmid: 2249
INSERT_REGION: (988, 2249)
SpliceJunctionCount: 1
ExonCount: 2
ExonLengths: 630,632
IntronLengths: 145
Strand: +
Chromosome: chr12
GenomicStart: 47350
GenomicEnd: 48757

# Parser handles both:
@dataclass
class InsertRegionData:
    start: int
    end: int
    junction_count: int | None = None  # None for old files
    exon_count: int | None = None
    # ... more optional fields
```

### Pattern 2: Dataclass Parser Return (Type-Safe API Evolution)
**What:** Replace tuple return with dataclass to enable adding fields without breaking callers
**When to use:** Functions that return multiple values that may grow over time
**Example:**
```python
# Old API (tuple):
def parse_insert_region(file_path: str) -> tuple[int, int]:
    # ... parsing logic
    return (start, end)

# New API (dataclass, backward compatible via unpacking):
def parse_insert_region(file_path: str) -> InsertRegionData:
    # ... parsing logic
    return InsertRegionData(start=start, end=end, junction_count=junc_count)

# Old code still works:
start, end = parse_insert_region(file_path)  # tuple unpacking

# New code can access all fields:
data = parse_insert_region(file_path)
if data.junction_count is not None:
    print(f"Junctions: {data.junction_count}")
```

### Pattern 3: Primary-Only Filtering for Metadata
**What:** Use primary alignment (flags 0 or 16) for structural metadata, all alignments for reference extraction
**When to use:** When alignment metadata serves different purposes (accuracy vs completeness)
**Example:**
```python
# Filter for primary alignment only:
for read in bamfile.fetch():
    if read.is_secondary or read.is_supplementary:
        continue  # Skip non-primary alignments

    # Now process primary alignment's CIGAR
    junction_count = sum(1 for op, length in read.cigartuples if op == 3)
```

**Rationale:** PlasmiCheck already uses this pattern in compare_alignments.py (_best_read function picks primary alignment). Primary alignment represents the aligner's best judgment of the true mapping, while supplementary/secondary alignments are alternative mappings that shouldn't influence structural metadata.

### Pattern 4: Summary TSV as Category-Value Rows
**What:** Store metrics as two-column TSV with Category and Value/Count columns
**When to use:** Heterogeneous metrics (strings, floats, booleans) that need to be easily parseable
**Example (current format):**
```tsv
Category	Count
Plasmid	1500
Human	800
Tied	50
Verdict	Sample is contaminated with plasmid DNA
Ratio	1.875
MeanDepthInsert	45.2
```

**Extension for Phase 15:**
```tsv
SpliceJunctionCount	1
ExonCount	2
ExonLengths	630,632
IntronLengths	145
```

**How it's read:** generate_report.py uses pandas to read summary.tsv, filters by Category to extract specific metrics (see `_get_metric` function, line 615-617).

### Anti-Patterns to Avoid
- **Breaking tuple unpacking:** If parse_insert_region returns a dataclass, ensure it's iterable or provide __iter__ for backward compatibility
- **Reading supplementary alignments:** These are alternative mappings, not independent observations. Counting N operations in supplementary alignments would inflate junction counts
- **Forgetting NA values:** Unmapped plasmids should write "NA" for all new fields, not skip writing cDNA_positions.txt entirely
- **Hardcoding field positions:** When extending text file format, parse by key prefix ("SpliceJunctionCount: "), not by line number

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| BAM file parsing | Custom BGZF reader | pysam.AlignmentFile | Handles compression, indexes, malformed files, platform quirks |
| CIGAR parsing | String splitting | read.cigartuples | Pysam returns pre-parsed tuples, handles edge cases (empty CIGAR, None) |
| Dataclass iteration | Custom __iter__ | dataclasses.astuple() | Standard library function for tuple conversion |
| HTML table truncation | JavaScript string slicing | Jinja2 filters + CSS tooltips | Template-level logic keeps presentation separate from data |

**Key insight:** pysam already provides parsed CIGAR tuples as a list of (operation, length) pairs. Operation codes are integers (0=M, 3=N). Don't parse cigarstring manually.

## Common Pitfalls

### Pitfall 1: Forgetting Secondary/Supplementary Filtering
**What goes wrong:** Junction count includes alignments from alternative mappings, inflating the count
**Why it happens:** BAM files contain multiple alignment records per read (primary + secondary/supplementary)
**How to avoid:** Always check `read.is_secondary or read.is_supplementary` before processing alignment metadata
**Warning signs:** Junction count is 3x higher than expected; debug output shows multiple CIGAR strings per plasmid

**Code example:**
```python
# WRONG: processes all alignments
for read in bamfile.fetch():
    junction_count = sum(1 for op, length in read.cigartuples if op == 3)

# CORRECT: primary alignment only
for read in bamfile.fetch():
    if read.is_secondary or read.is_supplementary:
        continue
    junction_count = sum(1 for op, length in read.cigartuples if op == 3)
```

### Pitfall 2: Breaking Tuple Unpacking in Existing Code
**What goes wrong:** Changing parse_insert_region return type from tuple to dataclass breaks all existing callers
**Why it happens:** Python doesn't auto-convert dataclasses to tuples
**How to avoid:** Make dataclass iterable by defining `__iter__` that yields (start, end) in order, or use dataclasses.astuple() wrapper
**Warning signs:** TypeError: cannot unpack non-iterable InsertRegionData object

**Solution:**
```python
from dataclasses import astuple, dataclass

@dataclass
class InsertRegionData:
    start: int
    end: int
    junction_count: int | None = None

    def __iter__(self):
        """Support tuple unpacking for backward compatibility."""
        return iter((self.start, self.end))

# Now both work:
start, end = parse_insert_region(file)  # Old code
data = parse_insert_region(file)        # New code
```

### Pitfall 3: CIGAR None Handling
**What goes wrong:** read.cigartuples is None for unmapped reads, causing AttributeError when iterating
**Why it happens:** Unmapped reads have no CIGAR string
**How to avoid:** Check if read.is_unmapped or read.cigartuples is None before processing
**Warning signs:** AttributeError: 'NoneType' object is not iterable

**Code example:**
```python
# WRONG:
for op, length in read.cigartuples:  # Crashes if None
    pass

# CORRECT:
if read.cigartuples is None:
    return None  # Handle unmapped case
for op, length in read.cigartuples:
    pass
```

### Pitfall 4: Comma-Separated List Formatting
**What goes wrong:** Exon lengths like [630, 632, 544] written as "630632544" (no separators) or "[630, 632, 544]" (Python repr)
**Why it happens:** Forgetting to join with comma, or using str() on list
**How to avoid:** Use ','.join(map(str, lengths))
**Warning signs:** TSV columns contain Python list literals; Excel can't parse numbers

**Code example:**
```python
# WRONG:
exon_lengths = [630, 632, 544]
f.write(f"ExonLengths: {exon_lengths}\n")  # Outputs: [630, 632, 544]

# CORRECT:
f.write(f"ExonLengths: {','.join(map(str, exon_lengths))}\n")  # Outputs: 630,632,544
```

### Pitfall 5: Strand Inference from Flag
**What goes wrong:** Assuming read.is_reverse is the strand; this is the *read's* orientation, not the *transcript's* strand
**Why it happens:** Confusion between read mapping orientation and biological strand
**How to avoid:** For cDNA alignments, read.is_reverse (flag 16) indicates the strand of the alignment. For plus strand, flag 16 is unset (flag 0 or 16). Check flag 16 bit: strand = "-" if read.flag & 16 else "+"
**Warning signs:** All strands reported as plus when some should be minus

**Code example:**
```python
# CORRECT strand from SAM flag:
strand = "-" if read.is_reverse else "+"
# Or equivalently:
strand = "-" if read.flag & 16 else "+"
```

## Code Examples

Verified patterns from codebase inspection:

### Parsing Primary Alignment CIGAR
```python
# Source: spliced_alignment.py lines 106-161 (extract_plasmid_cdna_positions)
# Modified to add junction counting

import pysam

def extract_junction_metadata(bam_file: str) -> dict[str, any]:
    """Extract splice junction metadata from primary alignment only."""
    bamfile = pysam.AlignmentFile(bam_file, "rb")

    junction_count = None
    exon_count = None
    exon_lengths = []
    intron_lengths = []
    strand = None
    chromosome = None
    genomic_start = None
    genomic_end = None

    for read in bamfile.fetch():
        # Skip secondary and supplementary alignments
        if read.is_secondary or read.is_supplementary:
            continue

        if read.is_unmapped or read.cigartuples is None:
            continue

        # Process primary alignment
        junction_count = sum(1 for op, length in read.cigartuples if op == 3)
        exon_count = junction_count + 1

        # Extract exon and intron lengths
        for op, length in read.cigartuples:
            if op == 0:  # Match/mismatch (exon)
                exon_lengths.append(length)
            elif op == 3:  # Skipped region (intron)
                intron_lengths.append(length)

        # Extract genomic coordinates
        strand = "-" if read.is_reverse else "+"
        chromosome = read.reference_name
        genomic_start = read.reference_start
        genomic_end = read.reference_end

        # Only process first primary alignment
        break

    bamfile.close()

    if junction_count is None:
        # No mapped primary alignment found
        return {
            "junction_count": None,
            "exon_count": None,
            "exon_lengths": [],
            "intron_lengths": [],
            "strand": None,
            "chromosome": None,
            "genomic_start": None,
            "genomic_end": None,
        }

    return {
        "junction_count": junction_count,
        "exon_count": exon_count,
        "exon_lengths": exon_lengths,
        "intron_lengths": intron_lengths,
        "strand": strand,
        "chromosome": chromosome,
        "genomic_start": genomic_start,
        "genomic_end": genomic_end,
    }
```

### Writing Extended cDNA Positions File
```python
# Source: spliced_alignment.py lines 155-160 (current format)
# Extended with new fields

def write_cdna_positions(output_file: str, metadata: dict) -> None:
    """Write cDNA positions and junction metadata to text file.

    Always generates file, even for unmapped plasmids (writes NA values).
    """
    with open(output_file, "w") as f:
        # Existing fields (always present)
        start = metadata.get("cdna_start")
        end = metadata.get("cdna_end")

        if start is not None and end is not None:
            f.write(f"cDNA start position in plasmid: {start}\n")
            f.write(f"cDNA end position in plasmid: {end}\n")
            f.write(f"INSERT_REGION: ({start}, {end})\n")
        else:
            # Unmapped plasmid
            f.write("cDNA start position in plasmid: NA\n")
            f.write("cDNA end position in plasmid: NA\n")
            f.write("INSERT_REGION: NA\n")

        # New fields (Phase 15)
        junc = metadata.get("junction_count")
        exon = metadata.get("exon_count")
        exon_lens = metadata.get("exon_lengths", [])
        intron_lens = metadata.get("intron_lengths", [])
        strand = metadata.get("strand")
        chrom = metadata.get("chromosome")
        g_start = metadata.get("genomic_start")
        g_end = metadata.get("genomic_end")

        f.write(f"SpliceJunctionCount: {junc if junc is not None else 'NA'}\n")
        f.write(f"ExonCount: {exon if exon is not None else 'NA'}\n")
        f.write(f"ExonLengths: {','.join(map(str, exon_lens)) if exon_lens else 'NA'}\n")
        f.write(f"IntronLengths: {','.join(map(str, intron_lens)) if intron_lens else 'NA'}\n")
        f.write(f"Strand: {strand if strand is not None else 'NA'}\n")
        f.write(f"Chromosome: {chrom if chrom is not None else 'NA'}\n")
        f.write(f"GenomicStart: {g_start if g_start is not None else 'NA'}\n")
        f.write(f"GenomicEnd: {g_end if g_end is not None else 'NA'}\n")
```

### Backward-Compatible Parser with Dataclass
```python
# Source: compare_alignments.py lines 58-72 (current parse_insert_region)
# Extended to return dataclass with optional fields

from dataclasses import dataclass

@dataclass
class InsertRegionData:
    """Parsed cDNA position and junction metadata.

    Supports tuple unpacking for backward compatibility:
        start, end = parse_insert_region(file)

    New code can access all fields:
        data = parse_insert_region(file)
        if data.junction_count is not None:
            print(f"Junctions: {data.junction_count}")
    """
    start: int
    end: int
    junction_count: int | None = None
    exon_count: int | None = None
    exon_lengths: list[int] | None = None
    intron_lengths: list[int] | None = None
    strand: str | None = None
    chromosome: str | None = None
    genomic_start: int | None = None
    genomic_end: int | None = None

    def __iter__(self):
        """Support tuple unpacking: start, end = parse_insert_region(file)"""
        return iter((self.start, self.end))

def parse_insert_region(cdna_positions_file: str) -> InsertRegionData:
    """Parse cDNA_positions.txt file to extract insert region and metadata.

    Returns:
        InsertRegionData with start/end and optional junction metadata.
        New fields are None if file format is old (pre-Phase 15).
    """
    if not os.path.exists(cdna_positions_file):
        raise FileNotFoundError(
            f"Expected cDNA_positions.txt at {cdna_positions_file}"
        )

    with open(cdna_positions_file) as f:
        lines = f.readlines()

    # Parse required fields (lines 0-2)
    cdna_start_str = lines[0].split(": ")[1].strip()
    cdna_end_str = lines[1].split(": ")[1].strip()

    # Handle NA values for unmapped plasmids
    if cdna_start_str == "NA" or cdna_end_str == "NA":
        start = 0
        end = 0
    else:
        start = int(cdna_start_str)
        end = int(cdna_end_str)

    # Parse optional fields (lines 3+) with defaults
    def get_field(prefix: str) -> str | None:
        for line in lines[3:]:
            if line.startswith(prefix):
                value = line.split(": ", 1)[1].strip()
                return None if value == "NA" else value
        return None

    junction_count = None
    exon_count = None
    exon_lengths = None
    intron_lengths = None
    strand = None
    chromosome = None
    genomic_start = None
    genomic_end = None

    if junc_str := get_field("SpliceJunctionCount"):
        junction_count = int(junc_str)

    if exon_str := get_field("ExonCount"):
        exon_count = int(exon_str)

    if exon_lens_str := get_field("ExonLengths"):
        exon_lengths = [int(x) for x in exon_lens_str.split(",")]

    if intron_lens_str := get_field("IntronLengths"):
        intron_lengths = [int(x) for x in intron_lens_str.split(",")]

    strand = get_field("Strand")
    chromosome = get_field("Chromosome")

    if g_start_str := get_field("GenomicStart"):
        genomic_start = int(g_start_str)

    if g_end_str := get_field("GenomicEnd"):
        genomic_end = int(g_end_str)

    return InsertRegionData(
        start=start,
        end=end,
        junction_count=junction_count,
        exon_count=exon_count,
        exon_lengths=exon_lengths,
        intron_lengths=intron_lengths,
        strand=strand,
        chromosome=chromosome,
        genomic_start=genomic_start,
        genomic_end=genomic_end,
    )
```

### Adding Fields to Summary TSV
```python
# Source: compare_alignments.py lines 471-527 (current summary writing)
# Extended with junction fields

def write_summary_with_junctions(
    summary_file_path: str,
    assigned_counts: dict,
    verdict: str,
    ratio: float,
    junction_metadata: dict,
) -> None:
    """Write summary.tsv with junction metadata fields."""
    with open(summary_file_path, "w") as f:
        f.write("Category\tCount\n")

        # Existing fields
        for category, count in assigned_counts.items():
            f.write(f"{category}\t{count}\n")
        f.write(f"Verdict\t{verdict}\n")
        f.write(f"Ratio\t{ratio}\n")

        # Junction metadata fields (new in Phase 15)
        junc = junction_metadata.get("junction_count")
        exon = junction_metadata.get("exon_count")
        exon_lens = junction_metadata.get("exon_lengths", [])
        intron_lens = junction_metadata.get("intron_lengths", [])

        f.write(f"SpliceJunctionCount\t{junc if junc is not None else 'NA'}\n")
        f.write(f"ExonCount\t{exon if exon is not None else 'NA'}\n")
        f.write(f"ExonLengths\t{','.join(map(str, exon_lens)) if exon_lens else 'NA'}\n")
        f.write(f"IntronLengths\t{','.join(map(str, intron_lens)) if intron_lens else 'NA'}\n")
```

### HTML Table Truncation with Tooltip
```python
# Jinja2 template filter for truncating long comma-separated lists
# Source: Based on existing template patterns in report_template.html

def truncate_list(value: str, max_items: int = 5) -> tuple[str, str]:
    """Truncate comma-separated list for display.

    Returns:
        (display_text, tooltip_text) tuple for HTML rendering
    """
    if not value or value == "NA":
        return (value, value)

    items = value.split(",")
    if len(items) <= max_items:
        return (value, value)

    display = ",".join(items[:max_items]) + ",..."
    tooltip = value  # Full list for hover
    return (display, tooltip)

# Usage in Jinja2 template:
# <td title="{{ full_list }}">{{ truncated_list }}</td>
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Tuple returns | Dataclass returns | Python 3.7+ | Type safety, named fields, optional fields without breaking API |
| Manual CIGAR parsing | pysam cigartuples | Always in project | Native (operation, length) tuples, no string parsing |
| All alignments | Primary only | Phase 15 decision | Accurate structural metadata (1 junction count per plasmid, not 3x) |
| Hardcoded line parsing | Key-based parsing | Phase 15 design | Forward-compatible file format extension |

**Deprecated/outdated:**
- Manual SAM flag bit manipulation: Use pysam properties (read.is_reverse, read.is_secondary, read.is_supplementary) instead of flag & 0x10 checks
- String-based CIGAR parsing: Use read.cigartuples, don't parse read.cigarstring
- Returning tuples from parsers: Use dataclasses for extensibility

## Open Questions

None. All technical questions resolved:

1. **How to filter primary alignment?** → `if read.is_secondary or read.is_supplementary: continue`
2. **How to parse CIGAR?** → `read.cigartuples` provides pre-parsed list
3. **How to maintain backward compatibility?** → Dataclass with `__iter__` for tuple unpacking
4. **How to handle unmapped plasmids?** → Write cDNA_positions.txt with "NA" values for all fields
5. **How to store in summary.tsv?** → Add as Category/Count rows, parsed by existing `_get_metric` pattern
6. **How to store comma-separated lists?** → `','.join(map(str, lengths))`
7. **How to truncate in HTML?** → Jinja2 filter + tooltip (show first 5, tooltip shows all)
8. **Does cache store cDNA_positions.txt?** → Yes, cache_manager.py puts/gets it as "cdna_path" (line 186, 220, 238)

## Sources

### Primary (HIGH confidence)
- Codebase inspection: spliced_alignment.py (lines 106-161) — extract_plasmid_cdna_positions implementation
- Codebase inspection: compare_alignments.py (lines 58-72) — parse_insert_region implementation
- Codebase inspection: generate_report.py (lines 608-630) — summary.tsv parsing with _get_metric
- Codebase inspection: generate_summary_reports.py (lines 91-138) — multi-sample table aggregation
- Codebase inspection: cache_manager.py (lines 166-249) — spliced alignment caching (includes cdna_positions)
- Codebase inspection: output_paths.py (lines 140-143) — cdna_positions.txt path definition
- pysam documentation: AlignmentFile.fetch, AlignedSegment properties (cigartuples, is_reverse, reference_name, reference_start, reference_end)

### Secondary (MEDIUM confidence)
- SAM/BAM specification: CIGAR operation codes (0=M, 3=N), flag bits (0x100=secondary, 0x800=supplementary, 0x10=reverse)
- Python dataclasses documentation: __iter__ for unpacking, optional fields with defaults

### Tertiary (LOW confidence)
None — all findings verified with codebase inspection

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All libraries already in pyproject.toml, verified with Read tool
- Architecture: HIGH - Inspected all relevant modules (spliced_alignment, compare_alignments, generate_report, cache_manager)
- Pitfalls: HIGH - Derived from pysam API inspection and existing codebase patterns
- File formats: HIGH - Examined real cDNA_positions.txt files in output directory

**Research date:** 2026-02-16
**Valid until:** 60 days (stable domain — pysam API and project architecture unlikely to change)
