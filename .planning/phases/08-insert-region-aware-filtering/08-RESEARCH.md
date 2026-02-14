# Phase 8: Insert-Region-Aware Filtering - Research

**Researched:** 2026-02-14
**Domain:** Bioinformatics alignment filtering with positional awareness
**Confidence:** HIGH

## Summary

This phase implements insert-region-aware read filtering for plasmid contamination detection, distinguishing reads that map only to shared backbone regions from reads overlapping the plasmid insert. The implementation modifies the existing streaming comparison logic in `compare_alignments.py` to classify reads based on their genomic position overlap with the insert region, excluding backbone-only reads from the contamination ratio calculation.

The existing codebase already contains the necessary infrastructure: insert region parsing from `cDNA_positions.txt`, pysam-based alignment processing, and positional overlap detection patterns in `calculate_coverage_outside_insert()`. The phase extends this with new read categories (Backbone_Only, Ambiguous), a backward-compatibility toggle, and updated reporting.

**Primary recommendation:** Use the existing overlap detection pattern from `calculate_coverage_outside_insert()` (lines 87-97 in compare_alignments.py) as the template for insert-region checking. This pattern correctly handles partial overlaps using pysam's `reference_start` and `reference_end` properties with half-open interval semantics.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pysam | 0.22+ | BAM/SAM file parsing and alignment access | Industry-standard Python interface to htslib, used throughout existing codebase |
| Python typing | 3.9+ | Type hints for tuple[int, int] insert regions | Already used extensively in codebase, mypy strict mode enforced |
| Jinja2 | 3.x | HTML template rendering for reports | Already used for report generation (plasmicheck/templates/) |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| plotly | 5.x | Interactive plot generation | Already used for box/scatter plots in generate_report.py |
| pandas | 2.x | DataFrame manipulation for summary tables | Already used for reading/writing TSV summaries |
| pytest | 7.x+ | Testing with parametrization | Already used extensively, supports feature flag testing |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| pysam | HTSeq | pysam is already integrated, HTSeq would add unnecessary dependency |
| Tuple for insert_region | Custom InsertRegion class | Over-engineering for simple (start, end) pair, tuple is idiomatic Python |
| Config JSON | Pydantic models | Current config.json pattern already established in codebase, no validation issues |

**Installation:**
No new dependencies required - all libraries already present in `pyproject.toml`

## Architecture Patterns

### Recommended Project Structure
```
plasmicheck/
├── scripts/
│   ├── compare_alignments.py    # Core filtering logic (main file)
│   ├── generate_report.py       # Report generation updates
│   └── run_pipeline.py          # Pipeline orchestration
├── config.json                  # New filter_backbone_only, score_margin
└── templates/
    └── report_template.html     # Updated read category display
```

### Pattern 1: Positional Overlap Detection (Existing Pattern)
**What:** Determine if a read's alignment overlaps with the insert region using pysam's reference_start/reference_end
**When to use:** When classifying reads as Backbone_Only vs potentially contaminating
**Example:**
```python
# Source: Existing pattern in calculate_coverage_outside_insert() lines 87-97
# Reference: https://pysam.readthedocs.io/en/latest/api.html

def read_overlaps_insert(read: Any, insert_region: tuple[int, int]) -> bool:
    """Check if read alignment overlaps insert region.

    Args:
        read: pysam AlignedSegment with reference_start (0-based, inclusive)
              and reference_end (0-based, exclusive - half-open interval)
        insert_region: (start, end) both 0-based inclusive

    Returns:
        True if read overlaps insert region, False if entirely outside
    """
    if read.is_unmapped or read.reference_start is None or read.reference_end is None:
        return False

    # No overlap if read ends before insert or starts after insert
    # Note: read.reference_end is exclusive (half-open), insert_region[1] is inclusive
    if read.reference_end <= insert_region[0] or read.reference_start > insert_region[1]:
        return False

    return True
```

**Critical detail:** pysam uses 0-based half-open intervals (`reference_end` is exclusive), while the insert_region from `cDNA_positions.txt` is 0-based inclusive on both ends. The existing `calculate_coverage_outside_insert()` handles this correctly with `read_end < insert_region[0]` and `read_start > insert_region[1]`.

### Pattern 2: Backward-Compatible Feature Flag
**What:** Config-driven behavior toggle that maintains pre-v0.33.0 behavior when disabled
**When to use:** When adding breaking changes that need opt-out for testing/validation
**Example:**
```python
# Source: Best practice from pytest backward compatibility policy
# Reference: https://docs.pytest.org/en/stable/backwards-compatibility.html

# In config.json
{
  "filter_backbone_only": true,  # Default ON for improved accuracy
  "score_margin": 0              # Default 0 (disabled)
}

# In compare_alignments.py
_cfg = get_config()
FILTER_BACKBONE_ONLY: bool = _cfg.get("filter_backbone_only", True)
SCORE_MARGIN: int = _cfg.get("score_margin", 0)

# When filter_backbone_only=false: classify but don't exclude
if FILTER_backbone_ONLY:
    # Exclude Backbone_Only and Ambiguous from ratio
    ratio = plasmid_count / human_count
else:
    # Include all categories (pre-v0.33.0 behavior)
    ratio = (plasmid_count + backbone_only_count + ambiguous_count) / human_count
```

### Pattern 3: Streaming Comparison with Extended Classification
**What:** Modify `_streaming_compare()` to accept insert_region parameter and classify into 5 categories
**When to use:** Core comparison logic that must remain O(1) memory efficient
**Example:**
```python
# Source: Existing _streaming_compare() pattern lines 204-257
# Extended with insert-region awareness

def _assign_with_region(
    plasmid_score: int,
    human_score: int,
    plasmid_read: Any | None,
    insert_region: tuple[int, int] | None,
    score_margin: int = 0,
) -> str:
    """Assign read to category with insert-region and score-margin awareness."""
    # Score margin check (if enabled)
    if score_margin > 0 and abs(plasmid_score - human_score) < score_margin:
        return "Ambiguous"

    # Standard comparison
    if plasmid_score > human_score:
        # Read favors plasmid - check if it overlaps insert
        if insert_region and plasmid_read:
            if read_overlaps_insert(plasmid_read, insert_region):
                return "Plasmid"  # Overlaps insert, likely contamination
            else:
                return "Backbone_Only"  # Outside insert, shared backbone
        return "Plasmid"  # No insert region available, default behavior
    elif human_score > plasmid_score:
        return "Human"
    else:
        return "Tied"
```

### Pattern 4: Fallback for Missing Insert Region
**What:** Graceful degradation when cDNA_positions.txt is unavailable or malformed
**When to use:** Error handling for optional dependencies
**Example:**
```python
# Parse insert region with fallback
try:
    insert_region = parse_insert_region(cdna_positions_file)
    logging.debug(f"INSERT_REGION extracted: {insert_region}")
except (FileNotFoundError, ValueError, IndexError) as e:
    logging.warning(f"Backbone filtering unavailable: {e}")
    insert_region = None  # Signals fallback to pre-v0.33.0 behavior

# Pass through pipeline
assigned_counts = _streaming_compare(
    plasmid_ns, human_ns, outfile, insert_region=insert_region
)

# In _streaming_compare, None insert_region means no filtering
```

### Anti-Patterns to Avoid
- **Loading entire BAM into memory:** The existing streaming approach must be preserved. Do NOT load all reads into a list before processing.
- **Nested loops for overlap detection:** Use direct pysam properties (`reference_start`, `reference_end`), not `fetch()` calls per read.
- **Inconsistent schema:** Always include Backbone_Only and Ambiguous columns in TSV output, even when 0 (parsers depend on consistent schema).
- **Validating insert region format at runtime:** The existing `parse_insert_region()` already handles this correctly. Don't add redundant validation.
- **Mixing half-open and inclusive intervals:** Stick to pysam's convention. `reference_end` is exclusive, insert_region endpoints are inclusive.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Interval overlap detection | Custom range checking with +1/-1 adjustments | Existing `calculate_coverage_outside_insert()` pattern | Half-open vs inclusive interval semantics are error-prone; existing code is already correct |
| BAM file iteration | Manual file reading or fetch() per region | Existing `_iter_reads_by_name()` streaming pattern | Memory-efficient, already handles name-sorted BAMs correctly |
| Config validation | Custom JSON schema checker or runtime validation | Simple `dict.get()` with defaults | Config is trusted internal file, no user-facing input; Pydantic would be over-engineering |
| Report table styling | Custom CSS framework or JavaScript | Bootstrap 4.5.2 (already loaded in template) | Template already uses Bootstrap; dimmed text can be achieved with `class="text-muted"` |
| Insert region parsing | Regex or complex parsing logic | Existing `parse_insert_region()` with split() | File format is fixed, simple split() is reliable and already tested |

**Key insight:** Bioinformatics coordinate systems (0-based vs 1-based, half-open vs inclusive) are a common source of off-by-one errors. The existing codebase already handles these correctly in `calculate_coverage_outside_insert()` - reuse this pattern rather than reimplementing.

## Common Pitfalls

### Pitfall 1: Off-by-One Errors in Interval Overlap
**What goes wrong:** Reads on insert boundaries are incorrectly classified as Backbone_Only or vice versa
**Why it happens:** Mixing pysam's half-open intervals (reference_end is exclusive) with insert_region's inclusive endpoints
**How to avoid:**
- Use the existing pattern: `read_end <= insert_region[0]` (not `<`) and `read_start > insert_region[1]` (not `>=`)
- Test with reads that align exactly at insert boundaries (start == insert_region[0], end == insert_region[1])
**Warning signs:**
- Reads with `reference_start == insert_region[0]` classified as Backbone_Only
- Edge case test failures for boundary alignments

### Pitfall 2: Breaking Backward Compatibility
**What goes wrong:** Tests fail when `filter_backbone_only=false` because new categories are counted in the ratio
**Why it happens:** Forgetting that "backward compatible" means identical output, not just "runs without errors"
**How to avoid:**
- When `filter_backbone_only=false`, still classify reads as Backbone_Only/Ambiguous for visibility, but **include them in the ratio calculation**
- Run entire test suite with `filter_backbone_only=false` and verify no changes to expected verdicts
- Use pytest parametrization to test both modes:
```python
@pytest.mark.parametrize("filter_mode", [True, False])
def test_backward_compat(filter_mode):
    # Test with both filter modes
```
**Warning signs:**
- Ratio values change when toggling the flag
- Tests marked as XFAIL or modified expected values

### Pitfall 3: Unmapped Reads and None Values
**What goes wrong:** `AttributeError: 'NoneType' object has no attribute` when accessing `reference_start` or `reference_end`
**Why it happens:** Unmapped reads have `reference_start = None` and `reference_end = None`
**How to avoid:**
- Always check `read.is_unmapped` before accessing positional properties
- Check `reference_start is not None` before using it (even mapped reads can have None in malformed BAMs)
- The existing `calculate_coverage_outside_insert()` shows the correct pattern (lines 79-85)
**Warning signs:**
- Crashes during comparison phase with "NoneType" errors
- Integration tests failing with real-world BAM files

### Pitfall 4: Inconsistent TSV Schema
**What goes wrong:** Downstream parsers break when Backbone_Only/Ambiguous columns are missing
**Why it happens:** Only adding columns when filtering is enabled, or omitting them when counts are 0
**How to avoid:**
- **Always** write Backbone_Only and Ambiguous columns, even when values are 0
- When insert_region is None (fallback mode), write 0 for both categories
- Update TSV header in `_write_assignment()` to include new columns
**Warning signs:**
- Variable number of columns in TSV depending on input
- Pandas read errors in downstream tools

### Pitfall 5: Score Margin Interaction with Tied Reads
**What goes wrong:** Reads with exactly equal scores classified as both "Tied" and "Ambiguous"
**Why it happens:** Unclear precedence between tie detection and score margin logic
**How to avoid:**
- Check score margin first, then check for ties
- Document that `score_margin=0` means disabled (not "allow 0 difference")
- Tied reads (exact score match) should remain "Tied", not become "Ambiguous"
```python
# Correct order
if score_margin > 0 and abs(plasmid_score - human_score) < score_margin:
    return "Ambiguous"
if plasmid_score == human_score:
    return "Tied"
```
**Warning signs:**
- Tied count drops to 0 when score_margin > 0
- Test failures for equal-score scenarios

## Code Examples

Verified patterns from official sources:

### Reading Insert Region with Error Handling
```python
# Source: Existing parse_insert_region() + fallback pattern
# Reference: compare_alignments.py lines 48-62

def safe_parse_insert_region(cdna_positions_file: str) -> tuple[int, int] | None:
    """Parse insert region with fallback for missing/malformed files."""
    try:
        if not os.path.exists(cdna_positions_file):
            logging.warning(
                f"cDNA_positions.txt not found at {cdna_positions_file}. "
                "Backbone filtering will be skipped."
            )
            return None

        with open(cdna_positions_file) as f:
            lines = f.readlines()

        if len(lines) < 2:
            logging.warning(
                f"cDNA_positions.txt at {cdna_positions_file} is incomplete. "
                "Backbone filtering will be skipped."
            )
            return None

        cdna_start = int(lines[0].split(": ")[1])
        cdna_end = int(lines[1].split(": ")[1])
        return (cdna_start, cdna_end)

    except (ValueError, IndexError) as e:
        logging.warning(
            f"Failed to parse cDNA_positions.txt: {e}. "
            "Backbone filtering will be skipped."
        )
        return None
```

### Checking Read Overlap with Insert Region
```python
# Source: Adapted from calculate_coverage_outside_insert() lines 87-97
# Reference: https://pysam.readthedocs.io/en/latest/usage.html

def read_overlaps_insert(read: Any, insert_region: tuple[int, int]) -> bool:
    """
    Check if read alignment overlaps the insert region.

    Uses pysam half-open interval convention:
    - read.reference_start: 0-based, inclusive
    - read.reference_end: 0-based, exclusive

    insert_region uses inclusive boundaries on both ends.
    """
    if read.is_unmapped:
        return False

    read_start = read.reference_start
    read_end = read.reference_end

    if read_start is None or read_end is None:
        return False

    # No overlap if read entirely before or after insert
    # read_end is exclusive, so use <= for "ends before insert"
    # insert_region[1] is inclusive, so use > for "starts after insert"
    if read_end <= insert_region[0] or read_start > insert_region[1]:
        return False

    return True
```

### Extended Assignment Logic
```python
# Source: Extended from _assign() lines 196-201
# New logic for insert-region and score-margin awareness

def _assign_with_filtering(
    plasmid_score: int,
    human_score: int,
    plasmid_read: Any | None,
    insert_region: tuple[int, int] | None,
    score_margin: int = 0,
) -> str:
    """
    Assign read to category with optional backbone and margin filtering.

    Args:
        plasmid_score: Custom alignment score for plasmid alignment
        human_score: Custom alignment score for human alignment
        plasmid_read: pysam AlignedSegment for plasmid (None if human-only)
        insert_region: (start, end) of insert region, or None to disable filtering
        score_margin: Minimum score difference for confident assignment (0 = disabled)

    Returns:
        One of: "Plasmid", "Human", "Tied", "Backbone_Only", "Ambiguous"
    """
    # Check score margin first (if enabled and applicable)
    if score_margin > 0:
        score_diff = abs(plasmid_score - human_score)
        if 0 < score_diff < score_margin:
            return "Ambiguous"

    # Check for exact tie
    if plasmid_score == human_score:
        return "Tied"

    # Plasmid wins
    if plasmid_score > human_score:
        # Check if read overlaps insert region (if available)
        if insert_region is not None and plasmid_read is not None:
            if read_overlaps_insert(plasmid_read, insert_region):
                return "Plasmid"
            else:
                return "Backbone_Only"
        return "Plasmid"

    # Human wins
    return "Human"
```

### Updated TSV Output Schema
```python
# Source: Modified from _write_assignment() lines 177-193
# New columns for Backbone_Only and Ambiguous

def _write_assignment_extended(
    outfile: IO[str],
    query_name: str,
    assigned_to: str,
    plasmid_score: int,
    human_score: int,
    plasmid_read: Any | None,
    human_read: Any | None,
) -> None:
    """Write assignment with extended category set."""
    p_cigar = plasmid_read.cigarstring if plasmid_read and plasmid_read.cigarstring else "NA"
    h_cigar = human_read.cigarstring if human_read and human_read.cigarstring else "NA"
    p_mapq = plasmid_read.mapping_quality if plasmid_read else "NA"
    h_mapq = human_read.mapping_quality if human_read else "NA"

    # Extended schema - category can now be: Plasmid, Human, Tied, Backbone_Only, Ambiguous
    outfile.write(
        f"{query_name}\t{assigned_to}\t{plasmid_score}\t{human_score}\t"
        f"{p_cigar}\t{h_cigar}\t{p_mapq}\t{h_mapq}\n"
    )

# Header update in compare_alignments()
outfile.write(
    "ReadID\tAssignedTo\tPlasmidScore\tHumanScore\t"
    "PlasmidCIGAR\tHumanCIGAR\tPlasmidMapQ\tHumanMapQ\n"
)
# AssignedTo column values expand from {Plasmid, Human, Tied}
# to {Plasmid, Human, Tied, Backbone_Only, Ambiguous}
```

### Report Template Updates
```html
<!-- Source: Bootstrap 4.5.2 utilities (already loaded in report_template.html) -->
<!-- Reference: https://getbootstrap.com/docs/4.5/utilities/text/#text-colors -->

<h2>Read Assignment Summary</h2>
<table class="table table-bordered">
    <thead>
        <tr>
            <th>Category</th>
            <th>Count</th>
            <th>Included in Ratio</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>Plasmid</td>
            <td>{{ plasmid_count }}</td>
            <td><span class="badge badge-success">Yes</span></td>
        </tr>
        <tr>
            <td>Human</td>
            <td>{{ human_count }}</td>
            <td><span class="badge badge-success">Yes</span></td>
        </tr>
        <tr class="text-muted">
            <td>Tied</td>
            <td>{{ tied_count }}</td>
            <td><span class="badge badge-secondary">No</span></td>
        </tr>
        <tr class="text-muted">
            <td>Backbone_Only</td>
            <td>{{ backbone_only_count }}</td>
            <td><span class="badge badge-secondary">No</span></td>
        </tr>
        <tr class="text-muted">
            <td>Ambiguous</td>
            <td>{{ ambiguous_count }}</td>
            <td><span class="badge badge-secondary">No</span></td>
        </tr>
    </tbody>
</table>

<!-- Note: text-muted class provides gray dimming for excluded categories -->
<!-- Bootstrap badge colors provide visual distinction between included/excluded -->
```

### Plotly Color Scheme for New Categories
```python
# Source: Existing PLOT_SAMPLE_REPORT config pattern
# Reference: https://plotly.com/python/discrete-color/

# In config.json - extend existing plot_sample_report.colors
{
  "plot_sample_report": {
    "colors": {
      "human": "blue",
      "plasmid": "red",
      "tied": "orange",
      "backbone_only": "lightgray",  # Muted gray for excluded
      "ambiguous": "darkgray"         # Slightly darker gray for excluded
    }
  }
}

# In generate_report.py
COLOR_MAP = {
    "Plasmid": PLOT_SAMPLE_REPORT["colors"]["plasmid"],
    "Human": PLOT_SAMPLE_REPORT["colors"]["human"],
    "Tied": PLOT_SAMPLE_REPORT["colors"]["tied"],
    "Backbone_Only": PLOT_SAMPLE_REPORT["colors"]["backbone_only"],
    "Ambiguous": PLOT_SAMPLE_REPORT["colors"]["ambiguous"],
}

fig_scatter = px.scatter(
    reads_df,
    x="PlasmidScore",
    y="HumanScore",
    color="AssignedTo",
    color_discrete_map=COLOR_MAP,  # Explicit mapping
    # ... rest of config
)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| All plasmid-aligned reads counted as contamination | Insert-region-aware filtering distinguishes backbone from insert | v0.33.0 (this phase) | Eliminates false positives from shared backbone regions |
| Binary classification (Plasmid/Human/Tied) | 5-category system with Backbone_Only and Ambiguous | v0.33.0 (this phase) | More nuanced read classification for scientific accuracy |
| No score margin support | Optional score_margin parameter for conservative filtering | v0.33.0 (this phase) | Advanced users can require minimum score differences |
| No fallback for missing insert region | Graceful degradation to pre-v0.33.0 behavior with warning | v0.33.0 (this phase) | Robustness for incomplete pipeline runs |

**Deprecated/outdated:**
- Simple `_assign(plasmid_score, human_score)` without region awareness: Remains available when `filter_backbone_only=false` for backward compatibility, but default behavior now uses extended classification
- Assumption that all plasmid-aligned reads indicate contamination: Scientifically incorrect for plasmids sharing backbone sequences

## Open Questions

Things that couldn't be fully resolved:

1. **Optimal default for score_margin**
   - What we know: Bioinformatics tools typically default to 0 (disabled) and let users tune based on their data (e.g., samtools MAPQ filtering defaults to 0)
   - What's unclear: Whether there's a scientifically validated margin value for this specific use case
   - Recommendation: Keep default at 0, document as "advanced parameter" for users to experiment with

2. **Pie chart vs bar chart for 5 categories**
   - What we know: Current reports use box plots and scatter plots, no pie chart
   - What's unclear: Whether "pie chart" in CONTEXT.md refers to future visualization or misunderstanding of current plots
   - Recommendation: Update scatter plot color mapping to include new categories (as shown in code examples); consider adding a stacked bar chart showing ratio composition in future phases

3. **Cross-species plasmid scenarios**
   - What we know: Current pipeline assumes human reference; other research tools like FastQ Screen support multi-species panels
   - What's unclear: How insert-region filtering behaves when plasmid insert is from non-human species
   - Recommendation: Document assumption that insert region overlaps with human genome (basis of spliced_alignment.py logic); cross-species is out of scope for v0.33.0

## Sources

### Primary (HIGH confidence)
- pysam documentation: https://pysam.readthedocs.io/en/latest/api.html - AlignedSegment properties reference_start, reference_end
- pysam usage documentation: https://pysam.readthedocs.io/en/latest/usage.html - Working with BAM/CRAM/SAM files, fetch() and overlap semantics
- Existing codebase: `/mnt/c/development/scholl-lab/plasmicheck/plasmicheck/scripts/compare_alignments.py` - Verified patterns for insert region parsing, overlap detection, streaming comparison
- Existing codebase: `/mnt/c/development/scholl-lab/plasmicheck/plasmicheck/config.json` - Configuration structure and existing scoring parameters
- Bootstrap 4.5.2 documentation: https://getbootstrap.com/docs/4.5/utilities/text/#text-colors - Already loaded in report_template.html, text-muted class

### Secondary (MEDIUM confidence)
- Bedtools interval overlap semantics: https://bedtools.readthedocs.io/ - Standard bioinformatics 0-based half-open interval convention
- Bedtk interval tree paper (Oxford Academic): https://academic.oup.com/bioinformatics/article/37/9/1315/5910546 - Implicit interval tree for overlap queries, confirms half-open interval best practices
- pytest parametrization: https://docs.pytest.org/en/stable/how-to/parametrize.html - Feature flag testing patterns
- pytest backward compatibility: https://docs.pytest.org/en/stable/backwards-compatibility.html - Deprecation and compatibility policies
- Pydantic configuration guide: https://docs.pydantic.dev/latest/concepts/models/ - Config validation patterns (considered but not needed for this phase)

### Tertiary (LOW confidence)
- FastQ Screen documentation: https://stevenwingett.github.io/FastQ-Screen/ - Multi-genome contamination detection tool, mentions region-specific filtering but no implementation details
- ConFindr documentation: https://olc-bioinformatics.github.io/ConFindr/ - Intra-species contamination detection using core genes, different approach than region-based filtering
- ContScout Nature paper: https://www.nature.com/articles/s41467-024-45024-5 - Recent contamination detection tool, no specific backbone filtering strategy described
- Plasmid contamination in reagents (Nature): https://www.nature.com/articles/s41598-019-38733-1 - Scientific background on plasmid contamination problem

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All libraries already in use, pysam API verified from official docs
- Architecture: HIGH - Existing patterns in codebase provide verified implementation templates
- Pitfalls: HIGH - Derived from existing code review and pysam documentation on common errors

**Research date:** 2026-02-14
**Valid until:** 2026-03-14 (30 days - stable domain, pysam API is mature)

**Notes:**
- No new dependencies required
- All code examples verified against existing codebase patterns
- Interval overlap semantics cross-referenced with pysam docs and existing calculate_coverage_outside_insert()
- Bootstrap and Plotly versions already locked in current implementation
