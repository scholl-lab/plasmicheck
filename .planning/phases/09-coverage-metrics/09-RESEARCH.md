# Phase 9: Coverage Metrics - Research

**Researched:** 2026-02-15
**Domain:** Per-base depth calculation, coverage statistics (breadth, mean/median depth, uniformity)
**Confidence:** HIGH

## Summary

Coverage metrics quantify how well sequencing reads cover different regions of a reference sequence. This phase requires computing per-region (insert vs backbone) depth, breadth, and uniformity metrics from BAM alignments and writing them to summary.tsv.

The standard approach uses **pysam.AlignmentFile.count_coverage()** for per-base depth arrays, then applies numpy statistics (mean, median, std) to compute region-specific metrics. Breadth measures the fraction of bases with ≥N reads (typically 1x or 5x), while uniformity uses coefficient of variation (CV = std/mean) to detect uneven coverage.

**Critical finding:** The current pipeline filters out unmapped reads during alignment (`samtools view -F 4`), meaning `plasmid_alignment.bam` contains ONLY reads that mapped to the plasmid reference. This is correct for coverage calculation - we measure depth of aligned reads, not total input reads. Coverage metrics answer "how well do the reads that DID align cover the reference?" not "what fraction of input reads aligned?"

**Primary recommendation:** Use pysam.count_coverage() with quality_threshold=0 (to match existing alignment filtering), compute region-specific statistics using numpy, and add new rows to summary.tsv in CamelCase format matching existing conventions.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pysam | >=0.22.0 | BAM I/O and per-base depth calculation | Industry standard for SAM/BAM manipulation, already used in project |
| numpy | >=1.20 | Statistical computations (mean, median, std, CV) | De facto standard for numerical computing in bioinformatics |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| scipy | >=1.10 | scipy.stats.variation() for CV | Optional - can compute CV manually with numpy |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| pysam.count_coverage() | samtools depth + parsing | samtools depth is faster but requires subprocess + file I/O, less Pythonic |
| pysam | mosdepth (external tool) | mosdepth is 2x faster than samtools but adds external dependency, no programmatic API |
| Manual CV calculation | scipy.stats.variation() | scipy adds dependency for one-line function, manual is clearer |

**Installation:**
```bash
# Already in project dependencies
pip install pysam numpy
```

## Architecture Patterns

### Recommended Function Structure
```python
# Add to plasmicheck/scripts/compare_alignments.py
def compute_coverage_metrics(
    bam_path: str,
    insert_region: tuple[int, int] | None = None,
    breadth_thresholds: list[int] | None = None
) -> dict[str, dict[str, float]]:
    """Compute coverage metrics per region.

    Returns:
        {
            "insert": {"mean_depth": X, "median_depth": Y, ...},
            "backbone": {"mean_depth": X, "median_depth": Y, ...}
        }
    """
```

### Pattern 1: Per-Base Depth Extraction
**What:** Use pysam.count_coverage() to get per-base depth arrays for specific regions
**When to use:** When computing any coverage metric (breadth, mean depth, uniformity)
**Example:**
```python
# Source: pysam documentation + existing project patterns
import pysam
import numpy as np

def get_region_depth_array(bam_path: str, contig: str, start: int, end: int) -> np.ndarray:
    """Extract per-base depth for a genomic region."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # count_coverage returns tuple of 4 arrays (A, C, G, T counts)
        coverage_arrays = bam.count_coverage(
            contig=contig,
            start=start,
            stop=end,
            quality_threshold=0,  # Match existing alignment filtering
            read_callback="all"
        )
        # Sum across nucleotides to get total depth per position
        depth_array = np.array(coverage_arrays[0]) + np.array(coverage_arrays[1]) + \
                      np.array(coverage_arrays[2]) + np.array(coverage_arrays[3])
    return depth_array
```

### Pattern 2: Breadth Calculation
**What:** Compute fraction of bases with depth ≥ threshold
**When to use:** For COV-02 (breadth ≥1x) and COV-03 (breadth ≥5x)
**Example:**
```python
def compute_breadth(depth_array: np.ndarray, threshold: int = 1) -> float:
    """Compute fraction of bases with depth >= threshold."""
    if len(depth_array) == 0:
        return 0.0
    bases_above_threshold = np.sum(depth_array >= threshold)
    return bases_above_threshold / len(depth_array)
```

### Pattern 3: Coefficient of Variation
**What:** CV = std / mean, measures coverage uniformity
**When to use:** For COV-04 (uniformity metric)
**Example:**
```python
def compute_cv(depth_array: np.ndarray) -> float:
    """Compute coefficient of variation for coverage uniformity."""
    if len(depth_array) == 0 or np.mean(depth_array) == 0:
        return 0.0
    return np.std(depth_array, ddof=1) / np.mean(depth_array)
```

### Pattern 4: Region-Specific Metrics
**What:** Split plasmid into insert vs backbone regions, compute metrics separately
**When to use:** Always (per Phase 8 context - regions are fundamental)
**Example:**
```python
def split_by_region(
    bam_path: str,
    insert_region: tuple[int, int] | None
) -> tuple[np.ndarray, np.ndarray]:
    """Split coverage into insert and backbone arrays."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        contig = bam.get_reference_name(0)
        plasmid_length = bam.lengths[0]

        if insert_region is None:
            # Fallback: whole plasmid as single region
            whole_depth = get_region_depth_array(bam_path, contig, 0, plasmid_length)
            return whole_depth, np.array([])

        # Insert region
        insert_depth = get_region_depth_array(
            bam_path, contig, insert_region[0], insert_region[1] + 1
        )

        # Backbone regions (before + after insert)
        backbone_before = get_region_depth_array(
            bam_path, contig, 0, insert_region[0]
        ) if insert_region[0] > 0 else np.array([])

        backbone_after = get_region_depth_array(
            bam_path, contig, insert_region[1] + 1, plasmid_length
        ) if insert_region[1] + 1 < plasmid_length else np.array([])

        backbone_depth = np.concatenate([backbone_before, backbone_after])

        return insert_depth, backbone_depth
```

### Pattern 5: Summary.tsv Integration
**What:** Append new coverage rows after existing metrics, maintain backward compatibility
**When to use:** Writing metrics to summary.tsv (COV-05)
**Example:**
```python
# In compare_alignments.py, extend summary.tsv writing:
with open(f"{output_basename}.summary.tsv", "w") as summary_file:
    summary_file.write("Category\tCount\n")
    # Existing rows (unchanged)
    for category, count in assigned_counts.items():
        summary_file.write(f"{category}\t{count}\n")
    summary_file.write(f"Verdict\t{verdict}\n")
    summary_file.write(f"Ratio\t{ratio}\n")
    summary_file.write(f"CoverageOutsideINSERT\t{coverage_outside_insert:.4f}\n")
    summary_file.write(f"MismatchesNearINSERT\t{mismatches_near_insert}\n")

    # New coverage rows (Phase 9)
    summary_file.write(f"MeanDepthInsert\t{metrics['insert']['mean_depth']:.2f}\n")
    summary_file.write(f"MedianDepthInsert\t{metrics['insert']['median_depth']:.2f}\n")
    summary_file.write(f"BreadthInsert\t{metrics['insert']['breadth_1x']:.2f}\n")
    summary_file.write(f"CoverageCV_Insert\t{metrics['insert']['cv']:.2f}\n")
    summary_file.write(f"MeanDepthBackbone\t{metrics['backbone']['mean_depth']:.2f}\n")
    summary_file.write(f"MedianDepthBackbone\t{metrics['backbone']['median_depth']:.2f}\n")
    summary_file.write(f"BreadthBackbone\t{metrics['backbone']['breadth_1x']:.2f}\n")
    summary_file.write(f"BreadthBackbone_5x\t{metrics['backbone']['breadth_5x']:.2f}\n")
    summary_file.write(f"CoverageCV_Backbone\t{metrics['backbone']['cv']:.2f}\n")
```

### Anti-Patterns to Avoid
- **Don't iterate reads to compute coverage manually:** Use pysam.count_coverage() - it's optimized and handles edge cases
- **Don't use pysam.pileup():** Slower than count_coverage(), designed for variant calling not coverage stats
- **Don't call samtools depth subprocess:** Adds I/O overhead, requires parsing, no benefit over pysam API
- **Don't compute whole-plasmid aggregate metrics:** Violates Phase 8 region-aware design, less informative than per-region

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Per-base depth from BAM | Manual CIGAR parsing + position tracking | pysam.count_coverage() | Handles CIGAR operations (M/I/D/S/H), quality filtering, edge cases correctly |
| Reference length lookup | Parse FASTA headers | pysam.AlignmentFile.lengths | Already loaded in BAM header, no file I/O |
| Contig name resolution | Hardcode "plasmid" or guess | pysam.AlignmentFile.get_reference_name(0) | Handles arbitrary reference names |
| Coverage array merging | Manual list concatenation + conversion | numpy.concatenate() | Efficient, type-safe, handles empty arrays |
| Coefficient of variation | Manual std/mean with edge case checks | numpy std/mean (with explicit zero-check) | Standard, tested, efficient - but manually check for zero mean |

**Key insight:** BAM files embed reference metadata (lengths, contig names) - use pysam accessors, don't re-parse FASTAs.

## Common Pitfalls

### Pitfall 1: Off-by-One Errors in Region Boundaries
**What goes wrong:** Insert region from cDNA_positions.txt uses inclusive boundaries [start, end], but pysam uses half-open [start, end). Passing insert_region[1] directly to count_coverage() misses the last base.
**Why it happens:** Mixing biological coordinates (1-based or inclusive) with pysam's 0-based half-open intervals.
**How to avoid:** Always use `end + 1` when passing inclusive end to pysam: `bam.count_coverage(start=insert_region[0], stop=insert_region[1] + 1)`
**Warning signs:** Insert breadth always slightly less than expected, last base of insert shows zero depth in manual checks

### Pitfall 2: Empty Regions Causing Division by Zero
**What goes wrong:** CV computation divides by mean depth - if a region has zero depth (clean sample, no backbone coverage), mean=0 causes runtime error or NaN.
**Why it happens:** Contaminated samples are common in testing, clean samples (where backbone_depth array is all zeros) expose the edge case.
**How to avoid:** Explicit zero-checks before division: `if len(arr) == 0 or np.mean(arr) == 0: return 0.0`
**Warning signs:** Tests pass on contaminated data, fail on clean samples with "RuntimeWarning: invalid value encountered in double_scalars"

### Pitfall 3: Quality Threshold Mismatch
**What goes wrong:** Using pysam.count_coverage() default quality_threshold=15 excludes bases that were included in alignments (minimap2 doesn't filter by base quality), causing coverage to appear lower than actual.
**Why it happens:** count_coverage() defaults to Phred≥15 for variant calling use cases, but contamination detection cares about all aligned bases.
**How to avoid:** Always set `quality_threshold=0` to match alignment behavior: `bam.count_coverage(..., quality_threshold=0)`
**Warning signs:** Coverage metrics seem too low compared to read counts, contaminated samples show unexpectedly sparse backbone coverage

### Pitfall 4: Double-Counting Overlapping Paired Reads
**What goes wrong:** pysam.count_coverage() counts both mates in overlapping paired-end reads, inflating depth by ~2x in overlap regions.
**Why it happens:** count_coverage() processes each alignment record independently, doesn't know about mate pairs.
**How to avoid:** Accept this behavior - it's consistent with samtools depth default, and contamination detection compares plasmid vs human depth (both double-counted equally). Document in code comments.
**Warning signs:** Mean depth higher than expected from read count / region length calculation

### Pitfall 5: Missing Insert Region (Fallback Behavior)
**What goes wrong:** When cDNA_positions.txt is missing or malformed, coverage metrics fail instead of falling back gracefully.
**Why it happens:** Aggressive error handling in parse_insert_region() raises exceptions, but coverage computation happens later.
**How to avoid:** Check `if insert_region is None:` and compute whole-plasmid metrics, add warning to HTML report: "Insert region not defined — metrics computed for whole plasmid"
**Warning signs:** Pipeline crashes on samples where spliced alignment failed, but comparison still ran

### Pitfall 6: count_coverage() Returns Tuple of 4 Arrays (Not Total Depth)
**What goes wrong:** Directly using count_coverage() return value as depth array - it's actually (A_array, C_array, G_array, T_array), not total depth.
**Why it happens:** count_coverage() is designed for nucleotide-specific analysis (variant calling), not just depth.
**How to avoid:** Always sum the 4 arrays: `depth = sum(bam.count_coverage(...))` or explicitly `np.array(A) + np.array(C) + ...`
**Warning signs:** Type errors ("unsupported operand type for +: 'tuple'"), or using only first array and getting 4x lower depth

## Code Examples

Verified patterns from official sources:

### Extract Per-Base Depth for a Region
```python
# Source: pysam documentation https://pysam.readthedocs.io/en/latest/api.html
import pysam
import numpy as np

def get_depth_array(bam_path: str, contig: str, start: int, end: int) -> np.ndarray:
    """Get per-base depth array for a region [start, end)."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        cov = bam.count_coverage(
            contig=contig,
            start=start,
            stop=end,
            quality_threshold=0,  # Include all bases (match alignment behavior)
            read_callback="all"
        )
        # Sum A, C, G, T counts to get total depth per position
        depth = np.array(cov[0]) + np.array(cov[1]) + np.array(cov[2]) + np.array(cov[3])
    return depth
```

### Compute Coverage Statistics
```python
# Source: numpy documentation + bioinformatics conventions
import numpy as np

def compute_coverage_stats(depth_array: np.ndarray, breadth_thresholds: list[int]) -> dict[str, float]:
    """Compute mean, median, breadth, and CV from depth array."""
    if len(depth_array) == 0:
        return {
            "mean_depth": 0.0,
            "median_depth": 0.0,
            **{f"breadth_{t}x": 0.0 for t in breadth_thresholds},
            "cv": 0.0
        }

    stats = {
        "mean_depth": float(np.mean(depth_array)),
        "median_depth": float(np.median(depth_array)),
    }

    # Breadth at each threshold
    for threshold in breadth_thresholds:
        breadth = np.sum(depth_array >= threshold) / len(depth_array)
        stats[f"breadth_{threshold}x"] = breadth

    # Coefficient of variation (check for zero mean)
    mean_val = stats["mean_depth"]
    if mean_val == 0:
        stats["cv"] = 0.0
    else:
        stats["cv"] = float(np.std(depth_array, ddof=1) / mean_val)

    return stats
```

### Full Integration Pattern
```python
# Source: Existing compare_alignments.py patterns
def compute_region_coverage_metrics(
    plasmid_bam: str,
    insert_region: tuple[int, int] | None,
    breadth_thresholds: list[int]
) -> dict[str, dict[str, float]]:
    """Compute coverage metrics for insert and backbone regions."""
    with pysam.AlignmentFile(plasmid_bam, "rb") as bam:
        contig = bam.get_reference_name(0)
        plasmid_length = bam.lengths[0]

    if insert_region is None:
        # Fallback: whole plasmid as insert, empty backbone
        whole_depth = get_depth_array(plasmid_bam, contig, 0, plasmid_length)
        return {
            "insert": compute_coverage_stats(whole_depth, breadth_thresholds),
            "backbone": {k: 0.0 for k in ["mean_depth", "median_depth", "cv"] +
                         [f"breadth_{t}x" for t in breadth_thresholds]}
        }

    # Insert region (inclusive boundaries -> half-open for pysam)
    insert_depth = get_depth_array(plasmid_bam, contig, insert_region[0], insert_region[1] + 1)

    # Backbone regions
    backbone_parts = []
    if insert_region[0] > 0:
        backbone_parts.append(get_depth_array(plasmid_bam, contig, 0, insert_region[0]))
    if insert_region[1] + 1 < plasmid_length:
        backbone_parts.append(get_depth_array(plasmid_bam, contig, insert_region[1] + 1, plasmid_length))

    backbone_depth = np.concatenate(backbone_parts) if backbone_parts else np.array([])

    return {
        "insert": compute_coverage_stats(insert_depth, breadth_thresholds),
        "backbone": compute_coverage_stats(backbone_depth, breadth_thresholds)
    }
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| samtools depth subprocess | pysam.count_coverage() API | ~2018 (pysam 0.15+) | Faster, Pythonic, no I/O |
| Whole-genome metrics | Region-specific metrics | Phase 8 (v0.33.0) | Insert vs backbone insight |
| Read-count approximations | Per-base depth arrays | Phase 9 (this phase) | Accurate breadth, uniformity |
| Percentage format (95%) | Fraction format (0.95) | Phase 9 context decision | Consistent with scientific norms |

**Deprecated/outdated:**
- **samtools depth -a**: Still works but slower than pysam API for programmatic use
- **pysamstats**: More powerful than needed, adds dependency, slower for simple depth stats
- **BEDTools genomecov**: External tool, requires BED conversion, no advantage for single-region stats

## Critical Investigation: Read Mapping Flow

### Question: What Reads Are Available for Coverage Computation?

**User's concern:** Does the pipeline use unmapped reads to map against the plasmid, or only reads extracted from the gene region mapped to the plasmid?

**Finding:** The pipeline uses **all input reads** for both alignments, but filters out unmapped reads after alignment.

### Current Read Flow

```
1. Input: FASTQ files (all reads, R1 + R2)
   ↓
2. align_reads.py (plasmid):
   minimap2 -ax sr plasmid.mmi input.fastq \
   | samtools view -F 4    # ← FILTERS OUT UNMAPPED
   | samtools sort
   → plasmid_alignment.bam (only reads that mapped to plasmid)

3. align_reads.py (human):
   minimap2 -ax sr human.mmi input.fastq \
   | samtools view -F 4    # ← FILTERS OUT UNMAPPED
   | samtools sort
   → spliced_human_alignment.bam (only reads that mapped to human)

4. compare_alignments.py:
   - Compares reads present in BOTH BAMs (mapped to both)
   - Compares reads present in ONLY plasmid BAM (plasmid-specific)
   - Compares reads present in ONLY human BAM (human-specific)
   - Missing: reads that didn't map to either reference
```

### Key Findings

**`samtools view -F 4` means:** Exclude reads with flag 4 (UNMAP). This removes unmapped reads from the BAM file.

**Implication for coverage:**
- `plasmid_alignment.bam` contains ONLY reads that successfully aligned to the plasmid reference
- Coverage metrics computed from this BAM measure "depth of aligned reads per base"
- This is **correct** for coverage calculation - we want to know how well the aligned reads cover the reference

**What about unmapped reads?**
- Reads that don't align to either reference are lost (never enter the comparison)
- This is intentional - contamination detection compares reads that aligned to plasmid vs human
- Coverage metrics answer: "For reads that DID align to plasmid, how thoroughly do they cover it?"

### User's Suggestion: `samtools view -o unmap.bam input.bam '*'`

**What it does:** The `'*'` in RNAME column indicates unmapped reads in SAM format. However, FLAG-based filtering (`-f 4`) is more reliable.

**Relevance to coverage metrics:** Not applicable. Coverage metrics measure depth of **mapped** reads. Unmapped reads by definition contribute zero coverage.

**When this IS relevant:** If we wanted to quantify "how many input reads aligned to plasmid vs didn't align at all" - but that's alignment rate, not coverage depth.

### Conclusion: Current Approach is Correct

The pipeline correctly:
1. Aligns all input reads to both references (plasmid and human)
2. Retains only mapped reads in output BAMs
3. Computes coverage metrics from mapped reads

Coverage breadth answers: "What fraction of the plasmid reference is covered by aligned reads?"

If breadth is low (e.g., 0.10), it means only 10% of plasmid bases have ≥1 aligned read - indicating either:
- Very low contamination (few reads align)
- Uneven coverage (reads cluster in specific regions)

This is exactly what Phase 9 coverage metrics should detect.

## Open Questions

None - research domain is well-established.

## Sources

### Primary (HIGH confidence)
- [pysam API documentation - AlignmentFile.count_coverage()](https://pysam.readthedocs.io/en/latest/api.html)
- [samtools-coverage manual](http://www.htslib.org/doc/samtools-coverage.html) - Coverage metric definitions
- [samtools-depth manual](http://www.htslib.org/doc/samtools-depth.html) - Per-base depth behavior
- [Coefficient of Variation - Wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_variation) - Statistical definition
- [NumPy statistics documentation](https://numpy.org/doc/stable/reference/routines.statistics.html)

### Secondary (MEDIUM confidence)
- [Sequencing coverage and breadth of coverage](https://www.reneshbedre.com/blog/sequencing-coverage.html) - Bioinformatics definitions
- [Coverage depth - Metagenomics Wiki](https://www.metagenomics.wiki/pdf/qc/coverage-read-depth) - Formula explanations
- [mosdepth GitHub](https://github.com/brentp/mosdepth) - Alternative tool comparison
- [Coefficient of Variation - Statistics By Jim](https://statisticsbyjim.com/basics/coefficient-variation/) - CV interpretation
- [Samtools: How to Filter Mapped and Unmapped Reads](https://www.reneshbedre.com/blog/filter-reads-samtools.html) - Flag filtering

### Tertiary (LOW confidence)
- [pysam count_coverage GitHub issues](https://github.com/pysam-developers/pysam/issues/598) - Community discussion on edge cases

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - pysam and numpy are established, already in project dependencies
- Architecture: HIGH - count_coverage() API is well-documented, patterns verified in existing codebase
- Pitfalls: MEDIUM - Some edge cases inferred from issue discussions, not all tested in PlasmiCheck context

**Research date:** 2026-02-15
**Valid until:** 60 days (stable bioinformatics domain, slow-changing APIs)
