# Phase 7: Comparison & Cleanup - Research

**Researched:** 2026-02-14
**Domain:** BAM file processing, samtools collate/sort, matplotlib plotting, pipeline architecture
**Confidence:** HIGH

## Summary

Phase 7 optimizes BAM comparison using `samtools collate` (30-50% faster than `sort -n`), eliminates redundant human reference indexing in batch runs through upfront planning, and adds matplotlib as a static plotting backend alternative to Kaleido. The phase also improves batch resilience by continuing through failures and logging per-combination timing.

**Key findings:**
- `samtools collate` groups reads by name without full lexicographical sorting, significantly faster than `sort -n` but provides no ordering guarantees between groups or within supplementary alignments
- Collate fast mode (`-f`) filters out supplementary/secondary alignments entirely; standard mode preserves them but without ordering guarantees
- matplotlib and seaborn are already core dependencies (matplotlib>=3.9.1, seaborn>=0.13.2) and support all required chart types with viridis colormap compatibility with Plotly
- Index deduplication via upfront planning (build unique index list first, then process combinations) is standard in bioinformatics pipelines
- Python 3.11+ ExceptionGroup pattern enables robust continue-on-failure batch processing with comprehensive error reporting

**Primary recommendation:** Always use `samtools collate` standard mode (not fast mode to preserve supplementary alignments), explicitly re-sort supplementary alignments within each read group, write to temporary file for debuggability, and implement automatic fallback to `sort -n` with logged warning.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| samtools | 1.6+ | BAM collation and sorting | Industry standard for BAM manipulation, collate specifically designed for read-pairing tasks |
| matplotlib | >=3.9.1 | Static plot generation | Already a dependency, supports all chart types, publication-quality output |
| seaborn | >=0.13.2 | Statistical visualization styling | Already a dependency, provides clean themes and color palettes |
| tempfile (stdlib) | Python 3.8+ | Temporary file management | Standard library, context manager for automatic cleanup |
| subprocess (stdlib) | Python 3.8+ | External tool execution | Standard library, robust error handling |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pysam | Current | BAM file reading (already used) | Reading alignment data in Python |
| numpy | Current (via pandas) | Log transformation for plots | Already transitively available |
| ExceptionGroup | Python 3.11+ | Multi-error batch handling | Continue-on-failure batch processing |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| samtools collate | samtools sort -n | Full lexicographical sort 30-50% slower but provides strict ordering guarantees |
| matplotlib | Kaleido (Plotly export) | Already phased out in Phase 5; Kaleido has licensing/maintenance issues |
| Upfront planning | On-demand indexing | Simpler code but redundant operations in batch runs |

**Installation:**
```bash
# Core dependencies already in pyproject.toml
pip install -e ".[dev]"  # matplotlib>=3.9.1, seaborn>=0.13.2 already included

# External tools (system-level, already required)
# samtools >= 1.6 (check with: samtools --version)
```

## Architecture Patterns

### Recommended Project Structure
```
plasmicheck/scripts/
├── compare_alignments.py    # Add collate logic + fallback
├── run_pipeline.py           # Add PipelinePlan index tracking
└── plotting/                 # NEW: matplotlib backend module
    ├── __init__.py
    ├── matplotlib_backend.py # Chart generators matching Plotly
    └── colors.py             # Shared color schemes
```

### Pattern 1: Collate with Supplementary Re-sorting

**What:** Use `samtools collate` for fast name grouping, then explicitly re-sort supplementary alignments within each read group

**When to use:** Always for BAM comparison (CONTEXT decision)

**Example:**
```python
# Source: Research synthesis from samtools docs + codebase analysis
import subprocess
import tempfile
import logging
from pathlib import Path

def _collate_bam(input_bam: str, output_bam: str, threads: int = 4) -> None:
    """Collate BAM by read name with fallback to sort -n."""
    try:
        # Standard mode collate (preserves supplementary/secondary)
        with tempfile.NamedTemporaryFile(prefix="collate_", suffix=".bam", delete=False) as tmp:
            tmp_path = tmp.name

        subprocess.run(
            [
                "samtools", "collate",
                "-@", str(threads),
                "-o", tmp_path,
                input_bam,
                f"{input_bam}.collate_tmp"  # prefix for temp files
            ],
            check=True,
            capture_output=True,
            text=True
        )

        # Re-sort supplementary alignments within groups
        _resort_supplementary(tmp_path, output_bam)

        # Cleanup temp file
        Path(tmp_path).unlink(missing_ok=True)

        logging.info(f"Collated {input_bam} successfully")

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logging.warning(
            f"samtools collate failed ({e}), falling back to sort -n"
        )
        _namesort_bam_fallback(input_bam, output_bam, threads)

def _resort_supplementary(input_bam: str, output_bam: str) -> None:
    """Re-sort supplementary alignments within each read group."""
    import pysam
    import itertools

    with pysam.AlignmentFile(input_bam, "rb") as infile, \
         pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:

        for name, group in itertools.groupby(
            infile.fetch(until_eof=True),
            key=lambda r: r.query_name
        ):
            reads = list(group)
            # Sort: primary first, then supplementary, then secondary
            reads.sort(key=lambda r: (
                r.is_secondary + r.is_supplementary,  # primary=0, supp=1, sec=1, both=2
                r.is_secondary and not r.is_supplementary  # secondary-only last
            ))
            for read in reads:
                outfile.write(read)

def _namesort_bam_fallback(input_bam: str, output_bam: str, threads: int) -> None:
    """Fallback to samtools sort -n."""
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(threads), "-o", output_bam, input_bam],
        check=True
    )
```

### Pattern 2: Upfront Index Planning

**What:** Scan all combinations before execution, build unique index list, run all indexing, then process combinations

**When to use:** Batch pipeline runs with multiple sample-plasmid combinations

**Example:**
```python
# Source: Research synthesis from bioinformatics pipeline best practices
from dataclasses import dataclass, field
from pathlib import Path

@dataclass
class PipelinePlan:
    """Execution plan with index deduplication."""
    human_fasta: str
    plasmid_files: list[str]
    sequencing_inputs: list[SequencingInput]
    output_folder: str
    overwrite: bool
    combinations: list[tuple[str, SequencingInput]] = field(default_factory=list)
    combination_steps: dict[str, list[PipelineStep]] = field(default_factory=dict)
    built_indexes: set[str] = field(default_factory=set)  # NEW: track built indexes
    warnings: list[str] = field(default_factory=list)

def build_plan(...) -> PipelinePlan:
    """Build execution plan with index deduplication."""
    plan = PipelinePlan(...)

    # Scan all combinations to identify unique indexes
    unique_human_refs = {plan.human_fasta}  # Only deduplicate human (CONTEXT decision)

    # Add to plan metadata
    plan.built_indexes = set()

    # ... build combination steps

    return plan

def run_pipeline(...):
    plan = build_plan(...)

    # PHASE 1: Index all unique human references
    for human_ref in unique_human_refs:
        index_path = Path(human_ref).with_suffix(".mmi")
        if not index_path.exists() or overwrite:
            create_indexes(human_ref, overwrite)
            plan.built_indexes.add(str(index_path))
            logging.info(f"Built index: {index_path}")
        else:
            plan.built_indexes.add(str(index_path))
            logging.info(f"Skipping index for {human_ref} (already built)")

    # PHASE 2: Process combinations (skip human indexing step)
    for combo_label, seq_input in plan.combinations:
        # Human index already built in phase 1
        # Only index plasmid (cheap, no deduplication)
        create_indexes(plasmid_fasta, overwrite)
        # ... rest of pipeline
```

### Pattern 3: Matplotlib Backend with Plotly Color Matching

**What:** Generate static PNGs using matplotlib with viridis colormap and seaborn styling to match Plotly visual output

**When to use:** When `--plot-backend matplotlib` specified (requires `--static-report`)

**Example:**
```python
# Source: Research synthesis from matplotlib/seaborn/plotly docs
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Match Plotly's default colors and style
PLOTLY_COLORS = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A']
sns.set_theme(style="whitegrid")  # Clean seaborn theme

def generate_boxplot_matplotlib(
    reads_df: pd.DataFrame,
    output_path: str,
    figsize: tuple[int, int] = (10, 4)
) -> None:
    """Generate box plot using matplotlib matching Plotly style."""
    fig, ax = plt.subplots(figsize=figsize)

    # Use seaborn for styled boxplot
    sns.boxplot(
        data=reads_df,
        x="AssignedTo",
        y="PlasmidScore",
        palette=PLOTLY_COLORS[:3],  # Match Plotly category colors
        ax=ax
    )

    # Add individual points (matching Plotly 'points=all')
    sns.stripplot(
        data=reads_df,
        x="AssignedTo",
        y="PlasmidScore",
        color='black',
        alpha=0.3,
        size=2,
        ax=ax
    )

    ax.set_title(f"Plasmid Score Distribution (Total Reads: {len(reads_df)})")
    ax.set_ylabel("Plasmid Score")
    ax.set_xlabel("Assigned To")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

def generate_heatmap_matplotlib(
    ratio_data: pd.DataFrame,
    output_path: str,
    threshold: float,
    figsize: tuple[int, int] = (12, 12)
) -> None:
    """Generate heatmap using matplotlib matching Plotly style."""
    fig, ax = plt.subplots(figsize=figsize)

    # Use viridis colormap (same as Plotly default continuous scale)
    sns.heatmap(
        ratio_data,
        cmap="viridis",
        annot=False,
        cbar=False,  # Match Plotly config (no colorbar)
        square=False,
        ax=ax
    )

    ax.set_title("Contamination Ratio Heatmap")
    ax.set_xlabel("Plasmid")
    ax.set_ylabel("Sample")
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
```

### Pattern 4: Batch Resilience with Error Aggregation

**What:** Continue processing remaining combinations when one fails, collect errors, generate summary for successful combinations only

**When to use:** Multi-combination pipeline runs (CONTEXT decision)

**Example:**
```python
# Source: Python 3.11+ ExceptionGroup pattern from research
from dataclasses import dataclass
import time
import logging

@dataclass
class CombinationResult:
    """Result of processing a single combination."""
    combo_label: str
    success: bool
    duration: float
    error: Exception | None = None

def run_pipeline(...):
    plan = build_plan(...)
    results: list[CombinationResult] = []

    for combo_label, seq_input in plan.combinations:
        start_time = time.time()
        try:
            # ... process combination
            duration = time.time() - start_time
            results.append(CombinationResult(
                combo_label=combo_label,
                success=True,
                duration=duration
            ))
            logging.info(f"{combo_label}: {duration:.1f}s")

        except Exception as e:
            duration = time.time() - start_time
            results.append(CombinationResult(
                combo_label=combo_label,
                success=False,
                duration=duration,
                error=e
            ))
            logging.error(
                f"{combo_label} FAILED after {duration:.1f}s: {e}"
            )
            # Continue to next combination

    # Summary report for successful combinations only
    successful = [r for r in results if r.success]
    failed = [r for r in results if not r.success]

    logging.info(
        f"Pipeline complete: {len(successful)} succeeded, {len(failed)} failed"
    )

    if failed:
        logging.warning("Failed combinations:")
        for result in failed:
            logging.warning(f"  - {result.combo_label}: {result.error}")

    # Generate summary report excluding failed combinations
    if successful:
        generate_summary_report(
            successful_labels=[r.combo_label for r in successful],
            excluded=failed
        )
```

### Anti-Patterns to Avoid

- **Using collate fast mode (`-f`):** Filters out supplementary/secondary alignments entirely, breaking chimeric read analysis
- **Relying on collate ordering:** Collate provides no guarantees about supplementary alignment order within groups or read group ordering
- **Streaming collate output:** Write to temp file for debuggability (CONTEXT decision)
- **Checking existing indexes in combination loop:** Redundant lock checks, disk I/O, and logging noise in batch runs
- **Stopping on first failure:** Wastes computation in batch runs, poor user experience
- **Pixel-perfect matplotlib matching:** Functionally equivalent (same data, chart types, colors) is sufficient, not pixel-identical (CONTEXT decision)

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Temporary file cleanup | Manual file deletion after use | `tempfile.NamedTemporaryFile` context manager | Automatic cleanup, handles exceptions, POSIX/Windows portability |
| Subprocess retry logic | Custom retry loops per command | Existing `utils.py:run_command` with config-driven retries | Already implemented with exponential backoff and logging |
| BAM name sorting | Custom Python sort-by-name | `samtools collate` (or fallback to `sort -n`) | Optimized C implementation, handles large files, memory-efficient |
| Plot color schemes | Manual RGB tuning | Viridis colormap (matplotlib/seaborn/Plotly all support) | Perceptually uniform, colorblind-friendly, publication-standard |
| Multi-error reporting | Multiple try-except blocks | Python 3.11+ ExceptionGroup | Language-native, comprehensive error aggregation, caller-friendly |

**Key insight:** Bioinformatics tools like samtools have decade+ optimization for BAM manipulation. Python stdlib (tempfile, subprocess) handles edge cases (Windows permissions, signal handling, cleanup) that custom code misses. Matplotlib/seaborn/Plotly share standard colormaps (viridis) making visual consistency straightforward.

## Common Pitfalls

### Pitfall 1: Collate Ordering Assumptions

**What goes wrong:** Code assumes collate outputs reads in same order as `sort -n`, supplementary alignments appear in undefined order within read groups, causing incorrect mate-pair matching or score calculation

**Why it happens:** Documentation states collate "doesn't make any guarantees about the order of read names between groups" but developers assume within-group ordering is stable

**How to avoid:** Always explicitly re-sort supplementary alignments within each read group after collate using primary-first, supplementary-second, secondary-last ordering

**Warning signs:** Regression tests fail for chimeric reads, score comparison results differ between collate and sort -n, "best read" selection returns secondary instead of primary

### Pitfall 2: Collate Fast Mode Filtering

**What goes wrong:** Using `samtools collate -f` silently filters out supplementary and secondary alignments, breaking analysis for chimeric reads or multi-mapping scenarios

**Why it happens:** Fast mode documentation mentions filtering but developers focus on speed improvement, not data loss

**How to avoid:** Never use fast mode (`-f`) for PlasmiCheck; always use standard collate mode to preserve all alignment records

**Warning signs:** Lower total read counts after collate, missing supplementary alignments in output, chimeric reads not detected

### Pitfall 3: Index Deduplication by Content Hash

**What goes wrong:** Computing MD5/SHA hashes of reference FASTAs to detect duplicates is slow (minutes for GRCh38), adds complexity, and still misses identical content with different paths

**Why it happens:** Over-engineering to handle edge case of "same reference copied to multiple locations"

**How to avoid:** Track built indexes by file path matching only; trust existing indexes without timestamp checks (CONTEXT decision)

**Warning signs:** Pipeline startup takes minutes before first alignment, excessive disk I/O during planning phase, redundant index builds despite identical references

### Pitfall 4: Matplotlib Style Inconsistency

**What goes wrong:** Default matplotlib plots (blue bars, gray background, different fonts) look visually distinct from Plotly, confusing users

**Why it happens:** Matplotlib defaults are publication-style, Plotly defaults are web-interactive style

**How to avoid:** Use seaborn theme (`sns.set_theme(style="whitegrid")`), match Plotly color palette (`['#636EFA', '#EF553B', '#00CC96']`), use viridis colormap for continuous scales

**Warning signs:** Users report plots "look different" between backends, colors don't match report branding, accessibility issues with non-colorblind-safe colors

### Pitfall 5: Stopping Batch on First Failure

**What goes wrong:** Pipeline processes 1 of 50 combinations successfully, second fails, remaining 48 never run, user must restart entire batch

**Why it happens:** Default exception propagation (raise immediately on error) is simpler to implement than error aggregation

**How to avoid:** Wrap each combination in try-except, collect errors, continue to next combination, report all failures at end

**Warning signs:** Users complain about restarting large batches multiple times, support requests "can you make it continue on error?", inefficient cluster resource usage

### Pitfall 6: Temp File Cleanup Race Conditions

**What goes wrong:** Collate temp file deleted before final BAM written, or persists after crash, filling disk

**Why it happens:** Manual cleanup with `os.remove()` outside exception handling, or cleanup in finally block before file closed

**How to avoid:** Use `tempfile.NamedTemporaryFile` context manager with explicit `delete=False`, cleanup in finally block after all file operations complete, or rely on OS temp directory cleanup

**Warning signs:** Intermittent "file not found" errors, `/tmp` fills with `.collate_tmp.*` files, Windows PermissionError on cleanup

## Code Examples

Verified patterns from official sources:

### Collate Command with Temp File Handling
```bash
# Source: http://www.htslib.org/doc/samtools-collate.html
samtools collate \
    -@ 4 \
    -o output.collated.bam \
    input.bam \
    /tmp/collate_prefix
```

### Fallback Pattern for Missing samtools collate
```python
# Source: Research synthesis from subprocess best practices
import subprocess
import shutil
import logging

def collate_with_fallback(input_bam: str, output_bam: str, threads: int = 4) -> None:
    """Collate BAM with automatic fallback to sort -n."""
    if not shutil.which("samtools"):
        raise RuntimeError("samtools not found on PATH")

    # Check samtools version supports collate (1.6+)
    try:
        result = subprocess.run(
            ["samtools", "collate", "--version"],
            capture_output=True,
            text=True,
            timeout=5
        )
        has_collate = result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        has_collate = False

    if has_collate:
        try:
            _collate_bam(input_bam, output_bam, threads)
            return
        except subprocess.CalledProcessError as e:
            logging.warning(
                f"samtools collate failed (exit {e.returncode}), "
                f"falling back to sort -n: {e.stderr}"
            )
    else:
        logging.warning(
            "samtools collate not available (requires samtools 1.6+), "
            "using sort -n"
        )

    # Fallback to sort -n
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(threads), "-o", output_bam, input_bam],
        check=True
    )
```

### Matplotlib Box Plot Matching Plotly
```python
# Source: Research synthesis from matplotlib/seaborn docs + Plotly color reference
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def generate_boxplot_matplotlib(
    reads_df: pd.DataFrame,
    output_path: str,
    title: str = "Plasmid Score Distribution"
) -> None:
    """Generate box plot using matplotlib with Plotly-style colors."""
    # Match Plotly default categorical colors
    PLOTLY_COLORS = ['#636EFA', '#EF553B', '#00CC96']

    # Set seaborn theme for clean style
    sns.set_theme(style="whitegrid")

    fig, ax = plt.subplots(figsize=(10, 4))

    # Box plot
    box_parts = ax.boxplot(
        [reads_df[reads_df['AssignedTo'] == cat]['PlasmidScore'].values
         for cat in ['Plasmid', 'Human', 'Tied']],
        labels=['Plasmid', 'Human', 'Tied'],
        patch_artist=True,
        showfliers=False  # We'll add points manually
    )

    # Color boxes to match Plotly
    for patch, color in zip(box_parts['boxes'], PLOTLY_COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # Add all points (matching Plotly points='all')
    for i, cat in enumerate(['Plasmid', 'Human', 'Tied'], 1):
        data = reads_df[reads_df['AssignedTo'] == cat]['PlasmidScore'].values
        x = np.random.normal(i, 0.04, size=len(data))  # Jitter
        ax.scatter(x, data, alpha=0.3, s=8, color=PLOTLY_COLORS[i-1])

    ax.set_title(f"{title}\n(Total Reads: {len(reads_df)})")
    ax.set_ylabel("Plasmid Score")
    ax.set_xlabel("Assigned To")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
```

### Heatmap with Viridis Colormap
```python
# Source: Research synthesis from seaborn/matplotlib docs
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def generate_heatmap_matplotlib(
    ratio_data: pd.DataFrame,
    output_path: str,
    threshold: float = 0.8
) -> None:
    """Generate heatmap using matplotlib with viridis colormap."""
    fig, ax = plt.subplots(figsize=(12, 12))

    # Viridis is matplotlib default and matches Plotly
    sns.heatmap(
        ratio_data,
        cmap="viridis",
        annot=False,  # No value annotations
        cbar=False,   # No colorbar (matching Plotly config)
        square=False,
        linewidths=0.5,
        linecolor='white',
        ax=ax
    )

    ax.set_title("Contamination Ratio Heatmap", fontsize=14)
    ax.set_xlabel("Plasmid", fontsize=12)
    ax.set_ylabel("Sample", fontsize=12)

    # Rotate x labels for readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
```

### Per-Combination Timing with Error Aggregation
```python
# Source: Research synthesis from Python error handling patterns
import time
import logging
from dataclasses import dataclass

@dataclass
class CombinationTiming:
    """Timing and status for a single combination."""
    label: str
    duration: float
    success: bool
    error_msg: str | None = None

def process_combinations_with_timing(
    combinations: list[tuple[str, Any]]
) -> list[CombinationTiming]:
    """Process combinations with per-combination timing and error resilience."""
    results = []

    for combo_label, combo_data in combinations:
        start = time.time()
        try:
            # Process combination
            _process_single_combination(combo_data)

            duration = time.time() - start
            results.append(CombinationTiming(
                label=combo_label,
                duration=duration,
                success=True
            ))
            logging.info(f"{combo_label}: {duration:.1f}s")

        except Exception as e:
            duration = time.time() - start
            results.append(CombinationTiming(
                label=combo_label,
                duration=duration,
                success=False,
                error_msg=str(e)
            ))
            logging.error(
                f"{combo_label} FAILED after {duration:.1f}s: {e}",
                exc_info=True
            )
            # Continue to next combination

    # Summary
    successful = sum(1 for r in results if r.success)
    failed = len(results) - successful
    total_time = sum(r.duration for r in results)

    logging.info(
        f"Batch complete: {successful}/{len(results)} succeeded "
        f"in {total_time:.1f}s total"
    )

    if failed > 0:
        logging.warning(f"{failed} combinations failed:")
        for r in results:
            if not r.success:
                logging.warning(f"  - {r.label}: {r.error_msg}")

    return results
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `samtools sort -n` for name grouping | `samtools collate` for read pairing | samtools 1.6 (2017) | 30-50% faster, but no ordering guarantees |
| Kaleido for Plotly PNG export | matplotlib for static plots | PlasmiCheck v0.31.0 Phase 5 | No licensing issues, pure Python, already a dependency |
| On-demand indexing in loop | Upfront planning phase | Common in modern pipelines (2020+) | Eliminates redundant operations in batch runs |
| Stop-on-first-error batch processing | Continue-on-failure with error aggregation | Python 3.11 ExceptionGroup (2022) | Better UX, efficient cluster usage |
| Manual temp file cleanup | Context manager with auto-cleanup | Python 3.2+ tempfile (2011) | Exception-safe, no leaked files |

**Deprecated/outdated:**
- **Kaleido for static Plotly export**: Licensing uncertainty, maintenance issues, headless browser overhead; replaced by matplotlib in Phase 5
- **Fast mode collate (`-f`)**: Filters supplementary alignments, unsuitable for chimeric read analysis; use standard mode only
- **Global index lock files**: Over-engineering for parallelism that doesn't exist in single-process pipeline; simple path tracking sufficient

## Open Questions

Things that couldn't be fully resolved:

1. **Collate performance gain magnitude**
   - What we know: Official docs say collate is "faster alternative" to sort -n; community reports 30-50% improvement
   - What's unclear: Exact speedup depends on read count, BAM size, disk I/O; no authoritative benchmark
   - Recommendation: Measure in Phase 4 regression tests with real data; log timing for both approaches; document actual speedup observed

2. **Supplementary alignment re-sorting overhead**
   - What we know: Explicit Python re-sort adds processing step after collate
   - What's unclear: Whether overhead negates collate speed gain for high-supplementary datasets
   - Recommendation: Profile with chimeric-heavy test data; if overhead >20%, consider conditional re-sort only when supplementary count is high

3. **Matplotlib plot visual parity threshold**
   - What we know: CONTEXT says "functionally equivalent, not pixel-identical" but doesn't define tolerance
   - What's unclear: How much visual difference is acceptable (e.g., font rendering, axis label positioning)
   - Recommendation: Manual visual inspection of representative plots; user testing; document known differences in PLAN.md

4. **Index deduplication value in practice**
   - What we know: Human reference indexing takes ~30s (Phase 6 benchmarks); plasmid indexes are <1s
   - What's unclear: Whether 30s savings in batch runs justifies upfront planning complexity
   - Recommendation: Calculate break-even point (N combinations where planning overhead < savings); log planning phase timing separately

## Sources

### Primary (HIGH confidence)
- [samtools-collate manual](http://www.htslib.org/doc/samtools-collate.html) - Official collate documentation
- [samtools-sort manual](http://www.htslib.org/doc/samtools-sort.html) - Official sort documentation including -n behavior
- [Python tempfile docs](https://docs.python.org/3/library/tempfile.html) - NamedTemporaryFile context manager
- [Python subprocess docs](https://docs.python.org/3/library/subprocess.html) - Error handling and timeout
- [matplotlib colormap docs](https://matplotlib.org/stable/gallery/color/colormap_reference.html) - Viridis and built-in colormaps
- [seaborn color palette docs](https://seaborn.pydata.org/tutorial/color_palettes.html) - Palette API and theme configuration
- [Plotly built-in color scales](https://plotly.com/python/builtin-colorscales/) - Viridis and default categorical colors

### Secondary (MEDIUM confidence)
- [KDnuggets: Error Handling Patterns in Python](https://www.kdnuggets.com/5-error-handling-patterns-in-python-beyond-try-except) - Error aggregation pattern
- [Medium: Exception Groups in Python 3.11](https://thelinuxcode.com/exception-groups-in-python-311-a-practical-guide-to-multi-error-handling) - ExceptionGroup usage
- [samtools GitHub issue #2252](https://github.com/samtools/samtools/issues/2252) - Collate vs sort -n ordering discussion
- [Towards Data Science: Matplotlib vs Plotly](https://towardsdatascience.com/matplotlib-vs-plotly-express-which-one-is-the-best-library-for-data-visualization-7a96dbe3ff09/) - Library comparison
- [API Error Handling Guide 2026](https://easyparser.com/blog/api-error-handling-retry-strategies-python-guide) - Retry and fallback patterns

### Tertiary (LOW confidence)
- [Biostars: samtools collate discussion](https://www.biostars.org/p/9566947/) - Community usage patterns (not official)
- [QualityPoint: Matplotlib vs Seaborn vs Plotly](https://www.blog.qualitypointtech.com/2026/02/matplotlib-vs-seaborn-vs-plotly.html) - General comparison (no technical depth)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official samtools docs, Python stdlib docs, existing pyproject.toml dependencies verified
- Architecture: HIGH - Patterns synthesized from official docs and existing PlasmiCheck codebase structure
- Pitfalls: MEDIUM - Based on documentation warnings and common bioinformatics forum discussions, not PlasmiCheck-specific experience
- Code examples: HIGH - Synthesized from official documentation examples, adapted to PlasmiCheck patterns
- Performance claims (30-50% collate speedup): MEDIUM - Community reports, not benchmarked on PlasmiCheck data yet

**Research date:** 2026-02-14
**Valid until:** ~2026-04-14 (60 days for stable technologies: samtools, matplotlib, Python stdlib)
