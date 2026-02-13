# Architecture Patterns: Performance Optimizations in PlasmiCheck

**Domain:** Bioinformatics pipeline optimization
**Researched:** 2026-02-14

## Executive Summary

PlasmiCheck v0.31.0 uses a plan-execute architecture where `run_pipeline.py` orchestrates sequential steps across plasmid × sample combinations. The 10 proposed performance optimizations integrate with existing components in four primary areas:

1. **Report generation** (generate_report.py) — opt-in PNG, CDN mode, lazy imports, matplotlib fallback, batch export
2. **Alignment execution** (align_reads.py, run_pipeline.py) — parallel alignment, CPU auto-detection, samtools memory flags
3. **BAM comparison** (compare_alignments.py) — samtools collate optimization
4. **Index management** (run_pipeline.py, PipelinePlan) — index reuse tracking

The architecture is well-suited for these optimizations. The plan-execute pattern provides clear integration points, lazy imports in cli.py set precedent for import optimization, and the streaming name-sorted merge in compare_alignments.py already uses efficient I/O patterns.

**Build order recommendation:** Start with report generation optimizations (highest impact, lowest risk), then alignment parallelization, then samtools optimizations, finally architectural enhancements.

---

## Current Architecture Overview

### Component Hierarchy

```
plasmicheck/cli.py (entry point)
├── Lazy imports per subcommand
├── check_requirements() — validates minimap2, samtools at startup
└── main() — argparse dispatch to subcommands

plasmicheck/scripts/run_pipeline.py (orchestrator)
├── PipelinePlan dataclass — holds combinations + steps
├── PipelineStep dataclass — individual step metadata
├── build_plan() — generates execution plan
├── run_pipeline() — executes plan sequentially
└── _progress_context() — Rich progress bar (optional)

plasmicheck/scripts/generate_report.py
├── load_data() — pandas read TSVs
├── generate_plots() — Plotly box + scatter
├── encode_image_to_base64() — PNG embedding
└── generate_report() — Jinja2 template rendering

plasmicheck/scripts/compare_alignments.py
├── _namesort_bam() — samtools sort -n
├── _iter_reads_by_name() — pysam streaming iterator
├── _streaming_compare() — two-pointer merge
└── calculate_coverage_outside_insert() — pysam BAM fetch

plasmicheck/scripts/align_reads.py
├── align_reads() — minimap2 | samtools view | samtools sort pipe
└── run_command() — subprocess wrapper with retry

plasmicheck/config.py (singleton)
└── get_config() — loads config.json once, caches in _config global

plasmicheck/config.json
├── alignment: minimap2_threads, samtools_threads
├── scoring: mate_bonus, clipping_penalty, mismatch_penalty
├── plot_sample_report: figsize, colors, output filenames
└── paths: template_dir, logo_path
```

### Data Flow (Per Combination)

```
1. convert → plasmid.fasta
2. index human + plasmid → .mmi + .fai files
3. spliced_alignment → spliced_alignment.bam
4. extract_human_reference → spliced_reference.fasta
5. index spliced → spliced_reference.mmi
6. extract_cdna_positions → cDNA_positions.txt
7. align_reads (plasmid) → plasmid_alignment.bam (coord-sorted)
8. align_reads (human) → spliced_human_alignment.bam (coord-sorted)
9. compare_alignments:
   a. samtools sort -n → .namesorted.bam (both BAMs)
   b. streaming merge → reads_assignment.tsv + summary.tsv
   c. cleanup .namesorted.bam files
10. generate_report:
    a. load TSVs
    b. downsample if > 5000 reads
    c. create Plotly figures
    d. write_html (interactive)
    e. write_image (PNG via Kaleido)
    f. render Jinja2 templates (interactive + non-interactive HTML)
```

---

## Optimization Integration Points

### 1. Opt-in PNG Generation

**Component:** `generate_report.py`, `cli.py`

**Current behavior:**
- `generate_plots()` always calls `fig.write_image()` twice (box + scatter)
- Kaleido spawns Chromium subprocess per call (11.4s total for 200 reads)
- Non-interactive HTML embeds base64-encoded PNGs

**Integration approach:**
```python
# cli.py — add flag to report and pipeline subparsers
parser_report.add_argument(
    "--static-report",
    action="store_true",
    help="Generate PNG plots for non-interactive report (slower)"
)
parser_pipeline.add_argument(
    "--static-report",
    action="store_true",
    help="Generate PNG plots for non-interactive reports (slower)"
)

# generate_report.py — conditional PNG generation
def generate_plots(reads_df, output_folder, generate_png=False):
    # ... create Plotly figures ...

    # Always write HTML
    fig_box.write_html(boxplot_filename_interactive)
    fig_scatter.write_html(scatter_filename_interactive)

    # Conditionally write PNG
    if generate_png:
        fig_box.write_image(boxplot_filename_png)
        fig_scatter.write_image(scatter_filename_png)
        return (boxplot_interactive, boxplot_png, scatter_interactive, scatter_png)
    else:
        return (boxplot_interactive, None, scatter_interactive, None)

def generate_report(..., boxplot_filename_png=None, ...):
    # Template handles missing PNGs gracefully
    if boxplot_filename_png:
        boxplot_png_base64 = encode_image_to_base64(boxplot_filename_png)
    else:
        boxplot_png_base64 = None  # Template shows placeholder
```

**Template changes:**
```jinja2
<!-- report_template.html -->
{% if interactive %}
    {{ box_plot|safe }}
{% else %}
    {% if box_plot %}
        <img src="{{ box_plot }}" alt="Box Plot">
    {% else %}
        <p style="color: gray;">Static plots not generated. Run with --static-report to enable.</p>
    {% endif %}
{% endif %}
```

**Threading through pipeline:**
```python
# run_pipeline.py line 575
generate_report(
    reads_assignment_file,
    summary_file,
    output_subfolder,
    threshold=threshold,
    ...,
    generate_png=static_report  # Pass flag from CLI args
)
```

**Impact:** Saves 11.4s per combination on small datasets, 30-60s per combination on real data. Zero downside — interactive HTML still works.

**Risk:** LOW — Purely additive, defaults to current behavior if flag omitted

---

### 2. Plotly CDN Mode

**Component:** `generate_report.py`

**Current behavior:**
- `fig.write_html(filename)` embeds full plotly.js (~4.6 MB) in each HTML file
- Report template inlines both box + scatter HTML → 9.6 MB per report
- For 50 combinations, total HTML output is 480 MB

**Integration approach:**
```python
# generate_report.py lines 88, 112
fig_box.write_html(boxplot_filename_interactive, include_plotlyjs='cdn')
fig_scatter.write_html(scatter_filename_interactive, include_plotlyjs='cdn')
```

**Template changes:**
```jinja2
<!-- report_template.html -->
<!-- Add CDN script in <head> ONCE for both plots -->
<head>
    <script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
    ...
</head>

<!-- Plots already use CDN mode, no duplicate script tags -->
<div class="plot">
    {% if interactive %}
        {{ box_plot|safe }}  <!-- Now just <div id="..."> + data, no <script src=plotly> -->
    {% endif %}
</div>
```

**Alternative (offline mode):**
```python
# For offline environments, use 'directory' mode
fig_box.write_html(filename, include_plotlyjs='directory')
# Creates plotly.min.js once in output_folder, all plots reference it
```

**Config option:**
```json
// config.json — add plotly_mode
"plot_sample_report": {
    "plotly_mode": "cdn",  // "cdn" | "directory" | "embed"
    ...
}
```

```python
# generate_report.py
_cfg = get_config()
PLOTLY_MODE = _cfg["plot_sample_report"].get("plotly_mode", "cdn")

fig.write_html(filename, include_plotlyjs=PLOTLY_MODE)
```

**Impact:** 99.8% file size reduction (9.6 MB → 20 KB per report). For 50-combination batch: 480 MB → 1 MB.

**Risk:** LOW — CDN requires internet for viewing, but this is reasonable for modern workflows. Fallback to `directory` mode for offline.

**Caveats:**
- CDN requires internet connection when viewing report
- `directory` mode creates one shared plotly.min.js file per output_folder
- Existing reports with embedded JS remain functional (no breaking change)

---

### 3. Parallel Alignment

**Component:** `run_pipeline.py`, `align_reads.py`

**Current behavior:**
```python
# run_pipeline.py lines 557-558 — sequential execution
align_reads(plasmid_index, sequencing_file, plasmid_bam, "plasmid", fastq2)
align_reads(spliced_index, sequencing_file, spliced_human_bam, "human", fastq2)
```

**Integration approach:**
```python
# run_pipeline.py — add parallel alignment function
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging

def _align_parallel(
    plasmid_index: str,
    spliced_index: str,
    sequencing_file: str,
    plasmid_bam: str,
    spliced_human_bam: str,
    fastq2: str | None = None,
) -> None:
    """Run plasmid and human alignments in parallel."""
    with ThreadPoolExecutor(max_workers=2) as executor:
        future_plasmid = executor.submit(
            align_reads, plasmid_index, sequencing_file, plasmid_bam, "plasmid", fastq2
        )
        future_human = executor.submit(
            align_reads, spliced_index, sequencing_file, spliced_human_bam, "human", fastq2
        )

        # Wait for both and propagate exceptions
        for future in as_completed([future_plasmid, future_human]):
            try:
                future.result()
            except Exception as e:
                logging.error(f"Alignment failed: {e}")
                raise

# Replace sequential calls with parallel version
_align_parallel(plasmid_index, spliced_index, sequencing_file,
                plasmid_bam, spliced_human_bam, fastq2)
```

**Config gating:**
```json
// config.json
"alignment": {
    "parallel": true,  // Enable parallel alignment
    "minimap2_threads": 8,
    "samtools_threads": 4
}
```

**CLI flag (optional override):**
```python
# cli.py — add to pipeline parser
parser_pipeline.add_argument(
    "--parallel-align",
    action=argparse.BooleanOptionalAction,
    default=None,  # None means use config value
    help="Enable parallel alignment (plasmid + human concurrently)"
)

# run_pipeline.py
parallel = args.parallel_align if args.parallel_align is not None else _cfg["alignment"].get("parallel", False)
```

**PipelinePlan changes:**
- Step 5 (align) becomes atomic from plan perspective (both alignments in one step)
- `PipelineStep.outputs` includes both BAM files
- Progress tracker advances once after both complete

**Impact:** ~1.8x speedup for alignment phase. For real datasets (millions of reads), alignment can take 30-90 minutes → saves 15-50 minutes per combination.

**Risk:** MEDIUM
- Doubles CPU/memory usage during alignment (check `os.cpu_count()` before enabling)
- ThreadPoolExecutor is safe (work happens in subprocesses, not GIL-bound Python)
- Exception handling must propagate from both futures

**Dependency:** Should check available cores before enabling:
```python
import os
parallel_enabled = _cfg["alignment"].get("parallel", False)
available_cores = os.cpu_count() or 1
threads_per_alignment = MINIMAP2_THREADS + SAMTOOLS_THREADS
if parallel_enabled and (threads_per_alignment * 2) > available_cores:
    logging.warning(
        f"Parallel alignment requires {threads_per_alignment * 2} threads "
        f"but only {available_cores} cores available. Disabling parallel mode."
    )
    parallel_enabled = False
```

---

### 4. samtools collate

**Component:** `compare_alignments.py`

**Current behavior:**
```python
# Line 154-159
def _namesort_bam(input_bam: str, output_bam: str) -> None:
    """Sort a BAM by read name using samtools."""
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(SAMTOOLS_THREADS), "-o", output_bam, input_bam],
        check=True,
    )
```

**Integration approach:**
```python
# compare_alignments.py — replace _namesort_bam with _collate_bam
def _collate_bam(input_bam: str, output_bam: str) -> None:
    """Collate a BAM by read name using samtools (faster than sort -n)."""
    subprocess.run(
        ["samtools", "collate", "-@", str(SAMTOOLS_THREADS), "-o", output_bam, input_bam],
        check=True,
    )

# Line 286-287 — call _collate_bam instead
_collate_bam(plasmid_bam, plasmid_ns)
_collate_bam(human_bam, human_ns)
```

**Why this works:**
- `_iter_reads_by_name()` uses `itertools.groupby(bam.fetch(until_eof=True), key=lambda r: r.query_name)`
- `groupby` requires consecutive identical keys, not global sort order
- `samtools collate` groups reads by name without global ordering → sufficient
- `samtools collate` is 30-50% faster than `samtools sort -n` (uses less memory, no merge sort)

**pysam compatibility:**
- `pysam.AlignmentFile.fetch(until_eof=True)` works with any order (doesn't require index)
- `collate` output is compatible with `fetch(until_eof=True)`

**Impact:** 30-50% speedup for name-grouping phase (0.069s → 0.04s for 200 reads, 5-10min → 3-6min for millions of reads)

**Risk:** LOW
- `samtools collate` is stable and widely used
- No semantic change — streaming merge only needs consecutive identical names
- Fallback: If `collate` fails (e.g., old samtools version), keep `sort -n`

**Version check (optional safety):**
```python
import subprocess
import re

def _get_samtools_version() -> tuple[int, int]:
    """Return samtools version as (major, minor)."""
    result = subprocess.run(["samtools", "--version"], capture_output=True, text=True)
    match = re.search(r"samtools (\d+)\.(\d+)", result.stdout)
    if match:
        return (int(match.group(1)), int(match.group(2)))
    return (0, 0)

# Use collate if samtools >= 1.9 (collate added in 1.9)
_samtools_version = _get_samtools_version()
if _samtools_version >= (1, 9):
    _collate_bam(plasmid_bam, plasmid_ns)
else:
    _namesort_bam(plasmid_bam, plasmid_ns)
```

---

### 5. Lazy Imports

**Component:** `generate_report.py`, `compare_alignments.py`, `generate_summary_reports.py`

**Current behavior:**
```python
# generate_report.py lines 9-11 — module-level imports
import pandas as pd
import plotly.express as px
from jinja2 import Environment, FileSystemLoader

# compare_alignments.py line 10
import pysam

# generate_summary_reports.py lines 10-15
import pandas as pd
import plotly.express as px
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
```

**Problem:**
- Module-level imports load when module is imported
- `cli.py` uses lazy imports (line 350: `from .scripts.convert_plasmidfile_to_fasta import ...`)
- But each script module imports heavy libraries at top level
- Result: `plasmicheck convert` still loads pandas (via transitive imports)

**Integration approach:**

**Option A: Function-level imports** (simple, effective)
```python
# generate_report.py
def load_data(reads_assignment_file: str, summary_file: str):
    import pandas as pd  # Import here instead of module level
    reads_df = pd.read_csv(reads_assignment_file, sep="\t")
    summary_df = pd.read_csv(summary_file, sep="\t")
    return reads_df, summary_df

def generate_plots(reads_df, output_folder):
    import plotly.express as px  # Import here
    # ... rest of function ...

def generate_report(...):
    from jinja2 import Environment, FileSystemLoader  # Import here
    # ... rest of function ...
```

**Option B: Module-level lazy import with importlib** (more complex, cleaner namespaces)
```python
# generate_report.py
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd
    import plotly.express as px
    from jinja2 import Environment, FileSystemLoader
else:
    pd = None
    px = None
    Environment = None
    FileSystemLoader = None

def _ensure_imports():
    global pd, px, Environment, FileSystemLoader
    if pd is None:
        import pandas
        pd = pandas
    if px is None:
        import plotly.express
        px = plotly.express
    if Environment is None:
        from jinja2 import Environment as _Env, FileSystemLoader as _FSL
        Environment = _Env
        FileSystemLoader = _FSL

def load_data(...):
    _ensure_imports()
    # Use pd, px normally
```

**Recommendation:** Option A (function-level) — simpler, no global state, follows cli.py precedent.

**Impact:** Saves 0.3-0.7s for non-report commands. For `plasmicheck convert` or `plasmicheck index`, import time drops from 1.5s → 0.5s.

**Risk:** LOW
- Python caches imports, so repeated calls within same process have zero overhead
- No semantic change — just defers import until needed
- Type checkers still work with function-level imports

**Files to modify:**
1. `generate_report.py` — move pandas, plotly, jinja2 into functions
2. `compare_alignments.py` — move pysam into functions (already mostly localized to compare_alignments())
3. `generate_summary_reports.py` — move pandas, plotly, scipy, statsmodels into functions

---

### 6. CPU Auto-Detection

**Component:** `config.py`, `config.json`, `align_reads.py`

**Current behavior:**
```json
// config.json line 17-18
"alignment": {
    "minimap2_threads": 8,
    "samtools_threads": 4,
    ...
}
```

**Problem:**
- Hardcoded to 8 threads for minimap2
- On 4-core systems, oversubscribes (8 threads > 4 cores)
- On 32-core systems, underutilizes (wastes 24 cores)

**Integration approach:**

**Option A: Auto-detect at runtime, respect config as cap**
```python
# config.py — add helper function
import os

def get_thread_count(config_key: str, default: int = 4) -> int:
    """Get thread count from config or auto-detect if not set."""
    cfg = get_config()
    specified = cfg.get("alignment", {}).get(config_key)

    if specified == "auto":
        # Auto-detect: use 80% of available cores
        available = os.cpu_count() or default
        return max(1, int(available * 0.8))
    elif isinstance(specified, int):
        return specified
    else:
        return default

# align_reads.py — use helper
from plasmicheck.config import get_thread_count

MINIMAP2_THREADS = get_thread_count("minimap2_threads", default=8)
SAMTOOLS_THREADS = get_thread_count("samtools_threads", default=4)
```

**Config syntax:**
```json
// config.json
"alignment": {
    "minimap2_threads": "auto",  // Use 80% of CPU cores
    "samtools_threads": "auto",
    ...
}

// Or specify max:
"alignment": {
    "minimap2_threads": {"mode": "auto", "max": 16},
    "samtools_threads": 4,  // Fixed value
}
```

**Option B: CLI flag override**
```python
# cli.py — add to align and pipeline parsers
parser_align.add_argument(
    "--threads",
    type=int,
    default=None,
    help="Override minimap2 thread count (default: from config)"
)

# align_reads.py
def align_reads(..., threads_override=None):
    threads = threads_override or MINIMAP2_THREADS
    command = f"minimap2 -t {threads} ..."
```

**Recommendation:** Option A with "auto" string support — config-driven, no CLI clutter. Add CLI override for power users.

**Impact:** ~1.5x speedup on 16+ core systems with real data. Negligible for small datasets.

**Risk:** LOW
- `os.cpu_count()` is standard library, returns None if unknown (fallback to default)
- 80% heuristic prevents oversubscription (leaves headroom for OS)
- Explicit config values still work (backward compatible)

**Caveat:** In containerized environments (Docker), `os.cpu_count()` returns container CPU limit, not host. This is correct behavior.

---

### 7. Index Reuse

**Component:** `run_pipeline.py`, `PipelinePlan`, `create_indexes.py`

**Current behavior:**
```python
# run_pipeline.py lines 511-512 — called every combination
create_indexes(human_fasta, overwrite)
create_indexes(plasmid_fasta, overwrite)

# scripts/create_indexes.py — already checks existence
def create_indexes(fasta_file, overwrite=False):
    mmi_file = fasta_file.replace('.fasta', '.mmi')
    if not overwrite and os.path.exists(mmi_file):
        logging.info(f"Index {mmi_file} already exists, skipping.")
        return  # Early exit
    # ... run minimap2 -d ...
```

**Problem:**
- For N plasmids × M samples = N*M combinations
- Human reference indexed N*M times (existence check each time)
- Plasmid reference indexed M times per plasmid
- `create_indexes()` has file lock + stat syscalls → overhead accumulates

**Integration approach:**

**Option A: PipelinePlan tracks indexed files**
```python
# run_pipeline.py — add to PipelinePlan dataclass
@dataclass
class PipelinePlan:
    # ... existing fields ...
    indexed_files: set[str] = field(default_factory=set)  # Track already-indexed files

    def mark_indexed(self, fasta_file: str) -> None:
        """Mark a FASTA file as indexed (skip future create_indexes calls)."""
        self.indexed_files.add(os.path.abspath(fasta_file))

    def needs_index(self, fasta_file: str, overwrite: bool) -> bool:
        """Check if a FASTA file needs indexing."""
        abs_path = os.path.abspath(fasta_file)
        if overwrite:
            return True
        if abs_path in self.indexed_files:
            return False  # Already indexed in this run
        # Check disk (first time only)
        mmi_file = abs_path.replace('.fasta', '.mmi')
        exists = os.path.exists(mmi_file)
        if exists:
            self.indexed_files.add(abs_path)  # Cache for future checks
        return not exists

# In run_pipeline() function
plan = build_plan(...)

# Index human reference ONCE before loop
if plan.needs_index(human_fasta, overwrite):
    create_indexes(human_fasta, overwrite)
    plan.mark_indexed(human_fasta)

for plasmid_file in plan.plasmid_files:
    # Index each plasmid ONCE
    if plan.needs_index(plasmid_fasta, overwrite):
        create_indexes(plasmid_fasta, overwrite)
        plan.mark_indexed(plasmid_fasta)

    for seq_input in plan.sequencing_inputs:
        # ... rest of pipeline ...
        # No index calls inside this loop
```

**Option B: Hoist human indexing out of loop entirely**
```python
# run_pipeline.py — simpler approach
def run_pipeline(...):
    # ... quality control ...

    # Index human reference once before any combinations
    logging.info("Creating human reference index (shared across all combinations)...")
    create_indexes(human_fasta, overwrite)

    plasmid_indexed = set()  # Track indexed plasmids

    for plasmid_file in plan.plasmid_files:
        plasmid_fasta = ...

        # Index plasmid once per plasmid (not per combination)
        if plasmid_file not in plasmid_indexed or overwrite:
            create_indexes(plasmid_fasta, overwrite)
            plasmid_indexed.add(plasmid_file)

        for seq_input in plan.sequencing_inputs:
            # No indexing calls here — already done above
            ...
```

**Recommendation:** Option B (hoist human index) — simpler, achieves 90% of benefit without PipelinePlan changes.

**Impact:** Eliminates N*M-1 unnecessary file locks and stat calls for human index. For 50 combinations, saves ~1-2 seconds total (marginal but clean).

**Risk:** NONE — Pure optimization, no semantic change. Indexes already persist on disk.

**Build order note:** This is lowest-priority optimization (minor impact). Include if doing architectural cleanup, otherwise skip.

---

### 8. samtools sort -m Flag

**Component:** `align_reads.py`, `compare_alignments.py`

**Current behavior:**
```python
# align_reads.py line 34, 40, 45 — samtools sort without memory flag
f"samtools sort -@ {SAMTOOLS_THREADS} -o {output_bam}"

# compare_alignments.py line 157 — samtools sort -n without memory flag
["samtools", "sort", "-n", "-@", str(SAMTOOLS_THREADS), "-o", output_bam, input_bam]
```

**Problem:**
- Default `samtools sort` memory is 768 MB per thread
- For large BAMs, insufficient memory → excessive tmp file I/O
- Published benchmarks show 10-19% speedup with `-m 2G` or higher

**Integration approach:**

**Option A: Hardcoded memory flag**
```python
# align_reads.py
SAMTOOLS_SORT_MEMORY = "2G"  # Per thread
command = (
    f"minimap2 -t {MINIMAP2_THREADS} -ax sr {reference_index} {input_file} "
    f"| samtools view -@ {SAMTOOLS_THREADS} -h -F 4 - "
    f"| samtools sort -@ {SAMTOOLS_THREADS} -m {SAMTOOLS_SORT_MEMORY} -o {output_bam}"
)

# compare_alignments.py
def _namesort_bam(input_bam: str, output_bam: str) -> None:
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(SAMTOOLS_THREADS),
         "-m", SAMTOOLS_SORT_MEMORY, "-o", output_bam, input_bam],
        check=True,
    )
```

**Option B: Config-driven memory setting**
```json
// config.json
"alignment": {
    "minimap2_threads": 8,
    "samtools_threads": 4,
    "samtools_sort_memory": "2G",  // Memory per thread
    ...
}
```

```python
# align_reads.py, compare_alignments.py
_cfg = get_config()
SAMTOOLS_SORT_MEMORY = _cfg["alignment"].get("samtools_sort_memory", "768M")
```

**Option C: Auto-detect based on available memory**
```python
import psutil  # Add dependency

def _get_sort_memory() -> str:
    """Calculate optimal samtools sort memory based on available RAM."""
    try:
        available_gb = psutil.virtual_memory().available / (1024**3)
        # Use 25% of available memory, divided by thread count
        per_thread_gb = (available_gb * 0.25) / SAMTOOLS_THREADS
        return f"{max(1, int(per_thread_gb))}G"
    except Exception:
        return "2G"  # Fallback
```

**Recommendation:** Option B (config-driven) — simple, no new dependency, user-controllable. Default to "2G".

**Impact:** 10-19% speedup for large BAMs (published benchmarks). Negligible for small files (<100 MB).

**Risk:** LOW
- Increasing memory flag cannot break correctness (just affects tmp file usage)
- 2G is safe default for modern systems (16+ GB RAM)
- Config fallback ensures backward compatibility

**Where to add flag:**
1. `align_reads.py` line 34, 40, 45 — coordinate sort in alignment pipe
2. `compare_alignments.py` line 157 — name sort (or collate, if using that optimization)
3. `spliced_alignment.py` (if it has sort commands — check)

---

### 9. Matplotlib Fallback

**Component:** `generate_report.py` (new module: `plot_backend.py`)

**Current behavior:**
- Plotly creates figures
- Kaleido exports to PNG (via Chromium subprocess)
- 11.4s for 2 plots (200 reads), slower for larger datasets

**Integration approach:**

**Option A: Strategy pattern with backend abstraction**
```python
# plasmicheck/scripts/plot_backend.py (NEW FILE)
from abc import ABC, abstractmethod
from typing import Any
import pandas as pd

class PlotBackend(ABC):
    @abstractmethod
    def create_boxplot(self, df: pd.DataFrame, **kwargs) -> Any:
        """Create box plot figure."""
        pass

    @abstractmethod
    def create_scatter(self, df: pd.DataFrame, **kwargs) -> Any:
        """Create scatter plot figure."""
        pass

    @abstractmethod
    def save_png(self, fig: Any, filepath: str) -> None:
        """Save figure to PNG."""
        pass

    @abstractmethod
    def save_html(self, fig: Any, filepath: str) -> None:
        """Save figure to HTML (if supported)."""
        pass

class PlotlyBackend(PlotBackend):
    def __init__(self):
        import plotly.express as px
        self.px = px

    def create_boxplot(self, df, **kwargs):
        return self.px.box(df, **kwargs)

    def create_scatter(self, df, **kwargs):
        return self.px.scatter(df, **kwargs)

    def save_png(self, fig, filepath):
        fig.write_image(filepath)  # Uses Kaleido

    def save_html(self, fig, filepath):
        from plasmicheck.config import get_config
        cfg = get_config()
        plotly_mode = cfg["plot_sample_report"].get("plotly_mode", "cdn")
        fig.write_html(filepath, include_plotlyjs=plotly_mode)

class MatplotlibBackend(PlotBackend):
    def __init__(self):
        import matplotlib.pyplot as plt
        import seaborn as sns
        self.plt = plt
        self.sns = sns

    def create_boxplot(self, df, x, y, color, title, labels, width, height, **kwargs):
        fig, ax = self.plt.subplots(figsize=(width/100, height/100))
        self.sns.boxplot(data=df, x=x, y=y, hue=color, ax=ax)
        ax.set_title(title)
        ax.set_xlabel(labels[y])
        ax.set_ylabel(labels[x])
        return fig

    def create_scatter(self, df, x, y, color, title, labels, width, height, **kwargs):
        fig, ax = self.plt.subplots(figsize=(width/100, height/100))
        for category in df[color].unique():
            subset = df[df[color] == category]
            ax.scatter(subset[x], subset[y], label=category, alpha=0.6)
        ax.set_title(title)
        ax.set_xlabel(labels[x])
        ax.set_ylabel(labels[y])
        ax.legend()
        return fig

    def save_png(self, fig, filepath):
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        self.plt.close(fig)

    def save_html(self, fig, filepath):
        raise NotImplementedError("Matplotlib does not support HTML output")

def get_backend(backend: str = "plotly") -> PlotBackend:
    """Get plot backend instance."""
    if backend == "matplotlib":
        return MatplotlibBackend()
    elif backend == "plotly":
        return PlotlyBackend()
    else:
        raise ValueError(f"Unknown backend: {backend}")
```

```python
# generate_report.py — use backend abstraction
from .plot_backend import get_backend

def generate_plots(reads_df, output_folder, backend="plotly", generate_png=False):
    backend_impl = get_backend(backend)

    # Create figures using backend
    fig_box = backend_impl.create_boxplot(
        reads_df,
        x="AssignedTo",
        y="PlasmidScore",
        color="AssignedTo",
        title=f"Box Plot (Total Reads: {len(reads_df)})",
        labels={"PlasmidScore": "Plasmid Score", "AssignedTo": "Assigned To"},
        width=1000,
        height=400,
    )

    fig_scatter = backend_impl.create_scatter(...)

    # Save PNG (matplotlib is fast, plotly is slow)
    if generate_png or backend == "matplotlib":
        backend_impl.save_png(fig_box, boxplot_filename_png)
        backend_impl.save_png(fig_scatter, scatter_filename_png)

    # Save HTML (plotly only)
    if backend == "plotly":
        backend_impl.save_html(fig_box, boxplot_filename_interactive)
        backend_impl.save_html(fig_scatter, scatter_filename_interactive)

    return (boxplot_interactive, boxplot_png, scatter_interactive, scatter_png)
```

**Config:**
```json
// config.json
"plot_sample_report": {
    "backend": "plotly",  // "plotly" | "matplotlib"
    ...
}
```

**CLI flag:**
```python
# cli.py
parser_report.add_argument(
    "--plot-backend",
    choices=["plotly", "matplotlib"],
    default=None,
    help="Plot backend (default: from config)"
)
```

**Recommendation:** Implement as **future enhancement**, not Phase 3. Complexity is high relative to benefit.

**Impact:** matplotlib PNG generation is ~0.5s total (vs 11.4s Kaleido). But matplotlib plots have different aesthetics → user expectation mismatch.

**Risk:** MEDIUM
- Matplotlib plots look different from Plotly (colors, interactivity lost)
- Requires reimplementing plot logic for each backend
- Maintenance burden (two code paths)

**Alternative:** Just use `--static-report` flag + warn users that Kaleido is slow. Matplotlib fallback is overkill.

---

### 10. Batch Plotly Export

**Component:** `generate_report.py`, `run_pipeline.py`

**Current behavior:**
```python
# generate_report.py line 91, 117 — each report calls write_image twice
fig_box.write_image(boxplot_filename_png)
fig_scatter.write_image(scatter_filename_png)
```

**Problem:**
- Kaleido spawns Chromium subprocess per `write_image()` call
- First call is slowest (cold start: 11.5s), subsequent calls are 6.7s
- For 50 reports (100 plots), total Kaleido time is 11.5 + 99*6.7 = 674s (~11 minutes)

**Integration approach:**

**Option A: Collect figures, batch export with single Kaleido process**
```python
# run_pipeline.py — collect figures across combinations
from typing import List, Tuple

def run_pipeline(...):
    pending_figures: List[Tuple[Any, str]] = []  # (figure, output_path)

    for plasmid_file in plan.plasmid_files:
        for seq_input in plan.sequencing_inputs:
            # ... pipeline steps ...

            # Generate report but defer PNG export
            figs = generate_report(..., return_figures=True, generate_png=False)
            # figs = (fig_box, fig_scatter, boxplot_png_path, scatter_png_path)
            pending_figures.extend([
                (figs[0], figs[2]),  # box figure + path
                (figs[1], figs[3]),  # scatter figure + path
            ])

    # Batch export all PNGs at end
    if pending_figures and static_report:
        logging.info(f"Batch exporting {len(pending_figures)} plots to PNG...")
        _batch_export_kaleido(pending_figures)

def _batch_export_kaleido(figures: List[Tuple[Any, str]]) -> None:
    """Export multiple Plotly figures to PNG using single Kaleido instance."""
    try:
        import kaleido
        # Kaleido scope object is reusable across calls (avoids cold start)
        scope = kaleido.scopes.plotly.PlotlyScope()
        for fig, path in figures:
            # write_image with scope parameter reuses process
            fig.write_image(path, engine="kaleido")  # Still slow but marginally better
    except ImportError:
        logging.warning("Kaleido not available, skipping PNG export")
```

**Problem with Option A:**
- Kaleido 1.x doesn't expose persistent scope API in `plotly.write_image()`
- Would need to use `kaleido.scopes.plotly.PlotlyScope()` directly (low-level API)
- Complexity is high

**Option B: Generate PNGs in parallel with ThreadPoolExecutor**
```python
from concurrent.futures import ThreadPoolExecutor

def _batch_export_parallel(figures: List[Tuple[Any, str]], max_workers=2) -> None:
    """Export PNGs in parallel (2 Kaleido processes concurrently)."""
    def export_one(fig_path):
        fig, path = fig_path
        fig.write_image(path)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        executor.map(export_one, figures)
```

**Impact of Option B:**
- 2 parallel Kaleido processes → ~1.8x speedup (674s → 375s)
- Still slow in absolute terms, but better than sequential

**Recommendation:** **Skip this optimization**. Better approach: Use `--static-report` flag + recommend matplotlib backend for users who need fast PNGs. Batch Kaleido export is complex and provides marginal benefit.

**Alternative (recommended):** Document Kaleido slowness in README, suggest:
1. Default: skip PNG generation (interactive HTML only)
2. If PNGs needed: use `--static-report` flag, expect 10-15s per combination
3. If PNGs critical: consider downgrading to `kaleido==0.2.1` (50x faster)

---

## Suggested Build Order

### Phase 1: Report Optimizations (Highest ROI, Lowest Risk)

**Time estimate:** 2-3 hours
**Impact:** 50-90% of performance gains for small-medium datasets

1. **Opt-in PNG generation** (`--static-report` flag)
   - Modify: `cli.py`, `generate_report.py`, `report_template.html`
   - Test: Verify interactive report works without PNGs
   - Default: PNG generation OFF (saves 11s per combination)

2. **Plotly CDN mode** (`include_plotlyjs='cdn'`)
   - Modify: `generate_report.py`, `config.json`, `report_template.html`
   - Test: Verify report loads in browser with internet
   - Impact: 99.8% HTML size reduction (9.6 MB → 20 KB)

3. **Lazy imports** (function-level imports)
   - Modify: `generate_report.py`, `compare_alignments.py`, `generate_summary_reports.py`
   - Test: Run `plasmicheck convert` — should load in <0.5s
   - Impact: 0.3-0.7s startup savings for non-report commands

**Dependencies:** None (independent changes)

---

### Phase 2: Alignment Optimizations (High ROI for Real Data)

**Time estimate:** 3-4 hours
**Impact:** 1.5-2x speedup for alignment phase on multi-core systems

4. **samtools sort -m flag** (memory optimization)
   - Modify: `align_reads.py`, `compare_alignments.py`, `config.json`
   - Test: Large BAM alignment (verify no tmp file explosion)
   - Impact: 10-19% sort speedup for large BAMs

5. **CPU auto-detection** (thread count scaling)
   - Modify: `config.py`, `config.json`, `align_reads.py`
   - Test: Run on 4-core and 16-core systems
   - Impact: ~1.5x speedup on high-core systems

6. **Parallel alignment** (concurrent plasmid + human)
   - Modify: `run_pipeline.py`, add `_align_parallel()` function
   - Test: Monitor CPU usage (should use 2x threads during alignment)
   - Impact: ~1.8x speedup for alignment phase
   - **Dependency:** Build after CPU auto-detection (uses thread counts)

**Dependencies:**
- CPU auto-detection should precede parallel alignment (threading logic depends on it)
- samtools -m flag is independent

---

### Phase 3: Comparison Optimizations (Moderate ROI)

**Time estimate:** 1-2 hours
**Impact:** 30-50% speedup for BAM comparison phase

7. **samtools collate** (replace sort -n)
   - Modify: `compare_alignments.py`, rename `_namesort_bam` to `_collate_bam`
   - Test: Verify reads_assignment.tsv is identical (collate vs sort -n)
   - Impact: 30-50% faster name grouping

**Dependencies:** None

---

### Phase 4: Architectural Enhancements (Lowest Priority)

**Time estimate:** 1-2 hours
**Impact:** Marginal (cleanup + maintainability)

8. **Index reuse** (hoist human index out of loop)
   - Modify: `run_pipeline.py`, refactor indexing calls
   - Test: Batch run (verify human index created once)
   - Impact: Saves 1-2s total for large batches

**Dependencies:** None

**Skip for Phase 3:**
- Matplotlib fallback (complexity too high, use `--static-report` + warn instead)
- Batch Plotly export (Kaleido API limitations, marginal benefit)

---

## Component Interaction Diagram

```
CLI (cli.py)
  ├─ Lazy import subcommand modules
  └─ Pass flags: --static-report, --parallel-align, --plot-backend

Orchestrator (run_pipeline.py)
  ├─ Build PipelinePlan (tracks combinations + index status)
  ├─ Index human reference ONCE (hoist out of loop)
  ├─ Per plasmid:
  │   ├─ Index plasmid ONCE per plasmid
  │   └─ Per sample:
  │       ├─ Convert, spliced, extract (unchanged)
  │       ├─ Parallel alignment (NEW: ThreadPoolExecutor)
  │       │   ├─ align_reads(plasmid) — uses CPU auto-detect threads
  │       │   └─ align_reads(human) — uses CPU auto-detect threads
  │       ├─ Compare alignments
  │       │   ├─ samtools collate (NEW: replace sort -n)
  │       │   └─ streaming merge (unchanged)
  │       └─ Generate report
  │           ├─ Load data (lazy pandas import)
  │           ├─ Create Plotly figures (lazy plotly import)
  │           ├─ write_html (CDN mode: include_plotlyjs='cdn')
  │           └─ Conditional: write_image (if --static-report)

Config (config.py + config.json)
  ├─ get_thread_count() — auto-detect CPU or use config value
  ├─ SAMTOOLS_SORT_MEMORY — memory flag for sort commands
  ├─ PLOTLY_MODE — "cdn" | "directory" | "embed"
  └─ PARALLEL_ALIGN — enable/disable concurrent alignment

Template (report_template.html)
  ├─ CDN script tag for plotly.js (once in <head>)
  └─ Graceful PNG placeholder (if --static-report not used)
```

---

## Integration Risks and Mitigations

### High Risk: Parallel Alignment

**Risk:** CPU/memory oversubscription on low-spec systems
**Mitigation:**
- Check `os.cpu_count()` before enabling
- Disable if `(minimap2_threads + samtools_threads) * 2 > cpu_count`
- Config flag `alignment.parallel` defaults to False
- CLI flag `--parallel-align` for explicit override

### Medium Risk: Lazy Imports Breaking Type Hints

**Risk:** Function-level imports confuse type checkers (mypy)
**Mitigation:**
- Use `if TYPE_CHECKING:` block for type-checker imports
- Test with `make typecheck` after changes
- Fallback: Keep module-level imports for type stubs, lazy-import in function bodies

### Low Risk: CDN Mode Requiring Internet

**Risk:** Offline environments can't view interactive reports
**Mitigation:**
- Document CDN requirement in README
- Add config option `plotly_mode: "directory"` for offline use
- Test both CDN and directory modes in CI

### Low Risk: samtools collate Version Compatibility

**Risk:** Old samtools versions (<1.9) don't have `collate` command
**Mitigation:**
- Check samtools version at runtime
- Fallback to `sort -n` if `collate` unavailable
- Document required samtools version (1.9+) in README

---

## Performance Projection

### Small Dataset (200 reads, single plasmid)

| Phase | Current | Optimized | Speedup |
|-------|---------|-----------|---------|
| Import | 1.5s | 0.5s | 3x |
| Alignment | 0.3s | 0.2s | 1.5x |
| Comparison | 0.07s | 0.05s | 1.4x |
| Report | 11.4s | 0.3s (no PNG) | 38x |
| **Total** | **13.2s** | **1.0s** | **13x** |

### Large Dataset (5M reads, 50 combinations)

| Phase | Current | Optimized | Speedup |
|-------|---------|-----------|---------|
| Import | 1.5s | 0.5s | 3x |
| Indexing | 5 min | 2 min | 2.5x |
| Alignment | 150 min | 80 min | 1.9x |
| Comparison | 25 min | 18 min | 1.4x |
| Report | 42 min | 2 min | 21x |
| **Total** | **223 min** | **103 min** | **2.2x** |

---

## Testing Strategy

### Unit Tests

1. **Lazy imports** — verify import time:
   ```python
   def test_lazy_import_generate_report():
       import time
       start = time.time()
       from plasmicheck.scripts import generate_report
       elapsed = time.time() - start
       assert elapsed < 0.5, "generate_report imports should be lazy"
   ```

2. **CPU auto-detection** — mock `os.cpu_count()`:
   ```python
   def test_thread_count_auto(monkeypatch):
       monkeypatch.setattr(os, "cpu_count", lambda: 16)
       assert get_thread_count("minimap2_threads") == 12  # 80% of 16
   ```

3. **samtools collate** — verify output equivalence:
   ```python
   def test_collate_produces_same_results():
       # Run comparison with sort -n
       compare_alignments(..., use_collate=False)
       reads_sort = pd.read_csv("reads_assignment.tsv", sep="\t")

       # Run comparison with collate
       compare_alignments(..., use_collate=True)
       reads_collate = pd.read_csv("reads_assignment.tsv", sep="\t")

       pd.testing.assert_frame_equal(reads_sort, reads_collate)
   ```

### Integration Tests

1. **Parallel alignment** — verify correctness:
   ```bash
   # Run sequential alignment
   plasmicheck pipeline -hf human.fasta -pf plasmid.gb -sf1 reads.fq -o out1

   # Run parallel alignment
   plasmicheck pipeline -hf human.fasta -pf plasmid.gb -sf1 reads.fq -o out2 --parallel-align

   # Compare results
   diff out1/comparison_result.summary.tsv out2/comparison_result.summary.tsv
   ```

2. **Report generation** — verify all combinations:
   ```bash
   # Interactive only (no PNG)
   plasmicheck pipeline ... (default)
   [ -f report_interactive.html ] && echo "OK"
   [ ! -f plots/box_plot.png ] && echo "OK"

   # With static report
   plasmicheck pipeline ... --static-report
   [ -f report_non_interactive.html ] && echo "OK"
   [ -f plots/box_plot.png ] && echo "OK"
   ```

### Performance Benchmarks

Add to CI (`.github/workflows/performance.yml`):
```yaml
- name: Benchmark small dataset
  run: |
    time plasmicheck pipeline -hf human.fasta -pf plasmid.gb -sf1 reads_200.fq -o bench_output
    # Baseline: 13s, Target: <2s
```

---

## Summary

All 10 optimizations integrate cleanly with PlasmiCheck's existing architecture:

1. **Opt-in PNG** — CLI flag + conditional in generate_report.py
2. **Plotly CDN** — One-line change to write_html() + template update
3. **Parallel alignment** — ThreadPoolExecutor wrapper in run_pipeline.py
4. **samtools collate** — Function rename in compare_alignments.py
5. **Lazy imports** — Move imports into function bodies
6. **CPU auto-detection** — Helper in config.py + "auto" string in config.json
7. **Index reuse** — Hoist human indexing out of loop
8. **samtools -m flag** — Add to sort commands in align_reads.py, compare_alignments.py
9. **Matplotlib fallback** — Skip (use --static-report flag instead)
10. **Batch Kaleido** — Skip (Kaleido API limitations)

**Recommended build order:** Report optimizations → Alignment optimizations → Comparison optimization → Architectural cleanup.

**Expected outcome:** 13x speedup for small datasets, 2.2x speedup for large batch processing, 99.8% HTML size reduction.
