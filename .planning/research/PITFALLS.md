# Domain Pitfalls: Performance Optimization in Bioinformatics Python Tools

**Domain:** Adding performance optimizations to existing bioinformatics Python CLI tool
**Project:** PlasmiCheck v0.31.0 performance improvements
**Researched:** 2026-02-14
**Overall Confidence:** HIGH

---

## Critical Pitfalls

Mistakes that cause rewrites, broken functionality, or data corruption.

### Pitfall 1: Kaleido Downgrade Breaking Plotly 6.x Compatibility

**What goes wrong:**

Downgrading from Kaleido v1.2.0 to v0.2.1 for the 50-100x performance improvement breaks compatibility with Plotly 6.x. Kaleido v0.2.1 only works with Plotly 5.x, while PlasmiCheck currently uses Plotly 6.3.0. Attempting to use the incompatible versions results in runtime errors during `fig.write_image()` calls.

**Why it happens:**

Kaleido v1.0.0+ introduced breaking API changes and requires Plotly v6.1.1 or later. The v0.2.1 architecture (persistent background process) was completely rewritten in v1.x (Chrome subprocess per invocation), making them mutually incompatible. Developers see the dramatic performance difference (20ms vs 2000ms per plot) and downgrade Kaleido without checking Plotly version constraints.

**Consequences:**

- Report generation fails completely in the pipeline
- Users with existing Plotly 6.x installations get cryptic errors
- CI/CD pipelines break when dependency resolution picks incompatible versions
- Conda/pip may silently install incompatible combinations if constraints aren't explicit

**Prevention:**

```python
# In pyproject.toml, make the constraint explicit:
dependencies = [
    # Option A: Stay on Plotly 6.x, accept Kaleido slowness
    "plotly>=6.1.1",
    "kaleido>=1.0.0",

    # Option B: Downgrade to Plotly 5.x for fast Kaleido
    "plotly>=5.23.0,<6.0.0",
    "kaleido==0.2.1",

    # Option C: Make Kaleido optional (recommended)
    # Move kaleido to optional-dependencies
]

[project.optional-dependencies]
static-reports = ["kaleido>=1.0.0"]  # Only install if --static-report needed
```

**Detection:**

- Import test at pipeline start: `import kaleido; import plotly; check_version_compatibility()`
- Integration test that exercises `fig.write_image()` with actual Plotly figures
- CI matrix testing Plotly 5.x + Kaleido 0.2.1 vs Plotly 6.x + Kaleido 1.x

**Phase Assignment:** Phase 1 (Static Report Toggle) — must verify compatibility before making Kaleido optional

**Source Confidence:** HIGH — [Plotly GitHub issue #5241](https://github.com/plotly/plotly.py/issues/5241), [Kaleido issue #400](https://github.com/plotly/Kaleido/issues/400), [Plotly Community discussion](https://community.plotly.com/t/kaleido1-2-0-is-much-slower-than-the-lower-version-0-2-1/95537)

---

### Pitfall 2: CDN plotly.js Breaks in Air-Gapped Environments

**What goes wrong:**

Switching from `include_plotlyjs='cdn'` to reduce HTML file size from 9.6 MB to 20 KB renders reports completely broken in offline/air-gapped environments. The HTML loads but shows blank spaces where plots should be, with browser console errors: `Failed to load resource: net::ERR_INTERNET_DISCONNECTED https://cdn.plot.ly/plotly-latest.min.js`.

PlasmiCheck is used in:
- Hospital research networks (strict firewall policies)
- HPC clusters without internet access
- Secure genomics facilities (air-gapped for HIPAA/data protection)
- Docker containers in CI environments with no external network

**Why it happens:**

Developers optimize for the common case (connected laptop) and test reports by opening them in their browser with internet access. The CDN approach works perfectly in development, so the offline failure mode is never discovered until deployed.

**Consequences:**

- Reports generated on HPC clusters cannot be viewed
- Users copy HTML to their laptop to view → CDN works → they don't report the issue
- Batch processing generates hundreds of "broken" HTML files
- Users lose trust in the tool and revert to older versions

**Prevention:**

**Multi-tier fallback strategy:**

```python
# generate_report.py — implement configurable plotly.js inclusion

from plasmicheck.config import get_config

_cfg = get_config()
PLOTLYJS_MODE = _cfg.get("plotly_mode", "directory")  # cdn | directory | embedded

def generate_plots(reads_df, output_folder):
    # ... create figures ...

    if PLOTLYJS_MODE == "cdn":
        include_plotlyjs = "cdn"
    elif PLOTLYJS_MODE == "directory":
        # Copy plotly.min.js to output folder once
        plotlyjs_path = Path(output_folder) / "plotly.min.js"
        if not plotlyjs_path.exists():
            copy_plotly_bundle(plotlyjs_path)
        include_plotlyjs = "plotly.min.js"  # Relative path
    else:  # embedded
        include_plotlyjs = True

    fig_box.write_html(box_html, include_plotlyjs=include_plotlyjs)
    fig_scatter.write_html(scatter_html, include_plotlyjs=include_plotlyjs)
```

**Config.json:**

```json
{
  "plotly_mode": "directory",  // Safe default for all environments
  "plotly_fallback_modes": ["directory", "embedded", "cdn"]  // Priority order
}
```

**CLI flag for override:**

```bash
plasmicheck pipeline --plotly-mode cdn    # User explicitly enables CDN
plasmicheck pipeline --plotly-mode offline  # Alias for 'directory'
```

**Detection:**

- Integration test that generates report and validates it in a network-restricted container
- Documentation warning in README about offline usage
- Runtime check: attempt to fetch CDN plotly.js and warn if unreachable

**Phase Assignment:** Phase 2 (HTML Optimization) — must address when implementing CDN mode

**Source Confidence:** MEDIUM — [Medium article on air-gapped Plotly](https://foongminwong.medium.com/plotting-data-with-plotly-offline-mode-in-an-air-gapped-environment-5844df874537), [Plotly offline issue #1742](https://github.com/plotly/plotly.py/issues/1742), [Plotly.js CDN discussion](https://community.plotly.com/t/plotly-js-cdn-down/76253)

---

### Pitfall 3: Parallel Alignment Subprocess File Handle Exhaustion

**What goes wrong:**

Running plasmid and human alignments concurrently with `ThreadPoolExecutor` causes sporadic failures:
```
OSError: [Errno 24] Too many open files
minimap2: failed to open file '/tmp/minimap2.12345.tmp'
samtools sort: truncated file. Aborting
```

The pipeline creates ~17 subprocesses per sample-plasmid combination. For batch processing (10 samples × 5 plasmids = 50 combinations), running 2 concurrent alignments per combination could spike to 200+ simultaneous processes.

**Why it happens:**

Each `minimap2 | samtools view | samtools sort` pipe spawns 3 processes. With plasmid + human alignments in parallel, that's 6 processes per combination. Python's `ThreadPoolExecutor` doesn't limit subprocesses — only threads. In batch mode, if 10 combinations run concurrently (via outer parallelism), you hit 60 processes × 3-4 open files each = 240+ file handles.

Linux default `ulimit -n` is often 1024, but on some HPC systems it's as low as 256. Docker containers may inherit even lower limits.

**Consequences:**

- Pipeline fails intermittently (race condition depending on process timing)
- Batch processing hangs or crashes
- Temp files (`/tmp/minimap2.*.mmi.tmp`) left behind, filling disk
- samtools fails with truncated BAM errors, corrupting output

**Prevention:**

**1. Limit concurrent subprocess-heavy operations:**

```python
# run_pipeline.py

from concurrent.futures import ThreadPoolExecutor
import os

# Don't blindly parallelize everything
def should_parallelize_alignments():
    """Only parallelize if resources allow."""
    cpu_count = get_cpu_count()  # cgroup-aware
    ulimit_nofile = get_file_handle_limit()

    # Each concurrent alignment uses ~6 processes × 4 file handles = 24
    if ulimit_nofile < 512:
        logging.warning(
            f"File handle limit ({ulimit_nofile}) too low for parallel alignment. "
            "Running sequentially."
        )
        return False

    if cpu_count < 16:
        # Not enough CPUs to benefit from parallel alignment
        return False

    return True

# In pipeline execution:
if should_parallelize_alignments():
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [
            pool.submit(align_reads, plasmid_index, seq_file, plasmid_bam, "plasmid", fq2),
            pool.submit(align_reads, spliced_index, seq_file, human_bam, "human", fq2)
        ]
        wait(futures)
        for f in futures:
            f.result()  # Raise exceptions
else:
    # Sequential fallback
    align_reads(plasmid_index, seq_file, plasmid_bam, "plasmid", fq2)
    align_reads(spliced_index, seq_file, human_bam, "human", fq2)
```

**2. Increase file handle limits in documentation:**

```bash
# README.md section on HPC usage
# Increase file descriptor limit before running batch jobs:
ulimit -n 4096

# For Docker:
docker run --ulimit nofile=4096:4096 plasmicheck:latest ...
```

**3. Use file handle context managers rigorously:**

```python
# Ensure pysam files are closed promptly
def _namesort_bam(input_bam, output_bam):
    # Close input handles before spawning samtools
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(SAMTOOLS_THREADS), "-o", output_bam, input_bam],
        check=True
    )
    # Don't hold handles open during subprocess execution
```

**Detection:**

- Integration test that runs batch pipeline (5+ combinations) in parallel
- Monitor file handle usage: `lsof | grep minimap2 | wc -l` during test runs
- CI test with artificially low `ulimit -n 256` to trigger failure mode
- Load test: 50 combinations × 2 parallel alignments

**Phase Assignment:** Phase 3 (Parallel Alignment) — critical before implementing concurrency

**Source Confidence:** MEDIUM — Based on general subprocess/threading patterns; no specific 2026 source found for minimap2 file handle issues, but this is a well-known class of problems in concurrent bioinformatics pipelines.

---

### Pitfall 4: `samtools collate` Non-Equivalence to `sort -n`

**What goes wrong:**

Replacing `samtools sort -n` with `samtools collate` to speed up read-name grouping causes subtle correctness bugs in `compare_alignments.py`. The `_iter_reads_by_name()` merger assumes reads with the same name appear in a specific order: primary alignments before supplementary/secondary. `collate` doesn't guarantee this, leading to:

- Tied reads being misassigned
- Incorrect mate bonuses (mate appears "unmapped" because primary alignment wasn't seen yet)
- Supplementary alignments processed before primary → wrong alignment scores

**Why it happens:**

The samtools documentation states: "collate ensures that reads of the same name are grouped together in contiguous groups, but doesn't make any guarantees about the order of read names between groups." Critically, within a group (same QNAME), `collate` doesn't order by flags.

`samtools sort -n` has secondary sort: "records with the same name will be ordered according to the values of the READ1 and READ2 flags. When that flag is also equal, ties are resolved with primary alignments first, then SUPPLEMENTARY, SECONDARY."

PlasmiCheck's comparison logic has an **implicit dependency** on this ordering:

```python
def _iter_reads_by_name(plasmid_bam, human_bam):
    # Assumes reads from same template are in primary→supplementary order
    for read_name in merged_names:
        plasmid_reads = list(plasmid_iter)  # May get supplementary BEFORE primary
        human_reads = list(human_iter)
        # Scoring logic assumes primary is first...
```

**Consequences:**

- Contamination ratios shift by 5-15% (not catastrophic, but wrong)
- Integration tests may still pass if they don't verify read assignment order
- Results differ between `collate` and `sort -n` runs → reproducibility failure
- Published results become non-reproducible

**Prevention:**

**Option A: Use collate correctly (recommended for performance):**

```python
# compare_alignments.py

def _iter_reads_by_name(plasmid_bam_path, human_bam_path):
    """
    Iterate over reads grouped by name.

    IMPORTANT: This function handles reads in any order within a QNAME group.
    It separates primary from supplementary/secondary and processes primaries first.
    """
    plasmid_bam = pysam.AlignmentFile(plasmid_bam_path, "rb")
    human_bam = pysam.AlignmentFile(human_bam_path, "rb")

    for qname, plasmid_group, human_group in _merge_collated_reads(plasmid_bam, human_bam):
        # Sort reads within group: primary first, then supplementary
        plasmid_primary = [r for r in plasmid_group if not r.is_supplementary and not r.is_secondary]
        plasmid_other = [r for r in plasmid_group if r.is_supplementary or r.is_secondary]

        human_primary = [r for r in human_group if not r.is_supplementary and not r.is_secondary]
        human_other = [r for r in human_group if r.is_supplementary or r.is_secondary]

        # Process in correct order
        yield qname, plasmid_primary + plasmid_other, human_primary + human_other
```

**Option B: Keep `sort -n` and accept the performance cost:**

If the speedup from `collate` is <20% and the code complexity increase is significant, keep `sort -n`. The correctness guarantee is worth the cost.

**Detection:**

- Regression test with known BAM files containing supplementary alignments
- Compare output of `sort -n` vs `collate` on same input → must be identical
- Integration test with aligned data containing:
  - Chimeric reads (supplementary alignments)
  - Multiple mappings per read
  - Secondary alignments
- Assert: contamination ratio within 0.1% tolerance

**Phase Assignment:** Phase 4 (Sort Optimization) — requires extensive testing before deployment

**Source Confidence:** HIGH — [samtools manual](http://www.htslib.org/doc/samtools-collate.html), [samtools issue #2010](https://github.com/samtools/samtools/issues/2010) on supplementary ordering, [Ubuntu manpage](https://manpages.ubuntu.com/manpages/jammy/man1/samtools-collate.1.html)

---

### Pitfall 5: Lazy Imports Breaking Type Checking and Module Constants

**What goes wrong:**

Moving `import pandas`, `import plotly.express`, `import pysam` into function bodies to reduce startup time from 1.5s to 0.5s breaks:

1. **Type checking with mypy:**
   ```python
   # generate_report.py
   def load_data(file: str) -> pd.DataFrame:  # NameError: pd not defined
       import pandas as pd  # mypy can't see this
       return pd.read_csv(file)
   ```

2. **Module-level constants that depend on imports:**
   ```python
   # compare_alignments.py (current code)
   import pysam  # Module-level

   _cfg = get_config()
   MATE_BONUS: int = _cfg["scoring"]["mate_bonus"]  # Runs at import time

   # After lazy import:
   def calculate_alignment_score(read):
       import pysam  # pysam not available when MATE_BONUS was initialized
       score = MATE_BONUS  # This works, but...
   ```

3. **Type annotations in function signatures:**
   ```python
   def align_reads(ref: str, input: str) -> pysam.AlignmentFile:  # Type not defined
       import pysam
       return pysam.AlignmentFile(...)
   ```

**Why it happens:**

Python's type checkers (mypy, pyright) perform static analysis and can't track runtime imports inside function bodies. The `if TYPE_CHECKING:` pattern requires duplicating imports, which is error-prone.

**Consequences:**

- `make typecheck` fails with hundreds of errors
- Type hints become useless → IDE autocomplete breaks
- Pre-commit hooks fail
- CI fails
- Developers disable mypy to "fix" it → type safety lost

**Prevention:**

**Use `if TYPE_CHECKING:` pattern correctly:**

```python
# generate_report.py
from __future__ import annotations
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd
    import plotly.express as px
    import pysam

def load_data(file: str) -> pd.DataFrame:
    """Load data from TSV file."""
    import pandas as pd  # Runtime import
    return pd.read_csv(file, sep="\t")

def create_plot(df: pd.DataFrame) -> Any:  # Can't use px.Figure without import
    import plotly.express as px
    return px.scatter(df, x="x", y="y")
```

**For module constants, use lazy initialization:**

```python
# compare_alignments.py
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pysam

# Don't initialize constants at module level
_mate_bonus: int | None = None
_clipping_penalty: int | None = None

def _get_scoring_params() -> dict[str, int]:
    """Lazy load scoring parameters."""
    global _mate_bonus, _clipping_penalty
    if _mate_bonus is None:
        _cfg = get_config()
        _mate_bonus = _cfg["scoring"]["mate_bonus"]
        _clipping_penalty = _cfg["scoring"]["clipping_penalty"]
    return {"mate_bonus": _mate_bonus, "clipping_penalty": _clipping_penalty}

def calculate_alignment_score(read: pysam.AlignedSegment) -> int:
    import pysam  # Runtime import
    params = _get_scoring_params()
    score = read.mapping_quality
    # ... use params["mate_bonus"] ...
    return score
```

**Detection:**

- `make typecheck` must pass after lazy import changes
- CI runs mypy in strict mode
- Pre-commit hook runs mypy

**Phase Assignment:** Phase 5 (Import Optimization) — high risk, defer to later phase

**Source Confidence:** HIGH — [PEP 690 Lazy Imports](https://peps.python.org/pep-0690/), [Type checking with lazy imports](https://github.com/scientific-python/lazy-loader/issues/28), [TYPE_CHECKING explanation](https://vickiboykis.com/2023/12/11/why-if-type_checking/)

---

## Moderate Pitfalls

Mistakes that cause delays, confusion, or technical debt.

### Pitfall 6: Performance Regression Without Correctness Verification

**What goes wrong:**

Optimizations like `samtools collate`, parallel alignment, or lazy imports improve performance but subtly break correctness. Without regression testing, these bugs ship to production:

- Contamination ratios differ by 10% between v0.30.0 and v0.31.0
- Users report: "Why do I get different results with the new version?"
- Trust in the tool erodes
- Papers citing PlasmiCheck results become questionable

Bioinformatics tools have a **reproducibility crisis**: a 2025 study found only 5.9% of Jupyter notebooks in biomedical articles could be reproduced. Performance optimizations often introduce subtle bugs that standard unit tests don't catch.

**Why it happens:**

1. **Unit tests test isolated functions**, not end-to-end correctness
2. **Integration tests use tiny synthetic data**, where bugs may not manifest
3. **Developers trust that if tests pass, correctness is preserved**
4. **No regression suite comparing outputs before/after optimization**

**Consequences:**

- Silent data corruption
- Published scientific results become non-reproducible
- Users file GitHub issues: "Results changed after update"
- Developers spend weeks debugging, find bug was introduced in "simple" optimization

**Prevention:**

**Establish a regression test suite (golden outputs):**

```python
# tests/test_regression.py

import pytest
from pathlib import Path
import pandas as pd

@pytest.mark.slow
@pytest.mark.integration
class TestRegressionSuite:
    """
    Regression tests comparing v0.31.0 outputs against v0.30.0 golden outputs.

    Golden outputs generated with v0.30.0 are stored in tests/data/golden/.
    Any deviation >0.5% in contamination ratio is considered a regression.
    """

    def test_contaminated_sample_ratio_unchanged(self, golden_outputs):
        """Contamination ratio must match v0.30.0 within tolerance."""
        from plasmicheck.scripts.run_pipeline import run_pipeline

        # Run pipeline with same inputs as v0.30.0
        result_dir = run_pipeline(...)

        # Load new output
        new_summary = pd.read_csv(result_dir / "contamination_summary.tsv", sep="\t")
        new_ratio = new_summary.loc[0, "contamination_ratio"]

        # Load golden output from v0.30.0
        golden_summary = pd.read_csv(golden_outputs / "contamination_summary.tsv", sep="\t")
        golden_ratio = golden_summary.loc[0, "contamination_ratio"]

        # Assert within tolerance
        assert abs(new_ratio - golden_ratio) < 0.005, (
            f"Contamination ratio changed: {golden_ratio:.4f} → {new_ratio:.4f}. "
            f"Regression introduced by optimization!"
        )

    def test_read_assignments_unchanged(self, golden_outputs):
        """Read-level assignments must be identical to v0.30.0."""
        # Load read assignments
        new_reads = pd.read_csv("reads_assignment.tsv", sep="\t")
        golden_reads = pd.read_csv(golden_outputs / "reads_assignment.tsv", sep="\t")

        # Compare row-by-row
        discrepancies = new_reads.compare(golden_reads)
        assert len(discrepancies) == 0, f"Read assignments differ:\n{discrepancies}"
```

**Makefile target:**

```makefile
# Makefile
test-regression:
	# Compare outputs against v0.30.0 golden files
	pytest tests/test_regression.py -v --tb=short

# Run before merging any optimization PR
test-pre-merge: test typecheck lint test-regression
```

**CI gate:**

```yaml
# .github/workflows/ci.yml
- name: Regression tests
  run: make test-regression
```

**Detection:**

- Automated: CI fails if regression tests fail
- Manual: Code review checklist includes "Did you run regression tests?"

**Phase Assignment:** Phase 0 (Pre-work) — establish regression suite BEFORE implementing optimizations

**Source Confidence:** HIGH — [Bioinformatics reproducibility study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7575271/), [Essential benchmarking guidelines](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8), [Five pillars of computational reproducibility](https://pmc.ncbi.nlm.nih.gov/articles/PMC10591307/)

---

### Pitfall 7: CPU Count Auto-Detection in Containers/HPC

**What goes wrong:**

Replacing the hardcoded `minimap2_threads: 8` with auto-detected `os.cpu_count()` causes:

1. **Docker containers over-subscribe CPUs:**
   - Container limited to 4 CPUs via cgroups
   - `os.cpu_count()` returns 64 (host machine CPU count)
   - minimap2 spawns 64 threads → thrashing, slowdown instead of speedup

2. **HPC SLURM jobs violate resource limits:**
   - Job allocated 8 cores via `#SBATCH --cpus-per-task=8`
   - `os.cpu_count()` returns 128 (entire node)
   - Job uses 128 threads → scheduler kills it for exceeding allocation

3. **Shared HPC nodes become unresponsive:**
   - Multiple users on same node
   - Each PlasmiCheck job spawns 128 threads
   - Node crashes, affecting all users

**Why it happens:**

Python's `os.cpu_count()` returns the number of CPUs visible to the OS, not the number allocated to the process. Container orchestrators (Docker, Kubernetes) and HPC schedulers (SLURM, PBS) use cgroups to limit CPU access, but `os.cpu_count()` doesn't read cgroup limits.

**Consequences:**

- Performance degradation (thrashing)
- HPC jobs killed by scheduler
- Angry HPC admins
- Tool banned from HPC clusters

**Prevention:**

**Implement cgroup-aware CPU detection:**

```python
# plasmicheck/utils.py

import os
from pathlib import Path

def get_cpu_count() -> int:
    """
    Get available CPU count, respecting container/cgroup limits.

    Checks (in order):
    1. SLURM_CPUS_PER_TASK environment variable
    2. Cgroup v2 CPU quota
    3. Cgroup v1 CPU quota
    4. os.cpu_count() fallback
    """
    # 1. SLURM job allocation
    slurm_cpus = os.getenv("SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        logging.debug(f"Using SLURM allocation: {slurm_cpus} CPUs")
        return int(slurm_cpus)

    # 2. Cgroup v2 (unified hierarchy)
    cgroup_v2_cpu = Path("/sys/fs/cgroup/cpu.max")
    if cgroup_v2_cpu.exists():
        quota, period = cgroup_v2_cpu.read_text().strip().split()
        if quota != "max":
            cpus = int(quota) // int(period)
            logging.debug(f"Using cgroup v2 limit: {cpus} CPUs")
            return max(1, cpus)

    # 3. Cgroup v1 (legacy)
    cgroup_v1_quota = Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
    cgroup_v1_period = Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
    if cgroup_v1_quota.exists() and cgroup_v1_period.exists():
        quota = int(cgroup_v1_quota.read_text().strip())
        period = int(cgroup_v1_period.read_text().strip())
        if quota > 0:
            cpus = quota // period
            logging.debug(f"Using cgroup v1 limit: {cpus} CPUs")
            return max(1, cpus)

    # 4. Fallback to system CPU count
    cpus = os.cpu_count() or 8
    logging.debug(f"Using system CPU count: {cpus}")
    return cpus
```

**Use in config:**

```python
# config.py
from plasmicheck.utils import get_cpu_count

def get_config():
    with open(config_path) as f:
        cfg = json.load(f)

    # Override thread counts with detected values if set to "auto"
    if cfg["alignment"]["minimap2_threads"] == "auto":
        cfg["alignment"]["minimap2_threads"] = get_cpu_count()
    if cfg["alignment"]["samtools_threads"] == "auto":
        cfg["alignment"]["samtools_threads"] = max(1, get_cpu_count() // 2)

    return cfg
```

**Config.json:**

```json
{
  "alignment": {
    "minimap2_threads": "auto",  // Or explicit integer
    "samtools_threads": "auto"
  }
}
```

**Detection:**

- Integration test in Docker with `--cpus=4` limit, assert thread count ≤ 4
- Integration test with `SLURM_CPUS_PER_TASK=8`, assert thread count = 8
- Documentation section on HPC usage

**Phase Assignment:** Phase 6 (Thread Auto-Detection) — defer until basic optimizations proven

**Source Confidence:** HIGH — [Python issue #36054](https://bugs.python.org/issue36054) on cgroup CPU detection, [SLURM cgroups docs](https://slurm.schedmd.com/cgroups.html), [container CPU detection](https://github.com/agile6v/container_cpu_detection)

---

### Pitfall 8: matplotlib vs Plotly Visual Inconsistency

**What goes wrong:**

Implementing matplotlib as a fallback for static PNG generation to avoid Kaleido creates visually inconsistent reports:

- Plotly default colors: `["#636EFA", "#EF553B", "#00CC96"]` (blue, red, teal)
- matplotlib default colors: `["#1f77b4", "#ff7f0e", "#2ca02c"]` (different blue, orange, green)
- Font sizes, axis labels, legend positions all differ
- Users complain: "Why do some reports look different?"

**Why it happens:**

Plotly and matplotlib have different default themes. Developers implement matplotlib version focusing on functionality (does it produce a PNG?) without checking visual consistency.

**Consequences:**

- Reports look unprofessional
- Users think there's a bug
- Scientific presentations have inconsistent figures
- Manual post-processing required to fix aesthetics

**Prevention:**

**Create a shared style configuration:**

```python
# plasmicheck/plotting.py

from plasmicheck.config import get_config

_cfg = get_config()
COLORS = _cfg["plot_sample_report"]["colors"]  # {"human": "blue", "tied": "orange", "plasmid": "red"}

# Map to specific hex codes for consistency
COLOR_MAP = {
    "human": "#636EFA",     # Plotly blue
    "tied": "#FFA500",      # Orange
    "plasmid": "#EF553B"    # Plotly red
}

def create_boxplot_plotly(df):
    import plotly.express as px
    fig = px.box(
        df, x="AssignedTo", y="PlasmidScore",
        color="AssignedTo",
        color_discrete_map=COLOR_MAP,
        template="plotly_white"  # Clean white background
    )
    return fig

def create_boxplot_matplotlib(df):
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Use same colors
    palette = [COLOR_MAP["human"], COLOR_MAP["tied"], COLOR_MAP["plasmid"]]

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=df, x="AssignedTo", y="PlasmidScore",
        hue="AssignedTo", palette=palette, ax=ax
    )
    ax.set_xlabel("Assigned To", fontsize=12)
    ax.set_ylabel("Plasmid Score", fontsize=12)
    ax.set_title("Box Plot of Plasmid Scores by Assignment", fontsize=14)

    # Match Plotly's white background
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    return fig
```

**Validation test:**

```python
# tests/test_visual_consistency.py

def test_plotly_matplotlib_color_consistency():
    """Plotly and matplotlib versions use identical colors."""
    from plasmicheck.plotting import COLOR_MAP, create_boxplot_plotly, create_boxplot_matplotlib

    # Create both plots
    df = pd.DataFrame({"AssignedTo": ["Human", "Tied", "Plasmid"], "PlasmidScore": [10, 20, 30]})
    plotly_fig = create_boxplot_plotly(df)
    mpl_fig = create_boxplot_matplotlib(df)

    # Extract colors from Plotly
    plotly_colors = plotly_fig.layout.colorway or [COLOR_MAP[k] for k in COLOR_MAP]

    # Extract colors from matplotlib
    mpl_colors = [patch.get_facecolor() for patch in mpl_fig.axes[0].patches]

    # Compare (requires hex conversion)
    # ... assert colors match ...
```

**Documentation:**

Include visual comparison in docs showing Plotly vs matplotlib output side-by-side.

**Detection:**

- Visual regression test comparing screenshots
- Manual review during PR: "Do the plots look identical?"

**Phase Assignment:** Phase 7 (matplotlib fallback) — only if Kaleido alternative needed

**Source Confidence:** MEDIUM — [Plotly color scales](https://plotly.com/python/colorscales/), [matplotlib colormaps](https://www.analyticsvidhya.com/blog/2020/09/colormaps-matplotlib/), general best practices

---

## Minor Pitfalls

Mistakes that cause annoyance but are easily fixable.

### Pitfall 9: CLI Flag Changes Breaking Snakemake Workflows

**What goes wrong:**

Adding a `--static-report` flag and making PNG generation opt-in breaks existing Snakemake workflows:

**Old behavior (v0.30.0):**
```bash
plasmicheck pipeline --output results/
# Always generates both interactive and non-interactive HTML
```

**New behavior (v0.31.0):**
```bash
plasmicheck pipeline --output results/
# Only generates interactive HTML (no PNGs)

plasmicheck pipeline --output results/ --static-report
# Generates both
```

**User's Snakefile:**
```python
rule plasmicheck:
    output:
        "results/{sample}/report_non_interactive.html"
    shell:
        "plasmicheck pipeline --output results/{wildcards.sample}"
```

After upgrading to v0.31.0, `report_non_interactive.html` is never created → Snakemake fails.

**Why it happens:**

Developers change default behavior without considering downstream automation. CLI flags are part of the API — changing defaults is a breaking change.

**Consequences:**

- Snakemake workflows break silently (output files missing)
- Users' production pipelines fail
- GitHub issues flood in: "Workflow stopped working after update"
- Users pin to v0.30.0 and never upgrade

**Prevention:**

**Maintain backward compatibility:**

```bash
# Option A: Keep old default, add --skip-static-report flag
plasmicheck pipeline --output results/  # Still generates both (v0.30.0 behavior)
plasmicheck pipeline --output results/ --skip-static-report  # New: opt-out

# Option B: Add --report-mode flag with explicit values
plasmicheck pipeline --output results/ --report-mode both        # Default
plasmicheck pipeline --output results/ --report-mode interactive  # New: skip PNGs
plasmicheck pipeline --output results/ --report-mode static       # Generate only PNGs

# Option C: Deprecation cycle
# v0.31.0: Add --skip-static-report flag, keep both as default, warn about future change
# v0.32.0: Change default to interactive-only, deprecation warning if flag not specified
# v0.33.0: Remove warning, new default is interactive-only
```

**Recommended: Option A (least disruptive)**

**Changelog communication:**

```markdown
## v0.31.0 (2026-02-20)

### Added
- `--skip-static-report` flag to skip PNG generation (saves 11s per sample)

### Performance
- Reduced HTML report size from 9.6 MB to 20 KB (99.8% reduction)
- Added parallel alignment support (1.8x speedup)

### Migration Notes
- **No breaking changes**: Default behavior unchanged (both report types generated)
- To skip PNGs for faster pipeline, add `--skip-static-report` flag
- Snakemake workflows require no changes
```

**Detection:**

- Integration test that runs CLI without flags, asserts both reports exist
- Backward compatibility test suite
- Documentation audit: search for CLI examples and update

**Phase Assignment:** Phase 1 (Static Report Toggle) — must preserve backward compatibility

**Source Confidence:** HIGH — [Snakemake migration guide](https://snakemake.readthedocs.io/en/stable/getting_started/migration.html), general software engineering best practices

---

### Pitfall 10: Benchmark Cold Start Cache Contamination

**What goes wrong:**

Benchmarking performance optimizations gives misleading results:

**Scenario 1: Cache contamination**
```python
# Benchmark: "Lazy imports save 1.2s"
import time

# Cold start
start = time.time()
from plasmicheck.scripts.generate_report import main
cold_time = time.time() - start  # 1.5s

# Warm start (re-import)
start = time.time()
from plasmicheck.scripts.generate_report import main  # Already imported!
warm_time = time.time() - start  # 0.0s

print(f"Speedup: {cold_time - warm_time}s")  # 1.5s — misleading!
```

**Scenario 2: Measuring the wrong thing**
```python
# Benchmark: "Parallel alignment is 1.8x faster"
# But benchmark only tests 200-read synthetic file where alignment takes 0.1s
# Real datasets: 10M reads, 30 min alignment → 1.8x speedup = 17 min saved (real benefit)
# Developer optimizes based on synthetic benchmark → wastes time on wrong bottlenecks
```

**Scenario 3: System state variation**
```bash
# First run: minimap2 index cached in RAM (fast)
# Second run: system under load, index evicted (slow)
# Conclusion: optimization made it slower (wrong!)
```

**Why it happens:**

- Not isolating benchmark runs (process state leaks)
- Using too-small datasets where bottlenecks differ
- Not accounting for filesystem cache
- Single-run measurements (no statistical significance)

**Consequences:**

- Wrong optimizations prioritized
- Real bottlenecks ignored
- Performance claims in README/docs are false
- Users don't see promised speedups

**Prevention:**

**Rigorous benchmarking protocol:**

```python
# tests/benchmark_pipeline.py

import subprocess
import time
import statistics
from pathlib import Path

def benchmark_pipeline(input_data, iterations=5):
    """
    Benchmark pipeline with cold-start isolation.

    Each iteration runs in a fresh subprocess to avoid import caching.
    """
    times = []

    for i in range(iterations):
        # Clear filesystem cache (requires sudo)
        # subprocess.run(["sudo", "sync"])
        # subprocess.run(["sudo", "sh", "-c", "echo 3 > /proc/sys/vm/drop_caches"])

        # Run in fresh subprocess
        start = time.time()
        result = subprocess.run(
            ["plasmicheck", "pipeline", "--output", f"/tmp/bench_{i}", ...],
            capture_output=True,
            check=True
        )
        elapsed = time.time() - start
        times.append(elapsed)

        # Cleanup
        shutil.rmtree(f"/tmp/bench_{i}")

    return {
        "mean": statistics.mean(times),
        "median": statistics.median(times),
        "stdev": statistics.stdev(times),
        "min": min(times),
        "max": max(times)
    }

# Use realistic dataset
BENCHMARK_DATA = {
    "small": "200 reads (synthetic)",       # For quick iteration
    "medium": "100K reads (real data)",     # For realistic bottlenecks
    "large": "10M reads (production scale)" # For production validation
}

def test_optimization_speedup():
    """Validate that optimization actually speeds up realistic workloads."""
    baseline = benchmark_pipeline(BENCHMARK_DATA["medium"], iterations=3)
    optimized = benchmark_pipeline(BENCHMARK_DATA["medium"], iterations=3)

    speedup = baseline["mean"] / optimized["mean"]
    assert speedup > 1.5, f"Expected 1.5x speedup, got {speedup:.2f}x"
```

**Benchmark report format:**

```markdown
## Benchmark Results: Parallel Alignment

**Dataset:** 100K paired-end reads (real Illumina data)
**Iterations:** 5
**System:** 16-core AMD EPYC, 64 GB RAM, NVMe SSD

| Metric | Baseline (v0.30.0) | Optimized (v0.31.0) | Speedup |
|--------|-------------------|---------------------|---------|
| Mean   | 1847s (30.8 min)  | 1024s (17.1 min)    | 1.80x   |
| Median | 1839s             | 1019s               | 1.80x   |
| Stdev  | 23s               | 18s                 | -       |

**Bottleneck breakdown:**
- Alignment: 1680s → 920s (1.83x — expected from parallelism)
- Report: 120s → 62s (1.94x — from skipping PNGs)
- Other: 47s → 42s (no change)

**Validation:** Contamination ratios identical to v0.30.0 (regression test passed)
```

**Detection:**

- CI runs benchmark suite on realistic data
- PR description includes benchmark results with methodology
- Code review checklist: "Benchmark uses realistic dataset?"

**Phase Assignment:** Phase 0 (Pre-work) — establish benchmarking protocol before optimizing

**Source Confidence:** HIGH — [Bioinformatics benchmarking guidelines](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8), [Performance profiling best practices](https://pmc.ncbi.nlm.nih.gov/articles/PMC7575271/)

---

## Phase-Specific Warnings

Mapping pitfalls to roadmap phases (from PERFORMANCE_ANALYSIS.md recommended order):

| Phase | Topic | Critical Pitfalls | Prevention Priority |
|-------|-------|-------------------|---------------------|
| **Phase 0: Pre-work** | Establish regression tests & benchmarks | Pitfall 6 (correctness regression), Pitfall 10 (benchmark misleading) | **MUST complete before any optimization** |
| **Phase 1: Static Report Toggle** | Make PNG export opt-in | Pitfall 1 (Kaleido compatibility), Pitfall 9 (CLI breaking changes) | Test with Plotly 5.x and 6.x |
| **Phase 2: HTML Optimization** | Use `include_plotlyjs='cdn'` or `'directory'` | Pitfall 2 (CDN offline failure) | Implement multi-tier fallback |
| **Phase 3: Parallel Alignment** | `ThreadPoolExecutor` for plasmid + human | Pitfall 3 (file handle exhaustion) | Test with low `ulimit`, monitor file handles |
| **Phase 4: Sort Optimization** | Replace `sort -n` with `collate` | Pitfall 4 (collate non-equivalence) | Extensive regression testing on real BAMs |
| **Phase 5: Import Optimization** | Lazy imports in function bodies | Pitfall 5 (type checking breaks) | Use `if TYPE_CHECKING:`, verify mypy passes |
| **Phase 6: Thread Auto-Detection** | Replace hardcoded threads with `get_cpu_count()` | Pitfall 7 (container over-subscription) | Test in Docker, SLURM, bare metal |
| **Phase 7: matplotlib Fallback** | Optional matplotlib for static plots | Pitfall 8 (visual inconsistency) | Visual regression tests |

---

## Cross-Cutting Concerns

### Testing Strategy for Performance Optimizations

All optimization PRs must include:

1. **Regression tests** (Pitfall 6):
   - Output identical to v0.30.0 within tolerance
   - Run on realistic datasets (>100K reads)

2. **Integration tests** (Pitfall 3, 4):
   - Full pipeline end-to-end
   - Real minimap2/samtools execution
   - Resource-constrained environments

3. **Type checking** (Pitfall 5):
   - `make typecheck` passes
   - No new mypy errors

4. **Backward compatibility** (Pitfall 9):
   - CLI interface unchanged OR explicitly versioned
   - Deprecated flags have warning messages

5. **Benchmark validation** (Pitfall 10):
   - Realistic dataset sizes
   - Multiple iterations
   - Statistical analysis (mean, stdev)

### Documentation Requirements

For each optimization:

1. **Changelog entry** explaining behavior changes
2. **Migration guide** if backward compatibility affected
3. **README update** with performance numbers and caveats
4. **HPC usage section** for environment-specific issues (Pitfall 7)
5. **Offline usage section** for CDN concerns (Pitfall 2)

---

## Sources

### High Confidence (Official documentation, recent studies)

1. [Kaleido v1.2.0 performance regression](https://github.com/plotly/Kaleido/issues/400) — Kaleido #400
2. [Plotly 5.x vs 6.x compatibility](https://github.com/plotly/plotly.py/issues/5241) — Plotly.py #5241
3. [Plotly offline in air-gapped environments](https://foongminwong.medium.com/plotting-data-with-plotly-offline-mode-in-an-air-gapped-environment-5844df874537)
4. [samtools collate vs sort -n](http://www.htslib.org/doc/samtools-collate.html) — Official manual
5. [samtools supplementary alignment ordering](https://github.com/samtools/samtools/issues/2010)
6. [PEP 690 Lazy Imports](https://peps.python.org/pep-0690/)
7. [TYPE_CHECKING pattern for lazy imports](https://github.com/scientific-python/lazy-loader/issues/28)
8. [Python os.cpu_count() cgroup issue](https://bugs.python.org/issue36054)
9. [SLURM cgroups documentation](https://slurm.schedmd.com/cgroups.html)
10. [Bioinformatics benchmarking guidelines](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8)
11. [Computational reproducibility five pillars](https://pmc.ncbi.nlm.nih.gov/articles/PMC10591307/)
12. [Snakemake migration guide](https://snakemake.readthedocs.io/en/stable/getting_started/migration.html)

### Medium Confidence (Community forums, blog posts)

13. [Plotly CDN down incident](https://community.plotly.com/t/plotly-js-cdn-down/76253)
14. [Kaleido 0.2.1 performance discussion](https://community.plotly.com/t/kaleido1-2-0-is-much-slower-than-the-lower-version-0-2-1/95537)
15. [matplotlib vs Plotly visual consistency](https://www.analyticsvidhya.com/blog/2020/09/colormaps-matplotlib/)
16. [Container CPU detection](https://github.com/agile6v/container_cpu_detection)

### Context from PlasmiCheck Codebase

17. PERFORMANCE_ANALYSIS.md — Profiling data showing Kaleido bottleneck
18. config.json — Current thread settings, scoring parameters
19. compare_alignments.py — Module-level constants depending on imports
20. Test suite (1370 lines) — Current testing strategy
