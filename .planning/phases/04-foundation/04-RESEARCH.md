# Phase 4: Foundation - Research

**Researched:** 2026-02-14
**Domain:** Python regression testing and benchmarking for bioinformatics pipelines
**Confidence:** HIGH

## Summary

This phase establishes local developer tooling to validate that performance optimizations in Phases 5-7 preserve correctness. Based on user decisions, this is standalone scripts (not pytest suite), local-only (no CI integration), and focused on practical comparison rather than over-engineered infrastructure.

The standard approach for this problem domain combines:
1. **Regression testing**: Compare final pipeline outputs (contamination ratios, read assignments, verdicts) using pandas DataFrame comparison with floating-point tolerance
2. **Benchmarking**: Per-step timing using `time.perf_counter()` for multi-second operations, with statistical reporting (mean/std over 3 iterations)
3. **Baseline caching**: JSON files stored locally (not in repo) with atomic writes to prevent corruption

PlasmiCheck already has synthetic test fixtures (200 reads, 448K reads available), pytest infrastructure, and well-defined output formats (.summary.tsv, .reads_assignment.tsv). The regression script runs v0.31.0 pipeline once to establish baseline, then compares future runs against cached results. The benchmark script times individual pipeline steps and outputs markdown tables for PR review.

**Primary recommendation:** Use manual `time.perf_counter()` for benchmarking (not pytest-benchmark or pyperf) because PlasmiCheck pipeline steps are multi-second subprocess chains. Use pandas `assert_frame_equal()` with `check_exact=False` for regression comparison to handle floating-point variance.

## Standard Stack

The established libraries/tools for Python regression testing and benchmarking:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pandas | 2.2.2+ | DataFrame comparison for regression testing | Already in project, industry-standard for tabular data comparison with built-in tolerance handling |
| time (stdlib) | - | Benchmarking via `perf_counter()` | Python standard library, high-precision monotonic clock, no dependencies |
| pathlib (stdlib) | - | Cache directory and file management | Python 3.10+ standard for filesystem operations, cross-platform |
| json (stdlib) | - | Baseline cache serialization | Built-in, simple key-value storage for benchmark/baseline data |
| argparse (stdlib) | - | CLI interface for standalone scripts | Standard library argument parser, sufficient for simple tools |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pytest | 8.0+ | Approx for float comparison | Already in dev deps, use `pytest.approx()` for small floating-point tolerances |
| subprocess (stdlib) | - | Run v0.31.0 pipeline for baseline | Execute external commands, capture output |
| statistics (stdlib) | - | Calculate mean/stdev for benchmark results | Built-in, sufficient for basic statistical reporting |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Manual `time.perf_counter()` | pytest-benchmark | pytest-benchmark adds overhead for test suite integration; unnecessary for standalone scripts timing multi-second operations |
| Manual timing | pyperf | pyperf isolates benchmarks in subprocesses with CPU affinity—overkill for timing entire pipeline steps that already spawn subprocesses |
| JSON cache | pytest cache | pytest cache is for test metadata; user explicitly wants standalone scripts, not pytest integration |
| pandas comparison | Line-by-line diff | pandas handles floating-point tolerance, column ordering, and provides clear error messages |

**Installation:**
```bash
# All dependencies already in project
# No new packages needed (pandas, pytest in dev deps)
```

## Architecture Patterns

### Recommended Project Structure
```
scripts/
├── regression_test.py      # Compare v0.31.0 baseline with current
├── benchmark.py            # Time pipeline steps, output markdown
└── .regression_cache/      # Local cache (gitignored)
    ├── baseline.json       # Cached v0.31.0 results
    └── last_run.json       # Timestamp metadata
```

### Pattern 1: Baseline Generation and Caching
**What:** Run v0.31.0 pipeline once, cache results locally as JSON
**When to use:** First time regression test runs or when baseline needs regeneration
**Example:**
```python
# Source: Python atomicwrites pattern + pathlib best practices
from pathlib import Path
import json
import tempfile
import os

def save_baseline(data: dict, cache_path: Path) -> None:
    """Atomically write baseline to prevent corruption."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    # Write to temp file first
    fd, temp_path = tempfile.mkstemp(
        dir=cache_path.parent,
        prefix=cache_path.name,
        suffix='.tmp'
    )

    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(data, f, indent=2)
        # Atomic rename (POSIX guarantees atomicity)
        os.replace(temp_path, cache_path)
    except Exception:
        os.unlink(temp_path)
        raise

def load_baseline(cache_path: Path) -> dict | None:
    """Load cached baseline or return None if missing."""
    if not cache_path.exists():
        return None
    with cache_path.open() as f:
        return json.load(f)
```

### Pattern 2: Per-Step Pipeline Timing
**What:** Time each pipeline step separately using context manager
**When to use:** Benchmarking individual phases (convert, index, align, compare, report)
**Example:**
```python
# Source: time.perf_counter() best practices for multi-second operations
import time
from contextlib import contextmanager
from typing import Generator

@contextmanager
def timed_step(name: str, timings: dict) -> Generator[None, None, None]:
    """Context manager to time a pipeline step."""
    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        timings[name] = elapsed

# Usage
timings = {}
with timed_step("align_plasmid", timings):
    align_reads(...)  # Multi-second subprocess

with timed_step("align_human", timings):
    align_reads(...)
```

### Pattern 3: DataFrame Regression Comparison
**What:** Compare TSV outputs using pandas with floating-point tolerance
**When to use:** Validating pipeline outputs haven't changed
**Example:**
```python
# Source: pandas.testing.assert_frame_equal official docs
import pandas as pd
from pandas.testing import assert_frame_equal

def compare_summary_files(baseline_tsv: str, current_tsv: str, tol: float = 0.001) -> None:
    """Compare summary.tsv files with tolerance for float columns."""
    baseline_df = pd.read_csv(baseline_tsv, sep='\t')
    current_df = pd.read_csv(current_tsv, sep='\t')

    # Separate numeric and categorical columns
    numeric_cols = baseline_df.select_dtypes(include='number').columns

    # Check exact equality for non-numeric (verdicts, categories)
    if not baseline_df['Category'].equals(current_df['Category']):
        raise AssertionError("Categories don't match")

    # Check float columns with tolerance
    for col in numeric_cols:
        if col in baseline_df.columns and col in current_df.columns:
            assert_frame_equal(
                baseline_df[[col]],
                current_df[[col]],
                check_exact=False,  # Allow float tolerance
                rtol=tol,  # Relative tolerance
                check_dtype=False  # Allow int vs float
            )
```

### Pattern 4: Statistical Benchmark Reporting
**What:** Run multiple iterations, report mean/std in markdown
**When to use:** Generating BENCHMARK.md for PR review
**Example:**
```python
# Source: Python statistics stdlib + markdown table formatting
import statistics
from pathlib import Path

def run_benchmark_iterations(n: int = 3) -> dict[str, list[float]]:
    """Run benchmark n times, collect per-step timings."""
    all_timings = {step: [] for step in PIPELINE_STEPS}

    for i in range(n):
        print(f"Iteration {i+1}/{n}...")
        timings = run_single_benchmark()
        for step, elapsed in timings.items():
            all_timings[step].append(elapsed)

    return all_timings

def generate_benchmark_markdown(timings: dict[str, list[float]], output: Path) -> None:
    """Generate markdown table with mean/std."""
    lines = [
        "# Benchmark Results",
        "",
        f"**Iterations:** {len(next(iter(timings.values())))}",
        "",
        "| Step | Mean (s) | Std Dev (s) | Min (s) | Max (s) |",
        "|------|----------|-------------|---------|---------|"
    ]

    for step, times in timings.items():
        mean = statistics.mean(times)
        stdev = statistics.stdev(times) if len(times) > 1 else 0.0
        min_t = min(times)
        max_t = max(times)
        lines.append(f"| {step} | {mean:.3f} | {stdev:.3f} | {min_t:.3f} | {max_t:.3f} |")

    output.write_text('\n'.join(lines))
```

### Anti-Patterns to Avoid
- **Brittle exact float comparison:** Don't use `==` for contamination ratios; use `pytest.approx()` or pandas rtol
- **Caching without atomicity:** Direct `json.dump()` to final path can corrupt cache if process crashes
- **Timing inside loops:** For multi-second operations, wrap entire step, not individual iterations
- **Over-engineering:** Don't build CI integration, persistent databases, or web dashboards—user wants simple local scripts

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Float comparison tolerance | Custom epsilon logic | `pytest.approx(rel=1e-6)` or pandas `rtol` | Handles edge cases (0, inf, nan), well-tested |
| Atomic file writes | Direct write then rename | `tempfile.NamedTemporaryFile` + `os.replace()` | POSIX atomic rename guarantees, handles cleanup on error |
| DataFrame comparison | Row-by-row iteration | `pandas.testing.assert_frame_equal()` | Column type coercion, detailed error messages, flexible tolerance |
| CLI argument parsing | Manual `sys.argv` parsing | `argparse` with subparsers | Auto-generates help, type conversion, error messages |
| Statistical calculations | Manual mean/variance | `statistics.mean()`, `statistics.stdev()` | Handles edge cases (n=1, empty lists), numerically stable |

**Key insight:** PlasmiCheck outputs are structured TSV files. Parsing them line-by-line is error-prone; pandas handles type conversion, missing values, and provides tolerance-aware comparison built-in.

## Common Pitfalls

### Pitfall 1: Floating-Point Comparison Without Tolerance
**What goes wrong:** Contamination ratios like `0.666666666667` vs `0.666666666666` fail equality check
**Why it happens:** IEEE 754 floating-point arithmetic, different rounding in calculation order
**How to avoid:** Use `pytest.approx(expected, rel=0.001)` for ratios, `pandas rtol=0.001` for DataFrames
**Warning signs:** Regression tests fail with tiny differences (1e-10) in numeric columns

### Pitfall 2: Cache Corruption on Interrupted Baseline Run
**What goes wrong:** Baseline generation crashes mid-pipeline, writes partial JSON, future runs compare against corrupted baseline
**Why it happens:** `json.dump()` writes directly to final path, not atomic
**How to avoid:** Write to temp file first, then `os.replace()` (atomic on POSIX)
**Warning signs:** JSON parse errors, missing keys in baseline dict

### Pitfall 3: Benchmark Timing Includes Import Overhead
**What goes wrong:** First pipeline run takes 12s, subsequent runs take 9s—benchmark shows 30% variance
**Why it happens:** Python lazy imports, module caching across iterations
**How to avoid:** Run one warm-up iteration before benchmarking, or time steps individually (imports already done)
**Warning signs:** First iteration significantly slower than subsequent iterations

### Pitfall 4: Comparing Verdicts as Floats
**What goes wrong:** Summary.tsv has mixed types—"Verdict" row has string, "Ratio" row has float—pandas reads entire "Count" column as object type
**Why it happens:** TSV has inconsistent data types in single column
**How to avoid:** Parse summary.tsv manually or use `pd.read_csv(dtype={'Count': str})` then convert numeric rows
**Warning signs:** `assert_frame_equal()` fails with dtype mismatch (object vs float64)

### Pitfall 5: Benchmark Runs on Busy System
**What goes wrong:** High variance in timings (std dev > 20% of mean) because CPU is shared
**Why it happens:** Other processes, background tasks competing for resources
**How to avoid:** Document expected variance, run 5+ iterations to smooth noise, warn if stdev > threshold
**Warning signs:** Benchmark shows 2x difference between min and max run times

### Pitfall 6: Baseline Cache Invalidation Strategy
**What goes wrong:** Code changes but baseline still compares against old v0.31.0—false negatives (regressions missed)
**Why it happens:** No mechanism to detect when baseline is stale
**How to avoid:** Store version tag in baseline metadata, check on load, prompt user to regenerate if mismatch
**Warning signs:** Tests pass but manual verification shows different results

## Code Examples

Verified patterns from research:

### Baseline Generation Script Structure
```python
#!/usr/bin/env python3
"""Regression test comparing v0.31.0 baseline with current implementation."""
from pathlib import Path
import subprocess
import json
import sys

CACHE_DIR = Path(__file__).parent / ".regression_cache"
BASELINE_FILE = CACHE_DIR / "baseline.json"

def run_v0310_pipeline(output_dir: Path) -> dict:
    """Run v0.31.0 pipeline and extract results."""
    # Use subprocess to run pipeline command
    cmd = [
        "plasmicheck", "pipeline",
        "--human_fasta", "tests/data/synthetic/human_ref.fasta",
        "--plasmid_files", "tests/data/synthetic/plasmid.gb",
        "--sequencing_files_r1", "tests/data/synthetic/contaminated_R1.fastq",
        "--sequencing_files_r2", "tests/data/synthetic/contaminated_R2.fastq",
        "--output_folder", str(output_dir),
        "--overwrite"
    ]

    subprocess.run(cmd, check=True, capture_output=True, text=True)

    # Extract results from output
    summary_files = list(output_dir.rglob("*.summary.tsv"))
    assignment_files = list(output_dir.rglob("*.reads_assignment.tsv"))

    # Parse into dict for caching
    return {
        "summary": summary_files[0].read_text(),
        "assignments": assignment_files[0].read_text(),
        "version": "0.31.0"
    }

if __name__ == "__main__":
    if BASELINE_FILE.exists():
        print("Baseline already exists. Delete to regenerate.")
        sys.exit(0)

    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        baseline = run_v0310_pipeline(Path(tmpdir))
        save_baseline(baseline, BASELINE_FILE)  # From Pattern 1

    print(f"Baseline saved to {BASELINE_FILE}")
```

### Benchmark Script CLI Interface
```python
#!/usr/bin/env python3
"""Benchmark PlasmiCheck pipeline performance."""
import argparse
from pathlib import Path

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark PlasmiCheck pipeline and generate markdown report"
    )
    parser.add_argument(
        "--iterations", "-n",
        type=int,
        default=3,
        help="Number of benchmark iterations (default: 3)"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=Path("BENCHMARK.md"),
        help="Output markdown file (default: BENCHMARK.md)"
    )
    parser.add_argument(
        "--dataset",
        choices=["synthetic", "real"],
        default="synthetic",
        help="Dataset to benchmark (default: synthetic 200 reads)"
    )

    args = parser.parse_args()

    # Run benchmark
    timings = run_benchmark_iterations(args.iterations)
    generate_benchmark_markdown(timings, args.output)

    print(f"Benchmark complete: {args.output}")

if __name__ == "__main__":
    main()
```

### DataFrame Comparison with Mixed Types
```python
# Source: Handling PlasmiCheck summary.tsv mixed-type columns
import pandas as pd
from pytest import approx

def compare_summary_tsv(baseline_path: Path, current_path: Path) -> None:
    """Compare summary.tsv files with proper type handling."""
    # Read both as strings to preserve all data
    baseline = pd.read_csv(baseline_path, sep='\t', dtype=str)
    current = pd.read_csv(current_path, sep='\t', dtype=str)

    # Check row-by-row with appropriate comparison
    for idx, row in baseline.iterrows():
        category = row['Category']
        baseline_val = row['Count']
        current_val = current.loc[current['Category'] == category, 'Count'].values[0]

        if category in ['Plasmid', 'Human', 'Tied']:
            # Integer counts must match exactly
            assert int(baseline_val) == int(current_val), \
                f"{category} count mismatch: {baseline_val} vs {current_val}"
        elif category == 'Ratio':
            # Float ratio with tolerance
            assert float(current_val) == approx(float(baseline_val), rel=0.001), \
                f"Ratio mismatch: {baseline_val} vs {current_val}"
        elif category == 'Verdict':
            # String must match exactly
            assert baseline_val == current_val, \
                f"Verdict mismatch: {baseline_val} vs {current_val}"
        elif category == 'CoverageOutsideINSERT':
            # Float with tolerance
            assert float(current_val) == approx(float(baseline_val), rel=0.001)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `timeit` for all timing | `time.perf_counter()` for scripts | Python 3.7+ | perf_counter has nanosecond precision, monotonic clock |
| unittest assertions | pytest.approx() | pytest 3.0+ (2016) | Cleaner syntax, relative/absolute tolerance |
| os.path | pathlib | Python 3.4+ std, 3.10+ dominant | Object-oriented, cross-platform |
| Manual temp file + rename | tempfile.NamedTemporaryFile context | stdlib | Automatic cleanup, exception-safe |
| pandas 1.x | pandas 2.x | 2023 | Faster, better type handling, nullable dtypes |

**Deprecated/outdated:**
- `pytest-benchmark` for standalone scripts — overkill for multi-second subprocess chains, adds pytest dependency to production workflow (user wants standalone)
- `time.time()` for benchmarking — affected by system clock adjustments, use `time.perf_counter()` instead
- Manual epsilon comparison for floats — use `pytest.approx()` or pandas `rtol` parameter

## Open Questions

Things that couldn't be fully resolved:

1. **Baseline version tagging**
   - What we know: User wants to compare against v0.31.0 specifically
   - What's unclear: How to handle baseline when v0.31.0 is no longer installed (checkout old tag? keep binary?)
   - Recommendation: Document that regression test requires v0.31.0 to be installed, or generate baseline once and commit as JSON (user said "not committed to repo" but might reconsider for convenience)

2. **Real dataset location**
   - What we know: User wants to benchmark 448K reads dataset
   - What's unclear: Is this in tests/data/real/ (which is gitignored)? User will need to specify path
   - Recommendation: Add `--dataset-path` argument to benchmark script, default to synthetic

3. **Statistical significance threshold**
   - What we know: Report mean/std over 3 iterations
   - What's unclear: What coefficient of variation (CV = std/mean) should trigger a warning?
   - Recommendation: Document expected CV < 0.1 (10%) for stable steps, flag if higher

## Sources

### Primary (HIGH confidence)
- [Python timeit documentation](https://docs.python.org/3/library/timeit.html) - Official stdlib docs on timing
- [time.perf_counter() vs timeit](https://superfastpython.com/time-time-vs-timeit/) - Best practices for script vs function timing
- [pandas.testing.assert_frame_equal documentation](https://pandas.pydata.org/docs/reference/api/pandas.testing.assert_frame_equal.html) - Official API for DataFrame comparison
- [pytest.approx documentation](https://docs.pytest.org/en/stable/reference/reference.html) - Float comparison with tolerance
- [Python pathlib documentation](https://docs.python.org/3/library/pathlib.html) - Official stdlib docs
- [argparse documentation](https://docs.python.org/3/library/argparse.html) - CLI argument parsing

### Secondary (MEDIUM confidence)
- [Pytest Approx comprehensive guide](https://pytest-with-eric.com/pytest-advanced/pytest-approx/) - Practical examples
- [Atomic file writes in Python](https://gist.github.com/therightstuff/cbdcbef4010c20acc70d2175a91a321f) - Pattern verification
- [pytest-benchmark vs pyperf comparison](https://thelinuxcode.com/how-to-check-the-execution-time-of-a-python-script-practical-2026-guide/) - 2026 best practices
- [Markdown table formatting for benchmarks](https://www.tomarkdown.org/guides/markdown-table) - Common patterns
- [Performance baseline testing](https://oneuptime.com/blog/post/2026-01-30-performance-baseline-testing/view) - Statistical approaches

### Tertiary (LOW confidence)
- Community discussions on pytest-regressions vs custom scripts - marked for validation based on user's "standalone scripts" decision
- Blog posts on bioinformatics pipeline benchmarking - general patterns, not PlasmiCheck-specific

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All stdlib + pandas already in project, official docs verified
- Architecture: HIGH - Patterns verified from official docs and standard practices
- Pitfalls: MEDIUM - Based on common Python testing issues, PlasmiCheck-specific ones require validation

**Research date:** 2026-02-14
**Valid until:** ~90 days (Python stdlib stable, pandas 2.x mature)
