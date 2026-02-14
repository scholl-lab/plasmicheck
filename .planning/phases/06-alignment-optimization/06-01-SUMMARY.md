---
phase: 06-alignment-optimization
plan: 01
subsystem: alignment
status: complete
tags: [threading, cpu-detection, cgroups, slurm, docker, alignment, optimization]

# Dependencies
requires:
  - 05-04  # Test coverage and CI validation
provides:
  - thread_config module with CPU detection and allocation
  - parameterized align_reads() function
  - samtools sort memory configuration
affects:
  - 06-02  # CLI wiring for thread optimization
  - 06-03  # Integration testing

# Technical
tech-stack:
  added: []
  patterns:
    - CPU detection chain (SLURM → cgroup v2 → cgroup v1 → os.cpu_count → fallback)
    - Thread allocation strategy (80/20 minimap2/samtools split)
    - Optional parameters with config defaults pattern

key-files:
  created:
    - plasmicheck/thread_config.py
    - tests/test_thread_config.py
  modified:
    - plasmicheck/scripts/align_reads.py
    - plasmicheck/config.json

# Decisions
decisions:
  - id: CPU_DETECTION_CHAIN
    decision: "Implement 5-tier detection: SLURM → cgroup v2 → cgroup v1 → os.cpu_count → fallback(4)"
    rationale: "Handles Docker (cgroups), SLURM clusters, and bare metal environments"
    alternatives: "Single detection method would fail in some environments"
    impact: "Works across all deployment scenarios"

  - id: THREAD_ALLOCATION_BOUNDS
    decision: "Min 2, max 16 total threads; 80% to minimap2, remainder to samtools (max 4)"
    rationale: "Minimap2 is CPU-bound and benefits from parallelism; samtools sort scales poorly beyond 4 threads"
    alternatives: "Equal split wastes resources; no bounds risks resource exhaustion"
    impact: "Optimal resource utilization in 2-16 CPU environments"

  - id: SAMTOOLS_SORT_MEMORY
    decision: "Add -m 2G flag to samtools sort, configurable via config.json"
    rationale: "Controls memory usage, prevents OOM in memory-constrained environments"
    alternatives: "Default memory might be too high/low for different datasets"
    impact: "Predictable memory footprint, configurable per deployment"

  - id: BACKWARD_COMPATIBILITY
    decision: "Make thread/memory parameters optional with config defaults"
    rationale: "Existing callers don't need changes; new callers can override"
    alternatives: "Required parameters would break all existing code"
    impact: "Zero-impact deployment, gradual adoption"

# Metrics
duration: 5.7min
completed: 2026-02-14
---

# Phase 06 Plan 01: Thread Detection and Allocation Summary

**One-liner:** CPU auto-detection across SLURM/Docker/bare-metal with 80/20 minimap2/samtools thread allocation and configurable samtools sort memory.

## What Was Built

Created the foundation for alignment optimization through intelligent CPU detection and thread allocation:

1. **Thread Configuration Module** (`plasmicheck/thread_config.py`):
   - `detect_cpu_count()` — 5-tier detection chain with fallbacks
   - `allocate_threads()` — bounded 80/20 split between minimap2 and samtools

2. **Parameterized Alignment** (`plasmicheck/scripts/align_reads.py`):
   - Added optional `minimap2_threads`, `samtools_threads`, `samtools_sort_memory` parameters
   - Removed module-level thread count globals
   - Added `-m` flag to all `samtools sort` commands

3. **Configuration** (`plasmicheck/config.json`):
   - Added `samtools_sort_memory: "2G"` to alignment section

4. **Comprehensive Unit Tests** (`tests/test_thread_config.py`):
   - 14 tests covering all detection scenarios and allocation bounds
   - Mock-based testing (no actual system file reads)
   - All tests marked `@pytest.mark.unit`

## Technical Implementation

### CPU Detection Chain

**Priority order:**

1. **SLURM_CPUS_PER_TASK** — HPC cluster allocation
2. **cgroup v2** (`/sys/fs/cgroup/cpu.max`) — Modern Docker
3. **cgroup v1** (`/sys/fs/cgroup/cpu/cpu.cfs_quota_us` + `cpu.cfs_period_us`) — Legacy Docker
4. **os.cpu_count()** — Bare metal Python detection
5. **Fallback to 4** — Safe minimum when all detection fails

**Edge cases handled:**

- cgroup v2 quota = "max" (unlimited) → skip
- cgroup v1 quota = -1 (unlimited) → skip
- os.cpu_count() returns None → skip
- Invalid values, IOError, OSError → silent fallthrough (logged only at final fallback)

### Thread Allocation Logic

**Bounds:**
- Min: 2 CPUs (enforced)
- Max: 16 CPUs (clamped)

**Split:**
- Minimap2: ~80% of total (int(total * 0.8)), min 2
- Samtools: remainder, min 2, max 4

**Examples:**
- 2 CPUs → minimap2=2, samtools=2 (min enforcement may exceed total)
- 8 CPUs → minimap2=6, samtools=2
- 16 CPUs → minimap2=12, samtools=4
- 32 CPUs → clamped to 16 → minimap2=12, samtools=4

### Backward Compatibility

All existing calls to `align_reads()` work unchanged:

```python
# Old code (still works)
align_reads(ref_idx, input_file, output_bam, "human")

# New code (override defaults)
align_reads(ref_idx, input_file, output_bam, "human",
            minimap2_threads=12, samtools_threads=4, samtools_sort_memory="4G")
```

Module-level globals (`MINIMAP2_THREADS`, `SAMTOOLS_THREADS`) removed — config loading now happens inside `align_reads()`.

## Testing

**Coverage:**

- **SLURM detection:** monkeypatch `os.environ`
- **cgroup v2 detection:** mock Path to fake cgroup file, test both normal and "max" quota
- **cgroup v1 detection:** mock Path to fake quota/period files, test both normal and -1 quota
- **os.cpu_count() detection:** mock Path to skip cgroups
- **Fallback detection:** mock Path + `os.cpu_count() → None`
- **Allocation bounds:** parametrized tests for inputs 1-32, verify all constraints

**Results:**

- 144 tests pass (14 new, 130 existing)
- mypy strict mode: 0 errors (19 files checked)
- ruff: 0 issues
- Zero regressions in existing test suite

## Deviations from Plan

None — plan executed exactly as written.

## Code Quality

**Type safety:**
- Full type hints with `from __future__ import annotations`
- mypy strict mode passes
- Proper handling of `os.cpu_count() -> int | None`

**Error handling:**
- Try/except around each detection step
- Silent fallthrough (detection failures are normal)
- Warning logged only on final fallback

**Testability:**
- Monkeypatch and mock-based tests
- No actual system file reads in tests
- Each detection tier independently testable

## Next Phase Readiness

**Blockers:** None

**Ready for Plan 06-02:**
- ✅ `thread_config` module available for import
- ✅ `align_reads()` accepts thread/memory parameters
- ✅ Backward compatibility verified
- ✅ All tests passing

**Integration requirements for 06-02:**

1. CLI must call `detect_cpu_count()` and `allocate_threads()` at startup
2. Pass allocated thread counts to `align_reads()` calls
3. Add `--threads` CLI flag for manual override
4. Add `--samtools-sort-memory` CLI flag

**Known dependencies:**

- Plan 06-02 (CLI wiring) directly consumes this module
- Plan 06-03 (integration tests) will verify real-world behavior

## Files Changed

**Created:**

- `plasmicheck/thread_config.py` (110 lines)
- `tests/test_thread_config.py` (260 lines)

**Modified:**

- `plasmicheck/scripts/align_reads.py` (+21, -11 lines)
  - Added 3 optional parameters
  - Removed 2 module-level globals
  - Updated 3 command strings (added `-m` flag)
  - Enhanced logging message
- `plasmicheck/config.json` (+1 line)
  - Added `samtools_sort_memory: "2G"`

**Total impact:** +381 lines, 4 files touched

## Commits

1. **337adc5** `feat(06-01): add thread detection and allocation module`
   - Created `plasmicheck/thread_config.py`
   - Created `tests/test_thread_config.py`
   - 370 lines added

2. **5d8c1f3** `feat(06-01): parameterize align_reads with thread and memory args`
   - Updated `plasmicheck/scripts/align_reads.py`
   - Updated `plasmicheck/config.json`
   - 21 additions, 11 deletions

## Decisions Made

See frontmatter `decisions` section for full details.

**Key takeaways:**

1. **5-tier detection chain** handles all deployment environments
2. **80/20 thread split** optimizes for minimap2's CPU-bound nature
3. **2G samtools sort memory** prevents OOM while maintaining performance
4. **Optional parameters** preserve backward compatibility

## Performance Impact

**Expected improvements (to be verified in 06-03):**

- **SLURM clusters:** Respects job allocation (no oversubscription)
- **Docker:** Uses actual CPU quota (not host CPU count)
- **Bare metal:** Auto-scales to available cores (up to 16)
- **Memory:** Predictable samtools sort footprint

**Not yet integrated:** This plan provides the building blocks only. Performance gains realized after CLI wiring (06-02) and pipeline integration.

## Conclusion

Successfully created the thread detection and allocation foundation for alignment optimization. All detection scenarios tested and working. align_reads() now accepts thread/memory parameters while maintaining full backward compatibility. Zero regressions, all quality checks passing.

**Status:** ✅ Complete — Ready for 06-02 (CLI wiring)
