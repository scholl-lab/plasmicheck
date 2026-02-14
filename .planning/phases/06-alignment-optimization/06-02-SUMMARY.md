---
phase: 06-alignment-optimization
plan: 02
subsystem: performance
tags: [cli, threading, alignment, pipeline, optimization]

# Dependency graph
requires:
  - phase: 06-01
    provides: thread_config module with CPU detection and allocation functions
provides:
  - --threads CLI flag on pipeline subcommand with auto-detection fallback
  - Thread detection and allocation integrated into pipeline orchestration
  - Transparent logging of thread count, source, and allocation
  - Parameterized thread/memory values passed to align_reads()
affects: [documentation, deployment]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - CLI flag with None default enables optional override pattern
    - Thread detection at pipeline start with source logging
    - Sequential alignment with full thread allocation per step

key-files:
  created:
    - tests/test_thread_cli.py
  modified:
    - plasmicheck/cli.py
    - plasmicheck/scripts/run_pipeline.py

key-decisions:
  - "Pass threads=None by default to enable auto-detection"
  - "Log thread source (CLI vs detection method) for transparency"
  - "Maintain sequential alignment (no concurrent execution)"
  - "Load samtools_sort_memory from config.json (not hardcoded)"

patterns-established:
  - "Optional override pattern: CLI flag (default=None) → auto-detect if None → log source"
  - "Thread parameter propagation: CLI → run_pipeline() → allocate_threads() → align_reads()"
  - "Integration testing with comprehensive mocking to avoid external tool dependencies"

# Metrics
duration: 4min
completed: 2026-02-14
---

# Phase 06 Plan 02: CLI Thread Integration Summary

**Users can now control alignment thread allocation via --threads flag or rely on automatic CPU detection across SLURM, cgroup, and bare metal environments**

## Performance

- **Duration:** 4 minutes
- **Started:** 2026-02-14T12:46:27Z
- **Completed:** 2026-02-14T12:50:43Z
- **Tasks:** 2/2 completed
- **Files modified:** 2 (+ 1 created)

## Accomplishments

- Wired --threads CLI flag through pipeline subcommand with auto-detection fallback
- Integrated thread detection, allocation, and logging into run_pipeline() orchestration
- Passed minimap2_threads, samtools_threads, and samtools_sort_memory to align_reads()
- Added 5 integration tests covering CLI parsing, thread allocation, and logging (149 total tests)
- Maintained sequential alignment execution (no concurrent processing) to avoid memory pressure

## Task Commits

Each task was committed atomically:

1. **Task 1: Add --threads CLI flag and wire through pipeline** - `9cac552` (feat)
2. **Task 2: Add integration tests for thread CLI and pipeline wiring** - `0204964` (test)

## Files Created/Modified

**Created:**
- `tests/test_thread_cli.py` - Integration tests for --threads flag parsing, thread allocation propagation, auto-detection fallback, and logging verification

**Modified:**
- `plasmicheck/cli.py` - Added --threads argument to parser_pipeline, passed threads=args.threads to run_pipeline()
- `plasmicheck/scripts/run_pipeline.py` - Added threads parameter, imported thread_config and config, added thread detection/allocation/logging block, passed thread/memory params to align_reads() calls, updated __main__ block

## Decisions Made

**1. Thread source logging format**
- CLI override logs as "CLI --threads=N"
- Auto-detection logs as "[detection method]" (e.g., "SLURM_CPUS_PER_TASK", "cgroup v2", "os.cpu_count()")
- Rationale: Transparent debugging for users to understand where thread count came from

**2. Sequential alignment execution**
- Kept alignments sequential (plasmid, then human) with full thread allocation to each
- Did NOT add ThreadPoolExecutor or concurrent execution
- Rationale: User explicitly decided to avoid memory pressure from concurrent minimap2 processes; speedup comes from optimal thread allocation, not concurrency

**3. samtools_sort_memory from config**
- Load from `config.json` alignment.samtools_sort_memory (default "2G")
- Not hardcoded in run_pipeline.py
- Rationale: Keeps all tunable parameters in config.json for centralized control

**4. Optional --threads flag (default=None)**
- None triggers auto-detection, explicit value overrides
- Consistent with optional parameter pattern across codebase
- Rationale: Most users benefit from auto-detection; HPC users can override when needed

## Deviations from Plan

None - plan executed exactly as written.

## Verification

All success criteria met:
- ✓ --threads flag available on pipeline subcommand
- ✓ Pipeline logs thread count and detection source at startup
- ✓ Pipeline logs minimap2/samtools thread split
- ✓ align_reads() receives explicit thread counts from pipeline
- ✓ Alignments remain sequential (no concurrent execution)
- ✓ samtools sort memory loaded from config.json
- ✓ 149 tests passing (5 new)
- ✓ mypy strict type checking clean
- ✓ ruff linting clean

## Test Coverage

New integration tests in `test_thread_cli.py`:
1. `test_pipeline_help_shows_threads_flag` - Verifies --threads in help output
2. `test_pipeline_threads_override` - Verifies threads=8 override passes correct allocation to align_reads
3. `test_pipeline_threads_autodetect` - Verifies threads=None calls detect_cpu_count() and uses result
4. `test_pipeline_thread_logging` - Verifies CLI override logging format
5. `test_pipeline_thread_logging_autodetect` - Verifies auto-detection source logging

All tests use mocks to avoid requiring minimap2/samtools or real data files.

## Next Phase Readiness

**Ready for production testing:**
- Thread detection works across SLURM, Docker (cgroup v1/v2), and bare metal
- CLI integration complete with transparent logging
- Sequential alignment with optimized thread allocation
- All unit tests passing

**Future considerations:**
- Monitor actual performance improvement in production (alignment step time)
- Consider exposing minimap2/samtools thread split ratio in config.json if users need tuning
- Document --threads flag in user-facing documentation

## Integration Points

**CLI → Pipeline:**
- `plasmicheck/cli.py` parser_pipeline defines --threads argument (default=None)
- `args.threads` passed to `run_pipeline(threads=...)`

**Pipeline → Thread Config:**
- `run_pipeline()` imports `detect_cpu_count` and `allocate_threads` from thread_config
- Calls `detect_cpu_count()` if threads=None, else uses CLI value
- Calls `allocate_threads(total_threads)` to get mm2_threads, sam_threads

**Pipeline → Config:**
- `get_config().get("alignment", {}).get("samtools_sort_memory", "2G")` loads sort memory

**Pipeline → Align:**
- `align_reads()` called with `minimap2_threads=mm2_threads`, `samtools_threads=sam_threads`, `samtools_sort_memory=sort_memory`
- Alignment commands use `-t {minimap2_threads}` for minimap2, `--threads {samtools_threads}` for samtools view/sort, `-m {sort_memory}` for sort

## Implementation Notes

**Thread detection flow:**
1. User runs `plasmicheck pipeline --threads 8` OR `plasmicheck pipeline` (no --threads)
2. cli.py parses args.threads (8 or None)
3. run_pipeline() receives threads parameter
4. If threads is not None → use CLI value, set source = "CLI --threads=N"
5. If threads is None → call detect_cpu_count(), get (count, source) tuple
6. Call allocate_threads(total_threads) → get minimap2_threads, samtools_threads (80/20 split)
7. Load samtools_sort_memory from config.json
8. Log all thread parameters with source
9. Pass thread/memory params to both align_reads() calls (plasmid, human)

**Logging output example:**
```
INFO - Using 8 threads (CLI --threads=8)
INFO -   minimap2: 6 threads, samtools: 2 threads
INFO -   samtools sort memory: 2G per thread
```

**Error handling:**
- Invalid thread count (e.g., 0, negative) → argparse type=int prevents at CLI level
- detect_cpu_count() always returns valid value (fallback chain ends at 4)
- allocate_threads() enforces bounds (2-16 total, min 2 mm2, max 4 sam)
