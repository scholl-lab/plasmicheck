---
phase: 07-comparison-cleanup
plan: 02
subsystem: pipeline
tags: [pipeline, indexing, batch-processing, error-handling, resilience]

# Dependency graph
requires:
  - phase: 06-alignment-optimization
    provides: "Thread configuration and samtools optimization for alignment"
provides:
  - "Upfront human reference indexing (shared across combinations)"
  - "PipelinePlan.built_indexes tracking to avoid redundant index creation"
  - "Batch resilience with continue-on-failure processing"
  - "Per-combination timing and error logging"
  - "CombinationResult dataclass for success/failure tracking"
affects: [pipeline, batch-processing, error-handling]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Upfront resource initialization pattern (human index built once before combinations)"
    - "Batch resilience with try-except per iteration, all-fail detection"
    - "Result tracking with dataclass (CombinationResult)"

key-files:
  created: []
  modified:
    - plasmicheck/scripts/run_pipeline.py
    - tests/test_run_pipeline.py

key-decisions:
  - "Hoist human index creation to upfront phase (before combination loop) to eliminate redundant filelock overhead"
  - "Track built indexes in PipelinePlan.built_indexes set for defensive skip checks"
  - "No deduplication for plasmid indexes (per-combination, cheap to rebuild)"
  - "Continue processing remaining combinations after one fails (batch resilience)"
  - "Raise RuntimeError only if ALL combinations fail"
  - "Log per-combination timing at INFO level (format: 'label: NN.Ns')"

patterns-established:
  - "Upfront resource pattern: Create shared resources once before main loop, track in plan state"
  - "Batch resilience pattern: Wrap each iteration in try-except, track results, raise only if all fail"
  - "Result tracking pattern: Use dataclass with success/duration/error fields for structured logging"

# Metrics
duration: 8min
completed: 2026-02-14
---

# Phase 07 Plan 02: Index Dedup & Batch Resilience Summary

**Human reference indexing hoisted to upfront phase with built_indexes tracking, plus batch-resilient pipeline that continues on per-combination failures**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-14T16:22:38Z
- **Completed:** 2026-02-14T16:30:45Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Human reference index now created once upfront (before combination loop), eliminating redundant filelock overhead in batch runs with 10+ combinations
- PipelinePlan tracks built indexes in `built_indexes: set[str]` field to skip redundant creation attempts
- Pipeline continues processing remaining combinations when one fails, logging errors and timing for each
- RuntimeError raised only if ALL combinations fail, preventing one bad sample from halting entire batch
- Per-combination timing logged at INFO level for performance tracking

## Task Commits

Each task was committed atomically:

1. **Task 1: Add index tracking to PipelinePlan and upfront indexing phase** - `20727bc` (feat)
   - Added `built_indexes: set[str]` field to PipelinePlan dataclass
   - Added `CombinationResult` dataclass for tracking success/failure/timing
   - Hoisted human reference indexing out of combination loop into upfront phase
   - Wrapped each combination in try-except for batch resilience
   - Added per-combination timing and batch summary logging

2. **Task 2: Add unit tests for index deduplication and batch resilience** - `ac10fe2` (test)
   - 8 new unit tests covering all features
   - Verified human index built once, plasmid indexes per-combination
   - Verified continue-on-failure and all-fail-raises behavior
   - Verified timing logging format

**Formatting:** `8c861c2` (style: auto-format test file)

## Files Created/Modified
- `plasmicheck/scripts/run_pipeline.py` - Added upfront indexing phase, built_indexes tracking, batch resilience with CombinationResult tracking
- `tests/test_run_pipeline.py` - Added 8 new unit tests (163 total passing)

## Decisions Made
- **Upfront human indexing:** Moved human reference index creation to before combination loop to eliminate filelock overhead in batch runs. While filelock prevents redundant work, it still incurs lock acquisition cost and log noise. Upfront indexing is cleaner and faster.
- **No plasmid index dedup:** Plasmid indexes are cheap to rebuild and per-combination (each has unique output path), so no deduplication. This keeps code simple and avoids stale index issues.
- **Batch resilience:** Wrap each combination in try-except to continue processing after failures. Log each failure at ERROR level with timing. Only raise RuntimeError if ALL combinations fail, preventing one bad sample from halting entire batch.
- **Per-combination timing:** Log timing at INFO level (format: "label: NN.Ns") for performance monitoring and batch analysis.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - implementation proceeded smoothly. All 8 new tests passed on first run.

## Next Phase Readiness

- Index deduplication complete (ARCH-01, ARCH-02 from roadmap)
- Batch resilience implemented
- Ready for next comparison optimization: samtools collate (already completed in 07-01)
- Remaining Phase 7 items: matplotlib backend configuration (07-03)

---
*Phase: 07-comparison-cleanup*
*Completed: 2026-02-14*
