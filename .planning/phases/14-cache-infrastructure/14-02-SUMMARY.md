---
phase: 14-cache-infrastructure
plan: 02
subsystem: infra
tags: [cache, pipeline, performance, minimap2, spliced-alignment]

# Dependency graph
requires:
  - phase: 14-01
    provides: CacheManager class with index and spliced caching
provides:
  - Pipeline integration with --cache-dir and --no-cache CLI flags
  - Index caching for plasmid FASTAs with within-run deduplication
  - Spliced alignment caching (BAM, FASTA, cDNA positions)
  - Graceful degradation on cache errors
affects: [14-03-cli-cache-management, pipeline-performance-testing]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Cache-first pattern: check cache before expensive operations"
    - "Graceful degradation: cache errors log warnings, never fail pipeline"
    - "Within-run deduplication: first combo builds and caches, subsequent combos hit cache"

key-files:
  created: []
  modified:
    - plasmicheck/cli.py
    - plasmicheck/scripts/run_pipeline.py
    - tests/test_run_pipeline.py

key-decisions:
  - "Cache checks wrap in try/except (OSError, PermissionError) for resilience"
  - "Within-run dedup uses cache mechanism: same plasmid builds once, hits on subsequent combos"
  - "Spliced FASTA indexing NOT cached separately (fast operation, not worth complexity)"
  - "Cache initialization happens after plan build but before dry-run check"

patterns-established:
  - "Cache integration pattern: check cache (get), on miss build + store (put), on hit restore (copy2)"
  - "All cache operations resilient: errors logged as warnings, pipeline continues"

# Metrics
duration: 5min
completed: 2026-02-16
---

# Phase 14 Plan 02: Pipeline Cache Integration Summary

**Pipeline checks cache before indexing and spliced alignment, with --cache-dir/--no-cache CLI flags and graceful degradation**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-16T06:40:26Z
- **Completed:** 2026-02-16T06:45:21Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- CLI flags (--cache-dir, --no-cache) added to pipeline subcommand
- CacheManager initialized in run_pipeline() with graceful degradation on errors
- Index caching integrated for plasmid FASTAs with within-run deduplication
- Spliced alignment caching integrated (BAM, FASTA, cDNA positions)
- Cache hits restore both .mmi and .fai files via shutil.copy2
- All cache operations wrapped in try/except for resilience
- Five cache integration tests added (signature, no-cache bypass, cache-dir override, init failure, within-run dedup)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add CLI flags and wire CacheManager into pipeline** - `1c31d33` (feat)
2. **Task 2: Pipeline cache integration tests** - `36c1576` (test)

**Plan metadata:** Will be committed after SUMMARY.md creation

## Files Created/Modified
- `plasmicheck/cli.py` - Added --cache-dir and --no-cache flags to pipeline subcommand
- `plasmicheck/scripts/run_pipeline.py` - Integrated CacheManager with index and spliced caching, graceful degradation
- `tests/test_run_pipeline.py` - Added 5 cache integration tests (signature, bypass, override, failure, dedup)
- `tests/test_cache_manager.py` - Fixed lint issues from 14-01 (unused vars, strict zip)

## Decisions Made

1. **Cache initialization location:** After plan build but before dry-run check
   - Rationale: Cache path needs to be determined, but dry-run shouldn't initialize cache

2. **Within-run deduplication via cache:** Same plasmid builds once, subsequent combos hit cache
   - Rationale: Replaces and improves on previous `built_indexes` within-run dedup for plasmid indexes
   - Implementation: First combo checks cache (MISS), builds, stores. Second combo checks cache (HIT), copies, skips build.

3. **Spliced FASTA indexing NOT cached:** Always build after spliced alignment cache hit
   - Rationale: Fast operation (~1s), not worth caching complexity and storage

4. **Error handling:** All cache operations wrapped in try/except (OSError, PermissionError)
   - Rationale: Cache errors should never fail pipeline, only degrade performance

5. **Cache hit restoration:** Use shutil.copy2 instead of symlinks
   - Rationale: Ensures files exist at expected paths, no symlink handling complexity

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed lint errors in test_cache_manager.py from 14-01**
- **Found during:** Task 1 (make lint)
- **Issue:** Two tests had unused `cache` variable, one test used zip without strict=True
- **Fix:** Renamed unused vars to `_cache`, added strict=True to zip
- **Files modified:** tests/test_cache_manager.py
- **Verification:** make lint passes
- **Committed in:** 1c31d33 (Task 1 commit)

**2. [Rule 2 - Missing Critical] Moved cDNA extraction into cache miss path**
- **Found during:** Task 1 (spliced alignment caching implementation)
- **Issue:** Plan showed cDNA extraction as separate step, but it must be in cache miss path to be cached
- **Fix:** Moved extract_plasmid_cdna_positions() call into spliced cache miss block, adjusted prog.advance() calls
- **Files modified:** plasmicheck/scripts/run_pipeline.py
- **Verification:** make lint && make typecheck pass, pipeline flow correct
- **Committed in:** 1c31d33 (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (1 bug, 1 missing critical)
**Impact on plan:** Both auto-fixes necessary for correctness. Lint fix prevents CI failure. cDNA move ensures proper caching.

## Issues Encountered
None

## Next Phase Readiness
- Cache integration complete, ready for CLI cache management commands (14-03)
- All 379 tests passing (30 cache manager + 5 cache integration + 344 existing)
- Cache hit/miss logging at INFO level provides visibility
- Within-run dedup verified via test_cache_within_run_dedup

**Blockers:** None

**Concerns:** None

---
*Phase: 14-cache-infrastructure*
*Completed: 2026-02-16*
