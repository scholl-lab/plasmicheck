---
phase: 14-cache-infrastructure
plan: 01
subsystem: infra
tags: [caching, md5, content-addressing, lru-eviction, json-manifest, atomic-writes]

# Dependency graph
requires:
  - phase: 13-output-restructuring
    provides: pathlib.Path usage patterns established in OutputPaths
provides:
  - CacheManager class with content-addressed storage for pipeline artifacts
  - Index cache storing both .mmi and .fai files together
  - Spliced alignment cache storing BAM, FASTA, and cDNA positions
  - LRU eviction enforcing size limits
  - JSON manifest with atomic writes preventing corruption
affects: [14-02, 14-03, pipeline-integration]

# Tech tracking
tech-stack:
  added: []  # No new dependencies - stdlib only
  patterns:
    - Content-addressed caching using MD5 hashing
    - Atomic file writes with tempfile + os.replace
    - LRU eviction with collections.OrderedDict
    - JSON manifest for cache metadata

key-files:
  created:
    - plasmicheck/cache_manager.py
    - tests/test_cache_manager.py
  modified: []

key-decisions:
  - "Use MD5 for content addressing (already used in codebase, sufficient for non-adversarial data)"
  - "Store both .mmi and .fai files together in index cache (prevents partial index state)"
  - "Atomic writes via tempfile in same directory + os.replace (prevents corruption)"
  - "LRU eviction with OrderedDict move_to_end() (O(1) operations)"
  - "Lazy manifest creation (only created when first entry added)"
  - "Graceful manifest corruption handling (log warning, recreate empty)"

patterns-established:
  - "Content-addressed cache key: MD5 of input file for indexes"
  - "Combined cache key: MD5 of 'plasmid_md5:human_md5' for spliced artifacts"
  - "Manifest schema: {version, entries} with per-entry metadata (paths, sizes, timestamps, type)"
  - "Cache directory structure: cache_root/{indexes,spliced}/"
  - "Eviction logging: WARNING level for transparency"

# Metrics
duration: 5min
completed: 2026-02-16
---

# Phase 14 Plan 01: Cache Infrastructure Core Summary

**Content-addressed cache engine with MD5 lookup, atomic writes, LRU eviction, and comprehensive test coverage (30 unit tests)**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-16T06:31:34Z
- **Completed:** 2026-02-16T06:36:51Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- CacheManager class with full index and spliced artifact caching
- Index cache stores BOTH .mmi and .fai files together (critical for preventing partial index state)
- Content-addressed lookup using MD5 hashing (consistent with existing codebase)
- Atomic writes prevent cache corruption (tempfile + os.replace pattern)
- LRU eviction enforces size limits (OrderedDict with O(1) operations)
- Graceful manifest corruption handling (log warning, recreate empty)
- 30 comprehensive unit tests (100% pass rate)
- 373 total tests passing (up from 353 at Phase 13)

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement CacheManager class** - `c99943f` (feat)
2. **Task 2: Unit tests for CacheManager** - `d85e401` (test)

## Files Created/Modified
- `plasmicheck/cache_manager.py` - CacheManager class with get/put for indexes and spliced artifacts, atomic writes, LRU eviction, info/clean methods
- `tests/test_cache_manager.py` - 30 unit tests covering initialization, get/put round-trips, LRU eviction, manifest persistence, corruption recovery, info/clean

## Decisions Made

1. **Store both .mmi and .fai together in index cache**
   - Rationale: create_indexes() builds both files, and if only .mmi exists without .fai, it rebuilds BOTH from scratch (partial index state). Cache must store/restore both atomically.

2. **Use MD5 for content addressing**
   - Rationale: Already used throughout PlasmiCheck (calculate_md5 in utils.py), sufficient for non-adversarial data, faster than SHA256 for multi-GB references.

3. **Lazy manifest creation**
   - Rationale: Manifest only created when first entry added (not on initialization). Simpler, avoids creating empty files.

4. **Graceful manifest corruption handling**
   - Rationale: Log warning and recreate empty manifest rather than failing. Cache is optimization, not critical path. Better to lose cache history than break pipeline.

5. **LRU eviction with OrderedDict**
   - Rationale: O(1) access and move_to_end() operations. Simpler than custom linked list. Stdlib-based, well-tested.

6. **Combined cache key for spliced artifacts**
   - Rationale: Spliced alignment depends on BOTH plasmid and human references. Use `MD5(plasmid_md5:human_md5)` to ensure unique key for each combination. Delimiter prevents collision if boundaries shift.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - implementation followed research patterns from 14-RESEARCH.md. All stdlib patterns worked as expected (tempfile, os.replace, OrderedDict).

## Next Phase Readiness

**Ready for Plan 02 (Index Cache Integration):**
- CacheManager class complete with get_index/put_index methods
- Returns dict with "mmi" and "fai" keys (caller knows where both files are)
- Atomic writes ensure no partial cache entries
- LRU eviction prevents unbounded growth

**Ready for Plan 03 (CLI and Spliced Integration):**
- info() and clean_all() methods ready for CLI wiring
- get_spliced/put_spliced ready for spliced_alignment.py integration

**No blockers:**
- All tests passing
- mypy strict clean
- ruff linting clean
- No external dependencies added

---
*Phase: 14-cache-infrastructure*
*Completed: 2026-02-16*
