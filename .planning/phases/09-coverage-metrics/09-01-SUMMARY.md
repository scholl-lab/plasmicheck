---
phase: 09-coverage-metrics
plan: 01
subsystem: metrics
tags: [pysam, numpy, coverage-analysis, bioinformatics]

# Dependency graph
requires:
  - phase: 08-insert-region-aware-filtering
    provides: Insert region computation (cDNA_positions.txt) and read_overlaps_insert() function
provides:
  - Per-region coverage metrics (depth, breadth, uniformity) computed via pysam.count_coverage()
  - coverage_metrics.py module with compute_region_coverage_metrics() function
  - 10 new rows in summary.tsv: MeanDepthInsert, MedianDepthInsert, BreadthInsert, BreadthInsert_5x, CoverageCV_Insert, MeanDepthBackbone, MedianDepthBackbone, BreadthBackbone, BreadthBackbone_5x, CoverageCV_Backbone
  - breadth_thresholds configuration in config.json
affects: [10-resistance-gene-detection, 11-summary-report-integration, report-generation]

# Tech tracking
tech-stack:
  added: [numpy.typing for type annotations]
  patterns: [Per-base depth extraction via pysam.count_coverage(), Region-specific statistics computation, Graceful fallback for missing insert region]

key-files:
  created:
    - plasmicheck/scripts/coverage_metrics.py
    - tests/test_coverage_metrics.py
  modified:
    - plasmicheck/scripts/compare_alignments.py
    - plasmicheck/config.json
    - tests/conftest.py

key-decisions:
  - "Use pysam.count_coverage() with quality_threshold=0 to match existing alignment behavior"
  - "Compute breadth at 1x always, plus configurable thresholds (default: [5])"
  - "Handle single-element arrays in CV calculation (return 0.0, not NaN)"
  - "Fallback to whole-plasmid metrics when insert_region is None, write CoverageFallback row to summary.tsv"

patterns-established:
  - "Coverage metrics as separate module: coverage_metrics.py contains all coverage computation logic, imported by compare_alignments.py"
  - "Breadth thresholds configurable via config.json: breadth_thresholds array allows multiple thresholds (e.g., [1, 5, 10])"
  - "CamelCase row names in summary.tsv for backward compatibility: MeanDepthInsert matches existing CoverageOutsideINSERT convention"

# Metrics
duration: 8min
completed: 2026-02-15
---

# Phase 9 Plan 01: Coverage Metrics Summary

**Per-region coverage metrics (mean/median depth, breadth at 1x/5x, CV) computed via pysam and written to summary.tsv with backward compatibility**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-15T06:50:51Z
- **Completed:** 2026-02-15T06:58:46Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Created coverage_metrics.py module with per-base depth extraction and statistics computation
- Integrated coverage metrics into compare_alignments.py pipeline
- Added 10 new coverage rows to summary.tsv (insert and backbone regions)
- Implemented graceful fallback when insert region unavailable (whole-plasmid metrics)
- All 216 tests passing (up from 195) with 21 new coverage metric tests

## Task Commits

Each task was committed atomically:

1. **Task 1: Coverage metrics module + config + pipeline wiring** - `54bf60e` (feat)
   - Created coverage_metrics.py with get_depth_array(), compute_coverage_stats(), compute_region_coverage_metrics()
   - Added breadth_thresholds to config.json
   - Updated compare_alignments.py to call coverage computation and write new rows to summary.tsv
   - Handles edge cases: empty arrays, zero depth, off-by-one boundaries

2. **Task 2: Unit tests for coverage metrics** - `f643048` (test)
   - 21 new unit tests covering all coverage metric functions
   - Tests for get_depth_array: count_coverage calls, ACGT summing, empty regions
   - Tests for compute_coverage_stats: mean/median, breadth at 1x/5x, CV, edge cases (all-zero, empty, single-element)
   - Tests for compute_region_coverage_metrics: insert/backbone splitting, fallback mode, boundary handling
   - Tests for summary.tsv integration: row writing, existing rows unchanged, fallback indicator
   - Added sample_summary_df_with_coverage fixture to conftest.py

## Files Created/Modified
- `plasmicheck/scripts/coverage_metrics.py` - Coverage computation engine with per-region metrics (get_depth_array, compute_coverage_stats, compute_region_coverage_metrics)
- `plasmicheck/scripts/compare_alignments.py` - Calls coverage_metrics and writes 10 new rows to summary.tsv
- `plasmicheck/config.json` - Added breadth_thresholds: [5] configuration
- `tests/test_coverage_metrics.py` - 21 comprehensive unit tests for coverage metrics
- `tests/conftest.py` - Added sample_summary_df_with_coverage fixture for Plan 02 report tests

## Decisions Made

**1. pysam.count_coverage() with quality_threshold=0**
- Rationale: Matches existing alignment behavior (samtools view -F 4 filters unmapped reads but doesn't filter by base quality). Using default quality_threshold=15 would exclude bases that were included in alignments, causing coverage to appear artificially low.

**2. Always compute breadth_1x, plus configurable thresholds**
- Rationale: Breadth at 1x (fraction of bases with ≥1 read) is fundamental metric always needed. Additional thresholds (default: [5]) are configurable via breadth_thresholds in config.json. This supports multiple thresholds (e.g., [1, 5, 10]) without hardcoding.

**3. CV=0.0 for single-element arrays**
- Rationale: numpy.std(array, ddof=1) returns NaN for single-element arrays (degrees of freedom = 0). Explicitly check len(depth_array) == 1 and return cv=0.0 to avoid NaN propagation in summary.tsv.

**4. Fallback mode for missing insert region**
- Rationale: When cDNA_positions.txt is missing (spliced alignment failed), compute metrics for whole plasmid as "insert" and return empty backbone metrics (all 0.0). Write CoverageFallback row to summary.tsv so downstream report generation can detect and warn about fallback mode.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

**1. Single-element array CV calculation returned NaN**
- **Problem:** numpy.std([42], ddof=1) returns NaN (degrees of freedom = 0), causing test failure
- **Resolution:** Added explicit check `len(depth_array) == 1` to return cv=0.0 before calling np.std()
- **Verification:** test_compute_coverage_stats_single_element now passes

**2. Linting errors for unused mock variables**
- **Problem:** Helper function _mock_alignment_file() returns (mock_ctx, mock_bam) but some tests only used mock_ctx
- **Resolution:** Prefixed unused variables with underscore: `mock_ctx, _mock_bam = ...`
- **Verification:** make lint passes without warnings

## Next Phase Readiness

- Coverage metrics infrastructure complete and tested
- summary.tsv contains all 10 new coverage rows with proper formatting (2 decimal places)
- Fallback mode documented via CoverageFallback row
- Ready for Plan 02 (report display) - conftest.py has sample_summary_df_with_coverage fixture
- Ready for Phase 10 (resistance gene detection) - coverage metrics provide foundation for gene-level coverage analysis
- Ready for Phase 11 (summary report integration) - coverage metrics can be aggregated across samples

**Blockers:** None

**Concerns:** None - all success criteria (COV-01 through COV-05) met:
- COV-01: Per-region mean and median depth ✓
- COV-02: Breadth at 1x ✓
- COV-03: Breadth at configurable thresholds (5x) ✓
- COV-04: Coverage uniformity (CV) for insert and backbone ✓
- COV-05: All metrics written to summary.tsv with backward compatibility ✓

---
*Phase: 09-coverage-metrics*
*Completed: 2026-02-15*
