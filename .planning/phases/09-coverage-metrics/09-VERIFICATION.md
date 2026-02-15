---
phase: 09-coverage-metrics
verified: 2026-02-15T12:30:00Z
status: passed
score: 11/11 must-haves verified
---

# Phase 9: Coverage Metrics Verification Report

**Phase Goal:** Every sample-plasmid pair has quantitative depth, breadth, and uniformity metrics computed per-region (insert vs backbone) and written to summary.tsv.

**Verified:** 2026-02-15T12:30:00Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Running the pipeline produces summary.tsv with MeanDepthInsert, MeanDepthBackbone, MedianDepthInsert, MedianDepthBackbone rows | ✓ VERIFIED | compare_alignments.py lines 432-433, 443-444 write these rows. Test test_summary_tsv_coverage_rows_written validates output contains all rows. |
| 2 | summary.tsv contains BreadthInsert, BreadthBackbone, and BreadthBackbone_5x rows with fraction values 0.00-1.00 | ✓ VERIFIED | compare_alignments.py lines 434-439, 445-450 write breadth rows. Format is .2f (2 decimal places). Test validates "BreadthInsert\t0.95" format. |
| 3 | summary.tsv contains CoverageCV_Insert and CoverageCV_Backbone rows | ✓ VERIFIED | compare_alignments.py lines 440, 451 write CV rows. compute_coverage_stats returns cv key. |
| 4 | Existing summary.tsv rows (Plasmid, Human, Tied, Backbone_Only, Ambiguous, Verdict, Ratio, CoverageOutsideINSERT, MismatchesNearINSERT) are unchanged | ✓ VERIFIED | Coverage rows written AFTER existing rows (line 431 onwards, after line 424 MismatchesNearINSERT). Test test_summary_tsv_existing_rows_unchanged validates all old rows present. |
| 5 | When cDNA_positions.txt is missing, coverage metrics are computed for whole plasmid with fallback labels | ✓ VERIFIED | coverage_metrics.py lines 141-158 handle fallback mode (insert_region=None → whole plasmid as insert, backbone all 0.0). Fallback bool returned and written to summary.tsv line 454. Test test_summary_tsv_fallback_mode validates. |
| 6 | When no reads align, all coverage metrics are 0.00 | ✓ VERIFIED | compute_coverage_stats lines 77-87 handle empty array → all metrics 0.0. Test test_compute_coverage_stats_empty_array validates. Test test_compute_region_coverage_metrics_zero_depth validates BAM with zero depth. |
| 7 | Single-sample HTML report contains a Coverage Analysis card with depth, breadth, and uniformity metrics | ✓ VERIFIED | report_template.html lines 604-648 contains Coverage Analysis card with table. Test test_coverage_metrics_in_report validates card presence and metric display. |
| 8 | Coverage Analysis card displays insert and backbone metrics side by side | ✓ VERIFIED | Template has table with Insert Region and Backbone columns (lines 616-617), 5 metric rows (Mean Depth, Median Depth, Breadth 1x, Breadth 5x, Uniformity CV). |
| 9 | When all coverage values are zero, the card still renders with zero values displayed | ✓ VERIFIED | Template uses Jinja2 `default('0.00')` filter (lines 623-644) for graceful zero display. Test test_coverage_metrics_missing_graceful validates backward compatibility with pre-Phase-9 summary.tsv (no coverage rows). |
| 10 | When using whole-plasmid fallback, the card shows a warning: Insert region not defined | ✓ VERIFIED | Template lines 607-611 show fallback warning when coverage_fallback is true. Test test_coverage_fallback_warning validates warning appears. |
| 11 | Existing report sections (verdict, stats, gauge, plots, read assignment, insert region analysis) are unchanged | ✓ VERIFIED | Coverage Analysis card inserted at line 604, BEFORE Run Details section (line 650), AFTER Insert Region Analysis. No modifications to existing sections. Git diff shows only additions, no changes to existing report structure. |

**Score:** 11/11 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `plasmicheck/scripts/coverage_metrics.py` | Coverage computation engine with per-region metrics | ✓ VERIFIED | 185 lines. Exports get_depth_array, compute_coverage_stats, compute_region_coverage_metrics. Uses pysam.count_coverage with quality_threshold=0. Handles edge cases (empty arrays, single-element CV). |
| `plasmicheck/config.json` | breadth_thresholds configuration | ✓ VERIFIED | Line 6: "breadth_thresholds": [5]. Loaded in coverage_metrics.py via get_config(). |
| `tests/test_coverage_metrics.py` | Unit tests for coverage metric computation | ✓ VERIFIED | 598 lines, 21 unit tests. Covers get_depth_array (3 tests), compute_coverage_stats (8 tests), compute_region_coverage_metrics (7 tests), summary.tsv integration (3 tests). All tests pass. Min_lines requirement: 100 (actual: 598). |
| `plasmicheck/templates/report_template.html` | Coverage Analysis card in HTML report | ✓ VERIFIED | Lines 604-648 contain Coverage Analysis card with pc-table layout. Uses existing Bootstrap 5.3 CSS classes (pc-card, pc-table, pc-num, pc-notice). Card always visible, uses Jinja2 default() filters for backward compatibility. |
| `plasmicheck/scripts/generate_report.py` | Coverage metric parsing from summary.tsv and template variable passing | ✓ VERIFIED | Lines 591-599 extract 9 coverage metrics using _get_metric helper. Lines 601-602 detect coverage_fallback from CoverageFallback row. Lines 390, 461, 521, 707 pass coverage_metrics dict and coverage_fallback bool to template. |
| `tests/test_generate_report.py` | Tests for coverage metric display in reports | ✓ VERIFIED | Lines 103-260 contain 3 new tests: test_coverage_metrics_in_report, test_coverage_metrics_missing_graceful, test_coverage_fallback_warning. Tests validate Coverage Analysis card rendering, backward compatibility, and fallback warning. Min_lines requirement: 10 (actual: 158 lines across 3 tests). |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| compare_alignments.py | coverage_metrics.py | import and call compute_region_coverage_metrics() | ✓ WIRED | Line 14: `from .coverage_metrics import compute_region_coverage_metrics`. Line 360: function called with plasmid_bam, insert_region, BREADTH_THRESHOLDS. Result unpacked as (coverage_metrics, coverage_fallback). |
| compare_alignments.py | summary.tsv | summary_file.write() calls for new coverage rows | ✓ WIRED | Lines 432-454 write 10 coverage rows (MeanDepthInsert, MedianDepthInsert, BreadthInsert, BreadthInsert_5x, CoverageCV_Insert, MeanDepthBackbone, MedianDepthBackbone, BreadthBackbone, BreadthBackbone_5x, CoverageCV_Backbone) plus CoverageFallback indicator. All formatted with .2f. |
| coverage_metrics.py | pysam.count_coverage() | per-base depth extraction | ✓ WIRED | Line 43: `bam.count_coverage(contig=contig, start=start, stop=end, quality_threshold=0, read_callback="all")`. Lines 51-52: sum A+C+G+T arrays to get total depth. Used in get_depth_array function. |
| generate_report.py | summary.tsv | parsing MeanDepthInsert, BreadthBackbone etc. from summary_df | ✓ WIRED | Lines 591-599 extract coverage metrics from summary_df using _get_metric helper. Keys: MeanDepthInsert, MedianDepthInsert, BreadthInsert, BreadthInsert_5x, CoverageCV_Insert, MeanDepthBackbone, MedianDepthBackbone, BreadthBackbone, BreadthBackbone_5x, CoverageCV_Backbone. |
| generate_report.py | report_template.html | template.render() with coverage_metrics dict | ✓ WIRED | Lines 461, 521: template.render() called with coverage_metrics=coverage_metrics or {}, coverage_fallback=coverage_fallback. Both interactive and non-interactive templates receive variables. |
| report_template.html | coverage_metrics template variable | Jinja2 template rendering | ✓ WIRED | Lines 623-644: template accesses coverage_metrics.mean_depth_insert, coverage_metrics.mean_depth_backbone, etc. All use `| default('0.00')` filter for backward compatibility. Lines 607-611: coverage_fallback controls warning display. |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| COV-01: Per-region mean and median depth | ✓ SATISFIED | compute_coverage_stats returns mean_depth and median_depth. Written to summary.tsv as MeanDepthInsert, MedianDepthInsert, MeanDepthBackbone, MedianDepthBackbone. Displayed in report template. |
| COV-02: Breadth of coverage (fraction >=1 read) per region | ✓ SATISFIED | compute_coverage_stats returns breadth_1x (always computed, line 98-99). Written to summary.tsv as BreadthInsert, BreadthBackbone. Displayed in report template. |
| COV-03: Breadth at 5x threshold per region | ✓ SATISFIED | Config breadth_thresholds: [5]. compute_coverage_stats returns breadth_5x (lines 102-105). Written to summary.tsv as BreadthInsert_5x, BreadthBackbone_5x. Displayed in report template. |
| COV-04: Coverage uniformity (CV) for backbone | ✓ SATISFIED | compute_coverage_stats returns cv (coefficient of variation, line 108-111). Handles edge cases (zero mean, single element). Written to summary.tsv as CoverageCV_Insert, CoverageCV_Backbone. Displayed in report template. |
| COV-05: Metrics as additional rows in summary.tsv (backward compatible) | ✓ SATISFIED | Coverage rows written AFTER existing rows (line 431+). Test test_summary_tsv_existing_rows_unchanged validates all old rows unchanged. Template uses default('0.00') for missing rows (backward compat with pre-Phase-9 data). |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| None | - | - | - | No anti-patterns detected. Code follows established patterns: separate module for coverage logic, config-driven thresholds, comprehensive error handling, 219 passing tests including 24 new coverage tests (21 in test_coverage_metrics.py, 3 in test_generate_report.py). |

### Human Verification Required

None. All verifications completed programmatically:
- Code structure verified via file inspection
- Wiring verified via grep pattern matching
- Tests verified via pytest execution (219 passed)
- CI checks verified via make ci-check (all passed)
- Template structure verified via content inspection

### Gaps Summary

No gaps found. Phase 9 goal fully achieved:

**Core Deliverables:**
1. ✓ Coverage metrics computation engine (coverage_metrics.py)
2. ✓ Per-region statistics (insert vs backbone)
3. ✓ 10 new rows in summary.tsv with backward compatibility
4. ✓ Coverage Analysis card in HTML reports
5. ✓ Fallback mode for missing insert region
6. ✓ Comprehensive test coverage (24 new tests)

**Success Criteria (from ROADMAP.md):**
1. ✓ Running the pipeline produces summary.tsv with all 10 coverage metric rows
2. ✓ Existing summary.tsv rows unchanged (backward compatible)
3. ✓ Coverage metrics computed from actual per-base depth via pysam (not read count approximations)
4. ✓ Contaminated vs clean samples distinguishable via backbone coverage (backbone breadth/depth metrics enable this analysis)

**Requirements (COV-01 through COV-05):**
All 5 requirements satisfied (see Requirements Coverage table above).

**Quality:**
- All 219 tests pass (195 baseline + 24 new)
- make ci-check passes (lint, format, typecheck, tests)
- No stub patterns or anti-patterns detected
- Type annotations throughout (mypy strict mode)
- Edge cases handled (empty arrays, zero depth, single elements, missing insert region)

---

*Verified: 2026-02-15T12:30:00Z*
*Verifier: Claude (gsd-verifier)*
