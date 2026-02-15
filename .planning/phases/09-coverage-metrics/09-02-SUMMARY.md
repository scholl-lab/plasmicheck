---
phase: 09-coverage-metrics
plan: 02
subsystem: reporting
tags: [jinja2, html, bootstrap, plotly, coverage-analysis]

# Dependency graph
requires:
  - phase: 09-coverage-metrics
    plan: 01
    provides: Coverage metrics computation and 10 new rows in summary.tsv
provides:
  - Coverage Analysis card in single-sample HTML reports showing depth, breadth, and uniformity metrics for insert and backbone regions
  - Template variable coverage_metrics passed from generate_report.py to report_template.html
  - Fallback warning display when insert region not defined
  - 3 new report tests for coverage metrics display
affects: [11-summary-report-integration, reporting]

# Tech tracking
tech-stack:
  added: []
  patterns: [Jinja2 template variable passing with backward compatibility via default() filter, Bootstrap 5.3 table layout for metric display]

key-files:
  created: []
  modified:
    - plasmicheck/scripts/generate_report.py
    - plasmicheck/templates/report_template.html
    - tests/test_generate_report.py

key-decisions:
  - "Always show Coverage Analysis card (even with zero values) per CONTEXT decision - avoids confusion when metrics exist but card is hidden"
  - "Use Jinja2 default() filter for backward compatibility - old summary.tsv files without coverage rows show 0.00 instead of error"
  - "Display fallback warning inline in card (not as separate notice) - keeps coverage context together"

patterns-established:
  - "Coverage metrics passed as dict to template - enables clean template code with coverage_metrics.mean_depth_insert syntax"
  - "Fallback indicator passed as boolean (coverage_fallback) - template decides how to render warning"
  - "Template parameters default to empty dict/False - backward compatible with code that doesn't pass coverage_metrics"

# Metrics
duration: 5min
completed: 2026-02-15
---

# Phase 09 Plan 02: Coverage Metrics Display Summary

**Coverage Analysis card in per-sample HTML reports displaying depth, breadth, and uniformity metrics for insert and backbone regions with fallback warning**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-15T07:05:35Z
- **Completed:** 2026-02-15T07:08:28Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Coverage Analysis card added to report_template.html after Insert Region Analysis section
- Table displays 5 metrics (Mean Depth, Median Depth, Breadth ≥1x, Breadth ≥5x, Uniformity CV) for both Insert Region and Backbone
- Fallback warning shown when insert region not defined (whole-plasmid metrics computed)
- generate_report.py parses 10 coverage metrics from summary.tsv and passes to template
- Backward compatible with pre-Phase-9 summary.tsv files (missing coverage rows show 0.00)
- All 219 tests passing (up from 216) with 3 new report-specific coverage tests

## Task Commits

Each task was committed atomically:

1. **Task 1: Parse coverage metrics in generate_report.py and pass to template** - `c013759` (feat)
   - Extracted 10 coverage metrics from summary.tsv using _get_metric() helper
   - Added coverage_metrics dict and coverage_fallback bool to generate_report() signature
   - Passed coverage_metrics and coverage_fallback to both interactive and non-interactive template.render() calls
   - Backward compatible: coverage_metrics defaults to None, fallback defaults to False

2. **Task 2: Coverage Analysis card in HTML template + report tests** - `ac36805` (feat)
   - Added Coverage Analysis card in report_template.html with Bootstrap 5.3 pc-table styling
   - Table shows Insert Region and Backbone columns with 5 metric rows
   - Fallback warning displayed using existing pc-notice class when coverage_fallback=True
   - Card always visible, uses Jinja2 | default('0.00') for missing values
   - Added 3 comprehensive tests: test_coverage_metrics_in_report, test_coverage_metrics_missing_graceful, test_coverage_fallback_warning

## Files Created/Modified
- `plasmicheck/scripts/generate_report.py` - Parse coverage metrics from summary.tsv, pass coverage_metrics dict and coverage_fallback bool to template
- `plasmicheck/templates/report_template.html` - Coverage Analysis card with table showing insert/backbone metrics
- `tests/test_generate_report.py` - 3 new tests for coverage metrics display in reports

## Decisions Made

**1. Always show Coverage Analysis card (even with zero values)**
- Rationale: CONTEXT decision (09-CONTEXT.md) specified "always visible, even when all values are zero" to avoid confusion when metrics exist but card is hidden. Users should see the card structure to understand what coverage metrics are available.

**2. Use Jinja2 default() filter for backward compatibility**
- Rationale: Old summary.tsv files (pre-Phase-9) don't have coverage metric rows. Using `| default('0.00')` in template renders 0.00 instead of raising error, ensuring reports can still be generated from legacy data.

**3. Display fallback warning inline in card**
- Rationale: Keeps coverage context together. Warning appears at top of Coverage Analysis card (not as separate notice elsewhere) so users immediately understand that metrics are whole-plasmid, not insert/backbone split.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all tasks completed as planned, all tests passed on first run.

## Next Phase Readiness

- Coverage metrics now visible in single-sample HTML reports
- Ready for Phase 10 (Resistance Gene Detection) - can use same coverage computation infrastructure for gene-level analysis
- Ready for Phase 11 (Summary Report Integration) - coverage metrics can be aggregated across samples in summary reports (generate_summary_reports.py will need to parse coverage_metrics from individual summary.tsv files)

**Blockers:** None

**Concerns:** None - all success criteria met:
- Coverage Analysis card visible in HTML report ✓
- Table shows Insert Region and Backbone columns with depth, breadth, CV metrics ✓
- Card always visible (even with zero values) ✓
- Fallback warning displayed when insert region not defined ✓
- Existing report sections unchanged ✓
- All 219 tests passing including 3+ new report tests ✓

---
*Phase: 09-coverage-metrics*
*Completed: 2026-02-15*
