---
phase: 11-summary-report-integration
plan: 01
subsystem: reporting
tags: [plotly, matplotlib, pandas, summary-reports, heatmaps, boxplots]

# Dependency graph
requires:
  - phase: 10-resistance-gene-detection
    provides: ResGene_ metrics in summary.tsv and ResistanceGenesAvailable flag
  - phase: 09-coverage-metrics
    provides: Coverage metrics in summary.tsv (MeanDepthInsert, CoverageOutsideINSERT, etc.)
provides:
  - New plot generation functions (create_metric_heatmap_data, plot_metric_heatmap, plot_coverage_outside_insert_boxplot, create_resistance_gene_table)
  - Matplotlib fallback functions for new plots (generate_metric_heatmap_matplotlib, generate_coverage_boxplot_matplotlib)
  - depth_heatmap_cap configuration setting
affects: [11-02, summary-report-template-integration]

# Tech tracking
tech-stack:
  added: []
  patterns: [generic-metric-heatmap-pattern, resistance-gene-table-extraction]

key-files:
  created: []
  modified:
    - plasmicheck/scripts/generate_summary_reports.py
    - plasmicheck/scripts/plotting/matplotlib_backend.py
    - plasmicheck/config.json

key-decisions:
  - "Generic plot_metric_heatmap() function handles all metric heatmaps with consistent blue-to-red color scale"
  - "depth_heatmap_cap defaults to 100 to prevent extreme outliers from dominating color scale"
  - "create_resistance_gene_table() returns empty DataFrame with correct columns when no genes found (graceful degradation)"
  - "Matplotlib functions follow existing patterns: sns.set_theme, empty-check, 300 dpi, logger.info"

patterns-established:
  - "Metric extraction pattern: create_metric_heatmap_data(summary_df, metric_category) returns pivoted Sample x Plasmid DataFrame"
  - "Resistance gene table flattening: parse ResistanceGenes rows, split by comma, lookup ResGene_{sanitized}_MeanDepth and Breadth"
  - "Graceful degradation: empty DataFrames for missing/pre-v0.33.0 metrics prevent errors in report generation"

# Metrics
duration: 5min
completed: 2026-02-15
---

# Phase 11 Plan 01: Plot Building Blocks Summary

**Six new pure functions for Phase 11 summary report plots: metric heatmaps (Plotly + matplotlib), coverage boxplots, and resistance gene table extraction**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-15T08:48:38Z
- **Completed:** 2026-02-15T08:53:23Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Four new functions in generate_summary_reports.py for data extraction and Plotly plot generation
- Two new matplotlib fallback functions matching existing backend patterns
- Config setting for depth heatmap cap (100) to handle extreme outliers
- All functions handle missing data gracefully (empty DataFrame returns)
- Full type annotations and docstrings throughout

## Task Commits

Each task was committed atomically:

1. **Task 1: Add metric extraction and Plotly plot functions to generate_summary_reports.py** - `57d6e63` (feat)
2. **Task 2: Add matplotlib fallback functions to matplotlib_backend.py** - `616d1f7` (feat)

## Files Created/Modified
- `plasmicheck/scripts/generate_summary_reports.py` - Added create_metric_heatmap_data(), plot_metric_heatmap(), plot_coverage_outside_insert_boxplot(), create_resistance_gene_table() and DEPTH_HEATMAP_CAP config loading
- `plasmicheck/scripts/plotting/matplotlib_backend.py` - Added generate_metric_heatmap_matplotlib(), generate_coverage_boxplot_matplotlib()
- `plasmicheck/config.json` - Added depth_heatmap_cap: 100 to plot_summary section

## Decisions Made

**Generic metric heatmap function:**
- Single plot_metric_heatmap() function handles all metric types (depth, breadth, mismatches, coverage outside insert)
- Consistent blue-to-red color scale matching existing contamination ratio heatmap
- Optional vmax parameter allows capping extreme values

**Graceful degradation for pre-v0.33.0 data:**
- create_metric_heatmap_data() returns empty DataFrame when metric category not found
- Prevents errors when running on older summary.tsv files without new metrics
- Empty DataFrame check in plot functions logs warning and returns ("", None)

**Resistance gene table extraction:**
- Flattens ResistanceGenes comma-separated list into one row per gene
- Sanitizes gene names (replace "/" and " " with "_") for metric lookup
- Returns empty DataFrame with correct columns when ResistanceGenesAvailable is false

**Matplotlib backend patterns:**
- Both new functions follow existing patterns: sns.set_theme, empty-check with placeholder, 300 dpi save, logger.info
- Custom LinearSegmentedColormap from lightblue to red for metric heatmaps
- Proper handling of NaN values with mask parameter in sns.heatmap

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

**Minor lint/type fixes during development:**
- Removed unused pandas import from plot_metric_heatmap (only plotly needed)
- Fixed matplotlib boxplot call to avoid mypy type errors (separate labels list from values)
- Removed unused `bp` variable from boxplot function

All resolved immediately with no impact on plan execution.

## Next Phase Readiness

**Ready for Plan 02 (Wire new plots into report pipeline):**
- All building-block functions exist and pass typecheck/lint
- Functions are pure (no side effects) - can be called by Plan 02 without affecting existing flow
- Graceful empty DataFrame handling ensures backward compatibility
- 280 tests still passing (existing behavior unchanged)

**Key integration points for Plan 02:**
- Call create_metric_heatmap_data() for MeanDepthInsert, BreadthInsert, etc.
- Call plot_metric_heatmap() with vmax=DEPTH_HEATMAP_CAP for depth metrics
- Call plot_coverage_outside_insert_boxplot() for coverage distribution
- Call create_resistance_gene_table() for gene table in report template
- Pass HTML/PNG filenames to generate_report() function

**No blockers or concerns.**

---
*Phase: 11-summary-report-integration*
*Completed: 2026-02-15*
