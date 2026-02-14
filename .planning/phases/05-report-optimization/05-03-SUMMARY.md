---
phase: 05-report-optimization
plan: 03
subsystem: reporting
tags: [summary-reports, lazy-imports, plotly, kaleido, performance, png-export]

# Dependency graph
requires:
  - phase: 05-report-optimization
    plan: 01
    provides: CLI flags (--static-report, --plotly-mode) plumbed through to generate_summary_reports
provides:
  - Conditional PNG export for summary reports (default: interactive HTML only)
  - Plotly.js mode support (cdn/directory/embedded) in summary reports
  - Lazy imports for heavy dependencies (pandas, numpy, plotly, scipy, statsmodels, jinja2)
  - Kaleido performance optimization (single start_sync_server() call)
affects: [06-alignment-optimization]

# Tech tracking
tech-stack:
  added: []
  patterns: [lazy-imports-with-TYPE_CHECKING, conditional-png-export, kaleido-optimization]

key-files:
  created: []
  modified:
    - plasmicheck/scripts/generate_summary_reports.py
    - plasmicheck/templates/summary_template.html

key-decisions:
  - "Use TYPE_CHECKING pattern for type hints with lazy imports (pandas as pd in type hints, imported at function level)"
  - "Initialize kaleido once in main() before all plotting (not per-plot) to avoid repeated startup overhead"
  - "Generate PNGs and non-interactive HTML only when static_report=True (default: False)"
  - "Summary template defaults to embedded mode for backwards compatibility"

patterns-established:
  - "Lazy import pattern: TYPE_CHECKING block for types, function-level imports for runtime"
  - "Conditional resource generation: gate expensive operations behind boolean flag"
  - "Singleton resource initialization: kaleido.start_sync_server() once, not per-plot"
  - "Template mode fallback: {% if plotly_mode is not defined %}{% set plotly_mode = 'embedded' %}{% endif %}"

# Metrics
duration: 8min
completed: 2026-02-14
---

# Phase 05 Plan 03: Summary Reports PNG Export and Plotly.js Modes Summary

**Refactored generate_summary_reports.py with lazy imports, conditional PNG export, and plotly.js mode support for multi-sample summary reports**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-14T10:15:18Z
- **Completed:** 2026-02-14T10:23:09Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Moved heavy imports to function level: pandas (6 functions), numpy (2 functions), plotly.express (2 functions), scipy.stats (1 function), statsmodels.stats.multitest (1 function), jinja2 (1 function)
- Used TYPE_CHECKING pattern for type hints (pd.DataFrame annotations) without import-time overhead
- Added `static_report` parameter to `plot_boxplot()` and `plot_heatmap()` - PNG generated only when True
- Initialized kaleido once in `main()` before all plotting (not inside each plot function)
- Updated `generate_report()` to conditionally generate non-interactive HTML only when `static_report=True`
- Added `--static-report` and `--plotly-mode` CLI flags to __main__ block
- Updated summary_template.html with plotly.js mode conditionals (cdn/directory/embedded)
- Added CDN fallback for directory mode: `if (typeof Plotly === 'undefined') { document.write(...) }`
- Template defaults to embedded mode for backwards compatibility

## Task Commits

Each task was committed atomically:

1. **Task 1: Lazy imports and conditional PNG export in generate_summary_reports.py** - `26c5008` (feat)
2. **Task 2: Plotly.js mode support in summary template** - `9dd7c37` (feat)

## Files Created/Modified

- `plasmicheck/scripts/generate_summary_reports.py` - Lazy imports, conditional PNG, kaleido optimization, plotly.js asset management, new parameters
- `plasmicheck/templates/summary_template.html` - Plotly.js mode conditionals (cdn/directory/embedded), CDN fallback, backwards compatibility

## Decisions Made

- **TYPE_CHECKING pattern:** Import pandas under `if TYPE_CHECKING:` block for type hints, then import at function level where actually used - avoids import-time overhead while maintaining type safety
- **Single kaleido initialization:** Call `kaleido.start_sync_server()` once in `main()` before all plotting instead of inside each plot function - eliminates repeated startup overhead
- **Default to interactive-only:** `static_report` defaults to False, so default run skips PNG generation and non-interactive HTML - users must opt-in to the slow path
- **Template backwards compatibility:** Summary template checks `{% if plotly_mode is not defined %}` and defaults to embedded mode - ensures existing code continues to work

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

**Auto-revert by ruff:**
During initial implementation, ruff's auto-fix reverted file changes. Resolved by using `Write` tool to overwrite the entire file with final state instead of incremental `Edit` operations that could be auto-fixed.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for integration with cli.py:**
- summary_reports subcommand passes `--static-report` and `--plotly-mode` to `generate_summary_reports.main()` (requires Plan 05-04 to wire these through)
- Default behavior (no flags) generates only interactive HTML, skips PNG export
- Performance improvement: eliminates ~5s kaleido overhead and heavy import time for default runs

**Parallel with Plan 05-02:**
- Plan 05-02 applies same pattern to `generate_report.py` (single-sample reports)
- Plan 05-03 applies same pattern to `generate_summary_reports.py` (multi-sample reports)
- Both modify different files, no merge conflicts

**Blockers:** None

**Performance impact:**
- Default summary_reports run now skips PNG generation and kaleido initialization
- Heavy imports (pandas, numpy, plotly, scipy, statsmodels, jinja2) deferred to function call time
- Expected speedup: ~5-6s for small datasets (kaleido overhead + import overhead)
- Interactive HTML uses plotly.min.js from shared assets/ directory (Plan 05-01's output_root enables this)

---
*Phase: 05-report-optimization*
*Completed: 2026-02-14*
