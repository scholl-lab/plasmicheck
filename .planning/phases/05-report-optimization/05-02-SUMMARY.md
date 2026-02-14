---
phase: 05-report-optimization
plan: 02
subsystem: report-generation
tags: [plotly, kaleido, jinja2, pandas, performance, lazy-imports]

# Dependency graph
requires:
  - phase: 05-01
    provides: CLI flags for static-report, plotly-mode, output-root parameters
provides:
  - Conditional PNG export gated behind static_report flag
  - Plotly.js mode-aware template with directory/cdn/embedded support
  - Lazy imports for pandas, plotly, jinja2 (import at function level)
  - Kaleido optimization with start_sync_server() before write_image()
  - Shared plotly.min.js asset management in output_root/assets/
affects: [05-03, run-pipeline, report-generation]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "TYPE_CHECKING for lazy import type annotations"
    - "Conditional Kaleido initialization only when static_report=True"
    - "Shared asset directory pattern for plotly.min.js"

key-files:
  created: []
  modified:
    - plasmicheck/scripts/generate_report.py
    - plasmicheck/templates/report_template.html

key-decisions:
  - "Use TYPE_CHECKING for pd.DataFrame type hints with lazy imports"
  - "Default plotly_mode='directory' for optimal speed/offline balance"
  - "Always generate interactive HTML, conditionally generate non-interactive"
  - "write_html with include_plotlyjs=False, full_html=False for template embedding"

patterns-established:
  - "Lazy imports: Heavy libraries (pandas, plotly, jinja2) imported at function level"
  - "TYPE_CHECKING pattern: Import types for annotations without runtime import"
  - "Asset management: _ensure_plotly_assets() helper for shared plotly.min.js"
  - "Plotly.js modes: directory (shared), cdn (internet), embedded (self-contained)"

# Metrics
duration: 5min
completed: 2026-02-14
---

# Phase 05 Plan 02: Single-Sample Report Refactor Summary

**Conditional PNG export, plotly.js directory mode, and lazy imports eliminate 91.7% bottleneck from default pipeline**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-14T10:14:24Z
- **Completed:** 2026-02-14T10:19:27Z
- **Tasks:** 2 (committed together as tightly coupled)
- **Files modified:** 2

## Accomplishments
- Lazy imports for pandas, plotly, jinja2 using TYPE_CHECKING pattern
- Conditional PNG export: static_report parameter gates Kaleido usage
- Kaleido optimization: start_sync_server() called once before write_image()
- Plotly.js modes: directory (shared assets), cdn (CDN fallback), embedded (inline)
- Default behavior: Interactive HTML only (~19KB vs 9.6MB), no PNG, no Kaleido startup

## Task Commits

1. **Tasks 1-2: Lazy imports, conditional PNG, plotly.js modes** - `24f13b5` (feat)

_Note: Tasks 1 and 2 committed together as they're tightly coupled (report generation + template)_

## Files Created/Modified
- `plasmicheck/scripts/generate_report.py` - Lazy imports, conditional PNG export, plotly.js asset management
- `plasmicheck/templates/report_template.html` - Plotly.js mode-aware script loading with CDN fallback

## Decisions Made

**1. TYPE_CHECKING for lazy import type annotations**
- **Context:** When imports are moved to function level, pd.DataFrame isn't defined at module level
- **Solution:** Use `if TYPE_CHECKING: import pandas as pd` for type hints only
- **Benefit:** No runtime import cost, mypy still validates types

**2. write_html with include_plotlyjs=False, full_html=False**
- **Context:** Previously write_html generated full HTML documents with embedded plotly.js
- **Solution:** Generate div+script fragments, let template handle plotly.js loading
- **Benefit:** Template controls plotly.js mode (directory/cdn/embedded), plots are embeddable

**3. Always generate interactive HTML, conditionally generate non-interactive**
- **Context:** Original code always generated both reports
- **Solution:** Interactive HTML always generated (fast), non-interactive only when static_report=True
- **Benefit:** Default runs skip PNG export and non-interactive HTML rendering

**4. Kaleido start_sync_server() optimization**
- **Context:** Kaleido sync server startup is slow (~3s overhead)
- **Solution:** Call start_sync_server() once before both write_image() calls
- **Benefit:** Single startup cost instead of per-plot startup

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

**Issue: Mypy errors with lazy imports and pd.DataFrame in type hints**
- **Problem:** Moving pandas import to function level made pd.DataFrame undefined for type annotations
- **Resolution:** Added `if TYPE_CHECKING: import pandas as pd` pattern for type-only imports
- **Verification:** make typecheck passes (only expected kaleido untyped warning)

## Next Phase Readiness

**Ready for:**
- Plan 05-03: Multi-sample summary report refactor (parallel, no dependencies)
- Integration with run_pipeline.py (Plan 05-01 already wired the flags)

**Notes:**
- generate_summary_reports.py still needs same treatment (Plan 05-03 responsibility)
- Expected mypy errors in generate_summary_reports.py until Plan 05-03 completes
- Default pipeline behavior now: interactive HTML only, no PNG, no Kaleido startup

---
*Phase: 05-report-optimization*
*Completed: 2026-02-14*
