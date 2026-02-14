---
phase: 07-comparison-cleanup
plan: 03
subsystem: reporting
tags: [matplotlib, plotly, visualization, png, kaleido-alternative]

# Dependency graph
requires:
  - phase: 05-report-optimization
    provides: Plotly.js directory mode, static report flags, lazy imports
  - phase: 07-02
    provides: Human index upfront phase, batch resilience
provides:
  - matplotlib plotting backend for static PNG generation (ARCH-03)
  - --plot-backend CLI flag on pipeline, report, summary_reports
  - Kaleido-free PNG export option via matplotlib
  - 4 chart generators: boxplot, scatter, heatmap, summary boxplot
  - Plotly-matching color constants for visual consistency
affects: [users needing static reports without Kaleido, air-gapped deployments]

# Tech tracking
tech-stack:
  added: [matplotlib backend, seaborn for styling]
  patterns: [conditional backend dispatch, shared color constants, lazy imports for backends]

key-files:
  created:
    - plasmicheck/scripts/plotting/__init__.py
    - plasmicheck/scripts/plotting/colors.py
    - plasmicheck/scripts/plotting/matplotlib_backend.py
    - tests/test_matplotlib_backend.py
  modified:
    - plasmicheck/cli.py
    - plasmicheck/scripts/run_pipeline.py
    - plasmicheck/scripts/generate_report.py
    - plasmicheck/scripts/generate_summary_reports.py

key-decisions:
  - "matplotlib backend only for static PNG generation (interactive HTML always uses Plotly.js)"
  - "Default plot_backend='plotly' for backward compatibility"
  - "kaleido import only when static_report=True AND plot_backend='plotly'"
  - "Use seaborn whitegrid theme for matplotlib plots to match Plotly aesthetics"
  - "Plotly color constants shared via colors.py for consistent appearance"

patterns-established:
  - "Backend dispatch pattern: if plot_backend == 'matplotlib': ... else: (plotly path)"
  - "Lazy imports for plotting backends: import inside functions to avoid dependencies"
  - "Shared color constants module for cross-backend visual consistency"
  - "Empty DataFrame handling: create placeholder plots with warning, don't crash"

# Metrics
duration: 10min
completed: 2026-02-14
---

# Phase 07 Plan 03: Matplotlib Backend Summary

**matplotlib as alternative static PNG backend (no Kaleido), with --plot-backend CLI flag and 9 comprehensive tests**

## Performance

- **Duration:** 10 min
- **Started:** 2026-02-14T15:15:44Z
- **Completed:** 2026-02-14T15:25:50Z
- **Tasks:** 3/3 completed
- **Files modified:** 8 files (4 created, 4 modified)
- **Commits:** 4 (3 task commits + 1 style commit)
- **Tests:** 172 passing (9 new matplotlib tests)

## Accomplishments

- Created `plasmicheck/scripts/plotting/` subpackage with matplotlib backend
- 4 chart generators matching Plotly functionality: boxplot, scatter, heatmap, summary boxplot
- --plot-backend flag wired through CLI → pipeline → report generation (pipeline, report, summary_reports)
- Kaleido-free PNG export: `--plot-backend matplotlib --static-report` works without Kaleido installed
- Plotly-matching colors for visual consistency across backends
- 9 comprehensive tests covering PNG creation, color constants, empty DataFrame handling, CLI flags, kaleido-free operation

## Task Commits

Each task was committed atomically:

1. **Task 1: Create matplotlib plotting module** - `6c86108` (feat)
   - plasmicheck/scripts/plotting/__init__.py, colors.py, matplotlib_backend.py
   - PLOTLY_COLORS, ASSIGNMENT_COLORS, HEATMAP_COLORS constants
   - 4 chart generation functions with lazy imports and TYPE_CHECKING guards

2. **Task 2: Wire matplotlib backend through CLI and reports** - `fc40bf7` (feat)
   - CLI --plot-backend flag on _report_parser (shared by 3 subcommands)
   - run_pipeline.py: plot_backend parameter flow
   - generate_report.py: matplotlib dispatch for boxplot/scatter
   - generate_summary_reports.py: matplotlib dispatch for heatmap/summary boxplot
   - kaleido import only when plot_backend='plotly'

3. **Task 3: Add tests for matplotlib backend** - `9cdb223` (test)
   - tests/test_matplotlib_backend.py with 9 tests
   - PNG creation tests (4 chart types)
   - Color constant verification
   - Empty DataFrame handling (graceful degradation)
   - CLI flag existence and default value tests
   - Kaleido-free operation test (blocks kaleido import, verifies matplotlib works)

**Style formatting:** `8288e6e` (style: auto-format with ruff)

## Files Created/Modified

**Created:**
- `plasmicheck/scripts/plotting/__init__.py` - Plotting backends package marker
- `plasmicheck/scripts/plotting/colors.py` - PLOTLY_COLORS, ASSIGNMENT_COLORS, HEATMAP_COLORS
- `plasmicheck/scripts/plotting/matplotlib_backend.py` - 4 chart generation functions (boxplot, scatter, heatmap, summary boxplot)
- `tests/test_matplotlib_backend.py` - 9 comprehensive tests for matplotlib backend

**Modified:**
- `plasmicheck/cli.py` - Added --plot-backend flag to _report_parser, pass to subcommands
- `plasmicheck/scripts/run_pipeline.py` - plot_backend parameter, pass to generate_report
- `plasmicheck/scripts/generate_report.py` - plot_backend parameter, matplotlib dispatch in generate_plots()
- `plasmicheck/scripts/generate_summary_reports.py` - plot_backend parameter, matplotlib dispatch in plot_boxplot() and plot_heatmap()

## Decisions Made

**Interactive HTML always uses Plotly.js regardless of --plot-backend**
- Rationale: Interactive HTML requires Plotly.js for interactivity. matplotlib backend only affects static PNG generation (when --static-report is used).
- Impact: Users get consistent interactive experience. --plot-backend only matters for PNGs.

**Default plot_backend='plotly' for backward compatibility**
- Rationale: Existing users may rely on Kaleido behavior. Default to plotly to avoid breaking changes.
- Impact: Opt-in matplotlib backend. Users must explicitly request `--plot-backend matplotlib`.

**kaleido import only when static_report=True AND plot_backend='plotly'**
- Rationale: Avoid Kaleido dependency when not needed (default runs or matplotlib backend).
- Impact: Users can run without Kaleido installed if using matplotlib or not generating static reports.

**Use seaborn whitegrid theme for matplotlib plots**
- Rationale: Match Plotly's clean, professional aesthetic. Seaborn provides better defaults than matplotlib.
- Impact: Visual consistency between backends. Users can't easily tell which backend generated a plot.

**Shared color constants in colors.py**
- Rationale: ASSIGNMENT_COLORS (Plasmid/Human/Tied) must match exactly between backends for consistency.
- Impact: Single source of truth for colors. Future backends can import same constants.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

**mypy type errors with .values numpy arrays**
- Issue: matplotlib boxplot expects lists, but pandas .values returns numpy arrays. mypy complained about type mismatch.
- Resolution: Used .tolist() instead of .values to convert to Python lists.
- Impact: Cleaner type checking, no runtime difference.

**argparse doesn't show default in help text**
- Issue: Test expected "default: plotly" in `--help` output, but argparse doesn't show defaults for choices.
- Resolution: Changed test to inspect function signature instead of parsing help text.
- Impact: More robust test (checks actual default value, not help formatting).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for continuation:**
- matplotlib backend fully functional and tested
- Users can now choose between plotly/kaleido and matplotlib for static PNG reports
- No Kaleido dependency required when using matplotlib backend
- All 172 unit tests passing (9 new tests)
- CI checks clean (lint, format, typecheck, tests)

**No blockers identified.**

**Next steps (from 07-CONTEXT.md):**
- No additional matplotlib work planned in Phase 7
- Phase 7 focus shifts to comparison algorithm optimization (if planned)
- matplotlib backend is complete and ready for user adoption

---
*Phase: 07-comparison-cleanup*
*Completed: 2026-02-14*
