---
phase: 05-report-optimization
plan: 01
subsystem: cli
tags: [argparse, cli-flags, report-generation, plotly, performance]

# Dependency graph
requires:
  - phase: 04-foundation
    provides: baseline regression tests and benchmark script
provides:
  - CLI flags (--static-report, --plotly-mode) on pipeline, report, summary_reports subcommands
  - Flag pass-through from CLI to run_pipeline to generate_report
  - Foundation for conditional PNG export and plotly.js mode selection
affects: [05-02-report-png-export, 05-03-report-plotly-modes, 06-alignment-optimization, 07-comparison-optimization]

# Tech tracking
tech-stack:
  added: []
  patterns: [shared-parent-parser-pattern]

key-files:
  created: []
  modified:
    - plasmicheck/cli.py
    - plasmicheck/scripts/run_pipeline.py

key-decisions:
  - "Use shared _report_parser parent to avoid flag definition duplication across three subcommands"
  - "Pass output_root=output_folder to enable shared assets/ directory for directory-mode plotly.js"
  - "Default plotly-mode to 'directory' for optimal balance of speed and offline capability"

patterns-established:
  - "Shared parent parsers (_logging_parser, _threshold_parser, _report_parser) for flag reuse across subcommands"
  - "Flag pass-through: CLI → run_pipeline → generate_report with explicit parameter names"

# Metrics
duration: 3.5min
completed: 2026-02-14
---

# Phase 05 Plan 01: CLI Flags for Report Optimization Summary

**Added --static-report and --plotly-mode flags to all report-generating subcommands with full pass-through to generation functions**

## Performance

- **Duration:** 3.5 min
- **Started:** 2026-02-14T10:05:44Z
- **Completed:** 2026-02-14T10:09:13Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Created shared `_report_parser` with `--static-report` (boolean) and `--plotly-mode` (choices: cdn, directory, embedded) flags
- Added flags to `pipeline`, `report`, and `summary_reports` subcommands via parent parser inheritance
- Wired flags through `run_pipeline()` to `generate_report()` call including critical `output_root` parameter
- Updated standalone `run_pipeline.py` __main__ block with same flags for direct script execution

## Task Commits

Each task was committed atomically:

1. **Task 1: Add --static-report and --plotly-mode flags to CLI subcommands** - `63dd6dc` (feat)
2. **Task 2: Wire flags through run_pipeline to generate_report (including output_root)** - `ce5831e` (feat)

## Files Created/Modified

- `plasmicheck/cli.py` - Added shared _report_parser with two new flags, added parent to three subparsers, passed flags to three command handlers
- `plasmicheck/scripts/run_pipeline.py` - Added parameters to function signature, passed flags + output_root to generate_report(), added flags to __main__ argparse

## Decisions Made

- **Shared parent parser pattern:** Created `_report_parser` following existing `_logging_parser` and `_threshold_parser` patterns to avoid duplicating flag definitions across three subcommands
- **output_root parameter:** Passed `output_root=output_folder` to `generate_report()` to enable Plan 05-02's shared `output/assets/plotly.min.js` directory without fragile path walking
- **Default to directory mode:** Set `--plotly-mode` default to "directory" for optimal balance of speed (no 5s Kaleido overhead) and offline capability (self-contained assets)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None. Expected type errors from mypy documented in plan (generate_report.main() and generate_summary_reports.main() don't yet accept the new parameters - will be resolved in Plans 02 and 03).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Plan 05-02 (Report PNG Export):**
- `--static-report` flag plumbed through to generate_report
- Plan 02 can implement conditional PNG export logic in generate_report.main()

**Ready for Plan 05-03 (Plotly.js Modes):**
- `--plotly-mode` flag plumbed through to generate_report
- `output_root` parameter enables shared assets/ directory
- Plan 03 can implement three plotly.js inclusion modes

**Blockers:** None

**Type errors:** 7 expected mypy errors in cli.py and run_pipeline.py from calling generate_report.main() and generate_summary_reports.main() with new kwargs. These will be resolved when Plans 02 and 03 update the function signatures.

---
*Phase: 05-report-optimization*
*Completed: 2026-02-14*
