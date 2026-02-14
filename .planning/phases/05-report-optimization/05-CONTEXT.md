# Phase 5: Report Optimization - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Eliminate the report generation bottleneck (83.6% of pipeline time) by making static PNG export opt-in and reducing interactive HTML report size from 9.6 MB to ~19 KB. The interactive HTML report with Plotly graphs remains the default output. This phase covers the `pipeline`, `report`, and `summary_reports` subcommands.

</domain>

<decisions>
## Implementation Decisions

### Default behavior change
- Silent change: stop generating PNG reports by default (no deprecation warning)
- Interactive HTML reports (Plotly) are always generated — this is the primary output
- Non-interactive HTML reports (embedded PNGs) are also gated behind `--static-report` — only generated when user explicitly opts in
- When `--static-report` is used, PNGs are generated **alongside** interactive HTML (not replacing it)
- Same treatment for summary reports: skip static PNGs by default, opt-in with `--static-report`

### Plotly.js inclusion modes
- Default mode: **directory** — one shared `plotly.min.js` in `output/assets/`
- All directory-mode reports include CDN fallback: try local file first, fall back to CDN if not found
- Three modes available: `cdn`, `directory`, `embedded`
- Archives (`--archive_output`): do NOT include plotly.min.js — CDN fallback handles portability
- CDN version: pin to specific plotly.js version (not "latest") for reproducibility

### Air-gapped / offline behavior
- Air-gapped use is very rare — not a primary concern
- Drop TEST-03 (air-gapped Docker test) from this phase — not worth the complexity
- Embedded mode remains available for users who need fully self-contained reports
- If directory mode can't write plotly.min.js (permissions): **error and stop** — don't silently fall back

### CLI flag design
- `--static-report`: opt-in PNG generation. No short form (not used frequently enough)
- `--plotly-mode {cdn,directory,embedded}`: CLI flag only, no config.json setting
- Both flags available on all report-generating subcommands: `pipeline`, `report`, `summary_reports`
- Consistent behavior across all three subcommands

### Claude's Discretion
- Lazy import implementation details (which modules, import placement)
- Kaleido `start_sync_server()` initialization specifics
- Exact plotly.js version to pin to
- CDN fallback JavaScript implementation
- Output directory structure for `assets/` folder

</decisions>

<specifics>
## Specific Ideas

- CDN fallback pattern for directory-mode reports: `<script src="assets/plotly.min.js">` with inline JS that loads from CDN if `typeof Plotly === 'undefined'`
- Shared plotly.min.js goes in `output/assets/plotly.min.js` — one copy per batch run, all reports reference relative path

</specifics>

<deferred>
## Deferred Ideas

- TEST-03 air-gapped Docker test — dropped (air-gapped is very rare for this userbase)

</deferred>

---

*Phase: 05-report-optimization*
*Context gathered: 2026-02-14*
