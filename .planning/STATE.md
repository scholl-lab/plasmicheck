# PlasmiCheck — Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.32.0 Performance Optimization

## Current Position

**Phase:** 5 - Report Optimization (In Progress)
**Status:** Plan 03 of 03 complete
**Progress:** [████░░░░░░] 5/18 requirements (28%)

Last activity: 2026-02-14 — Completed 05-03-PLAN.md (multi-sample summary report refactor)

## Performance Metrics

**Baseline (v0.31.0):**
- Small dataset (200 reads): 5.5s (benchmark measurement)
- Large dataset (1M reads): ~4.5 min per combination
- Report generation: 91.7% of time (5.1s/5.5s) (benchmark measurement)
- Interactive HTML size: 9.6 MB (embedded plotly.js)

**Target (v0.32.0):**
- Small dataset: <2s (13x speedup)
- Large batch (50 combinations): 103 min (was 223 min, 2.2x speedup)
- Report size: 19 KB (99.8% reduction)

## Accumulated Context

### Decisions

- Performance analysis completed (PERFORMANCE_ANALYSIS.md) — profiled all pipeline steps
- Kaleido v1.2.0 is the root cause of report generation bottleneck (5.1s for report step)
- Streaming BAM comparison (v0.31.0) already reduced memory from O(n) to O(1)
- Four-phase roadmap prioritizes impact-per-effort: Foundation → Report → Alignment → Comparison
- Keep Kaleido v1.2.0 (not downgrade to 0.2.1) but make PNG export opt-in
- Use directory mode plotly.js with offline fallback for air-gapped environments
- ThreadPoolExecutor for parallel alignment (not ProcessPoolExecutor)
- samtools collate for name grouping (30-50% faster than sort -n)
- Use atomic write pattern (tempfile + os.replace) for baseline caching
- Store baselines as JSON for human readability and easy diff inspection
- Compare verdicts exactly, ratios with ±0.001 tolerance for regression tests
- Use shared _report_parser parent to avoid flag definition duplication (05-01)
- Pass output_root to enable shared assets/ directory for plotly.js (05-01)
- Default plotly-mode to 'directory' for optimal speed/offline balance (05-01)
- Use TYPE_CHECKING for pd.DataFrame type hints with lazy imports (05-02, 05-03)
- write_html with include_plotlyjs=False, full_html=False for template embedding (05-02, 05-03)
- Always generate interactive HTML, conditionally generate non-interactive (05-02, 05-03)
- Kaleido start_sync_server() called once before write_image() (05-02, 05-03)
- Summary reports apply same optimization pattern as single-sample reports (05-03)

### Todos

- [x] TEST-01: Regression tests (04-01) — Complete
- [x] TEST-02: Benchmark script (04-02) — Complete
- [ ] TEST-03: Air-gapped testing (deferred to Phase 5)
- [ ] Verify samtools version >=1.9 (collate requirement)
- [ ] Decide on matplotlib style config for visual consistency with Plotly

### Blockers

None currently identified.

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

## Session Continuity

**What we're building:** Performance optimization milestone (v0.32.0)

**What just happened:** Completed 05-03 (multi-sample summary refactor) — lazy imports, conditional PNG export, plotly.js mode support in generate_summary_reports.py

**Next step:** Execute 05-04 (integrate summary_reports flags) — wire CLI flags through to generate_summary_reports.main()

**Key context for next session:**
- Phase numbering starts at 4 (continues from v0.31.0 Phase 3)
- Regression test validates correctness: `python scripts/regression_test.py`
- Benchmark measures per-step timing: `python scripts/benchmark.py`
- Report generation is 91.7% of pipeline time — Phase 5 targets this bottleneck
- Default pipeline behavior now: interactive HTML only, no PNG, no Kaleido startup
- Both generate_report.py and generate_summary_reports.py refactored with lazy imports and conditional PNG
- Expected mypy errors resolved: generate_report.main() and generate_summary_reports.main() now accept static_report and plotly_mode parameters
- Phase 05 plans 02 and 03 ran in parallel, both complete

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-14 after 05-03 completion*
