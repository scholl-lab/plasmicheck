# PlasmiCheck — Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.32.0 Performance Optimization

## Current Position

**Phase:** 6 - Alignment Optimization (Complete)
**Plan:** 2 of 2 in current phase
**Status:** Phase 6 complete
**Progress:** [████████░░] 8/8 plans (100%)

Last activity: 2026-02-14 — Completed 06-02-PLAN.md (CLI thread integration)

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
- Use subprocess for lazy import verification tests (05-04)
- Add type: ignore[import-untyped] for kaleido imports (05-04)
- Parametrize CLI flag tests across subcommands (05-04)
- 5-tier CPU detection chain: SLURM → cgroup v2 → cgroup v1 → os.cpu_count → fallback(4) (06-01)
- Thread allocation: min 2, max 16 total; 80% minimap2, remainder samtools (max 4) (06-01)
- samtools sort -m 2G flag for predictable memory usage (06-01)
- align_reads() thread/memory parameters optional with config defaults (06-01)
- --threads CLI flag with None default enables optional override pattern (06-02)
- Thread detection at pipeline start with source logging for transparency (06-02)
- Sequential alignment with full thread allocation per step (no concurrency) (06-02)

### Todos

- [x] TEST-01: Regression tests (04-01) — Complete
- [x] TEST-02: Benchmark script (04-02) — Complete
- [x] REPT-01: Default runs WITHOUT PNGs (05-02, 05-03, 05-04) — Complete
- [x] REPT-02: --static-report flag (05-01, 05-02, 05-03, 05-04) — Complete
- [x] REPT-03: Directory mode shared plotly.min.js (05-02, 05-03, 05-04) — Complete
- [x] REPT-04: --plotly-mode CLI flag (05-01, 05-02, 05-03, 05-04) — Complete
- [x] REPT-05: Kaleido start_sync_server() optimization (05-02, 05-03, 05-04) — Complete
- [x] REPT-06: Lazy imports (05-02, 05-03, 05-04) — Complete
- [ ] TEST-03: Air-gapped testing (deferred)
- [ ] Verify samtools version >=1.9 (collate requirement)
- [ ] Decide on matplotlib style config for visual consistency with Plotly

### Blockers

None currently identified.

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

## Session Continuity

**What we're building:** Performance optimization milestone (v0.32.0)

**What just happened:** Phase 6 (Alignment Optimization) complete — Integrated thread detection/allocation into CLI and pipeline, 149 tests passing

**Next step:** Phase 6 complete. Performance optimization complete. Ready for production testing and benchmarking.

**Key context for next session:**
- Phase numbering starts at 4 (continues from v0.31.0 Phase 3)
- Regression test validates correctness: `python scripts/regression_test.py`
- Benchmark measures per-step timing: `python scripts/benchmark.py`
- Phase 5 complete: Report optimization (no Kaleido overhead by default)
- Phase 6 complete: Alignment optimization (automatic CPU detection + optimal thread allocation)
- 149 unit tests passing, mypy strict, ruff clean
- --threads CLI flag available on pipeline subcommand
- Thread detection: SLURM → cgroup v2 → cgroup v1 → os.cpu_count → fallback(4)
- Thread allocation: 80/20 minimap2/samtools split, 2-16 CPU bounds, 2G sort memory
- All thread parameters logged with source for transparency

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-14 after 06-02 completion*
