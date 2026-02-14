# PlasmiCheck — Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.32.0 Performance Optimization

## Current Position

**Phase:** 7 - Comparison & Cleanup ✓ Complete
**Status:** Phase 7 verified and complete (17/17 must-haves)
**Progress:** [██████████] 17/18 requirements (94%) — 1 dropped (TEST-03)

Last activity: 2026-02-14 — Phase 7 verified, all 4 phases complete

## Performance Metrics

**Baseline (v0.31.0):**
- Small dataset (200 reads): 5.5s (benchmark measurement)
- Large dataset (1M reads): ~4.5 min per combination
- Report generation: 91.7% of time (5.1s/5.5s) (benchmark measurement)
- Interactive HTML size: 9.6 MB (embedded plotly.js)

**Measured (v0.32.0-dev, after Phase 5+6):**
- Small dataset (200 reads): 0.58s (was 5.5s — 9.5x speedup)
- Real dataset (3M reads): 58.4s alignment (was 115.2s — 1.97x speedup)
- Report generation: 0.108s / 18.8% of total (was 91.7% — no longer the bottleneck)
- Key gain: samtools sort -m 2G reduced human alignment from 65.1s to 13.5s (4.8x)

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
- samtools sort -n for name grouping (collate benchmarked but 28-64x slower on filtered BAMs)
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
- Benchmarked: 1.97x alignment speedup on 3M read BAM (115.2s → 58.4s), -m 2G sort memory is the biggest win (06-benchmark)
- samtools collate benchmarked and reverted: 28-64x slower than sort -n on post-alignment filtered BAMs (<100K reads) due to hash-bucket fixed overhead (~2s) (07-01)
- BAM name grouping uses sort -n directly (collate code removed as dead code) (07-01)
- Hoist human index creation to upfront phase (before combination loop) to eliminate filelock overhead (07-02)
- Track built indexes in PipelinePlan.built_indexes set for defensive skip checks (07-02)
- No deduplication for plasmid indexes (per-combination, cheap to rebuild) (07-02)
- Batch resilience: Continue processing remaining combinations after one fails (07-02)
- Raise RuntimeError only if ALL combinations fail (07-02)
- Log per-combination timing at INFO level (format: 'label: NN.Ns') (07-02)
- matplotlib backend only for static PNG generation (interactive HTML always uses Plotly.js) (07-03)
- Default plot_backend='plotly' for backward compatibility (07-03)
- kaleido import only when static_report=True AND plot_backend='plotly' (07-03)
- Use seaborn whitegrid theme for matplotlib plots to match Plotly aesthetics (07-03)
- Plotly color constants shared via colors.py for consistent appearance (07-03)

### Todos

- [x] TEST-01: Regression tests (04-01) — Complete
- [x] TEST-02: Benchmark script (04-02) — Complete
- [x] REPT-01: Default runs WITHOUT PNGs (05-02, 05-03, 05-04) — Complete
- [x] REPT-02: --static-report flag (05-01, 05-02, 05-03, 05-04) — Complete
- [x] REPT-03: Directory mode shared plotly.min.js (05-02, 05-03, 05-04) — Complete
- [x] REPT-04: --plotly-mode CLI flag (05-01, 05-02, 05-03, 05-04) — Complete
- [x] REPT-05: Kaleido start_sync_server() optimization (05-02, 05-03, 05-04) — Complete
- [x] REPT-06: Lazy imports (05-02, 05-03, 05-04) — Complete
- [x] ALGN-01: Sequential alignment with full thread allocation (06-01, 06-02) — Complete
- [x] ALGN-02: CPU auto-detection with cgroup/SLURM awareness (06-01) — Complete
- [x] ALGN-03: --threads CLI flag (06-02) — Complete
- [x] ALGN-04: samtools sort -m 2G memory flag (06-01, 06-02) — Complete
- [x] COMP-01: BAM name grouping benchmarked, sort -n retained (07-01) — Complete
- [x] COMP-02: Supplementary ordering validated with sort -n (07-01) — Complete
- [x] ARCH-01: Hoist human index to upfront phase (07-02) — Complete
- [x] ARCH-02: PipelinePlan.built_indexes tracking (07-02) — Complete
- [x] ARCH-03: matplotlib backend for static PNG generation (07-03) — Complete
- [ ] TEST-03: Air-gapped testing (deferred)

### Blockers

None currently identified.

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

## Session Continuity

**What we're building:** Performance optimization milestone (v0.32.0)

**What just happened:** Phase 7 post-benchmark: collate reverted to sort -n (28-64x slower on filtered BAMs), dead code removed, all tests passing. All 4 phases (4-7) complete.

**Next step:** Audit milestone and complete v0.32.0

**Key context for next session:**
- Phase numbering starts at 4 (continues from v0.31.0 Phase 3)
- All 4 phases complete: Foundation (4), Report (5), Alignment (6), Comparison (7)
- 170 unit tests passing, mypy strict, ruff clean
- Regression test: `python scripts/regression_test.py`
- Benchmark: `python scripts/benchmark.py`
- Phase 5: Report optimization (no Kaleido overhead by default, 9.5x speedup)
- Phase 6: Alignment optimization (1.97x speedup, -m 2G sort memory)
- Phase 7: index dedup, batch resilience, matplotlib backend (collate reverted)
- New CLI flags: --static-report, --plotly-mode, --threads, --plot-backend
- Small dataset total: 0.704s (was 5.5s baseline, 8.2x speedup)

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-14 after collate revert and dead code cleanup*
