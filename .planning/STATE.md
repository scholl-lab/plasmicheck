# PlasmiCheck — Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.32.0 Performance Optimization

## Current Position

**Phase:** 4 - Foundation
**Plan:** 01 of 03
**Status:** In progress
**Progress:** [█░░░░░░░░░] 1/18 requirements (6%)

Last activity: 2026-02-14 — Completed plan 04-01 (regression testing infrastructure)

## Performance Metrics

**Baseline (v0.31.0):**
- Small dataset (200 reads): 6.5s (new benchmark measurement)
- Large dataset (1M reads): ~4.5 min per combination
- Report generation: 93.2% of time (5.1s/6.5s) (new benchmark measurement)
- Interactive HTML size: 9.6 MB (embedded plotly.js)

**Target (v0.32.0):**
- Small dataset: <2s (13x speedup)
- Large batch (50 combinations): 103 min (was 223 min, 2.2x speedup)
- Report size: 19 KB (99.8% reduction)

## Accumulated Context

### Decisions

- Performance analysis completed (PERFORMANCE_ANALYSIS.md) — profiled all pipeline steps
- Kaleido v1.2.0 is the root cause of report generation bottleneck (11s for 2 PNGs)
- Streaming BAM comparison (v0.31.0) already reduced memory from O(n) to O(1)
- Four-phase roadmap prioritizes impact-per-effort: Foundation → Report → Alignment → Comparison
- Keep Kaleido v1.2.0 (not downgrade to 0.2.1) but make PNG export opt-in
- Use directory mode plotly.js with offline fallback for air-gapped environments
- ThreadPoolExecutor for parallel alignment (not ProcessPoolExecutor)
- samtools collate for name grouping (30-50% faster than sort -n)
- **NEW (04-01):** Use atomic write pattern (tempfile + os.replace) for baseline caching to prevent corruption on crash
- **NEW (04-01):** Store baselines as JSON for human readability and easy diff inspection
- **NEW (04-01):** Compare verdicts exactly, ratios with ±0.001 tolerance for regression tests

### Todos

- [x] TEST-01: Regression tests (04-01) - Complete
- [ ] TEST-02: Benchmark script (04-02)
- [ ] TEST-03: Air-gapped testing
- [ ] Verify samtools version >=1.9 (collate requirement)
- [ ] Decide on matplotlib style config for visual consistency with Plotly

### Blockers

None currently identified.

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

## Session Continuity

**What we're building:** Performance optimization milestone (v0.32.0)

**What just happened:** Completed plan 04-01 (regression testing infrastructure)

**Last session:** 2026-02-14 07:06:06
**Stopped at:** Completed 04-01-SUMMARY.md
**Resume file:** None

**Next step:** Execute plan 04-02 (benchmark script)

**Key context for next session:**
- Regression test script complete and validated end-to-end
- Baseline captures v0.31.0 contamination ratios, verdicts, and read assignments
- All Phase 5-7 optimizations must pass: `python scripts/regression_test.py`
- Script uses atomic write pattern (tempfile + os.replace) to prevent baseline corruption

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-14 after plan 04-01 completion*
