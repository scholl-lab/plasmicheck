# PlasmiCheck — Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.32.0 Performance Optimization

## Current Position

**Phase:** 4 - Foundation
**Plan:** Not yet created
**Status:** Roadmap defined, awaiting Phase 4 planning
**Progress:** [░░░░░░░░░░] 0/18 requirements (0%)

Last activity: 2026-02-14 — Roadmap created for v0.32.0

## Performance Metrics

**Baseline (v0.31.0):**
- Small dataset (200 reads): 13.2s
- Large dataset (1M reads): ~4.5 min per combination
- Report generation: 83.6% of time (11.1s/13.2s)
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
- **NEW:** Four-phase roadmap prioritizes impact-per-effort: Foundation → Report → Alignment → Comparison
- **NEW:** Keep Kaleido v1.2.0 (not downgrade to 0.2.1) but make PNG export opt-in
- **NEW:** Use directory mode plotly.js with offline fallback for air-gapped environments
- **NEW:** ThreadPoolExecutor for parallel alignment (not ProcessPoolExecutor)
- **NEW:** samtools collate for name grouping (30-50% faster than sort -n)

### Todos

- [ ] Plan Phase 4 (Foundation) — regression tests, benchmarks, air-gapped testing
- [ ] Verify samtools version >=1.9 (collate requirement)
- [ ] Decide on matplotlib style config for visual consistency with Plotly

### Blockers

None currently identified.

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

## Session Continuity

**What we're building:** Performance optimization milestone (v0.32.0)

**What just happened:** Roadmap created with 4 phases (4-7) mapping all 18 requirements

**Next step:** Plan Phase 4 (Foundation) — establish regression testing infrastructure

**Key context for next session:**
- Phase numbering starts at 4 (continues from v0.31.0 Phase 3)
- Research identified critical pitfalls: CDN breaks air-gapped, collate != sort -n for supplementary alignments, parallel alignment file handle exhaustion
- All optimizations must pass regression tests (contamination ratios/verdicts unchanged)

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-14 after roadmap creation*
