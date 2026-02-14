# Phase 4: Foundation - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Establish regression testing and benchmarking infrastructure to validate that performance optimizations in Phases 5-7 preserve correctness. This is local-only developer tooling for the v0.32.0 optimization effort — no CI changes.

</domain>

<decisions>
## Implementation Decisions

### Regression baseline
- Compare final outputs only: contamination ratios, read assignments, and verdicts
- Small tolerance for floating-point ratios (e.g., +/-0.001), but verdicts must match exactly
- Baseline generated on first run (run v0.31.0 pipeline, cache results locally) — not committed to repo
- Use existing synthetic test fixtures (200 reads) as regression dataset

### Benchmark protocol
- Benchmark both synthetic (200 reads) and real (448K reads) datasets
- Per-step timing: time each pipeline step (convert, index, align, compare, report) separately
- 3 iterations per benchmark run, report mean/std
- Output as Markdown file (BENCHMARK.md) for PR review

### Air-gapped testing
- Deferred to Phase 5 — build air-gapped test when directory-mode reports are actually implemented
- Offline support is "nice to have", not critical for deployment

### Test form
- Standalone scripts (not pytest suite) — regression and benchmark are independent tools, not part of `make test`
- Local-only — no CI workflow changes needed for this phase

### Claude's Discretion
- Exact script structure and CLI interface for regression/benchmark scripts
- How baseline cache is stored and invalidated
- Markdown formatting for benchmark output
- Whether to use `time` module or `timeit` for measurements

</decisions>

<specifics>
## Specific Ideas

- "This is basically a one-time thing to optimize" — tooling should be practical, not over-engineered
- Regression tests validate that optimizations don't change results, not that results are "correct" in absolute terms

</specifics>

<deferred>
## Deferred Ideas

- Air-gapped/Docker testing — deferred to Phase 5 when directory-mode reports are implemented
- CI integration for regression tests — not needed for this local optimization effort

</deferred>

---

*Phase: 04-foundation*
*Context gathered: 2026-02-14*
