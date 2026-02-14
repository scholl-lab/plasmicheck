# PlasmiCheck v0.32.0 — Performance Optimization Roadmap

**Milestone:** v0.32.0 Performance Optimization
**Created:** 2026-02-14
**Depth:** Standard (from config)
**Phases:** 4 (Phase 4-7)
**Requirements:** 18 total

## Overview

This roadmap optimizes PlasmiCheck pipeline performance by targeting the three dominant bottlenecks: report generation (83.6% of time for small datasets), alignment execution (dominant for large datasets), and BAM comparison (11% of time). The approach delivers 13x speedup for small datasets and 2.2x for large batches through configuration changes and standard library features, with zero breaking changes.

## Phases

### Phase 4: Foundation

**Goal:** Establish regression testing and benchmarking infrastructure to validate that optimizations preserve correctness.

**Dependencies:** None (foundation phase)

**Plans:** 2 plans

Plans:
- [x] 04-01-PLAN.md -- Regression test script (TEST-01)
- [x] 04-02-PLAN.md -- Performance benchmark script (TEST-02)

**Requirements (delivered in Phase 4):**
- TEST-01: Regression test suite verifying optimization outputs match pre-optimization baseline
- TEST-02: Performance benchmark comparing v0.31.0 vs v0.32.0 on synthetic dataset

**Requirements (deferred from Phase 4 to Phase 5):**
- TEST-03: Air-gapped environment test for directory-mode reports -- deferred because directory-mode reports are implemented in Phase 5; test should be built alongside the feature it validates

**Success Criteria:**
1. Regression test suite compares v0.31.0 contamination ratios, read assignments, and verdicts byte-for-byte
2. Benchmark protocol runs 3 iterations on synthetic dataset and reports mean/std timing per pipeline step
3. ~~Docker container with `--network=none` successfully generates reports in offline mode~~ (deferred to Phase 5 with TEST-03)
4. ~~CI workflow runs regression tests automatically on every optimization PR~~ (not in scope -- local-only tooling per user decision)

### Phase 5: Report Optimization

**Goal:** Eliminate report generation bottleneck (83.6% of time) through opt-in PNG export and HTML size reduction.

**Dependencies:** Phase 4 (regression tests must pass)

**Plans:** 4 plans

Plans:
- [ ] 05-01-PLAN.md -- CLI flags and pipeline wiring (--static-report, --plotly-mode)
- [ ] 05-02-PLAN.md -- Single-sample report optimization (generate_report.py + template)
- [ ] 05-03-PLAN.md -- Summary report optimization (generate_summary_reports.py + template)
- [ ] 05-04-PLAN.md -- Test fixes, CI validation, and regression verification

**Requirements:**
- REPT-01: User can run pipeline without generating static PNG reports
- REPT-02: User can opt into static PNG report generation with `--static-report` CLI flag
- REPT-03: Interactive HTML reports use shared plotly.min.js file (directory mode)
- REPT-04: User can choose plotly.js inclusion mode via CLI flag (cdn, directory, embedded)
- REPT-05: Kaleido uses start_sync_server() initialization for faster PNG export
- REPT-06: Report-related imports are lazy-loaded inside functions
- ~~TEST-03: Air-gapped environment test for directory-mode reports~~ (DROPPED per user decision)

**Success Criteria:**
1. Default pipeline run (no `--static-report`) completes in <2s for small dataset (was 13.2s with Kaleido PNG export)
2. Interactive HTML reports are 19 KB (was 9.6 MB) when using directory mode with shared plotly.min.js
3. User can generate static PNG reports by adding `--static-report` flag, with no change to outputs
4. ~~Air-gapped Docker test passes with directory mode fallback to embedded mode (TEST-03)~~ (DROPPED)
5. CLI startup time reduced by 200-400ms through lazy imports of pandas/plotly/jinja2

### Phase 6: Alignment Optimization

**Goal:** Parallelize alignment execution and optimize threading for 1.5-2x speedup on real-world datasets.

**Dependencies:** Phase 4 (regression tests must pass)

**Requirements:**
- ALGN-01: Plasmid and human alignments run concurrently via ThreadPoolExecutor
- ALGN-02: Pipeline auto-detects CPU count with cgroup/SLURM awareness
- ALGN-03: User can override thread count with `--threads` CLI flag
- ALGN-04: All samtools sort commands use `-m 2G` memory flag

**Success Criteria:**
1. Pipeline runs plasmid and human alignments concurrently (observable via `ps aux` showing 2 minimap2 processes during alignment phase)
2. Auto-detection respects Docker CPU limits (e.g., `--cpus=4` results in 4 threads, not host's 16)
3. User can override auto-detected threads with `--threads 8` flag
4. Large dataset (1M reads) alignment completes 1.5-2x faster than v0.31.0 sequential execution
5. Regression tests pass with identical BAM outputs (same QNAME order, same alignment scores)

### Phase 7: Comparison & Cleanup

**Goal:** Optimize BAM comparison with faster name grouping and eliminate redundant index operations.

**Dependencies:** Phase 4 (regression tests must pass)

**Requirements:**
- COMP-01: BAM name grouping uses samtools collate instead of sort -n
- COMP-02: Supplementary alignment ordering handled explicitly after collate
- ARCH-01: Human reference indexing hoisted out of combination loop
- ARCH-02: PipelinePlan tracks which indexes are already built
- ARCH-03: User can generate static plots via matplotlib backend without Kaleido

**Success Criteria:**
1. BAM comparison completes 30-50% faster than v0.31.0 (samtools collate vs sort -n)
2. Supplementary alignments processed correctly (regression tests pass for chimeric reads)
3. Batch processing (10+ plasmid-sample combinations) skips redundant human index operations
4. User can run `--plot-backend matplotlib --static-report` without Kaleido installed (outputs visually consistent with Plotly)
5. All 125 existing tests pass plus new regression/performance tests from Phase 4

## Progress

| Phase | Requirements | Status | Completion |
|-------|--------------|--------|------------|
| 4 - Foundation | TEST-01, TEST-02 | Complete | 100% |
| 5 - Report Optimization | REPT-01 through REPT-06 | Planned | 0% |
| 6 - Alignment Optimization | ALGN-01 through ALGN-04 | Pending | 0% |
| 7 - Comparison & Cleanup | COMP-01, COMP-02, ARCH-01, ARCH-02, ARCH-03 | Pending | 0% |

**Overall:** 2/18 requirements completed (11%)

## Coverage

All 18 v1 requirements mapped to phases:
- Phase 4: 2 requirements (testing/validation) -- TEST-03 moved to Phase 5 then DROPPED
- Phase 5: 6 requirements (report optimization; TEST-03 DROPPED)
- Phase 6: 4 requirements (alignment optimization)
- Phase 7: 5 requirements (comparison & cleanup)

No orphaned requirements. 1 requirement dropped (TEST-03).

---
*Roadmap created: 2026-02-14*
*Last updated: 2026-02-14 after Phase 5 planning*
