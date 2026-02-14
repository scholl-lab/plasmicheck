# Research Summary: PlasmiCheck Performance Optimization

**Project:** PlasmiCheck v0.31.0 Performance Enhancements
**Research Date:** 2026-02-14
**Researcher:** GSD Research Synthesizer
**Overall Confidence:** HIGH

---

## Executive Summary

PlasmiCheck is a Python bioinformatics tool for detecting plasmid DNA contamination in sequencing data. The performance optimization research reveals that **report generation dominates execution time** (83.6% for small datasets), making it the highest-priority target. The recommended approach uses **existing stack components optimized through configuration changes** rather than major library replacements.

**Key findings:**
- Report generation bottleneck can be eliminated through opt-in PNG export (saves 11s per sample) and Plotly CDN mode (99.8% HTML size reduction)
- Alignment phase can achieve 1.8x speedup through parallel execution using standard library `ThreadPoolExecutor`
- BAM comparison can be accelerated 30-50% by replacing `samtools sort -n` with `samtools collate`
- **No new dependencies required** for core optimizations — all gains come from configuration changes and standard library features

**Key risks:**
- **CDN mode breaks air-gapped environments** (HPC clusters, secure facilities) — requires multi-tier fallback strategy
- **Parallel alignment risks file handle exhaustion** on systems with low `ulimit` values
- **`samtools collate` is not equivalent to `sort -n`** for reads with supplementary alignments — requires careful ordering logic
- **Kaleido downgrade incompatible with Plotly 6.x** — must stay on v1.x with `start_sync_server()` optimization

The research-backed approach delivers **13x speedup for small datasets and 2.2x for large batch processing** with minimal code changes and zero breaking changes to existing workflows.

---

## Key Findings

### From STACK.md

**Core technology recommendations:**

- **Plotly 6.5.2** (upgrade from 6.3.0) — better Kaleido compatibility, latest stable release
- **Kaleido 1.2.0** (keep current version) — use `start_sync_server()` for 57x speedup instead of downgrading to 0.2.1
- **matplotlib 3.10.8** (upgrade from 3.9.1) — optional fallback for static exports, already installed
- **Standard library ThreadPoolExecutor** — no new dependency needed for parallel alignment
- **Standard library os.cpu_count()** — auto-detect CPU cores with cgroup awareness

**Configuration changes (no version changes needed):**

- Plotly: `include_plotlyjs='directory'` mode creates one shared plotly.min.js (99.8% size reduction)
- Kaleido: Add `kaleido.start_sync_server()` initialization (8-10x speedup)
- samtools: Use `collate` instead of `sort -n` (30-50% faster), add `-m 2G` flag (10-19% faster)

**Dependencies to avoid:**

- Kaleido v0.2.1 downgrade (security vulnerabilities, Plotly 6.x incompatibility)
- orca (deprecated)
- Selenium/Playwright (massive overhead for screenshot task)
- lazy_loader for Python <3.11.9 (race condition)

**Stack confidence: HIGH** — All recommendations verified against official documentation.

---

### From FEATURES.md

**Table stakes (expected features users demand):**

- Opt-in report generation (`--no-report` flag standard in MultiQC, nf-core)
- Parallel alignment execution (Nextflow/Snakemake use dataflow parallelism by default)
- samtools threading flags (`-@ threads`, `-m memory` universally used)
- Index reuse across batch (already implemented via filelock in v0.31.0)
- CDN-based HTML reports (MultiQC v1.21+ uses CDN, embedding 4.6 MB plotly.js is obsolete)
- Lazy imports for CLI (PEP 810 makes this first-class in Python 3.15)
- Configurable threading (`--threads` or auto-detection, hardcoded values are inflexible)

**Differentiators (features that set tools apart):**

- `samtools collate` for name grouping (30-50% faster, rarely adopted)
- Batch-optimized Kaleido persistence (`start_sync_server()` fix from GitHub issue)
- Hybrid report mode (interactive HTML by default, opt-in static PNG)
- Report generation time profiling (`--profile` flag showing breakdown)
- Memory-tuned samtools operations (`-m 2G` flag)

**Anti-features (explicitly avoid):**

- Kaleido v0.2.1 pinning (security risk, unmaintained)
- Embedded Plotly.js in every report (9.6 MB files are archaic)
- ProcessPoolExecutor for subprocess-heavy work (ThreadPoolExecutor is lighter for I/O)
- Global sorting when only grouping needed (`sort -n` vs `collate`)
- Mandatory report generation (batch users want results only)
- Hardcoded resource limits (wastes cores on HPC, oversubscribes on laptops)

**Feature confidence: HIGH** — Based on MultiQC, nf-core patterns, Scientific Python SPEC 1.

---

### From ARCHITECTURE.md

**PlasmiCheck uses plan-execute architecture** where `run_pipeline.py` orchestrates sequential steps across plasmid × sample combinations. All 10 proposed optimizations integrate cleanly:

**Integration points:**

1. **Report generation** (generate_report.py) — opt-in PNG via CLI flag, CDN mode via config, lazy imports move pandas/plotly into functions
2. **Alignment execution** (align_reads.py, run_pipeline.py) — ThreadPoolExecutor wrapper for parallel plasmid+human alignment, CPU auto-detection in config.py
3. **BAM comparison** (compare_alignments.py) — replace `_namesort_bam()` with `_collate_bam()`, add `-m 2G` to sort commands
4. **Index management** (run_pipeline.py) — hoist human index out of loop, track indexed files in PipelinePlan

**Build order recommendation:**

- **Phase 1: Report Optimizations** (2-3 hours, 50-90% of gains, lowest risk)
  - Opt-in PNG generation, Plotly CDN mode, lazy imports
- **Phase 2: Alignment Optimizations** (3-4 hours, 1.5-2x speedup for real data)
  - samtools `-m` flag, CPU auto-detection, parallel alignment
- **Phase 3: Comparison Optimization** (1-2 hours, 30-50% speedup)
  - samtools collate
- **Phase 4: Architectural Cleanup** (1-2 hours, marginal gains)
  - Index reuse hoisting

**Skip for now:** matplotlib fallback (too complex), batch Kaleido export (API limitations)

**Performance projections:**

- Small dataset (200 reads): 13.2s → 1.0s = **13x speedup**
- Large batch (5M reads, 50 combos): 223 min → 103 min = **2.2x speedup**

**Architecture confidence: HIGH** — Existing patterns support all optimizations cleanly.

---

### From PITFALLS.md

**Critical pitfalls (cause rewrites/data corruption):**

1. **Kaleido downgrade breaking Plotly 6.x** — v0.2.1 only works with Plotly 5.x, must stay on v1.2.0 with `start_sync_server()`
2. **CDN plotly.js breaks in air-gapped environments** — HPC clusters, secure facilities have no internet; must implement multi-tier fallback (directory/embedded/CDN)
3. **Parallel alignment subprocess file handle exhaustion** — 6 processes × 2 alignments × N combos can spike to 200+ handles; need `ulimit` checks
4. **`samtools collate` non-equivalence to `sort -n`** — collate doesn't order by flags within QNAME group, breaks supplementary alignment logic
5. **Lazy imports breaking type checking** — mypy can't see function-level imports; must use `if TYPE_CHECKING:` pattern

**Moderate pitfalls (cause delays/confusion):**

6. **Performance regression without correctness verification** — optimization changes contamination ratios; need regression suite comparing to v0.30.0
7. **CPU count auto-detection in containers/HPC** — `os.cpu_count()` returns host cores, not container limit; need cgroup-aware detection
8. **matplotlib vs Plotly visual inconsistency** — different default colors/themes; need shared style config

**Minor pitfalls (annoying but fixable):**

9. **CLI flag changes breaking Snakemake workflows** — changing default behavior breaks downstream automation; maintain backward compatibility
10. **Benchmark cold start cache contamination** — import caching, wrong dataset sizes, single-run measurements give misleading speedup numbers

**Pitfall confidence: HIGH** — Based on official docs, GitHub issues, known bioinformatics patterns.

---

## Implications for Roadmap

### Suggested Phase Structure

Based on combined research, the roadmap should prioritize **impact-per-effort ratio** and **risk mitigation**:

#### Phase 0: Foundation (Pre-work, ~2 hours)

**Rationale:** Establish safety nets before making any optimizations.

**Deliverables:**
- Regression test suite comparing outputs to v0.30.0 (contamination ratios, read assignments)
- Benchmark protocol using realistic datasets (100K+ reads, multiple iterations)
- Check samtools version compatibility (collate requires >=1.9)

**Pitfalls addressed:** #6 (correctness regression), #10 (benchmark contamination)

**Research flag:** Standard testing patterns, no additional research needed

---

#### Phase 1: Report Optimization (Quick Wins, ~3 hours)

**Rationale:** Report generation is 83.6% of time for small datasets, 19% for large batches. Highest ROI.

**Features delivered:**
- Opt-in PNG export via `--static-report` flag (saves 11s per sample)
- Plotly CDN mode with offline fallback (`include_plotlyjs='directory'`)
- Lazy imports (pandas, plotly) in report generation scripts

**From FEATURES.md:** CDN HTML reports (table stakes), hybrid report mode (differentiator)

**Pitfalls to avoid:**
- #1: Verify Kaleido 1.2.0 + Plotly 6.5.2 compatibility
- #2: Implement directory/CDN/embedded fallback for air-gapped environments
- #9: Keep backward compatibility (default behavior unchanged)

**Research flag:** No additional research needed — patterns well-documented in Plotly docs

**Testing requirements:**
- Regression tests (outputs identical with/without PNGs)
- Air-gapped environment test (Docker with no network)
- Backward compatibility test (Snakemake workflow still works)

---

#### Phase 2: Alignment Optimization (Moderate Effort, ~4 hours)

**Rationale:** Alignment is dominant phase for real-world datasets (millions of reads). Parallel execution + threading optimization = 1.5-2x speedup.

**Features delivered:**
- samtools sort `-m 2G` flag (10-19% speedup for large BAMs)
- CPU auto-detection with cgroup awareness (scales to 16+ core systems)
- Parallel plasmid + human alignment via ThreadPoolExecutor (1.8x speedup)

**From FEATURES.md:** Parallel alignment execution (table stakes), configurable threading (table stakes)

**Pitfalls to avoid:**
- #3: Check file handle limits before enabling parallel mode (`ulimit -n`)
- #7: Implement cgroup-aware `get_cpu_count()` for Docker/SLURM
- #6: Validate alignment outputs identical to sequential version

**Research flag:** May need research on ThreadPoolExecutor exception handling patterns in pipeline context

**Testing requirements:**
- Low ulimit test (`ulimit -n 256`, verify graceful fallback)
- Docker container test (CPU limit respected)
- SLURM test (`SLURM_CPUS_PER_TASK` respected)
- Regression test (alignment outputs identical)

---

#### Phase 3: Comparison Optimization (~2 hours)

**Rationale:** BAM comparison is 11% of time for real datasets. `collate` optimization is straightforward but requires careful testing.

**Features delivered:**
- samtools collate replacing sort -n (30-50% faster name grouping)
- Explicit ordering of primary vs supplementary alignments

**From FEATURES.md:** samtools collate for name grouping (differentiator)

**Pitfalls to avoid:**
- #4: Handle supplementary alignment ordering explicitly (not implicit in collate)
- #6: Extensive regression testing on real BAMs with chimeric reads

**Research flag:** Needs research on pysam read ordering patterns with collate vs sort -n

**Testing requirements:**
- Regression test with supplementary alignments (outputs identical)
- Performance benchmark (verify 30-50% speedup)
- samtools version check (collate available in >=1.9)

---

#### Phase 4: Architectural Cleanup (Low Priority, ~1 hour)

**Rationale:** Minor cleanup, marginal performance gains. Only pursue if time allows.

**Features delivered:**
- Hoist human reference indexing out of loop
- PipelinePlan tracks indexed files

**From ARCHITECTURE.md:** Index reuse optimization (already mostly implemented via filelock)

**Pitfalls to avoid:** None critical

**Research flag:** No additional research needed

---

### Deferred to Post-v0.31.0

Based on complexity vs benefit analysis:

- **matplotlib fallback** (Phase 7 in ARCHITECTURE.md) — requires reimplementing plot styles, high maintenance burden, Kaleido optimization sufficient
- **Batch Kaleido export** (Phase 10 in ARCHITECTURE.md) — Kaleido API doesn't expose persistent scope cleanly
- **lazy_loader library** — requires Python 3.11.9+, PlasmiCheck supports 3.10+
- **Parallel report generation** — marginal benefit for typical batch sizes

---

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| **Stack** | HIGH | All recommendations from official docs (Plotly, samtools, Python stdlib) |
| **Features** | HIGH | Based on MultiQC, nf-core, Scientific Python SPEC 1 patterns |
| **Architecture** | HIGH | Existing PlasmiCheck architecture supports all optimizations cleanly |
| **Pitfalls** | HIGH | Critical issues verified via GitHub issues, official docs, academic papers |
| **Overall** | HIGH | Convergent evidence across all research files |

**Source quality breakdown:**

- Official documentation: Plotly, samtools, Python docs, PEPs (HIGH confidence)
- GitHub issues: Kaleido #400, samtools discussions (HIGH confidence)
- Community standards: SPEC 1, MultiQC patterns (HIGH confidence)
- Academic papers: Bioinformatics benchmarking, reproducibility studies (MEDIUM-HIGH confidence)
- Blog posts/forums: Kaleido performance, lazy imports (MEDIUM confidence)

**Gaps identified:**

1. **ThreadPoolExecutor exception handling in pipeline context** — STACK.md shows basic pattern, but PITFALLS.md flags file handle exhaustion as critical. May need research during Phase 2 on graceful degradation strategies.

2. **pysam read ordering with collate vs sort -n** — PITFALLS.md identifies supplementary alignment ordering as critical, but ARCHITECTURE.md assumes equivalence. Phase 3 needs research on explicit ordering logic.

3. **Air-gapped environment testing** — PITFALLS.md identifies CDN mode as critical issue for HPC, but no research file provides comprehensive offline fallback testing strategy. Phase 1 needs validation protocol.

All gaps are addressable during phase implementation — no blockers identified.

---

## Research Flags

**Phases that need `/gsd:research-phase` during planning:**

- **Phase 2 (Alignment):** May need research on:
  - ThreadPoolExecutor exception propagation patterns
  - File handle monitoring/limiting strategies
  - SLURM/Docker resource detection best practices

- **Phase 3 (Comparison):** Needs research on:
  - pysam read iteration order guarantees
  - Supplementary alignment flag ordering
  - Equivalence testing for collate vs sort -n

**Phases with well-documented patterns (skip additional research):**

- **Phase 0 (Foundation):** Regression testing is standard practice
- **Phase 1 (Report):** Plotly docs fully cover CDN/directory/embed modes
- **Phase 4 (Cleanup):** Index reuse is straightforward refactoring

---

## Ready for Requirements Definition

**SUMMARY.md complete.** Research synthesis identifies:

✓ Clear technology recommendations (Plotly 6.5.2, Kaleido 1.2.0, stdlib features)
✓ Four-phase roadmap structure (Foundation → Report → Alignment → Comparison)
✓ Critical pitfalls mapped to phases with prevention strategies
✓ Confidence assessment with identified gaps
✓ Research flags for phases needing deeper investigation

**Key recommendation for orchestrator:**

Proceed with requirements definition using the **four-phase structure**. Prioritize Phase 1 (Report Optimization) for immediate 13x speedup on small datasets with minimal risk. Defer matplotlib fallback and lazy_loader to post-v0.31.0.

**Expected outcome:** 13x speedup for small datasets, 2.2x for large batches, 99.8% HTML size reduction, zero breaking changes.

---

## Sources

### Official Documentation
- [Plotly write_html API](https://plotly.github.io/plotly.py-docs/generated/plotly.io.write_html.html)
- [Python concurrent.futures](https://docs.python.org/3/library/concurrent.futures.html)
- [samtools collate manual](http://www.htslib.org/doc/samtools-collate.html)
- [samtools sort manual](http://www.htslib.org/doc/samtools-sort.html)
- [PEP 810 — Explicit lazy imports](https://peps.python.org/pep-0810/)
- [SPEC 1 — Lazy Loading of Submodules](https://scientific-python.org/specs/spec-0001/)

### GitHub Issues & Discussions
- [Kaleido Performance Regression #400](https://github.com/plotly/Kaleido/issues/400)
- [Plotly 5.x vs 6.x compatibility #5241](https://github.com/plotly/plotly.py/issues/5241)
- [samtools collate vs sort -n #2252](https://github.com/samtools/samtools/issues/2252)
- [samtools supplementary alignment ordering #2010](https://github.com/samtools/samtools/issues/2010)
- [Python os.cpu_count() cgroup issue #36054](https://bugs.python.org/issue36054)

### Community Standards & Guides
- [MultiQC: A fresh coat of paint](https://seqera.io/blog/multiqc-plotly/)
- [Nextflow vs Snakemake 2026](https://www.tasrieit.com/blog/nextflow-vs-snakemake-2026)
- [ThreadPoolExecutor Best Practices](https://superfastpython.com/threadpoolexecutor-best-practices/)
- [Three times faster with lazy imports](https://hugovk.dev/blog/2025/lazy-imports/)

### Academic Papers & Best Practices
- [Bioinformatics benchmarking guidelines](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8)
- [Computational reproducibility five pillars](https://pmc.ncbi.nlm.nih.gov/articles/PMC10591307/)
- [Performance Analysis and Optimization of SAMtools Sorting](https://link.springer.com/chapter/10.1007/978-3-319-58943-5_33)
- [Bionitio: command-line bioinformatics best practices](https://academic.oup.com/gigascience/article/8/9/giz109/5572530)

### Environment-Specific Resources
- [Plotly offline in air-gapped environments](https://foongminwong.medium.com/plotting-data-with-plotly-offline-mode-in-an-air-gapped-environment-5844df874537)
- [SLURM cgroups documentation](https://slurm.schedmd.com/cgroups.html)
- [Container CPU detection](https://github.com/agile6v/container_cpu_detection)
