# Feature Landscape: Performance Optimization for Bioinformatics Python Tools

**Domain:** Python bioinformatics pipeline performance optimization
**Researched:** 2026-02-14
**Context:** PlasmiCheck v0.31.0 — adding performance optimizations to existing contamination detection pipeline

---

## Table Stakes

Features users expect in bioinformatics pipeline performance optimizations. Missing = incomplete optimization.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Opt-in report generation** | MultiQC (`--flat`), nf-core pipelines — users expect to skip expensive visualizations for batch processing | Low | `--no-report` or `--static-only` flags standard |
| **Parallel alignment execution** | Nextflow/Snakemake use dataflow parallelism by default; Python tools use `ThreadPoolExecutor` or `ProcessPoolExecutor` | Medium | Independent alignments (plasmid vs human) should run concurrently |
| **samtools threading flags** | Production pipelines universally use `-@ threads` and `-m memory` flags for samtools operations | Low | Doubling threads = ~1.8x speedup for sort/index |
| **Index reuse across batch** | Aligning N samples to same reference: index once, reuse N times — fundamental expectation | Low | Already implemented via filelock in v0.31.0 |
| **Static plot alternatives to Kaleido** | matplotlib is standard backend for non-interactive plots in bioinformatics (FastQC, MultiQC for large batches) | Medium | Kaleido v1.x has 50-250x regression; matplotlib = 0.5s vs Kaleido 11s |
| **CDN-based HTML reports** | MultiQC v1.21+ uses CDN for Plotly to keep reports <1 MB; embedding 4.6 MB plotly.js per report is obsolete | Low | `include_plotlyjs='cdn'` reduces 9.6 MB → 20 KB (99.8%) |
| **Lazy imports for CLI** | Click-based tools, modern Python CLIs — defer heavy imports (pandas, plotly) to subcommands | Low | PEP 810 (Python 3.15) makes this first-class; PlasmiCheck already has lazy CLI imports |
| **Configurable threading** | Users expect `--threads` or auto-detection (`os.cpu_count()`); hardcoded thread counts = inflexible | Low | Current: hardcoded 8 threads in config.json |
| **Progress indication for long operations** | Rich progress bars standard in modern Python tools (existing in v0.31.0) | Low | Already implemented |

---

## Differentiators

Features that set optimized bioinformatics tools apart. Not expected, but highly valued.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **`samtools collate` for name grouping** | 30-50% faster than `samtools sort -n` for paired-read iteration; rarely adopted outside cutting-edge pipelines | Low | PlasmiCheck only needs reads grouped by name, not globally sorted |
| **Report generation time profiling** | Expose where time is spent; users can make informed decisions about report trade-offs | Low | `--profile` flag showing Kaleido vs template vs data loading breakdown |
| **Hybrid report mode** | Interactive HTML by default, opt-in static PNG fallback for archival/PDF conversion | Medium | MultiQC supports both with `--pdf`; PlasmiCheck could offer `--static-report` |
| **Batch-optimized Kaleido persistence** | Kaleido v1.x launches fresh browser per call; calling `kaleido.start_sync_server()` restores v0.2.1 performance (8-10x speedup) | Low | Known fix from [Kaleido #400](https://github.com/plotly/Kaleido/issues/400) |
| **Downsampling for visualization** | For 10M+ reads, plot 5K representative points — matplotlib/Plotly both handle this; prevents bloated HTML | Low | Already implemented in v0.31.0 (`DOWNSAMPLE_LIMIT: 5000`) |
| **Memory-tuned samtools operations** | `-m 2G` flag for samtools sort reduces temp file thrashing on large BAMs (10-19% speedup) | Trivial | Documented in samtools performance studies |
| **Dry-run mode** | Show execution plan without running — valuable for debugging and understanding pipeline | Low | Already implemented in v0.31.0 |
| **Parallel report generation** | For N×M batch: generate reports concurrently after alignments finish | Medium | ThreadPoolExecutor with `max_workers=cpu_count()` |

---

## Anti-Features

Features to explicitly NOT build. Common mistakes or anti-patterns in this domain.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Kaleido v0.2.1 pinning** | Security risk; unmaintained; incompatible with Plotly 6.x+ | Use Kaleido v1.x with `start_sync_server()` OR switch to matplotlib |
| **Embedded Plotly.js in every report** | 9.6 MB HTML files are archaic; MultiQC migrated to CDN in 2024 | `include_plotlyjs='cdn'` (online) or `'directory'` (offline shared .js) |
| **ProcessPoolExecutor for subprocess-heavy work** | Python GIL doesn't affect subprocesses (minimap2/samtools); ThreadPoolExecutor is lighter | Use ThreadPoolExecutor for I/O-bound subprocess orchestration |
| **Global sorting when only grouping needed** | `samtools sort -n` is overkill for paired-read iteration; wastes 30-50% time | Use `samtools collate` for grouping without global ordering |
| **Module-level heavy imports in scripts** | pandas/plotly at top of `generate_report.py` = 2.3s load even for `plasmicheck align` | Move imports inside functions OR use PEP 810 `lazy import` (Python 3.15+) |
| **Mandatory report generation** | Batch users often want results only, view 1-2 reports manually | Make reports opt-in with `--report` OR opt-out with `--no-report` |
| **Hardcoded resource limits** | 8 threads works on developer laptop, wastes 8 cores on HPC node with 64 | Auto-detect `os.cpu_count()` or expose `--threads` flag |
| **Blocking parallel alignments on serial index** | Index human ref → align plasmid → align human = artificial dependency | Pre-index all references, then parallelize alignments |
| **Re-sorting coordinate-sorted BAMs** | Alignment produces coord-sorted + indexed BAM; compare step re-sorts by name = double work | Current PlasmiCheck workflow; mitigated by using `collate` instead of `sort -n` |

---

## Feature Dependencies

```
Performance Optimization Dependency Chain:

Report Generation:
├─ CDN Plotly.js (independent)
├─ Opt-in PNG export (independent)
│  ├─ matplotlib backend (alternative to Kaleido)
│  └─ Kaleido persistence fix (if using Kaleido)
└─ Parallel report generation (requires: alignments complete)

Alignment Optimization:
├─ Parallel plasmid+human alignment (requires: independent indexes)
├─ Auto-detected threading (independent)
└─ samtools memory tuning (independent)

Comparison Optimization:
├─ samtools collate (replaces sort -n, independent)
└─ Parallel comparison (requires: alignments complete)

Startup Optimization:
└─ Function-scoped imports (independent; PEP 810 in Python 3.15)

Batch Processing:
├─ Index reuse (already implemented)
└─ Parallel combination processing (requires: all per-combination optimizations)
```

**Critical path:** Report generation (83.6% of time for small datasets) → Kaleido + CDN optimization must come first.

---

## MVP Recommendation

For performance milestone, prioritize by impact/effort ratio:

### Priority 1: Quick Wins (1-2 hours, 95% of benefit)
1. **CDN Plotly.js** — 99.8% HTML size reduction (9.6 MB → 20 KB)
2. **Opt-in PNG export** — Skip Kaleido by default, add `--static-report` flag (saves 11s per combination)
3. **samtools collate** — Replace `sort -n` with `collate` in comparison (30-50% faster)
4. **samtools memory flag** — Add `-m 2G` to all sort commands (10-19% speedup on large BAMs)

### Priority 2: Moderate Impact (2-3 hours)
5. **Parallel alignment** — ThreadPoolExecutor for concurrent plasmid + human alignment (1.8x alignment phase)
6. **Auto-detect threads** — `os.cpu_count()` for minimap2/samtools (1.5x on 16+ core systems)
7. **Kaleido persistence** — Call `kaleido.start_sync_server()` if PNGs are generated (8-10x speedup)

### Defer to Post-MVP
- matplotlib fallback for static plots (requires reimplementing plot styles)
- Parallel report generation (complexity vs benefit for typical batch sizes)
- Function-scoped imports (marginal gain; wait for PEP 810 in Python 3.15)

---

## Expected Performance Impact

### Scenario: 10 samples × 5 plasmids = 50 combinations

| Optimization | Time Saved | Notes |
|--------------|------------|-------|
| Opt-in PNG export | **42 min** | 50 combos × 50s Kaleido overhead |
| Parallel alignment | **70 min** | 1.8x speedup on alignment phase |
| samtools collate | **7 min** | 30% faster than `sort -n` |
| CDN Plotly.js | 0 min (size only) | 480 MB → 1 MB total HTML output |
| samtools `-m` flag | **5 min** | 10% speedup on sort operations |
| Auto-detect threads (16 core) | **20 min** | 1.5x alignment speedup from 8→16 threads |
| **Total savings** | **~144 min** | 3.7 hours → 1.4 hours = **2.6x speedup** |

### Scenario: 1 sample × 1 plasmid (small dataset)

| Optimization | Time Saved | Notes |
|--------------|------------|-------|
| Opt-in PNG export | **11.1s** | 83% of total 13.2s pipeline time |
| samtools collate | **0.02s** | Negligible for 200 reads |
| **Total savings** | **11.1s** | 13.2s → 2.1s = **6.3x speedup** |

---

## Real-World Tool Comparisons

### MultiQC (Python, HTML reports)
- **Plotting:** Migrated from HighCharts + matplotlib to Plotly (2024)
- **Static plots:** Uses matplotlib with Agg backend for PDF mode
- **HTML size:** Uses CDN by default post-v1.21
- **Flags:** `--pdf` (flat plots), `-m` (module selection), `--no-data-dir` (skip TSV export)
- **Performance:** Switches to static matplotlib for 100+ samples to avoid HTML bloat

### FastQC (Java, HTML reports)
- **Plotting:** Custom Java plotting, embedded PNGs in HTML
- **Report size:** ~500 KB per sample (static images)
- **No interactive plots:** Trades interactivity for reliability and size

### nf-core Pipelines (Nextflow + MultiQC)
- **Parallelism:** Dataflow-based; all independent steps run concurrently
- **MultiQC integration:** Standard final step aggregating all QC outputs
- **Report flags:** Inherited from MultiQC (`--pdf`, module selection)
- **Index reuse:** Nextflow work directory caches all intermediate files automatically

### Snakemake Pipelines (Python workflow manager)
- **Parallelism:** DAG-based; `-j` flag for concurrent rule execution
- **Threading:** Per-rule `threads:` directive; auto-scaled to available cores
- **Batch processing:** Wildcards expand to N samples; index rules run once per reference

---

## Confidence Assessment

| Area | Confidence | Source |
|------|------------|--------|
| Report generation patterns | **HIGH** | MultiQC official docs, blog post on Plotly migration, Kaleido GitHub issues |
| Parallel alignment | **HIGH** | Nextflow/Snakemake docs, Python concurrent.futures best practices |
| samtools optimization | **MEDIUM** | Academic papers (Springer), Biostars discussions, samtools GitHub issues |
| CLI startup optimization | **HIGH** | PEP 810 (accepted), Hugo van Kemenade benchmarks (3x speedup), lazy import articles |
| matplotlib vs Kaleido | **HIGH** | Kaleido issue #400 (50-250x regression), MultiQC migration rationale |
| HTML report size | **HIGH** | MultiQC docs, Plotly documentation on `include_plotlyjs` modes |
| Batch processing patterns | **MEDIUM** | minimap2 man page, nf-core pipeline examples, general bioinformatics practice |
| Opt-in report flags | **MEDIUM** | MultiQC (`--pdf`, module flags), general CLI patterns — not universal standard |

---

## Sources

### Report Generation & Visualization
- [MultiQC: A fresh coat of paint](https://seqera.io/blog/multiqc-plotly/) — Migration from HighCharts + matplotlib to Plotly
- [Running MultiQC | Seqera Docs](https://docs.seqera.io/multiqc/getting_started/running_multiqc) — CLI flags for optimization
- [Performance Regression: v1.0.0 renders a lot slower than v0.2.1 · Issue #400 · plotly/Kaleido](https://github.com/plotly/Kaleido/issues/400) — 50x slowdown, `start_sync_server()` fix
- [Kaleido replacement? - Plotly Community Forum](https://community.plotly.com/t/kaleido-replacement/72599) — 250x regression reports
- [Plotly vs matplotlib: A quick comparison](https://www.fabi.ai/blog/plotly-vs-matplotlib-a-quick-comparison-with-visual-guides)
- [Matplotlib vs Seaborn vs Plotly](https://www.blog.qualitypointtech.com/2026/02/matplotlib-vs-seaborn-vs-plotly.html)

### Pipeline Parallelism
- [Nextflow vs Snakemake 2026: We've Run Both on 100+ Projects](https://www.tasrieit.com/blog/nextflow-vs-snakemake-2026) — Dataflow parallelism patterns
- [Bioinformatics Pipeline Frameworks (2025): Nextflow vs Flyte vs Airflow vs Snakemake](https://www.tracer.cloud/resources/bioinformatics-pipeline-frameworks-2025)
- [nf-core pipelines report generation MultiQC integration](https://nf-co.re/modules/multiqc/)
- [Using MultiQC in pipelines | Seqera Docs](https://docs.seqera.io/multiqc/usage/pipelines/)

### Python CLI Performance
- [PEP 810 – Explicit lazy imports](https://peps.python.org/pep-0810/) — Accepted for Python 3.15
- [Three times faster with lazy imports · Hugo van Kemenade](https://hugovk.dev/blog/2025/lazy-imports/) — 2.92x startup speedup
- [Python Steering Council greenlights "explicit lazy imports" (PEP 810)](https://medium.com/@virtualik/python-steering-council-greenlights-explicit-lazy-imports-pep-810-dc902231bc48)
- [Python lazy imports you can use today | PythonTest](https://pythontest.com/python-lazy-imports-now/)

### samtools Optimization
- [Performance Analysis and Optimization of SAMtools Sorting](https://link.springer.com/chapter/10.1007/978-3-319-58943-5_33) — Academic study
- [samtools multithreads - Biostars](https://www.biostars.org/p/9462175/) — Threading best practices
- [samtools sort number of threads in reading phase · Issue #891](https://github.com/samtools/samtools/issues/891)
- [Efficient And Fastest Way To Sort Large (>100Gb) Bam Files? - Biostars](https://www.biostars.org/p/18933/)

### minimap2 & Batch Processing
- [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html) — Official documentation on indexing and threading
- [GitHub - lh3/minimap2](https://github.com/lh3/minimap2)
- [Piping minimap2 result to samtools sort · Issue #667](https://github.com/lh3/minimap2/issues/667)

### Bioinformatics Best Practices
- [Bionitio: demonstrating and facilitating best practices for bioinformatics command-line software](https://academic.oup.com/gigascience/article/8/9/giz109/5572530)
- [Ten simple rules for getting started with command-line bioinformatics - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7891784/)
- [On best practices in the development of bioinformatics software - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4078907/)

### Python Concurrency
- [Python Parallelization: Deciding Between Subprocess, Multiprocessing, and Threading](https://www.pythontutorials.net/blog/deciding-among-subprocess-multiprocessing-and-thread-in-python/)
- [Multithreading vs. Multiprocessing in Python](https://thenewstack.io/python-threadpool-vs-multiprocessing/)
- [Managing Parallel Processing with Python's ThreadPoolExecutor](https://codezup.com/parallel-processing-with-threadpoolexecutor-python/)

---

## Notes for Roadmap Creation

### Phase Structure Recommendations

1. **Phase 1: Report Optimization** (highest impact/lowest effort)
   - CDN Plotly.js, opt-in PNG export, Kaleido persistence
   - Rationale: 83.6% of time for small datasets; 42 min savings for batch

2. **Phase 2: Alignment & Comparison Optimization**
   - Parallel alignment, samtools collate, memory tuning, auto-threading
   - Rationale: 1.8x speedup on alignment phase (dominant for real-world datasets)

3. **Phase 3: Polish & Observability** (optional)
   - Performance profiling flags, matplotlib fallback, parallel reports
   - Rationale: Nice-to-have; main gains already captured

### Research Flags for Phases

- **Phase 1:** Unlikely to need research — all patterns well-documented
- **Phase 2:** May need research on ThreadPoolExecutor exception handling in pipeline context
- **Phase 3:** May need research on matplotlib plot equivalents for PlasmiCheck's Plotly visuals

### Open Questions

1. **CDN reliability:** Should PlasmiCheck default to CDN (requires internet) or `'directory'` mode (one shared plotly.min.js)?
   - **Recommendation:** `'directory'` for offline support; users can override with config
2. **matplotlib migration scope:** Full replacement or fallback only?
   - **Recommendation:** Fallback only; Plotly interactive HTML is valuable UX
3. **Threading CLI surface:** `--threads` flag, config file, or auto-detect only?
   - **Recommendation:** Auto-detect with `--threads` override; config as fallback
