# PlasmiCheck

## What This Is

PlasmiCheck is a Python bioinformatics tool that detects plasmid DNA contamination in sequencing data. It aligns sequencing reads to both a plasmid-specific human reference and the plasmid sequence, then compares alignment scores to determine if contamination is present. Built for researchers and bioinformaticians working with sequencing data that may be affected by plasmid vector contamination.

## Core Value

Accurately detect plasmid contamination in sequencing data through comparative alignment scoring — the contamination verdict must be reliable.

## Requirements

### Validated

<!-- Shipped and confirmed valuable. -->

- ✓ **CORE-01**: Full pipeline: convert → index → spliced align → align → compare → report — v0.28.0
- ✓ **CORE-02**: Interactive (Plotly) and non-interactive (PNG) HTML reports — v0.28.0
- ✓ **CORE-03**: Multi-sample summary reports with heatmaps, boxplots, p-values — v0.28.0
- ✓ **CORE-04**: GenBank and xDNA plasmid format support — v0.28.0
- ✓ **CORE-05**: BAM and FASTQ input support — v0.28.0
- ✓ **INFRA-01**: Race-condition-safe index creation with filelock — v0.30.0
- ✓ **INFRA-02**: Config singleton (centralized config.json loading) — v0.30.0
- ✓ **INFRA-03**: Centralized CLI arguments with parent parsers — v0.30.0
- ✓ **INFRA-04**: Template/resource path resolution via importlib.resources — v0.30.0
- ✓ **INFRA-05**: MD5sum deduplication — v0.30.0
- ✓ **INFRA-06**: Excel boolean sanitization — v0.30.0
- ✓ **ENH-01**: Explicit -sf1/-sf2 paired-end file handling — v0.31.0
- ✓ **ENH-02**: O(1) memory streaming BAM comparison — v0.31.0
- ✓ **ENH-03**: Dry-run mode with plan-execute architecture — v0.31.0
- ✓ **ENH-04**: Rich progress bar for pipeline execution — v0.31.0
- ✓ **ENH-05**: Docker containerization with CI workflow — v0.31.0
- ✓ **ENH-06**: Synthetic test dataset + integration tests (125 tests) — v0.31.0

### Active

<!-- Current scope. Building toward these. -->

- [ ] Make PNG/static report generation opt-in (`--static-report` flag)
- [ ] Reduce interactive HTML report size (plotly.js CDN/directory mode)
- [ ] Replace `samtools sort -n` with `samtools collate` for faster name grouping
- [ ] Add `-m 2G` to samtools sort commands for large dataset performance
- [ ] Parallel plasmid + human alignment via ThreadPoolExecutor
- [ ] Move heavy imports into function bodies for faster CLI startup
- [ ] Auto-detect CPU count for minimap2/samtools thread settings
- [ ] Optional matplotlib fallback for static plots (eliminate Kaleido dependency)
- [ ] Batch Plotly export optimization
- [ ] Precomputed index registry in PipelinePlan to skip redundant lock checks

### Out of Scope

<!-- Explicit boundaries. Includes reasoning to prevent re-adding. -->

- Scientific enhancements (#82, #64, #68, #65, #50) — deferred to v0.33.0+
- ML classification (#44) — exploratory, current scoring works
- Interactive CLI mode (#25) — nice-to-have, not blocking
- Anonymize BAM outputs (#5) — privacy feature, lower priority
- IGV session generation (#49) — power user feature

## Current Milestone: v0.32.0 Performance Optimization

**Goal:** Optimize pipeline performance — eliminate Kaleido bottleneck, reduce HTML size, parallelize alignments, and improve startup time.

**Target features:**
- Opt-in static report generation (saves 83% of pipeline time for small datasets)
- 99.8% HTML report size reduction via CDN/directory plotly.js
- Parallel alignment execution (~1.8x speedup)
- samtools collate for 30-50% faster name grouping
- Lazy imports for faster CLI startup
- Auto CPU detection for thread tuning
- Optional matplotlib fallback for static plots
- Index reuse optimization for batch processing

## Context

- **Current state:** v0.31.0 on branch `feat/phase3-v0.31.0` — 125 tests passing, all Phase 3 features implemented
- **Performance profiling:** PERFORMANCE_ANALYSIS.md shows 83.6% of pipeline time (11.1s/13.2s) spent on Kaleido PNG export for small datasets. Interactive HTML is 9.6 MB due to embedded plotly.js. 17 subprocess calls per pipeline run.
- **30 open GitHub issues** remaining after 18 closures in Phases 1-3
- **Key bugs:** #75 (xDNA files), #78 (Snakemake — may be fixed by #86 race condition fix)
- **External tools:** minimap2 2.28, samtools 1.17, Python 3.10+, Kaleido 1.2.0, Plotly 6.3.0

## Constraints

- **Tech stack**: Python + minimap2 + samtools. No new external tool dependencies.
- **Backward compatibility**: CLI API changes must be backward-compatible (add flags, don't remove/rename existing ones)
- **Test coverage**: All optimizations must pass existing 125 tests + new performance-related tests
- **Kaleido dependency**: Cannot remove entirely (needed when `--static-report` is used), but can make optional

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Opt-in PNG instead of removing entirely | Non-interactive reports still valuable for email/PDF workflows | — Pending |
| CDN plotly.js as default (not directory) | Most users view reports with internet access | — Pending |
| ThreadPoolExecutor over multiprocessing | Actual work is in subprocesses (minimap2/samtools), not Python | — Pending |
| samtools collate over sort -n removal | Need coord-sorted BAMs for fetch()-based coverage functions | — Pending |

---
*Last updated: 2026-02-14 after milestone v0.32.0 initialization*
