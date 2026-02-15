# PlasmiCheck

## What This Is

PlasmiCheck is a Python bioinformatics tool that detects plasmid DNA contamination in sequencing data. It aligns sequencing reads to both a plasmid-specific human reference and the plasmid sequence, then compares alignment scores to determine if contamination is present. Built for researchers and bioinformaticians working with sequencing data that may be affected by plasmid vector contamination. Optimized for performance with opt-in static reports, auto-threaded alignment, and batch-resilient pipeline execution.

## Core Value

Accurately detect plasmid contamination in sequencing data through comparative alignment scoring — the contamination verdict must be reliable.

## Requirements

### Validated

<!-- Shipped and confirmed valuable. -->

- ✓ **CORE-01**: Full pipeline: convert -> index -> spliced align -> align -> compare -> report — v0.28.0
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
- ✓ **REPT-01**: Opt-in static PNG report generation (--static-report flag) — v0.32.0
- ✓ **REPT-02**: Directory-mode plotly.js for 99.8% HTML size reduction — v0.32.0
- ✓ **REPT-03**: Lazy-loaded report imports (pandas, plotly, jinja2) — v0.32.0
- ✓ **ALGN-01**: Auto CPU detection (SLURM/cgroup/bare metal) with --threads override — v0.32.0
- ✓ **ALGN-02**: samtools sort -m 2G memory flag for 4.8x human alignment speedup — v0.32.0
- ✓ **ARCH-01**: Matplotlib static plot backend (Kaleido-free PNG generation) — v0.32.0
- ✓ **ARCH-02**: Batch-resilient pipeline with hoisted human indexing — v0.32.0

### Active

<!-- Current scope. Building toward these. -->

**Current Milestone: v0.34.0 — Robustness & Infrastructure**

**Goal:** Harden format support, improve output organization, add caching for repeated computations, expand summary visualizations, and extract splice junction metadata.

- [ ] xDNA format hardening + SnapGene `.dna` support (#75)
- [ ] Verdict bar charts in summary reports (#51, plots 1-2)
- [ ] Cross-run index and plasmid caching (#39)
- [ ] Output folder restructuring with `OutputPaths` dataclass (#80)
- [ ] Splice junction count extraction and reporting (#68, Phase 1)

### Out of Scope

<!-- Explicit boundaries. Includes reasoning to prevent re-adding. -->

- Splice junction normalization + calibration (#68, Phases 2-3) — needs threshold recalibration
- Correct cutoffs for exon count/CDS length (#50) — deferred
- ML classification (#44) — exploratory, current scoring works
- Interactive CLI mode (#25) — nice-to-have, not blocking
- Anonymize BAM outputs (#5) — privacy feature, lower priority
- IGV session generation (#49) — power user feature
- Batch Kaleido export with write_images() — marginal gain over current approach
- Parallel report generation — I/O bound, minimal speedup expected
- lazy_loader library integration — requires Python 3.11.9+, current approach works

## Context

- **Current state:** v0.33.0 shipped — 297 tests passing, 21 requirements delivered
- **Performance:** Small dataset 0.577s, real 3M-read alignment 58.4s
- **Key issues targeted:** #75 (xDNA), #51 (plots), #39 (caching), #80 (output structure), #68 (splice junctions)
- **Key bugs remaining:** #78 (Snakemake — may be fixed by race condition fix in v0.30.0)
- **External tools:** minimap2 2.28, samtools 1.17, Python 3.10+, pysam, Biopython
- **v0.34.0 focus:** Format robustness, output organization, caching, visualization, scientific metadata

## Constraints

- **Tech stack**: Python + minimap2 + samtools. No new external tool dependencies.
- **Backward compatibility**: CLI API changes must be backward-compatible (add flags, don't remove/rename existing ones)
- **Test coverage**: All changes must pass existing 170 tests + new tests
- **Kaleido dependency**: Optional (needed only when `--static-report` with `--plot-backend plotly`)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Opt-in PNG instead of removing entirely | Non-interactive reports still valuable for email/PDF workflows | ✓ Good — eliminates 91.7% bottleneck from default path |
| Directory-mode plotly.js as default | Most users view reports with internet access; offline users get fallback | ✓ Good — 99.8% HTML size reduction |
| Sequential alignment (not concurrent) | Avoid memory pressure from concurrent minimap2 processes | ✓ Good — simpler, -m 2G sort memory was the real win |
| samtools sort -n retained (collate rejected) | Collate benchmarked 28-64x slower on filtered BAMs (<100K reads) | ✓ Good — evidence-based decision saved incorrect optimization |
| Matplotlib as alternative PNG backend | Kaleido-free option for air-gapped/minimal environments | ✓ Good — broadens deployment options |
| 5-tier CPU detection chain | SLURM/cgroup/os.cpu_count covers all deployment environments | ✓ Good — transparent with source logging |

---
*Last updated: 2026-02-15 after v0.34.0 milestone start*
