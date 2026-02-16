# PlasmiCheck Open Issues Audit Report

**Date:** 2026-02-16 (updated after v0.34.0 release)
**Codebase version:** v0.34.0 (branch `main`)
**Total open issues:** 16 (37 closed total)

---

## Executive Summary

Of the original 53 issues, **37 have been closed** across 6 milestones plus post-release triage. v0.34.0 (Robustness & Infrastructure, Phases 12-16) closed 6 issues: xDNA format hardening (#75), output restructuring (#80), index caching (#39), spliced alignment caching (#43), additional plasmid formats (#83), and summary visualization plots (#51). Post-release triage closed 3 more: scatter on boxplot (#61, already implemented), custom report template (#26, wontfix), favicon (#53, wontfix). **16 remain open**. The remaining issues span bugs (1), features (11), documentation (2), and partially resolved (2).

---

## 1. Issues Closed During Audit (2026-02-13)

These issues were confirmed resolved and closed on 2026-02-13 with comments referencing the resolving commits.

| # | Title | Resolved In | Status |
|---|-------|-------------|--------|
| **#76** | Bug: pip install -U kaleido & xlsxwriter | PR #90 (`dc32c22`) — added to `pyproject.toml` dependencies | **CLOSED** |
| **#77** | Bug: remove weasyprint completely | Commit `28c8311` — zero references remain | **CLOSED** |
| **#4** | Feat: Add semantic versioning | PR #90 + `f43bd54` — `version.py` + hatchling + tag v0.28.0 | **CLOSED** |
| **#69** | Feature Request: PEP8 Comments | PR #90 (`dc32c22`) — ruff + mypy strict + pre-commit + py.typed | **CLOSED** |
| **#28** | Feat: Implement shifted reference to fix MAPQ0 | Pre-existing — `convert_plasmidfile_to_fasta.py:81-102` | **CLOSED** |
| **#15** | Feat: Add an overwrite argument | Pre-existing — `-w`/`--overwrite` on convert, index, pipeline | **CLOSED** |

---

## 1b. Issues Closed in PR #91 (v0.30.0)

| # | Title | Fix Summary |
|---|-------|-------------|
| **#86** | Race condition in create_indexes.py | Added `filelock` dependency; `create_indexes()` uses `FileLock` with fast-path pre-check and partial-index detection |
| **#63** | MD5sum file duplicates | `write_md5sum()` reads existing entries into a set before appending |
| **#81** | Template files not found correctly | `plasmicheck/resources.py` uses `importlib.resources` for reliable path resolution |
| **#79** | "wahr" in Excel output | `_sanitize_for_excel()` converts boolean columns to strings before Excel export |
| **#87** | Singleton pattern for config management | `plasmicheck/config.py` lazy-loads `config.json` once; replaced duplicated `json.load()` in 10 source files |
| **#88** | Centralize common CLI arguments | Parent parsers in `cli.py`; shared helpers in `utils.py` |

---

## 1c. Issues Closed in PR #92 (v0.31.0 + v0.32.0)

| # | Title | Fix Summary |
|---|-------|-------------|
| **#89** | Improve paired-end file handling | Explicit `-sf1`/`-sf2` arguments. `SequencingInput` dataclass. |
| **#85** | Optimize memory in compare_alignments.py | Streaming name-sorted BAM merge. Memory O(1). |
| **#33** | Docker containerization | Multi-stage Dockerfile + `.dockerignore` + CI workflow. |
| **#32** | Dry-run mode | `--dry-run` / `-n` flag. Plan-execute refactor. |
| **#29** | Progress bar for long operations | Rich progress bar. `--no-progress` flag + auto-disable on non-TTY. |
| **#71** | Test dataset creation | Synthetic FASTA/FASTQ/GenBank test data (~50 KB). Integration tests. |

---

## 1d. Issues Closed in PR #93 (v0.33.0)

| # | Title | Fix Summary |
|---|-------|-------------|
| **#82** | Filter ambiguous reads / backbone-only reads | Insert-region-aware filtering. Reads mapping only to backbone classified as `Backbone_Only`. |
| **#64** | Coverage on antibiotic resistance genes | Resistance gene detection from GenBank annotations, per-gene depth/breadth. |
| **#58** | Integrate CoverageOutsideINSERT and MismatchesNearINSERT in summary report | Full coverage metrics integrated into summary with heatmaps, boxplots, TSV/Excel. |
| **#65** | Add coverage to summary report | Coverage for insert and backbone regions in single-sample and summary reports. |

---

## 1e. Issues Closed in PR #94 (v0.34.0)

These 6 issues were resolved in Phases 12-16 and closed when [PR #94](https://github.com/scholl-lab/plasmicheck/pull/94) merged to main on 2026-02-16, tagged `v0.34.0`.

| # | Title | Fix Summary |
|---|-------|-------------|
| **#75** | Bug: some xDNA files don't work | Phase 12: Magic byte format detection, sanitized filename as record ID, actionable error messages for unsupported xDNA versions |
| **#80** | Structure output folders with subfolders | Phase 13: `OutputPaths` frozen dataclass, subfolder layout (`alignments/`, `references/`, `data/`, `reports/`, `plots/`), depth-agnostic summary reader |
| **#39** | Reuse indices based on md5sum | Phase 14: Content-addressed `CacheManager` at `~/.plasmicheck/cache/`, cross-run index caching with manifest, CLI flags `--cache-dir`/`--no-cache`, `plasmicheck cache info\|clean` subcommand |
| **#43** | Reuse spliced alignment based on md5sum | Phase 14: Cross-run spliced alignment caching (BAM, spliced reference, cDNA positions) in `CacheManager` |
| **#83** | Other plasmid formats | Phase 12: SnapGene (`.dna`) + FASTA (`.fasta`/`.fa`) support via auto-detection with magic bytes |
| **#51** | Additional plots for summary report | Phase 16: Verdict bar charts by plasmid and by sample (Plotly + matplotlib). Violin + facet scatter deferred to v0.35.0+ (VIS-05/VIS-06) |

---

## 2. Issues Updated With Status Comments (Partially Resolved)

| # | Title | What's Done | What Remains |
|---|-------|------------|--------------|
| **#31** | Automated versioning with each commit | Manual semver works. Tags v0.28.0 through v0.34.0. GitHub releases created. | No automated bump (no `semantic-release` or CI release workflow). |
| **#68** | Normalize by splice junction count | Phase 15: Splice junction count + exon count extracted from CIGAR, stored in cDNA_positions.txt and summary.tsv/reports. | Actual normalization method not implemented. Threshold recalibration needed. Deferred to v0.35.0+. |
| **#84** | Extract splice positions (exon boundaries) | Phase 15: Junction count + exon count extracted. | Full exon boundary positions (start/end of each exon) not extracted. Needed for IGV screenshot annotations. |
| **#78** | Snakemake jobs randomly fail | Race condition (#86) fixed in v0.30.0. | Needs re-testing with v0.30.0+ to confirm full resolution. |

---

## 3. Active Bugs (Open, Valid, Should Fix)

### MEDIUM

| # | Title | Impact | Fix Complexity |
|---|-------|--------|----------------|
| **#78** | Snakemake jobs randomly fail | Pipeline unreliable under parallel execution | Race condition (#86) fixed in v0.30.0 — needs re-testing to confirm resolution |

> **Note:** #75 (xDNA files) was the only HIGH bug and is now closed in v0.34.0. No HIGH bugs remain.

---

## 4. Feature Requests — Scientific / Analytical

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#68** | Normalize by splice junction count | **MEDIUM** | Metadata extracted (Phase 15). Normalization method + threshold recalibration still needed. | Medium-High |
| **#50** | Correct cutoffs for exon count/CDS length | **MEDIUM** | Current thresholds may not be biologically justified. Research needed. | High (research-heavy) |
| **#22** | Synthetic plasmid standard for p-value | **LOW** | Research/exploration phase. Statistical framework for contamination detection. | High |
| **#21** | Control samples for p-value | **LOW** | Related to #22. Requires experimental design decisions. | High |
| **#44** | ML for plasmid/human classification | **LOW** | Exploratory. Current scoring method works. ML may improve edge cases. | Very High |
| **#84** | Extract splice positions (exon boundaries) | **LOW** | Junction/exon counts extracted in v0.34.0. Full boundary positions needed for IGV. | Medium |

---

## 5. Feature Requests — Reporting & Visualization

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#6** | Alignment screenshots in report | **LOW** | Would require IGV headless rendering or custom viz. Complex for small gain. | High |

> **Note:** #61 (scatter on boxplot) closed — already implemented via `points="all"`. #53 (favicon) and #26 (custom template) closed as wontfix. Violin + facet scatter plots (from #51) tracked as future VIS-05/VIS-06 requirements.

---

## 6. Feature Requests — Infrastructure & DevOps

| # | Title | Priority | Status | Notes |
|---|-------|----------|--------|-------|
| **#31** | Automated versioning with each commit | **LOW** | Open | Manual semver works. Automation (`python-semantic-release` in CI) is nice-to-have. |

---

## 7. Feature Requests — UX & CLI

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#49** | Generate IGV sessions | **LOW** | Auto-generate IGV XML session files. Nice for power users. | Low |
| **#25** | Interactive CLI mode | **LOW** | TUI wizard. Nice-to-have, not blocking. | Medium |
| **#5** | Anonymize BAM outputs | **LOW** | Privacy feature. Read name replacement in output BAMs. | Low |

---

## 8. Feature Requests — Testing & Data

| # | Title | Priority | Status | Notes |
|---|-------|----------|--------|-------|
| **#72** | Benchmarking dataset creation | **MEDIUM** | Open (partially addressed) | v0.32.0 added `scripts/benchmark.py` (timing tool) but still needs standardized multi-level contamination datasets. |
| **#24** | Synthetic spiked-in validation data | **MEDIUM** | Open | Ground-truth data at known contamination levels. Overlaps with #72. |

---

## 9. Documentation

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#35** | Describe tool logic visually with graphs | **MEDIUM** | Pipeline flow diagrams, scoring logic visualizations. | Low-Medium |
| **#34** | Detailed documentation with examples | **MEDIUM** | Comprehensive usage guide with real-world examples. | Medium |

---

## Priority Matrix Summary

### Closed (34 issues total)

| Tier | Issues | Milestone |
|------|--------|-----------|
| Tier 1 — Audit closures | `#4` `#15` `#28` `#69` `#76` `#77` | v0.29.0 |
| Tier 2 — Bug fixes & refactoring | `#86` `#63` `#81` `#79` `#87` `#88` | v0.30.0 |
| Tier 3 — Enhancements & infrastructure | `#89` `#85` `#33` `#32` `#29` `#71` | v0.31.0 |
| Tier 4 — Performance optimization | (no issues closed, internal requirements) | v0.32.0 |
| Tier 5 — Scientific & reporting | `#82` `#64` `#58` `#65` | v0.33.0 |
| Tier 6 — Robustness & infrastructure | `#75` `#80` `#39` `#43` `#83` `#51` | v0.34.0 |
| Post-release triage | `#61` `#26` `#53` | v0.34.0 (wontfix / already done) |

### Remaining Open: 16 issues

| Category | Count | Issues |
|----------|-------|--------|
| Bugs | 1 | #78 |
| Scientific/Analytical | 6 | #68, #50, #22, #21, #44, #84 |
| Reporting/Visualization | 1 | #6 |
| Infrastructure/DevOps | 1 | #31 |
| UX/CLI | 3 | #49, #25, #5 |
| Testing/Data | 2 | #72, #24 |
| Documentation | 2 | #35, #34 |

### Suggested for v0.35.0+

| Issue | Type | Effort | Why Next |
|-------|------|--------|----------|
| **#68** | Science: normalize by splice junction count | Medium-High | Phase 1 metadata done in v0.34.0, normalization is natural next step |
| **#50** | Science: correct cutoffs for exon count/CDS | High | Research-heavy but important for scientific rigor |
| **#72** | Testing: benchmarking dataset | Medium | Standardized datasets needed for validation |
| **#24** | Testing: spiked-in validation data | Medium | Ground-truth data, overlaps with #72 |
| **#34** | Docs: detailed documentation | Medium | Growing user base needs proper docs |
| **#35** | Docs: visual tool logic | Low-Medium | Helps onboarding and manuscript figures |

### Nice to Have (remaining 10 issues)

Everything else — valid but lower priority. See sections above for details.

---

## Completed Roadmap

### Phase 1: Bug Fixes & Housekeeping (v0.29.0) ~~done~~
### Phase 2: Code Quality (v0.30.0) ~~done~~
### Phase 3: Enhancements & Infrastructure (v0.31.0) ~~done~~
### Phases 4-7: Performance Optimization (v0.32.0) ~~done~~
### Phases 8-11: Scientific Enhancements (v0.33.0) ~~done~~
### Phases 12-16: Robustness & Infrastructure (v0.34.0) ~~done~~

---

## References

- [Bionitio: Best Practices for Bioinformatics CLI Software](https://academic.oup.com/gigascience/article/8/9/giz109/5572530)
- [PHA4GE: Public Health Pipeline Best Practices](https://github.com/pha4ge/public-health-pipeline-best-practices/blob/main/docs/pipeline-best-practices.md)
