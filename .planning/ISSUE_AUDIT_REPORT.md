# PlasmiCheck Open Issues Audit Report

**Date:** 2026-02-14 (updated after v0.32.0 milestone)
**Codebase version:** v0.32.0 (branch `main`)
**Total open issues:** 29 (24 closed total: 6 during audit + 6 in PR #91 + 6 in PR #92)

---

## Executive Summary

Of the original 48 open issues, **24 have been closed** (6 already resolved during audit + 6 fixed in PR #91 / v0.30.0 + 6 implemented in Phase 3 / v0.31.0 and merged in PR #92), and **29 remain open and valid**. v0.32.0 (Performance Optimization, Phases 4-7) did not close any GitHub issues but delivered significant internal improvements. The remaining issues span bugs (2), features (23), and documentation (4). Below is a full categorization with priority ratings and actionability assessment.

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

> **Note on #15:** Closed with comment that granular overwrite levels (if still desired) should be tracked in a new, more specific issue.

---

## 1b. Issues Closed in PR #91 (v0.30.0)

These 6 Tier 2 issues were fixed in [PR #91](https://github.com/scholl-lab/plasmicheck/pull/91), merged on 2026-02-13, tagged `v0.30.0`.

| # | Title | Fix Summary |
|---|-------|-------------|
| **#86** | Race condition in create_indexes.py | Added `filelock` dependency; `create_indexes()` uses `FileLock` with fast-path pre-check and partial-index detection |
| **#63** | MD5sum file duplicates | `write_md5sum()` reads existing entries into a set before appending, preventing duplicates |
| **#81** | Template files not found correctly | New `plasmicheck/resources.py` uses `importlib.resources` for reliable path resolution; all template/logo refs updated |
| **#79** | "wahr" in Excel output | `_sanitize_for_excel()` converts boolean columns to strings before Excel export |
| **#87** | Singleton pattern for config management | New `plasmicheck/config.py` lazy-loads `config.json` once; replaced duplicated `json.load()` in 10 source files |
| **#88** | Centralize common CLI arguments | Parent parsers `_logging_parser` / `_threshold_parser` in `cli.py`; shared helpers `add_logging_args()` / `configure_logging_from_args()` in `utils.py` |

---

## 1c. Issues Closed in PR #92 (v0.31.0 + v0.32.0)

These 6 issues were implemented in Phase 3 (v0.31.0) and closed when [PR #92](https://github.com/scholl-lab/plasmicheck/pull/92) merged to main on 2026-02-14.

| # | Title | Fix Summary |
|---|-------|-------------|
| **#89** | Improve paired-end file handling | Replaced comma-separated convention with explicit `-sf1`/`-sf2` arguments. `SequencingInput` dataclass. |
| **#85** | Optimize memory in compare_alignments.py | Streaming name-sorted BAM merge. Memory O(1) regardless of BAM size. |
| **#33** | Docker containerization | Multi-stage Dockerfile + `.dockerignore` + CI workflow. Builds minimap2/samtools from source. |
| **#32** | Dry-run mode | `--dry-run` / `-n` flag. Plan-execute refactor: `build_plan()` + `print_plan()` + execute. |
| **#29** | Progress bar for long operations | Rich progress bar wrapping pipeline steps. `--no-progress` flag + auto-disable on non-TTY. |
| **#71** | Test dataset creation | Synthetic FASTA/FASTQ/GenBank test data (~50 KB). Integration tests. 125 → 170 tests. |

---

## 1d. v0.32.0 Performance Optimization (Phases 4-7)

v0.32.0 did not close any GitHub issues (performance goals were tracked internally as requirements, not as issues). Key deliverables that affect future issue work:

| Improvement | Relevant Issues |
|-------------|-----------------|
| `PipelinePlan.built_indexes` — within-run index dedup | Partially addresses **#39** (cross-run md5sum caching still open) |
| `scripts/benchmark.py` — per-step timing tool | Partially addresses **#72** (standardized benchmark dataset still needed) |
| `--plot-backend matplotlib` — Kaleido-free PNG | Resolves Kaleido dependency concern (not a tracked issue) |
| Lazy imports, opt-in static reports | Performance improvements (not tracked issues) |

---

## 2. Issues Updated With Status Comments (Partially Resolved)

| # | Title | What's Done | What Remains |
|---|-------|------------|--------------|
| **#31** | Feat: Automated versioning with each commit | Manual semver works. Tags v0.28.0 through v0.32.0 exist. GitHub releases created. | No automated bump on commit (no `semantic-release` or CI release workflow). |
| **#78** | Bug: Jobs not finished in snakemake | Race condition (#86) fixed in v0.30.0. | Re-test with v0.30.0+ to confirm fix. If still failing, investigate other causes. |

---

## 3. Active Bugs (Open, Valid, Should Fix)

### HIGH

| # | Title | Impact | Root Cause | Fix Complexity |
|---|-------|--------|------------|----------------|
| **#75** | Some xDNA files don't work | Users can't process certain plasmid files | Unknown — needs debugging with specific failing xDNA files | **Unknown** — Needs reproduction data |

### MEDIUM

| # | Title | Impact | Fix Complexity |
|---|-------|--------|----------------|
| **#78** | Snakemake jobs randomly fail | Pipeline unreliable under parallel execution | Race condition (#86) fixed in v0.30.0 — needs re-testing to confirm resolution |

---

## 4. Feature Requests — Scientific / Analytical

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#82** | Filter ambiguous reads | **HIGH** | Highly contaminated samples produce elevated scores for unrelated plasmids. Core scientific accuracy issue. | Medium |
| **#64** | Coverage on antibiotic resistance genes | **HIGH** | AmpR/KanR coverage is strong contamination evidence (no human homolog). Scientifically sound enhancement. | Medium |
| **#68** | Normalize by splice junction count | **MEDIUM** | Standardize cross-sample comparisons accounting for junction usage differences. | Medium |
| **#65** | Add coverage to summary report | **MEDIUM** | Coverage depth/uniformity is standard QC metric missing from output. | Low-Medium |
| **#50** | Correct cutoffs for exon count/CDS length | **MEDIUM** | Current thresholds may not be biologically justified. Research needed. | High (research-heavy) |
| **#22** | Synthetic plasmid standard for p-value | **LOW** | Research/exploration phase. Statistical framework for contamination detection. | High |
| **#21** | Control samples for p-value | **LOW** | Related to #22. Requires experimental design decisions. | High |
| **#44** | ML for plasmid/human classification | **LOW** | Exploratory. Current scoring method works. ML may improve edge cases. | Very High |
| **#84** | Extract splice positions (exon boundaries) | **LOW** | Nice-to-have for visualization. Not blocking core functionality. | Medium |
| **#83** | Support other plasmid formats | **LOW** | Currently GenBank + xDNA. SnapGene, FASTA, APE formats could be added. | Medium per format |

---

## 5. Feature Requests — Reporting & Visualization

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#58** | Integrate CoverageOutsideINSERT and MismatchesNearINSERT in summary report | **HIGH** | These metrics are already computed but not shown in summary reports. Low-hanging fruit. | Low |
| **#51** | Additional plots for summary report | **MEDIUM** | 4 new plot types requested (verdict bar by plasmid, by sample, violin, facet scatter) | Medium |
| **#61** | Scatter plot on top of boxplot | **LOW** | Minor visualization improvement. Plotly supports this natively with `stripmode`. | Low |
| **#53** | Favicon for HTML output | **LOW** | Cosmetic. One-line `<link>` tag in template. | Trivial |
| **#6** | Alignment screenshots in report | **LOW** | Would require IGV headless rendering or a custom viz. Complex for small gain. | High |
| **#26** | Custom report template | **LOW** | Jinja2 already used. Just needs a `--template` CLI arg. | Low |

---

## 6. Feature Requests — Infrastructure & DevOps

| # | Title | Priority | Status | Notes |
|---|-------|----------|--------|-------|
| **#31** | Automated versioning with each commit | **LOW** | Open | Manual semver works. Automation (e.g., `python-semantic-release` in CI) is nice-to-have. |

---

## 7. Feature Requests — UX & CLI

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#80** | Structure output folders with subfolders | **MEDIUM** | Current flat structure in output subfolders mixes BAMs, TSVs, reports, plots. | Low-Medium |
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

## 9. Feature Requests — Caching & Performance

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#39** | Reuse indices based on md5sum | **MEDIUM** | Avoid redundant index creation across runs. v0.32.0 added within-run dedup (`PipelinePlan.built_indexes`), but cross-run md5sum caching is still needed. | Medium |
| **#43** | Reuse spliced alignment based on md5sum | **MEDIUM** | Extension of #39 to cache expensive spliced alignment step. Depends on #39. | Medium |

---

## 10. Documentation

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#35** | Describe tool logic visually with graphs | **MEDIUM** | Pipeline flow diagrams, scoring logic visualizations. Helps new users and reviewers. | Low-Medium |
| **#34** | Detailed documentation with examples | **MEDIUM** | Comprehensive usage guide with real-world examples. | Medium |

---

## Priority Matrix Summary

### Closed (24 issues total)

| Tier | Issues | Milestone |
|------|--------|-----------|
| Tier 1 — Audit closures | `#4` `#15` `#28` `#69` `#76` `#77` | v0.29.0 |
| Tier 2 — Bug fixes & refactoring | `#86` `#63` `#81` `#79` `#87` `#88` | v0.30.0 |
| Tier 3 — Enhancements & infrastructure | `#89` `#85` `#33` `#32` `#29` `#71` | v0.31.0 |
| Tier 4 — Performance optimization | (no issues closed, internal requirements) | v0.32.0 |

### Next Up — Tier 5: Scientific & Reporting (v0.33.0+)

| Issue | Type | Effort | Why Next |
|-------|------|--------|----------|
| **#82** | Science: filter ambiguous reads | Medium | Core scientific accuracy — highest impact |
| **#64** | Science: resistance gene coverage | Medium | Strong contamination signal, no human homolog |
| **#58** | Report: integrate existing metrics | Low | Already computed, just not displayed — quick win |
| **#65** | Report: coverage in summary | Low-Medium | Standard QC metric, natural complement to #58 |

### Tier 6 — Nice to Have (remaining 25 issues)

Everything else — valid but lower priority. See sections above for details.

---

## Recommended Roadmap

### Phase 1: Bug Fixes & Housekeeping (v0.29.0) ~~done~~
1. ~~Close 6 resolved issues~~ Done (2026-02-13)
2. ~~Update #78, #31 with status comments~~ Done (2026-02-13)
3. ~~Fix #86 (race condition)~~ Done (PR #91)
4. ~~Fix #63 (MD5 dedup)~~ Done (PR #91)
5. ~~Fix #81 (template path)~~ Done (PR #91)
6. ~~Fix #79 ("wahr" locale bug)~~ Done (PR #91)

### Phase 2: Code Quality (v0.30.0) ~~done~~
1. ~~Implement #87 (config singleton)~~ Done (PR #91)
2. ~~Implement #88 (centralize CLI args)~~ Done (PR #91)

### Phase 3: Enhancements & Infrastructure (v0.31.0) ~~done~~
1. ~~Implement #71 (synthetic test dataset)~~ Done
2. ~~Implement #89 (paired-end args)~~ Done
3. ~~Implement #85 (memory optimization)~~ Done
4. ~~Implement #32 (dry-run mode)~~ Done
5. ~~Implement #29 (progress bar)~~ Done
6. ~~Implement #33 (Docker)~~ Done

### Phases 4-7: Performance Optimization (v0.32.0) ~~done~~
1. ~~Regression testing & benchmarking infrastructure~~ Done
2. ~~Report optimization (opt-in PNG, directory plotly.js, lazy imports)~~ Done
3. ~~Alignment optimization (auto threads, -m 2G, SLURM/cgroup detection)~~ Done
4. ~~Comparison & cleanup (index dedup, batch resilience, matplotlib backend)~~ Done

### Next: Scientific Enhancements (v0.33.0+)
1. Implement #82 (filter ambiguous reads)
2. Implement #64 (resistance gene coverage)
3. Implement #58 (integrate metrics in summary)
4. Implement #65 (coverage in summary)
5. Remaining feature requests per user priority

---

## References

- [Bionitio: Best Practices for Bioinformatics CLI Software](https://academic.oup.com/gigascience/article/8/9/giz109/5572530)
- [PHA4GE: Public Health Pipeline Best Practices](https://github.com/pha4ge/public-health-pipeline-best-practices/blob/main/docs/pipeline-best-practices.md)
