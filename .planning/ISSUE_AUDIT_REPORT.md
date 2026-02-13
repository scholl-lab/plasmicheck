# PlasmiCheck Open Issues Audit Report

**Date:** 2026-02-13 (updated after Phase 3 implementation)
**Codebase version:** v0.31.0 (branch `feat/phase3-v0.31.0`)
**Branch:** `feat/phase3-v0.31.0`
**Total open issues:** 30 (18 closed total: 6 during audit + 6 in PR #91 + 6 in Phase 3)

---

## Executive Summary

Of the original 48 open issues, **18 have been closed** (6 already resolved during audit + 6 fixed in PR #91 / v0.30.0 + 6 implemented in Phase 3 / v0.31.0), **3 were updated** with status comments reflecting partial resolution, and **30 remain open and valid**. The remaining issues span bugs (2), features (24), and documentation (4). Below is a full categorization with priority ratings and actionability assessment.

---

## 1. Issues Closed During This Audit

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

## 2. Issues Updated With Status Comments (Partially Resolved)

These issues were updated on 2026-02-13 with comments clarifying remaining work.

| # | Title | What's Done | What Remains | Comment Link |
|---|-------|------------|--------------|--------------|
| **#31** | Feat: Automated versioning with each commit | Manual semver works. GitHub release v0.28.0 created. Basic semver (#4) closed. | No automated bump on commit (no `semantic-release` or CI release workflow). Issue now specifically tracks automation. | [comment](https://github.com/scholl-lab/plasmicheck/issues/31#issuecomment-3898193949) |
| **#78** | Bug: Jobs not finished in snakemake | Snakemake workflow exists in `snakemake/plasmicheck.smk` and is functional. Race condition (#86) now fixed in v0.30.0. | Re-test with v0.30.0 to confirm fix. If still failing, investigate other causes. | [comment](https://github.com/scholl-lab/plasmicheck/issues/78#issuecomment-3898194259) |

> **Note:** #81 was previously in this section but has been fully resolved in PR #91 (v0.30.0) — moved to "Issues Closed in PR #91" below.

---

## 3. Active Bugs (Open, Valid, Should Fix)

Sorted by priority (Critical > High > Medium > Low).

### HIGH

| # | Title | Impact | Root Cause | Fix Complexity |
|---|-------|--------|------------|----------------|
| **#75** | Some xDNA files don't work | Users can't process certain plasmid files | Unknown — needs debugging with specific failing xDNA files | **Unknown** — Needs reproduction data |

### MEDIUM

| # | Title | Impact | Fix Complexity |
|---|-------|--------|----------------|
| **#78** | Snakemake jobs randomly fail | Pipeline unreliable under parallel execution | Race condition (#86) fixed in v0.30.0 — needs re-testing to confirm resolution |

---

## 4. Code Quality & Refactoring — Closed in Phase 3

Both code quality issues from Tier 3 have been resolved in Phase 3 (v0.31.0):

| # | Title | Resolution |
|---|-------|------------|
| **#89** | Improve paired-end file handling | **CLOSED** — Replaced comma-separated convention with explicit `-sf1`/`-sf2` arguments. Legacy `-sf` removed entirely (alpha stage — no backward compatibility needed). `SequencingInput` dataclass normalizes all input styles. |
| **#85** | Optimize memory in compare_alignments.py | **CLOSED** — Replaced dict-based approach with streaming name-sorted BAM merge (`_iter_reads_by_name` + two-pointer merge). Memory usage now O(1) regardless of BAM size. |

---

## 5. Feature Requests — Infrastructure & DevOps

| # | Title | Priority | Status | Notes |
|---|-------|----------|--------|-------|
| **#33** | Docker containerization | **HIGH** | **CLOSED (Phase 3)** | Multi-stage Dockerfile + `.dockerignore` + CI workflow. Builds minimap2/samtools from source. |
| **#32** | Dry-run mode | **MEDIUM** | **CLOSED (Phase 3)** | `--dry-run` / `-n` flag. Plan-execute refactor: `build_plan()` + `print_plan()` + execute. |
| **#29** | Progress bar for long operations | **MEDIUM** | **CLOSED (Phase 3)** | Rich progress bar wrapping pipeline steps. `--no-progress` flag + auto-disable on non-TTY. |
| **#31** | Automated versioning with each commit | **LOW** | Open | Manual semver works. Automation (e.g., `python-semantic-release` in CI) is nice-to-have. |

---

## 6. Feature Requests — Scientific / Analytical

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

## 7. Feature Requests — Reporting & Visualization

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#58** | Integrate CoverageOutsideINSERT and MismatchesNearINSERT in summary report | **HIGH** | These metrics are already computed but not shown in summary reports. Low-hanging fruit. | Low |
| **#51** | Additional plots for summary report | **MEDIUM** | 4 new plot types requested (verdict bar by plasmid, by sample, violin, facet scatter) | Medium |
| **#61** | Scatter plot on top of boxplot | **LOW** | Minor visualization improvement. Plotly supports this natively with `stripmode`. | Low |
| **#53** | Favicon for HTML output | **LOW** | Cosmetic. One-line `<link>` tag in template. | Trivial |
| **#6** | Alignment screenshots in report | **LOW** | Would require IGV headless rendering or a custom viz. Complex for small gain. | High |
| **#26** | Custom report template | **LOW** | Jinja2 already used. Just needs a `--template` CLI arg. | Low |

---

## 8. Feature Requests — UX & CLI

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#80** | Structure output folders with subfolders | **MEDIUM** | Current flat structure in output subfolders mixes BAMs, TSVs, reports, plots. | Low-Medium |
| **#49** | Generate IGV sessions | **LOW** | Auto-generate IGV XML session files. Nice for power users. | Low |
| **#25** | Interactive CLI mode | **LOW** | TUI wizard. Nice-to-have, not blocking. | Medium |
| **#5** | Anonymize BAM outputs | **LOW** | Privacy feature. Read name replacement in output BAMs. | Low |

---

## 9. Feature Requests — Testing & Data

| # | Title | Priority | Status | Notes |
|---|-------|----------|--------|-------|
| **#71** | Test dataset creation | **HIGH** | **CLOSED (Phase 3)** | Synthetic FASTA/FASTQ/GenBank test data (~50 KB) in `tests/data/synthetic/`. Integration tests in `test_integration.py`. 125 total tests (unit + integration). |
| **#72** | Benchmarking dataset creation | **MEDIUM** | Open | Standardized contamination-level datasets for performance validation. |
| **#24** | Synthetic spiked-in validation data | **MEDIUM** | Open | Ground-truth data at known contamination levels. Overlaps with #72. |

---

## 10. Feature Requests — Caching & Performance

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#39** | Reuse indices based on md5sum | **MEDIUM** | Avoid redundant index creation for same input files. Pipeline already computes md5sums. | Medium |
| **#43** | Reuse spliced alignment based on md5sum | **MEDIUM** | Extension of #39 to cache expensive spliced alignment step. Depends on #39. | Medium |

---

## 11. Documentation

| # | Title | Priority | Rationale | Complexity |
|---|-------|----------|-----------|------------|
| **#35** | Describe tool logic visually with graphs | **MEDIUM** | Pipeline flow diagrams, scoring logic visualizations. Helps new users and reviewers. | Low-Medium |
| **#34** | Detailed documentation with examples | **MEDIUM** | Comprehensive usage guide with real-world examples. | Medium |

---

## Priority Matrix Summary

### Tier 1 — Closed (6 issues) ~~done~~
~~`#4` `#15` `#28` `#69` `#76` `#77`~~ — All closed on 2026-02-13.

### Tier 2 — Closed (6 issues) ~~done~~
~~`#86` `#63` `#81` `#79` `#87` `#88`~~ — All closed in PR #91 (v0.30.0) on 2026-02-13.

### Tier 3 — Closed (6 issues) ~~done~~
~~`#89` `#85` `#33` `#32` `#29` `#71`~~ — All implemented in Phase 3, v0.31.0 on branch `feat/phase3-v0.31.0`.

### Tier 4 — Important Enhancements (next up)
| Issue | Type | Effort |
|-------|------|--------|
| **#82** | Science: filter ambiguous reads | Medium |
| **#58** | Report: integrate existing metrics | Low |

### Tier 5 — Nice to Have (remaining 21 issues)
Everything else — valid but lower priority.

---

## Recommended Roadmap

### Phase 1: Bug Fixes & Housekeeping (v0.29.0) ~~done~~
1. ~~Close 6 resolved issues~~ Done (2026-02-13)
2. ~~Update #78, #31, #81 with status comments~~ Done (2026-02-13)
3. ~~Fix #86 (race condition — `filelock` library)~~ Done (PR #91)
4. ~~Fix #63 (MD5 dedup)~~ Done (PR #91)
5. ~~Fix #81 (template path — `importlib.resources`)~~ Done (PR #91)
6. ~~Fix #79 ("wahr" locale bug)~~ Done (PR #91)

### Phase 2: Code Quality (v0.30.0) ~~done~~
1. ~~Implement #87 (config singleton)~~ Done (PR #91)
2. ~~Implement #88 (centralize CLI args)~~ Done (PR #91)
3. Implement #89 (paired-end args) — deferred to Phase 3
4. Implement #85 (memory optimization) — deferred to Phase 3

### Phase 3: Enhancements & Infrastructure (v0.31.0) ~~done~~
1. ~~Implement #71 (synthetic test dataset)~~ Done
2. ~~Implement #89 (paired-end args — legacy `-sf` removed entirely)~~ Done
3. ~~Implement #85 (memory optimization — streaming name-sorted merge)~~ Done
4. ~~Implement #32 (dry-run mode — plan-execute refactor)~~ Done
5. ~~Implement #29 (progress bar — Rich)~~ Done
6. ~~Implement #33 (Docker — multi-stage build + CI)~~ Done
7. Bug fix: `.fastq.gz` extension handling in `align_reads.py` and `run_pipeline.py`
8. Real data validation: contaminated ratio 5.955, clean ratio 0.047 — correct verdicts

### Phase 4: Scientific Enhancements (v0.32.0+)
1. Implement #82 (filter ambiguous reads)
2. Implement #64 (resistance gene coverage)
3. Implement #58 (integrate metrics in summary)
4. Implement #65 (coverage in summary)
5. Remaining feature requests per user priority

---

## References

- [Bionitio: Best Practices for Bioinformatics CLI Software](https://academic.oup.com/gigascience/article/8/9/giz109/5572530)
- [PHA4GE: Public Health Pipeline Best Practices](https://github.com/pha4ge/public-health-pipeline-best-practices/blob/main/docs/pipeline-best-practices.md)
- [Python Race Conditions Best Practices 2025](https://medium.com/pythoneers/avoiding-race-conditions-in-python-in-2025-best-practices-for-async-and-threads-4e006579a622)
- [File Locking in Python: fcntl, msvcrt, portalocker](https://runebook.dev/en/docs/python/library/os/os.plock)
- [Bug Triage: How to Organize, Filter, and Prioritize](https://marker.io/blog/bug-triage)
