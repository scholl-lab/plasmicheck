# Roadmap: PlasmiCheck v0.33.0

**Milestone:** v0.33.0 -- Scientific & Reporting Enhancements
**Phases:** 4 (Phase 8 through Phase 11, continuing from v0.32.0)
**Depth:** Standard
**Coverage:** 21/21 requirements mapped

## Overview

This roadmap delivers multi-dimensional contamination evidence by correcting read assignment accuracy (insert-region-aware filtering), adding quantitative coverage metrics, detecting resistance gene markers, and surfacing all metrics in summary reports. Phases follow a strict dependency chain: filtering must come first because all downstream metrics depend on correct read assignments, coverage metrics establish per-region infrastructure, resistance gene detection builds on that infrastructure, and report integration renders all data sources.

## Phase 8: Insert-Region-Aware Filtering

**Goal:** Reads mapping only to shared backbone are excluded from contamination ratio, eliminating false positive calls for unrelated plasmids.

**Issue:** #82
**Dependencies:** None (first phase)
**Plans:** 2 plans

Plans:
- [ ] 08-01-PLAN.md -- Core filtering logic, config, pipeline wiring, and unit tests
- [ ] 08-02-PLAN.md -- Report display updates for new read categories

**Requirements:**
- FILT-01: Backbone_Only classification and exclusion from ratio
- FILT-02: Optional score margin parameter
- FILT-03: New read categories (Backbone_Only, Ambiguous) in comparison output
- FILT-04: filter_backbone_only config toggle with backward-compatible default
- FILT-05: New categories displayed in single-sample HTML reports

**Success Criteria:**
1. A sample contaminated with Plasmid A no longer shows elevated contamination scores for unrelated Plasmid B that shares only backbone sequence
2. Running with `filter_backbone_only=false` produces identical output to pre-v0.33.0 behavior (all existing tests pass without modification)
3. Comparison TSV output includes Backbone_Only and Ambiguous columns alongside existing Plasmid/Human/Tied counts
4. Single-sample HTML report displays the new read categories with counts and proportions
5. Setting a non-zero `score_margin` causes reads with small score differences to be classified as Ambiguous rather than Plasmid or Human

## Phase 9: Coverage Metrics

**Goal:** Every sample-plasmid pair has quantitative depth, breadth, and uniformity metrics computed per-region (insert vs backbone) and written to summary.tsv.

**Issue:** #65
**Dependencies:** Phase 8 (insert-region concept established, correct read assignments in place)

**Requirements:**
- COV-01: Per-region mean and median depth
- COV-02: Breadth of coverage (fraction >=1 read) per region
- COV-03: Breadth at 5x threshold per region
- COV-04: Coverage uniformity (CV) for backbone
- COV-05: Metrics as additional rows in summary.tsv (backward compatible)

**Success Criteria:**
1. Running the pipeline produces summary.tsv files containing new rows for mean_depth_insert, mean_depth_backbone, median_depth_backbone, breadth_insert, breadth_backbone, breadth_backbone_5x, and coverage_cv_backbone
2. Existing summary.tsv rows (Ratio, Plasmid, Human, Tied, Verdict, CoverageOutsideINSERT, MismatchesNearINSERT) are unchanged -- parsers consuming the old format continue to work
3. Coverage metrics are computed from actual per-base depth (via pysam), not from read count approximations
4. A contaminated sample shows measurable backbone breadth and depth, while a clean sample shows near-zero backbone coverage

## Phase 10: Resistance Gene Detection

**Goal:** Resistance genes are automatically detected from plasmid GenBank annotations and their per-gene coverage is reported as an independent contamination signal.

**Issue:** #64
**Dependencies:** Phase 9 (coverage infrastructure in place)

**Requirements:**
- RGENE-01: Resistance genes extracted from GenBank annotations via pattern matching
- RGENE-02: Per-gene coverage metrics (mean depth + breadth) via pysam
- RGENE-03: Resistance gene metrics in summary.tsv per gene
- RGENE-04: Works with existing plasmid GenBank files (>=90% recall)

**Success Criteria:**
1. Running the pipeline with a GenBank plasmid file containing a resistance gene (e.g., AmpR) produces per-gene mean depth and breadth values in summary.tsv
2. Testing against the 6 existing GenBank files in `reference/plasmid/` detects resistance genes in at least 90% of cases where annotations contain them
3. When no GenBank file is provided or no resistance genes are found, the pipeline completes without error and summary.tsv omits resistance gene rows gracefully
4. Per-gene coverage values are biologically coherent: a contaminated sample shows nonzero resistance gene depth, while a clean sample shows zero

## Phase 11: Summary Report Integration

**Goal:** All metrics -- existing (CoverageOutsideINSERT, MismatchesNearINSERT) and new (coverage depth/breadth, resistance genes) -- are visible in multi-sample summary reports with plots, tables, and exports.

**Issue:** #58 + visualization for #65 and #64
**Dependencies:** Phase 8, Phase 9, Phase 10 (all data sources must exist)

**Requirements:**
- REPT-04: CoverageOutsideINSERT in summary reports (heatmap + boxplot)
- REPT-05: MismatchesNearINSERT in summary reports (chart)
- REPT-06: Coverage depth heatmap (backbone mean depth)
- REPT-07: Coverage breadth heatmap (backbone breadth %)
- REPT-08: Resistance gene coverage table in summary reports
- REPT-09: Matplotlib fallback for all new summary plots
- REPT-10: TSV and Excel exports with new metric tables

**Success Criteria:**
1. Multi-sample summary report contains heatmaps for CoverageOutsideINSERT, backbone mean depth, and backbone breadth, plus a chart for MismatchesNearINSERT and a table for resistance gene coverage
2. All new plots render in both interactive (Plotly) and static (matplotlib) modes without error
3. When summary data lacks new metrics (e.g., reports generated from pre-v0.33.0 data), the report renders without error and omits the new sections gracefully
4. TSV and Excel exports include tables for all new metrics alongside existing tables
5. New plots follow the existing visual style (matching color scales, layout conventions, and dimensions of existing summary report plots)

## Progress

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 8 | Insert-Region-Aware Filtering | FILT-01..05 (5) | Planning Complete |
| 9 | Coverage Metrics | COV-01..05 (5) | Not Started |
| 10 | Resistance Gene Detection | RGENE-01..04 (4) | Not Started |
| 11 | Summary Report Integration | REPT-04..10 (7) | Not Started |

**Total:** 21 requirements across 4 phases. 0/21 complete.

## Coverage Map

| Requirement | Phase | Verified |
|-------------|-------|----------|
| FILT-01 | Phase 8 | Yes |
| FILT-02 | Phase 8 | Yes |
| FILT-03 | Phase 8 | Yes |
| FILT-04 | Phase 8 | Yes |
| FILT-05 | Phase 8 | Yes |
| COV-01 | Phase 9 | Yes |
| COV-02 | Phase 9 | Yes |
| COV-03 | Phase 9 | Yes |
| COV-04 | Phase 9 | Yes |
| COV-05 | Phase 9 | Yes |
| RGENE-01 | Phase 10 | Yes |
| RGENE-02 | Phase 10 | Yes |
| RGENE-03 | Phase 10 | Yes |
| RGENE-04 | Phase 10 | Yes |
| REPT-04 | Phase 11 | Yes |
| REPT-05 | Phase 11 | Yes |
| REPT-06 | Phase 11 | Yes |
| REPT-07 | Phase 11 | Yes |
| REPT-08 | Phase 11 | Yes |
| REPT-09 | Phase 11 | Yes |
| REPT-10 | Phase 11 | Yes |

**Orphaned requirements:** None
**Duplicate mappings:** None
**Coverage:** 21/21 (100%)

---
*Roadmap created: 2026-02-14*
*Last updated: 2026-02-14*
