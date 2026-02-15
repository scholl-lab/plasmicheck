# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.33.0 -- Scientific & Reporting Enhancements

## Current Position

**Phase:** 11 of 11 (Summary Report Integration) — COMPLETE ✓
**Plan:** 2 of 2 complete (11-01 ✓, 11-02 ✓)
**Status:** Phase 11 complete, v0.33.0 milestone ready for deployment
**Progress:** ████████░░ 100% (4/4 phases complete)

Last activity: 2026-02-15 -- Completed 11-02-PLAN.md (Summary report integration: wired all new plots/tables into report flow, added 17 tests)

## Milestones

- v0.29.0 Bug Fixes & Housekeeping (Phase 1) -- shipped 2026-02-13
- v0.30.0 Code Quality (Phase 2) -- shipped 2026-02-13
- v0.31.0 Enhancements & Infrastructure (Phase 3) -- shipped 2026-02-13
- v0.32.0 Performance Optimization (Phases 4-7) -- shipped 2026-02-14
- v0.33.0 Scientific & Reporting Enhancements (Phases 8-11) -- **in progress**

## Phase Overview (v0.33.0)

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 8 | Insert-Region-Aware Filtering | FILT-01..05 | Complete ✓ |
| 9 | Coverage Metrics | COV-01..05 | Complete ✓ |
| 10 | Resistance Gene Detection | RGENE-01..04 | Complete ✓ |
| 11 | Summary Report Integration | REPT-04..10 | Complete ✓ |

## Accumulated Context

### Key Outcomes

**From v0.32.0:**
- 9.5x faster pipeline for small datasets (5.5s -> 0.577s)
- 1.97x faster alignment on real datasets (115.2s -> 58.4s)
- 170 tests passing (up from 125 in v0.31.0)
- New CLI flags: --static-report, --plotly-mode, --threads, --plot-backend

**From Phase 8 (08-01, 08-02, report redesign):**
- 5-category read classification eliminates false positives from backbone contamination
- HTML reports display all 5 categories with visual distinction for excluded categories
- 195 tests passing (up from 170)
- Graceful fallback when cDNA_positions.txt unavailable
- Backward compatibility via filter_backbone_only config toggle
- Consistent color mapping across Plotly and matplotlib backends
- Complete report UI/UX redesign: Bootstrap 5.3 card layout, gradient verdict banner,
  contamination gauge, responsive Plotly charts (white bg, system fonts, consistent
  category ordering), sample/plasmid identity bar, comma-formatted numbers,
  parsed mismatches metrics, collapsible run details, print stylesheet

**From Phase 9 (Coverage Metrics - 09-01, 09-02):**
- Per-region coverage metrics (mean/median depth, breadth, CV) computed via pysam.count_coverage()
- 10 new rows in summary.tsv: MeanDepthInsert, MedianDepthInsert, BreadthInsert, BreadthInsert_5x, CoverageCV_Insert, MeanDepthBackbone, MedianDepthBackbone, BreadthBackbone, BreadthBackbone_5x, CoverageCV_Backbone
- 219 tests passing (up from 195)
- coverage_metrics.py module separates coverage computation logic
- breadth_thresholds configurable in config.json (default: [5])
- Graceful fallback for missing insert region (whole-plasmid metrics)
- Coverage Analysis card in HTML reports displays depth, breadth, and uniformity metrics for insert and backbone
- Fallback warning shown when insert region not defined
- Backward compatible template rendering (missing metrics show 0.00)

**From Phase 10 (Resistance Gene Detection - 10-01, 10-02):**
- resistance_gene_detection.py module with curated pattern matching for 14 bacterial, 7 mammalian, and 9 reporter gene patterns
- GenBank CDS feature parsing via Biopython SeqIO with CompoundLocation handling (origin-crossing features)
- Affinity tag exclusion (His-tag, FLAG, Myc, V5, HA) prevents false positives
- Conservative insert region overlap exclusion ensures backbone-only gene detection
- preserve_genbank parameter in convert_plasmidfile_to_fasta retains .gb files for gene detection
- resistance_gene_patterns extension point in config.json for custom patterns
- Per-gene coverage metrics in summary.tsv: ResGene_{name}_MeanDepth, Breadth, Start, End, Region
- ResistanceGenesAvailable flag distinguishes "no GenBank" from "GenBank but no genes"
- GenBank file preservation enabled by default in pipeline (preserve_genbank=True)
- 280 tests passing (274 from 10-01 + 6 new integration and recall validation tests)
- 90% recall validated against real plasmid GenBank files (pWPI-KCNJ5-WT.gb, pCMV6 files)
- Graceful fallback: gene detection failures don't crash pipeline

### Key Context for v0.33.0

- Phase 8 (filtering) must come first: corrects read assignments that all downstream metrics depend on
- compare_alignments.py is the core file for Phases 8-10
- generate_summary_reports.py is the core file for Phase 11
- insert_region already computed from cDNA_positions.txt in existing pipeline
- pysam.count_coverage() for coverage metrics (Phases 9-10)
- Biopython SeqIO for GenBank annotation parsing (Phase 10)
- Existing matplotlib backend at plasmicheck/scripts/plotting/matplotlib_backend.py
- Detailed milestone plan: .planning/v0.33.0-MILESTONE-PLAN.md

**Phase 11 outcomes (Plans 01 + 02):**
- Six new pure functions for summary report plots (4 in generate_summary_reports.py, 2 in matplotlib_backend.py)
- Generic metric extraction pattern: create_metric_heatmap_data(summary_df, metric_category)
- Resistance gene table flattening: create_resistance_gene_table() parses ResistanceGenes rows
- All functions handle empty DataFrames gracefully (backward compatibility with pre-v0.33.0 data)
- Restructured summary_template.html with 5 collapsible Bootstrap 4 card sections
- All new plots/tables wired into main() → generate_report() flow
- Extended TSV/Excel exports with extended_metrics_table.tsv and resistance_genes_table.tsv
- 297 tests passing (280 from Phase 10 + 17 new unit tests in 11-02)

### Decisions Log

| Phase | Decision | Rationale |
|-------|----------|-----------|
| 08-01 | filter_backbone_only defaults to true | New v0.33.0 behavior provides scientific accuracy by default; false restores pre-v0.33.0 |
| 08-01 | score_margin=0 means disabled | Preserves Tied category; margin > 0 creates Ambiguous category for borderline reads |
| 08-01 | Exact ties never Ambiguous | Tied (score_diff=0) is distinct from Ambiguous (0 < score_diff < margin) |
| 08-01 | Graceful fallback for missing insert region | Warning + pre-v0.33.0 behavior instead of crash enables partial functionality |
| 09-01 | pysam.count_coverage() with quality_threshold=0 | Matches existing alignment behavior (no base quality filtering) |
| 09-01 | Always compute breadth_1x plus configurable thresholds | Breadth ≥1 read is fundamental; additional thresholds configurable via breadth_thresholds |
| 09-01 | CV=0.0 for single-element arrays | Avoids NaN propagation from numpy.std(ddof=1) on single elements |
| 09-01 | CoverageFallback row in summary.tsv | Enables downstream report detection of whole-plasmid fallback mode |
| 09-02 | Always show Coverage Analysis card | Avoids confusion when metrics exist but card is hidden; users see card structure |
| 09-02 | Use Jinja2 default() filter for backward compatibility | Old summary.tsv files render 0.00 instead of error |
| 09-02 | Display fallback warning inline in card | Keeps coverage context together, users immediately understand whole-plasmid metrics |
| 10-01 | Extended bacterial patterns to include full antibiotic names | Broader coverage of GenBank annotation styles (ampicillin, kanamycin vs AmpR, KanR) |
| 10-01 | Curated exclusion list for affinity tags | Prevent false positives from His-tag, FLAG, Myc, V5, HA epitope tags |
| 10-01 | Conservative overlap exclusion for insert region | Any gene overlap (full or partial) with insert is excluded for backbone-only detection |
| 10-01 | Default preserve_genbank=False | Backward compatibility - no changes to existing convert_plasmidfile_to_fasta callers |
| 10-01 | GenBank file location priority (output_dir then plasmid_file) | Always finds most relevant .gb file (from conversion or original) |
| 10-02 | genbank_file parameter optional (None default) | Backward compatibility - existing compare_alignments callers work unchanged |
| 10-02 | Graceful exception handling for gene detection | Failures log warning, continue without gene data - pipeline never crashes |
| 10-02 | Per-gene coverage computed before streaming compare | Requires BAM file access - done early in flow |
| 10-02 | GenBank preservation enabled in pipeline | preserve_genbank=True in run_pipeline - gene detection needs annotations |
| 10-02 | ResistanceGenesAvailable flag | false = no GenBank file, true = GenBank present (may/may not find genes) |
| 11-01 | Generic metric heatmap function | Single plot_metric_heatmap() handles all metric types with blue-to-red color scale and optional vmax |
| 11-01 | depth_heatmap_cap defaults to 100 | Prevents extreme outliers from dominating color scale in depth heatmaps |
| 11-01 | Graceful degradation for missing metrics | Empty DataFrame returns when metric category not found (pre-v0.33.0 compatibility) |
| 11-01 | Matplotlib functions follow existing patterns | sns.set_theme, empty-check, 300 dpi, logger.info for consistency |
| 11-02 | Section ordering: Contamination → Coverage Outside Insert → Coverage Analysis → Mismatch Analysis → Resistance Genes | Per requirement prioritization (REPT-04 to REPT-08) |
| 11-02 | All sections default to expanded | Users see all data immediately without clicking (class="collapse show") |
| 11-02 | Bootstrap 4 syntax (data-toggle/data-target) | Matches existing template's Bootstrap 4.5.2, not Bootstrap 5 |
| 11-02 | Empty plot filenames return '' | Triggers {% if %} guards to skip sections when data unavailable |
| 11-02 | Extended metrics composition | Combine 4 categories into single DataFrame for export consistency |

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

### Open Blockers

None.

## Session Continuity

**Last session:** 2026-02-15
**Stopped at:** Completed Phase 11 Plan 02 (Summary Report Integration) — v0.33.0 milestone complete
**Resume file:** None

**What just happened:**
- Executed Phase 11 Plan 02 (Summary Report Integration) — 2 tasks, 8 minutes
- Task 1: Restructured summary_template.html with 5 collapsible sections, updated generate_report() with 10 new parameters, wired all plot functions in main(), extended save_tables_as_tsv_and_excel() with 2 optional DataFrames
- Task 2: Added 17 comprehensive unit tests (5 test classes: TestCreateMetricHeatmapData, TestCreateResistanceGeneTable, TestPlotMetricHeatmap, TestPlotCoverageOutsideInsertBoxplot, TestSaveTablesExtended)
- All success criteria met: 5 collapsible sections, dual-mode rendering (Plotly/matplotlib), TSV/Excel exports, 297 tests passing
- 3 commits: feat, test, style (ruff formatting)

**Phase 11 complete:**
- Total duration: 13 minutes (Plan 01: 5 min, Plan 02: 8 min)
- 4 tasks, 5 commits (2 feat, 1 test, 1 style, 1 docs from 11-01)
- 3 files modified (summary_template.html, generate_summary_reports.py, test_generate_summary_reports.py)
- No deviations from plan
- All 7 summary report requirements (REPT-04 to REPT-10) satisfied

**v0.33.0 milestone complete:**
- All 4 phases complete (Phases 8-11)
- 297 tests passing (125 → 280 → 297 across milestone)
- Ready for deployment and testing on real datasets

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-15 after Phase 11 Plan 02 execution (v0.33.0 milestone complete)*
