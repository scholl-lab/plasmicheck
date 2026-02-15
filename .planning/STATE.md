# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.33.0 -- Scientific & Reporting Enhancements

## Current Position

**Phase:** 10 of 11 (Resistance Gene Detection) — COMPLETE ✓
**Plan:** 2 of 2 complete (10-01 ✓, 10-02 ✓)
**Status:** Phase complete, ready for Phase 11
**Progress:** █████░░░░░ 45% (3/4 phases complete)

Last activity: 2026-02-15 -- Completed 10-02-PLAN.md (Gene detection integration with per-gene coverage, 280 tests passing)

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
| 11 | Summary Report Integration | REPT-04..10 | Not Started |

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

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

### Open Blockers

None.

## Session Continuity

**Last session:** 2026-02-15
**Stopped at:** Completed Phase 10 (Resistance Gene Detection)
**Resume file:** None

**What just happened:**
- Executed Phase 10 Plan 02 (Gene Detection Integration) — 2 tasks, 10 minutes
- Task 1: Integrated gene detection into compare_alignments.py with per-gene coverage; modified run_pipeline.py to preserve GenBank files and forward path
- Task 2: Added 6 integration and recall validation tests (3 for summary.tsv output, 3 for real GenBank file coverage)
- All success criteria met: Per-gene coverage in summary.tsv, GenBank preservation, graceful fallback, backward compatibility
- 280 tests passing (274 from 10-01 + 6 new) with all CI checks passing
- 2 commits: feat, test

**Phase 10 complete:**
- Plan 01: Resistance gene detection module (13 min, 3 commits, 55 tests added)
- Plan 02: Pipeline integration (10 min, 2 commits, 6 tests added)
- Total: 23 minutes, 5 commits, 61 tests added (219 → 280)

**Next step:** Execute Phase 11 (Summary Report Integration) to add resistance gene data to HTML reports.

**Key context for Phase 11:**
- ResGene_ rows in summary.tsv ready for parsing
- ResistanceGenesAvailable flag enables conditional report sections
- Gene data format: ResGene_{sanitized_name}_MeanDepth, Breadth, Start, End, Region
- HTML template at plasmicheck/templates/report_template.html needs new card
- generate_report.py parses summary.tsv and passes to template

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-15 after Phase 10 Plan 02 execution (Phase 10 complete)*
