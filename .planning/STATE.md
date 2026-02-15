# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.33.0 -- Scientific & Reporting Enhancements

## Current Position

**Phase:** 9 of 11 (Coverage Metrics) — COMPLETE ✓
**Plan:** 2 of 2 complete (09-01 ✓, 09-02 ✓)
**Status:** Phase 9 complete, ready for Phase 10
**Progress:** ███▓░░░░░░ 32% (2/4 phases complete)

Last activity: 2026-02-15 -- Completed 09-02-PLAN.md (Coverage metrics now visible in HTML reports, 219 tests passing)

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
| 10 | Resistance Gene Detection | RGENE-01..04 | Not Started |
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

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

### Open Blockers

None.

## Session Continuity

**Last session:** 2026-02-15
**Stopped at:** Completed 09-02-PLAN.md
**Resume file:** None

**What just happened:**
- Executed Phase 9 Plan 02 (Coverage Metrics Display) — 2 tasks, 5 minutes
- Task 1: Updated generate_report.py to parse coverage metrics from summary.tsv and pass to template
- Task 2: Added Coverage Analysis card to report_template.html with table showing insert/backbone metrics
- All success criteria met: Coverage metrics visible in HTML reports
- 219 tests passing (up from 216) with 3 new report-specific coverage tests
- All CI checks passing (lint, format, typecheck, test)
- Phase 9 complete ✓

**Next step:** Execute Phase 10 (Resistance Gene Detection) to parse plasmid GenBank annotations and compute gene-level coverage.

**Key context for Phase 10:**
- coverage_metrics.py infrastructure ready for gene-level coverage computation
- Biopython SeqIO for GenBank annotation parsing (CDS features)
- pysam.count_coverage() can be called for specific gene regions
- compare_alignments.py is the integration point (add gene detection after coverage metrics)
- New columns in summary.tsv for detected resistance genes with coverage stats

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-15 after Phase 9 Plan 02 execution (Phase 9 complete)*
