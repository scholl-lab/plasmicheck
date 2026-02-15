# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.33.0 -- Scientific & Reporting Enhancements

## Current Position

**Phase:** 10 of 11 (Resistance Gene Detection) — IN PROGRESS
**Plan:** 1 of 2 complete (10-01 ✓)
**Status:** Gene detection module ready, integration pending
**Progress:** ████░░░░░░ 38% (2.5/4 phases complete)

Last activity: 2026-02-15 -- Completed 10-01-PLAN.md (Resistance gene detection module with 90% recall, 274 tests passing)

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
| 10 | Resistance Gene Detection | RGENE-01..04 | Plan 01 Complete ✓ |
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

**From Phase 10 (Resistance Gene Detection - 10-01):**
- resistance_gene_detection.py module with curated pattern matching for 14 bacterial, 7 mammalian, and 9 reporter gene patterns
- GenBank CDS feature parsing via Biopython SeqIO with CompoundLocation handling (origin-crossing features)
- Affinity tag exclusion (His-tag, FLAG, Myc, V5, HA) prevents false positives
- Conservative insert region overlap exclusion ensures backbone-only gene detection
- preserve_genbank parameter in convert_plasmidfile_to_fasta retains .gb files for gene detection
- resistance_gene_patterns extension point in config.json for custom patterns
- 274 tests passing (219 existing + 55 new resistance gene detection tests)
- 90% recall validated against real plasmid GenBank files (pWPI-KCNJ5-WT.gb, pCMV6 files)

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

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

### Open Blockers

None.

## Session Continuity

**Last session:** 2026-02-15
**Stopped at:** Completed 10-01-PLAN.md
**Resume file:** None

**What just happened:**
- Executed Phase 10 Plan 01 (Resistance Gene Detection) — 2 tasks, 13 minutes
- Task 1: Created resistance_gene_detection.py module with curated pattern matching for bacterial, mammalian, and reporter genes; modified convert_plasmidfile_to_fasta to preserve GenBank files
- Task 2: Comprehensive unit tests (55 new tests) covering pattern matching, GenBank parsing, insert region exclusion, and real file validation
- All success criteria met: Gene detection module with 90% recall, GenBank preservation, config extension point
- 274 tests passing (219 existing + 55 new) with all CI checks passing
- 3 commits: feat, test, style

**Next step:** Execute Phase 10 Plan 02 (Gene Coverage Integration) to integrate gene detection into compare_alignments.py and compute per-gene coverage metrics.

**Key context for Plan 02:**
- resistance_gene_detection.py module ready for import
- DetectedGene dataclass provides structured data (name, sanitized_name, start, end, region)
- locate_genbank_file() helper finds .gb files in output_dir or original plasmid_file
- coverage_metrics.py get_depth_array() can compute coverage for specific gene regions
- Need to add gene detection call after insert region parsing in compare_alignments.py
- New columns in summary.tsv for each detected gene with coverage stats

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-15 after Phase 10 Plan 01 execution*
