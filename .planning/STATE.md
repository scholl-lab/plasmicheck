# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.33.0 -- Scientific & Reporting Enhancements

## Current Position

**Phase:** 8 of 11 (Insert-Region-Aware Filtering) — VERIFIED ✓
**Plan:** All plans complete, verified
**Status:** Phase 8 verified, ready for Phase 9
**Progress:** ██▓░░░░░░░ 25% (1/4 phases complete)

Last activity: 2026-02-14 -- Phase 8 verified (12/12 must-haves, 5/5 requirements)

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
| 9 | Coverage Metrics | COV-01..05 | Not Started |
| 10 | Resistance Gene Detection | RGENE-01..04 | Not Started |
| 11 | Summary Report Integration | REPT-04..10 | Not Started |

## Accumulated Context

### Key Outcomes

**From v0.32.0:**
- 9.5x faster pipeline for small datasets (5.5s -> 0.577s)
- 1.97x faster alignment on real datasets (115.2s -> 58.4s)
- 170 tests passing (up from 125 in v0.31.0)
- New CLI flags: --static-report, --plotly-mode, --threads, --plot-backend

**From Phase 8 (08-01, 08-02):**
- 5-category read classification eliminates false positives from backbone contamination
- HTML reports display all 5 categories with visual distinction for excluded categories
- 195 tests passing (up from 170)
- Graceful fallback when cDNA_positions.txt unavailable
- Backward compatibility via filter_backbone_only config toggle
- Consistent color mapping across Plotly and matplotlib backends

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

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

### Open Blockers

None.

## Session Continuity

**Last session:** 2026-02-14
**Stopped at:** Phase 8 verified and complete
**Resume file:** None

**What just happened:**
- Executed Phase 8 (Insert-Region-Aware Filtering) — 2 plans across 2 waves
- Plan 08-01: Core filtering logic with 5-category read classification, config toggle, 31 new tests
- Plan 08-02: Report display updates with structured 5-category table, Plotly/matplotlib color mapping
- Verification passed: 12/12 must-haves, 5/5 requirements (FILT-01..05) complete
- 195 tests passing (up from 170)

**Next step:** `/gsd:discuss-phase 9` or `/gsd:plan-phase 9` to plan Coverage Metrics.

**Key context for Phase 9:**
- Insert-region-aware classification established (08-01)
- Report display infrastructure updated (08-02)
- read_overlaps_insert() function available for coverage calculations
- Backbone_Only reads should be excluded from insert coverage metrics
- Requirements: COV-01 through COV-05
- Core files: compare_alignments.py (add coverage functions), plotting modules
- pysam.count_coverage() for per-base depth calculations

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-14 after Phase 8 verification*
