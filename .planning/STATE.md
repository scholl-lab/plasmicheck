# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.33.0 -- Scientific & Reporting Enhancements

## Current Position

**Phase:** 8 -- Insert-Region-Aware Filtering
**Plan:** Not yet created (run `/gsd:plan-phase 8`)
**Status:** Roadmap complete, ready to plan Phase 8
**Progress:** ░░░░░░░░░░ 0%

Last activity: 2026-02-14 -- Roadmap created for v0.33.0

## Milestones

- v0.29.0 Bug Fixes & Housekeeping (Phase 1) -- shipped 2026-02-13
- v0.30.0 Code Quality (Phase 2) -- shipped 2026-02-13
- v0.31.0 Enhancements & Infrastructure (Phase 3) -- shipped 2026-02-13
- v0.32.0 Performance Optimization (Phases 4-7) -- shipped 2026-02-14
- v0.33.0 Scientific & Reporting Enhancements (Phases 8-11) -- **in progress**

## Phase Overview (v0.33.0)

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 8 | Insert-Region-Aware Filtering | FILT-01..05 | Not Started |
| 9 | Coverage Metrics | COV-01..05 | Not Started |
| 10 | Resistance Gene Detection | RGENE-01..04 | Not Started |
| 11 | Summary Report Integration | REPT-04..10 | Not Started |

## Accumulated Context

### Key Outcomes (from v0.32.0)

- 9.5x faster pipeline for small datasets (5.5s -> 0.577s)
- 1.97x faster alignment on real datasets (115.2s -> 58.4s)
- 170 tests passing (up from 125 in v0.31.0)
- New CLI flags: --static-report, --plotly-mode, --threads, --plot-backend

### Key Context for v0.33.0

- Phase 8 (filtering) must come first: corrects read assignments that all downstream metrics depend on
- compare_alignments.py is the core file for Phases 8-10
- generate_summary_reports.py is the core file for Phase 11
- insert_region already computed from cDNA_positions.txt in existing pipeline
- pysam.count_coverage() for coverage metrics (Phases 9-10)
- Biopython SeqIO for GenBank annotation parsing (Phase 10)
- Existing matplotlib backend at plasmicheck/scripts/plotting/matplotlib_backend.py
- Detailed milestone plan: .planning/v0.33.0-MILESTONE-PLAN.md

### Known Issues

- #75: Some xDNA files don't work (needs reproduction data)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

### Open Blockers

None.

## Session Continuity

**What just happened:** Roadmap created for v0.33.0 with 4 phases (8-11) covering 21 requirements.

**Next step:** `/gsd:plan-phase 8` to create the execution plan for Insert-Region-Aware Filtering.

**Key context for next session:**
- Phase 8 targets issue #82 (ambiguous read filtering)
- 5 requirements: FILT-01 through FILT-05
- Core change: Add insert_region parameter to _streaming_compare() in compare_alignments.py
- Must maintain backward compatibility via filter_backbone_only config toggle
- Milestone plan has detailed implementation guidance: .planning/v0.33.0-MILESTONE-PLAN.md

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-14 after roadmap creation*
