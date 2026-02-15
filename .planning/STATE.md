# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.34.0 -- Robustness & Infrastructure

## Current Position

**Phase:** Not started (defining requirements)
**Plan:** —
**Status:** Defining requirements
**Progress:** ░░░░░░░░░░ 0%

Last activity: 2026-02-15 — Milestone v0.34.0 started

## Milestones

- v0.29.0 Bug Fixes & Housekeeping (Phase 1) -- shipped 2026-02-13
- v0.30.0 Code Quality (Phase 2) -- shipped 2026-02-13
- v0.31.0 Enhancements & Infrastructure (Phase 3) -- shipped 2026-02-13
- v0.32.0 Performance Optimization (Phases 4-7) -- shipped 2026-02-14
- v0.33.0 Scientific & Reporting Enhancements (Phases 8-11) -- shipped 2026-02-15
- v0.34.0 Robustness & Infrastructure (Phases 12+) -- **defining requirements**

## Accumulated Context

### Key Outcomes

**From v0.33.0:**
- 5-category read classification eliminates false positives from backbone contamination
- Per-region coverage metrics (mean/median depth, breadth, CV) via pysam
- Resistance gene detection from GenBank annotations (14 bacterial, 7 mammalian, 9 reporter patterns)
- Summary report integration with heatmaps, tables, and extended exports
- 297 tests passing (up from 170 at v0.32.0)
- 21/21 requirements complete across 4 phases

### Known Issues

- #75: Some xDNA files don't work (targeted for v0.34.0)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)

### Open Blockers

None.

## Session Continuity

**Last session:** 2026-02-15
**Stopped at:** Defining v0.34.0 requirements
**Resume file:** None

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-15 after v0.34.0 milestone start*
