# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.34.0 -- Robustness & Infrastructure

## Current Position

**Phase:** 12 of 16 (Format Hardening)
**Plan:** 2 of 2 complete
**Status:** Phase 12 complete
**Progress:** ██░░░░░░░░ 20% (2/10 plans in v0.34.0)

Last activity: 2026-02-15 — Completed 12-02-PLAN.md

## Milestones

- v0.29.0 Bug Fixes & Housekeeping (Phase 1) -- shipped 2026-02-13
- v0.30.0 Code Quality (Phase 2) -- shipped 2026-02-13
- v0.31.0 Enhancements & Infrastructure (Phase 3) -- shipped 2026-02-13
- v0.32.0 Performance Optimization (Phases 4-7) -- shipped 2026-02-14
- v0.33.0 Scientific & Reporting Enhancements (Phases 8-11) -- shipped 2026-02-15
- v0.34.0 Robustness & Infrastructure (Phases 12-16) -- **in progress (roadmap defined)**

## Phase Overview

| Phase | Name | Requirements | Plans | Status |
|-------|------|--------------|-------|--------|
| 12 | Format Hardening | FMT-01..04 (4) | 2/2 | Complete |
| 13 | Output Restructuring | OUT-01..04 (4) | 0/2 | Not started |
| 14 | Cache Infrastructure | CACHE-01..06 (6) | 0/3 | Not started |
| 15 | Splice Junction Metadata | SJ-01..02 (2) | 0/1 | Not started |
| 16 | Summary Visualization | VIS-01..04 (4) | 0/2 | Not started |

**Total v0.34.0:** 20 requirements, 10 plans (estimated)

## Accumulated Context

### Key Outcomes

**From v0.34.0 (Phase 12 - Complete):**
- Magic byte format detection for SnapGene and xDNA binary formats
- SnapGene .dna file support via BioPython native parser
- Auto-detection for all plasmid formats (GenBank, xDNA, SnapGene, FASTA)
- Record IDs now use sanitized filenames (fixes xDNA broken ID issue)
- Actionable error messages for conversion failures
- Pipeline auto-detects SnapGene .dna files via -pf flag (no user workaround needed)
- Comprehensive test suite: 18 format detection tests, 11 conversion tests
- 321 unit tests passing (up from 297 at v0.33.0)

**From v0.33.0:**
- 5-category read classification eliminates false positives from backbone contamination
- Per-region coverage metrics (mean/median depth, breadth, CV) via pysam
- Resistance gene detection from GenBank annotations (14 bacterial, 7 mammalian, 9 reporter patterns)
- Summary report integration with heatmaps, tables, and extended exports
- 297 tests passing (up from 170 at v0.32.0)
- 21/21 requirements complete across 4 phases

### Known Issues

- #75: Some xDNA files don't work — **partially addressed in 12-01** (actionable errors, record ID fix)
- #78: Snakemake jobs randomly fail (may be fixed by race condition fix in v0.30.0)
- #80: Output directory structure could be more organized (targeted for Phase 13)
- #39, #43: No caching for indexes and spliced alignments (targeted for Phase 14)

### Key Decisions

| ID | Phase | Context | Decision | Rationale |
|----|-------|---------|----------|-----------|
| FMT-D01 | 12-01 | xDNA record IDs unreliable | Override all record.id with sanitize_filename(stem) | Ensures consistent, valid IDs |
| FMT-D02 | 12-01 | BioPython errors cryptic | Wrap parse calls with actionable messages | Users know what to do on failure |
| FMT-D03 | 12-01 | sanitize_filename duplicated | Remove local copy, import from utils | Single source of truth |
| FMT-D04 | 12-01 | --file_type flag required | Add --format flag, auto-detect when omitted | Better UX, backward compatible |
| FMT-D05 | 12-02 | Pipeline format detection | Use detect_format() via inline import in loop | Follows lazy import pattern |
| FMT-D06 | 12-02 | Pipeline test mocking | Mock detect_format at module level | Avoids file I/O in unit tests |

### Open Blockers

None.

## Session Continuity

**Last session:** 2026-02-15
**Stopped at:** Completed 12-02-PLAN.md (Phase 12 complete)
**Resume file:** None

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-15 after completing 12-02-PLAN.md*
