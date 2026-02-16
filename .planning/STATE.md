# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.34.0 -- Robustness & Infrastructure

## Current Position

**Phase:** 14 of 16 (Cache Infrastructure)
**Plan:** 1 of 3 complete
**Status:** Phase 14 in progress
**Progress:** █████░░░░░ 50% (5/10 plans in v0.34.0)

Last activity: 2026-02-16 — Completed 14-01-PLAN.md

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
| 13 | Output Restructuring | OUT-01..04 (4) | 2/2 | Complete |
| 14 | Cache Infrastructure | CACHE-01..06 (6) | 1/3 | In progress |
| 15 | Splice Junction Metadata | SJ-01..02 (2) | 0/1 | Not started |
| 16 | Summary Visualization | VIS-01..04 (4) | 0/2 | Not started |

**Total v0.34.0:** 20 requirements, 10 plans (estimated)

## Accumulated Context

### Key Outcomes

**From v0.34.0 (Phase 14 - In Progress):**
- CacheManager class with content-addressed storage for pipeline artifacts
- Index cache stores both .mmi and .fai files together (prevents partial index state)
- Spliced alignment cache stores BAM, FASTA, and cDNA positions
- MD5-based content addressing (consistent with existing codebase)
- Atomic writes with tempfile + os.replace (prevents corruption)
- LRU eviction with OrderedDict (O(1) operations)
- Graceful manifest corruption handling
- 373 unit tests passing (up from 353 at Phase 13)

**From v0.34.0 (Phase 13 - Complete):**
- OutputPaths frozen dataclass centralizes all output path construction
- 5-subfolder output layout: alignments/, references/, data/, reports/, plots/
- Type-safe path management with pathlib.Path properties
- Backward-compatible optional parameters for downstream modules
- Depth-agnostic summary report reader supporting both flat and data/ subfolder layouts
- Backward compatibility for existing flat-structure outputs
- Forward compatibility for new 5-subfolder structure
- 353 unit tests passing (up from 321 at Phase 12)

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
| OUT-PATH-01 | 13-01 | Path construction approach | Frozen dataclass with pathlib.Path properties | Type safety, immutability, single source of truth |
| OUT-PATH-02 | 13-01 | comparison_prefix type | Returns str instead of Path | API compatibility with compare_alignments prefix usage |
| OUT-PATH-03 | 13-01 | Backward compatibility | Optional parameters for downstream modules | Zero breaking changes, gradual migration |
| OUT-PATH-04 | 13-01 | GenBank file location | references/ subdirectory | Reference data, matches FASTA/index location |
| OUT-PATH-05 | 13-01 | md5sum.txt location | Base level | Checksums files across all subdirectories |
| SUM-COMPAT-01 | 13-02 | Directory depth detection | Check parent name for 'data' | Simple, reliable, robust to arbitrary nesting |
| SUM-COMPAT-02 | 13-02 | Sample/plasmid extraction | Add _extract_sample_plasmid helper | Centralized logic, eliminates DRY violation, testable |
| CACHE-D01 | 14-01 | Index cache storage | Store both .mmi and .fai files together | create_indexes rebuilds both if only .mmi exists (partial state) |
| CACHE-D02 | 14-01 | Content addressing algorithm | Use MD5 for cache keys | Already used in codebase, sufficient for non-adversarial data |
| CACHE-D03 | 14-01 | Spliced cache key | Combined MD5 of 'plasmid_md5:human_md5' | Depends on both inputs, delimiter prevents collision |
| CACHE-D04 | 14-01 | Atomic writes | tempfile in same directory + os.replace | Same filesystem ensures atomic rename, prevents corruption |
| CACHE-D05 | 14-01 | Manifest corruption | Log warning and recreate empty | Cache is optimization, not critical path |
| CACHE-D06 | 14-01 | LRU eviction | OrderedDict with move_to_end() | O(1) operations, stdlib-based, simpler than custom linked list |

### Open Blockers

None.

## Session Continuity

**Last session:** 2026-02-16
**Stopped at:** Completed 14-01-PLAN.md (CacheManager core implementation)
**Resume file:** None

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-16 after completing 14-01-PLAN.md*
