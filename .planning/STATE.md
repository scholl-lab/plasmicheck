# PlasmiCheck -- Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Accurately detect plasmid contamination through comparative alignment scoring
**Current focus:** v0.34.0 -- Robustness & Infrastructure

## Current Position

**Phase:** 15 of 16 (Splice Junction Metadata)
**Plan:** 2 of 2 complete
**Status:** Phase 15 complete
**Progress:** █████████░ 90% (9/10 plans in v0.34.0)

Last activity: 2026-02-16 — Completed 15-02-PLAN.md

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
| 14 | Cache Infrastructure | CACHE-01..06 (6) | 3/3 | Complete |
| 15 | Splice Junction Metadata | SJ-01..02 (2) | 2/2 | Complete |
| 16 | Summary Visualization | VIS-01..04 (4) | 0/2 | Not started |

**Total v0.34.0:** 20 requirements, 10 plans (estimated)

## Accumulated Context

### Key Outcomes

**From v0.34.0 (Phase 15 - Complete):**
- Extended cDNA_positions.txt from 3 to 11 lines with junction metadata
- InsertRegionData dataclass with backward-compatible tuple unpacking (__iter__, __getitem__)
- Junction metadata extraction from primary alignment CIGAR only (skips secondary/supplementary)
- Summary.tsv includes SpliceJunctionCount, ExonCount, ExonLengths, IntronLengths rows
- Junction metadata display in per-sample HTML reports (Splice Junction Metadata card)
- Junction metadata table in multi-sample summary reports (Sample x Plasmid pivot)
- List truncation for long comma-separated values (5 items + "..." in HTML, full in TSV/Excel)
- junction_metadata_table.tsv and Junction Metadata Excel sheet in summary exports
- Key-based field parsing (not line-number based) for forward compatibility
- Graceful handling of unmapped plasmids (all-NA values, no ValueError)
- Backward compatibility: old 3-line cDNA_positions.txt files parse with None defaults
- 397 unit tests passing (up from 384 at Phase 14)

**From v0.34.0 (Phase 14 - Complete):**
- Pipeline integration with --cache-dir and --no-cache CLI flags
- Index caching for plasmid FASTAs with within-run deduplication
- Spliced alignment caching (BAM, FASTA, cDNA positions) integrated
- Cache hits restore both .mmi and .fai files via shutil.copy2
- Graceful degradation on cache errors (log warning, continue without cache)
- CacheManager class with content-addressed storage for pipeline artifacts
- Index cache stores both .mmi and .fai files together (prevents partial index state)
- Spliced alignment cache stores BAM, FASTA, and cDNA positions
- MD5-based content addressing (consistent with existing codebase)
- Atomic writes with tempfile + os.replace (prevents corruption)
- LRU eviction with OrderedDict (O(1) operations)
- Graceful manifest corruption handling
- Cache CLI with 'plasmicheck cache info' and 'plasmicheck cache clean --all' commands
- Cache commands work without minimap2/samtools installation (tool check skipped)
- 384 unit tests passing (up from 353 at Phase 13)

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
| CACHE-D07 | 14-02 | Cache initialization location | After plan build but before dry-run check | Cache path needs to be determined, but dry-run shouldn't initialize cache |
| CACHE-D08 | 14-02 | Within-run deduplication | Via cache mechanism | First combo builds and caches, subsequent combos get cache hits |
| CACHE-D09 | 14-02 | Spliced FASTA indexing | NOT cached separately | Fast operation (~1s), not worth caching complexity |
| CACHE-D10 | 14-02 | Cache error handling | Wrap in try/except (OSError, PermissionError) | Cache errors should never fail pipeline, only degrade performance |
| CACHE-D11 | 14-02 | Cache hit restoration | Use shutil.copy2 instead of symlinks | Ensures files exist at expected paths, no symlink complexity |
| CACHE-D12 | 14-03 | Cache CLI tool requirements | Skip check_requirements() for cache command | Cache commands don't need minimap2/samtools, improves portability |
| SJ-D01 | 15-01 | Junction metadata source | Extract from primary alignment only | Primary = canonical mapping, secondary/supplementary are alternatives |
| SJ-D02 | 15-01 | Unmapped plasmid handling | Write all-NA cDNA_positions.txt instead of raising | Graceful degradation, consistent with Phase 13 fallback behavior |
| SJ-D03 | 15-01 | Backward compatibility | InsertRegionData with __iter__ and __getitem__ | Mimics tuple for unpacking/subscripting, no breaking changes |
| SJ-D04 | 15-01 | Field parsing approach | Parse by prefix, not line number | Forward-compatible if we add more fields, handles missing fields gracefully |
| SJ-D05 | 15-01 | NA value mapping | start/end→0, junction fields→None | start/end=0 maintains existing logic, None signals "not available" for junction metadata |
| SJ-D06 | 15-02 | Junction card display | Hide when SpliceJunctionCount is N/A | Unmapped plasmids have no junction data, avoid showing empty sections |
| SJ-D07 | 15-02 | List truncation threshold | 5 items in HTML, full in exports | Maintains readability while preserving full data for analysis |
| SJ-D08 | 15-02 | Junction export structure | Add to extended_categories + dedicated file | Follows resistance genes pattern, easy access |

### Open Blockers

None.

## Session Continuity

**Last session:** 2026-02-16
**Stopped at:** Completed 15-02-PLAN.md (Splice junction metadata display layer)
**Resume file:** None

---
*State initialized: 2026-02-14*
*Last updated: 2026-02-16 after completing 15-02-PLAN.md — Phase 15 complete*
