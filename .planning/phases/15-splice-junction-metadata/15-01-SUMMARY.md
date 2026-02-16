---
phase: 15-splice-junction-metadata
plan: 01
subsystem: data-layer
tags: [pysam, BAM, CIGAR, splice-junctions, dataclass, backward-compatibility]

# Dependency graph
requires:
  - phase: 13-output-restructuring
    provides: OutputPaths structure for file organization
  - phase: 09-insert-backbone-coverage
    provides: parse_insert_region, cDNA_positions.txt format
provides:
  - Extended cDNA_positions.txt with junction metadata (11 fields)
  - InsertRegionData dataclass for structured junction data
  - Junction metadata in summary.tsv (SpliceJunctionCount, ExonCount, ExonLengths, IntronLengths)
  - Backward-compatible parser for old 3-line cDNA_positions.txt files
affects: [16-summary-visualization, future-splice-aware-features]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Dataclass with __iter__ and __getitem__ for backward-compatible tuple unpacking"
    - "Key-based parsing (not line-number based) for forward compatibility"
    - "Primary-only alignment filtering for junction counting"
    - "NA sentinel values for unmapped/missing data"

key-files:
  created: []
  modified:
    - plasmicheck/scripts/spliced_alignment.py
    - plasmicheck/scripts/compare_alignments.py
    - plasmicheck/scripts/coverage_metrics.py
    - plasmicheck/scripts/resistance_gene_detection.py
    - tests/test_spliced_alignment.py
    - tests/test_compare_alignments.py

key-decisions:
  - "Extract junction metadata from primary alignment only (secondary/supplementary skipped)"
  - "Write NA for unmapped plasmids instead of raising ValueError"
  - "Use InsertRegionData dataclass with __iter__/__getitem__ for tuple-like behavior"
  - "Parse by field prefix, not line number, for forward compatibility"
  - "NA values: 0 for start/end (maintains existing logic), None for junction fields"

patterns-established:
  - "Junction metadata extraction: primary-only filtering, CIGAR parsing"
  - "Backward-compatible dataclass: __iter__ for unpacking, __getitem__ for subscripting"
  - "Key-based parsing: field prefix matching for optional fields"
  - "NA handling: different semantics for coordinates (0) vs metadata (None)"

# Metrics
duration: 13min
completed: 2026-02-16
---

# Phase 15 Plan 01: Splice Junction Metadata Summary

**Extract splice junction metadata from BAM CIGAR, store in 11-line cDNA_positions.txt and summary.tsv, with backward-compatible InsertRegionData dataclass**

## Performance

- **Duration:** 13 min 26 sec
- **Started:** 2026-02-16T10:09:08Z
- **Completed:** 2026-02-16T10:22:35Z
- **Tasks:** 2/2
- **Files modified:** 6

## Accomplishments

- Extended cDNA_positions.txt from 3 to 11 lines with junction metadata (SpliceJunctionCount, ExonCount, ExonLengths, IntronLengths, Strand, Chromosome, GenomicStart, GenomicEnd)
- Created InsertRegionData dataclass with `__iter__` and `__getitem__` for backward-compatible tuple unpacking (`start, end = parse_insert_region(file)` still works)
- Primary-only alignment filtering ensures junction counts reflect the main alignment, not supplementary/secondary alignments
- Backward compatibility: old 3-line cDNA_positions.txt files parse with None defaults for new fields
- Summary.tsv includes junction metadata rows for downstream analysis

## Task Commits

Each task was committed atomically:

1. **Task 1: Junction metadata extraction and cDNA_positions.txt extension** - `8f35a18` (feat)
   - Extract junction metadata from primary alignment CIGAR
   - Extend cDNA_positions.txt from 3 to 11 lines
   - Write NA for unmapped plasmids (no ValueError)
   - Add 7 unit tests for junction extraction

2. **Task 2: InsertRegionData dataclass, parser update, and summary.tsv junction fields** - `1d1617e` (feat)
   - Create InsertRegionData with __iter__ and __getitem__
   - Update parse_insert_region to return InsertRegionData
   - Add junction fields to summary.tsv
   - Update type hints across 3 modules
   - Add 9 unit tests for dataclass and summary output

## Files Created/Modified

- **plasmicheck/scripts/spliced_alignment.py** - Extract junction metadata from primary alignment CIGAR, write 11-line cDNA_positions.txt
- **plasmicheck/scripts/compare_alignments.py** - InsertRegionData dataclass, updated parser, junction fields in summary.tsv
- **plasmicheck/scripts/coverage_metrics.py** - Updated type hints for InsertRegionData compatibility
- **plasmicheck/scripts/resistance_gene_detection.py** - Updated type hints for InsertRegionData compatibility
- **tests/test_spliced_alignment.py** - Tests for junction extraction, primary-only filtering, NA handling
- **tests/test_compare_alignments.py** - Tests for InsertRegionData, backward compatibility, summary.tsv junction fields

## Decisions Made

**1. Extract junction metadata from primary alignment only**
- **Context:** BAM files contain primary, secondary, and supplementary alignments. Secondary/supplementary represent alternative mappings or split reads.
- **Decision:** Count N operations (junctions) only from the primary alignment.
- **Rationale:** The primary alignment represents the "main" mapping. Including secondary/supplementary would inflate junction counts and create ambiguity about which alignment's structure is canonical.
- **Implementation:** Check `read.is_secondary` and `read.is_supplementary` before processing CIGAR for junction metadata.

**2. Write NA for unmapped plasmids instead of raising ValueError**
- **Context:** Original extract_plasmid_cdna_positions raised ValueError when no alignments found.
- **Decision:** Write all-NA cDNA_positions.txt instead of raising.
- **Rationale:** Graceful degradation allows pipeline to continue for unmapped plasmids. Downstream consumers (compare_alignments) already handle missing insert_region. Consistent with Phase 13 fallback behavior.
- **Implementation:** Check if cdna_start/cdna_end are None, write NA for all fields.

**3. Use InsertRegionData dataclass with __iter__ and __getitem__ for tuple-like behavior**
- **Context:** parse_insert_region previously returned `tuple[int, int]`. All callers use it as a tuple: `insert_region[0]`, `insert_region[1]`, `start, end = insert_region`.
- **Decision:** Return InsertRegionData dataclass but implement `__iter__` and `__getitem__` to mimic tuple.
- **Rationale:** Adds new fields without breaking existing code. Tuple unpacking still works (`start, end = result`). Subscripting still works (`result[0]`). Type hints updated but runtime behavior backward-compatible.
- **Tradeoff:** Slightly complex dataclass (custom methods) vs simpler return type change that would break all callers.

**4. Parse by field prefix, not line number**
- **Context:** Old parser used `lines[0]`, `lines[1]` for start/end. Need to support both old (3 lines) and new (11 lines) formats.
- **Decision:** Parse fields by prefix (`"SpliceJunctionCount:"`) not line number.
- **Rationale:** Forward-compatible if we add more fields in the future. Doesn't assume specific line ordering. Handles missing fields gracefully (returns None).
- **Implementation:** Helper function `_get_field(prefix)` searches all lines for matching prefix.

**5. NA handling: different semantics for coordinates (0) vs metadata (None)**
- **Context:** cDNA_positions.txt uses "NA" string for missing values. Need to map to Python types.
- **Decision:** NA for start/end becomes `0`. NA for junction fields becomes `None`.
- **Rationale:**
  - start=0, end=0 maintains existing behavior for calculate_coverage_outside_insert (it checks if region is too small).
  - None for junction fields signals "not available" to downstream code (can check `if junction_count is not None`).
  - Avoids mixing sentinel values and exceptions.
- **Implementation:** In parse_insert_region, check if value == "NA", return None; for start/end specifically, use 0 as fallback.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

**Type checking complexity for InsertRegionData compatibility**
- **Issue:** Functions expecting `tuple[int, int]` didn't accept `InsertRegionData` even with `__iter__` and `__getitem__`.
- **Cause:** mypy doesn't infer tuple-like behavior from magic methods alone.
- **Resolution:** Updated type hints to accept `tuple[int, int] | InsertRegionData` across 5 functions in 3 modules.
- **Files updated:** compare_alignments.py, coverage_metrics.py, resistance_gene_detection.py
- **Outcome:** All type checks pass, backward compatibility preserved, new functionality works.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Phase 15 Plan 02 (if exists) or Phase 16 (Summary Visualization):**
- Junction metadata extraction complete
- Data layer ready for visualization
- summary.tsv includes all junction fields

**Data structure established:**
- 11-line cDNA_positions.txt format documented
- InsertRegionData API stable
- Backward compatibility tested

**Testing infrastructure:**
- 397 unit tests passing (7 new for splice junction extraction, 9 new for InsertRegionData)
- Coverage for primary-only filtering, backward compatibility, NA handling
- Type checking strict mode passing

**Known limitations:**
- Junction metadata only available for plasmids that align to human genome
- Unmapped plasmids get all-NA values (expected behavior)
- Multi-exon insert detection relies on primary alignment quality

---
*Phase: 15-splice-junction-metadata*
*Completed: 2026-02-16*
