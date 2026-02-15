---
phase: 10-resistance-gene-detection
plan: 02
subsystem: contamination-detection
tags: [resistance-genes, genbank, coverage, pysam, biopython]

# Dependency graph
requires:
  - phase: 10-01
    provides: resistance_gene_detection.py module with curated pattern matching
  - phase: 09-coverage-metrics
    provides: get_depth_array() and compute_coverage_stats() for per-region coverage
provides:
  - Per-gene coverage metrics in summary.tsv (ResGene_ rows)
  - GenBank file preservation and forwarding through pipeline
  - ResistanceGenesAvailable flag distinguishing "no data" vs "looked but found nothing"
affects: [11-summary-report-integration]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Graceful fallback: gene detection failures don't crash pipeline"
    - "GenBank file path forwarding: pipeline determines path, passes to compare_alignments"
    - "Per-gene coverage via Phase 9 infrastructure (pysam-based)"

key-files:
  created: []
  modified:
    - plasmicheck/scripts/compare_alignments.py
    - plasmicheck/scripts/run_pipeline.py
    - tests/test_compare_alignments.py
    - tests/test_resistance_gene_detection.py

key-decisions:
  - "genbank_file parameter is optional (None default) for backward compatibility"
  - "Graceful exception handling: gene detection errors log warning, continue without gene data"
  - "Per-gene coverage computed before streaming compare (needs BAM file access)"
  - "GenBank preservation enabled by default in pipeline (preserve_genbank=True)"
  - "ResistanceGenesAvailable: false when no GenBank vs ResistanceGenes: None when GenBank but no genes"

patterns-established:
  - "Integration tests mock BAM files and coverage metrics for fast unit testing"
  - "Recall validation tests target specific known GenBank files for coverage"
  - "Summary.tsv writes gene data after coverage metrics, before end of file"

# Metrics
duration: 10 min
completed: 2026-02-15
---

# Phase 10 Plan 02: Resistance Gene Detection Integration Summary

**Per-gene coverage metrics in summary.tsv via integrated gene detection and Phase 9 pysam infrastructure**

## Performance

- **Duration:** 10 min
- **Started:** 2026-02-15T08:04:08Z
- **Completed:** 2026-02-15T08:13:45Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- Integrated resistance gene detection into compare_alignments.py with per-gene coverage
- GenBank file preservation in pipeline and path forwarding to compare_alignments
- ResGene_ rows in summary.tsv: MeanDepth, Breadth, Start, End, Region per detected gene
- ResistanceGenesAvailable flag distinguishes "no GenBank file" from "GenBank but no genes"
- 280 tests passing (up from 274): 3 integration tests + 3 recall validation tests
- Graceful fallback when GenBank unavailable or gene detection fails

## Task Commits

Each task was committed atomically:

1. **Task 1: Integrate gene detection into compare_alignments.py and run_pipeline.py** - `a80d027` (feat)
2. **Task 2: Integration tests for resistance gene pipeline flow** - `a0a1851` (test)

**Plan metadata:** (to be added after this summary)

## Files Created/Modified

- `plasmicheck/scripts/compare_alignments.py` - Added genbank_file parameter, gene detection call, per-gene coverage computation, summary.tsv ResGene_ rows
- `plasmicheck/scripts/run_pipeline.py` - Added preserve_genbank=True to conversion, GenBank path determination and forwarding
- `tests/test_compare_alignments.py` - Added 3 integration tests (ResGene_ rows, no GenBank, GenBank no genes)
- `tests/test_resistance_gene_detection.py` - Added 3 recall validation tests (general recall, pWPI-KCNJ5-WT AmpR, pCMV6 markers)

## Decisions Made

**genbank_file parameter design:**
- Optional parameter (None default) for backward compatibility
- Accepts explicit path or uses locate_genbank_file() fallback
- Graceful exception handling: failures log warning, continue without gene data
- Rationale: Existing pipeline runs without GenBank must continue working

**GenBank preservation in pipeline:**
- preserve_genbank=True passed to convert_plasmidfile_to_fasta
- GenBank path determined based on plasmid_file_type (genbank vs xdna)
- xDNA conversions check for .gb file in output_subfolder
- Rationale: Gene detection requires GenBank annotations, must preserve during conversion

**Per-gene coverage via Phase 9 infrastructure:**
- get_depth_array() extracts per-base depth for gene region
- compute_coverage_stats() computes mean, breadth for each gene
- Computed before streaming compare (requires BAM file access)
- Rationale: Reuse existing pysam-based coverage infrastructure, avoid duplication

**ResistanceGenesAvailable flag:**
- false when no GenBank file available (couldn't look)
- true when GenBank present (looked, may or may not have found genes)
- ResistanceGenes: None when looked but found nothing
- Rationale: Distinguish "no data source" from "looked but empty" for downstream analysis

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Phase 11 (Summary Report Integration):**
- Per-gene coverage metrics available in summary.tsv with ResGene_ prefix
- ResistanceGenesAvailable flag enables conditional report sections
- Gene detection fully integrated into pipeline flow
- Graceful fallback ensures existing reports continue working

**Verification:**
- 280 tests passing (6 new tests added)
- All CI checks pass (lint, typecheck, format, tests)
- Real GenBank file recall validated (pWPI-KCNJ5-WT.gb detects AmpR)

**Blockers:** None

**Recommendations for Phase 11:**
- Parse ResGene_ rows from summary.tsv for report display
- Create "Resistance Gene Detection" card in HTML report template
- Display gene coverage table with name, region, depth, breadth
- Show "Not Available" when ResistanceGenesAvailable: false
- Show "None Detected" when ResistanceGenes: None

---
*Phase: 10-resistance-gene-detection*
*Completed: 2026-02-15*
