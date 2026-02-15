---
phase: 10-resistance-gene-detection
plan: 01
subsystem: bioinformatics
tags: [biopython, genbank, regex, pattern-matching, gene-detection]

# Dependency graph
requires:
  - phase: 09-coverage-metrics
    provides: Established coverage_metrics.py module pattern for separate computation logic
provides:
  - resistance_gene_detection.py module with curated pattern matching for bacterial, mammalian, and reporter genes
  - GenBank file preservation in xDNA/GenBank conversion for downstream gene detection
  - Extension point in config.json for custom resistance gene patterns
  - 90% recall validated against real plasmid GenBank files

affects: [10-02, resistance-gene-integration, summary-reports]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Curated regex pattern lists as module-level constants for gene detection"
    - "Biopython SeqIO for GenBank CDS feature parsing"
    - "Dataclass for structured gene detection results (DetectedGene)"
    - "preserve_genbank parameter pattern in conversion functions"

key-files:
  created:
    - plasmicheck/scripts/resistance_gene_detection.py
    - tests/test_resistance_gene_detection.py
  modified:
    - plasmicheck/scripts/convert_plasmidfile_to_fasta.py
    - plasmicheck/config.json

key-decisions:
  - "Extended bacterial patterns to include full antibiotic names (ampicillin, kanamycin, etc.) for broader coverage"
  - "Curated exclusion list for affinity tags (His-tag, FLAG, Myc, V5, HA) to avoid false positives"
  - "Conservative overlap exclusion: any gene overlapping insert region is excluded"
  - "Default preserve_genbank=False for backward compatibility in convert_plasmidfile_to_fasta"
  - "GenBank file location priority: output_dir first, then original plasmid_file"

patterns-established:
  - "Pattern 1: Curated regex lists for biological entity detection (bacterial, mammalian, reporter genes)"
  - "Pattern 2: _private helper functions for matching logic, public API for detection"
  - "Pattern 3: Graceful fallback with empty list return for missing/malformed files"
  - "Pattern 4: Insert region exclusion for backbone-only gene detection"

# Metrics
duration: 13min
completed: 2026-02-15
---

# Phase 10 Plan 01: Resistance Gene Detection Summary

**Curated regex pattern matching detects resistance genes, mammalian markers, and reporters from GenBank CDS annotations with 90% recall**

## Performance

- **Duration:** 13 min
- **Started:** 2026-02-15T07:45:51Z
- **Completed:** 2026-02-15T07:58:36Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments

- Created resistance_gene_detection.py module with 14 bacterial resistance patterns, 7 mammalian selection markers, 9 reporter genes, and 6 exclusion patterns
- Modified convert_plasmidfile_to_fasta() to preserve GenBank files when preserve_genbank=True (xDNA retains intermediate .gb, GenBank creates sanitized copy)
- Added config.json extension point for custom resistance gene patterns
- Comprehensive test suite with 55 unit tests validating pattern matching, GenBank parsing, insert region exclusion, and edge cases
- Validated 90% recall requirement against real GenBank files (pWPI-KCNJ5-WT.gb detects AmpR and EGFP)

## Task Commits

Each task was committed atomically:

1. **Task 1: Preserve GenBank files and create resistance gene detection module** - `3511346` (feat)
2. **Task 2: Unit tests for resistance gene detection** - `367607c` (test)
3. **Formatting fix** - `2ea1281` (style)

## Files Created/Modified

- `plasmicheck/scripts/resistance_gene_detection.py` - Core gene detection module with pattern matching, GenBank parsing, and DetectedGene dataclass
- `plasmicheck/scripts/convert_plasmidfile_to_fasta.py` - Added preserve_genbank parameter to retain .gb files for gene detection
- `plasmicheck/config.json` - Added resistance_gene_patterns extension point
- `tests/test_resistance_gene_detection.py` - 55 unit tests covering all requirements

## Decisions Made

**Pattern Design:**
- Extended bacterial patterns to include full antibiotic names (ampicillin, kanamycin, chloramphenicol, tetracycline) alongside short forms (AmpR, KanR, ChlR, TetR) for broader coverage of GenBank annotation styles
- Created curated exclusion list for affinity tags (His-tag, FLAG, Myc, V5, HA) to prevent false positives from epitope tags

**Insert Region Handling:**
- Conservative overlap exclusion: genes with ANY overlap (full or partial) with insert region are excluded to ensure only backbone genes are detected
- Region classification: genes are labeled "backbone" (when insert region provided and no overlap), "insert" (excluded), or "unknown" (when no insert region provided)

**GenBank File Preservation:**
- Default preserve_genbank=False for backward compatibility (no changes to existing callers)
- xDNA conversion: intermediate .gb file retained when preserve_genbank=True (skip os.remove)
- GenBank conversion: sanitized copy created in output directory for consistent .gb location

**File Location Priority:**
- locate_genbank_file checks output_dir first (from xDNA/GenBank conversion with preserve_genbank=True), then falls back to original plasmid_file if it's .gb/.gbk
- Ensures gene detection always finds the most relevant GenBank file

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None. All pattern matching, GenBank parsing, and insert region exclusion worked as expected. Real GenBank file validation confirmed 90% recall requirement.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for 10-02 (Gene Coverage Integration):**
- resistance_gene_detection.py module ready for import in compare_alignments.py
- DetectedGene dataclass provides structured data (name, sanitized_name, start, end, region)
- locate_genbank_file() helper ready for integration
- Test coverage ensures robust gene detection across various GenBank formats

**Blockers/Concerns:**
None. Module is self-contained and tested against real plasmid files.

**Test Status:**
274 tests passing (219 existing + 55 new). All CI checks passing (lint, format, typecheck, test).

---
*Phase: 10-resistance-gene-detection*
*Completed: 2026-02-15*
