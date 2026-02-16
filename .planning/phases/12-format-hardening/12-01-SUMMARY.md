---
phase: 12
plan: 01
name: "Format Detection & SnapGene Support"
subsystem: conversion
status: complete
completed: 2026-02-15
duration: 6min

requires:
  - plasmicheck.scripts.utils.sanitize_filename
  - Bio.SeqIO (snapgene parser in BioPython 1.84+)

provides:
  - format_detection.py module with magic byte detection
  - SnapGene .dna file conversion
  - Auto-detection for all plasmid formats
  - Actionable error messages for conversion failures

affects:
  - 12-02 (pipeline format detection uses this module)
  - Future xDNA/SnapGene users (seamless conversion)

tech-stack:
  added:
    - plasmicheck.scripts.format_detection
  patterns:
    - Magic byte detection for binary formats
    - Text-based fallback for GenBank/FASTA
    - Sanitized filename as record ID override

key-files:
  created:
    - plasmicheck/scripts/format_detection.py
  modified:
    - plasmicheck/scripts/convert_plasmidfile_to_fasta.py
    - plasmicheck/cli.py
    - plasmicheck/config.json

decisions:
  - id: FMT-D01
    context: "Record IDs from xDNA files are unreliable (use comment field)"
    decision: "Override all record.id and record.name with sanitize_filename(Path(input_file).stem)"
    rationale: "Ensures consistent, valid IDs across all input formats"
    impact: "All converted FASTA files now use filename-based IDs"

  - id: FMT-D02
    context: "BioPython errors are cryptic for end users"
    decision: "Wrap BioPython parse calls with keyword detection and actionable error messages"
    rationale: "Users know exactly what to do when conversion fails"
    impact: "Error messages suggest GenBank export as fallback"

  - id: FMT-D03
    context: "sanitize_filename was duplicated in convert_plasmidfile_to_fasta.py and utils.py"
    decision: "Remove local copy, import from utils.py (single source of truth)"
    rationale: "Eliminates code duplication flagged in research"
    impact: "Cleaner codebase, single function to maintain"

  - id: FMT-D04
    context: "--file_type flag was required, no auto-detection"
    decision: "Add --format flag (optional), auto-detect when not specified, keep -t as deprecated alias"
    rationale: "Better UX, backward compatible, follows modern CLI patterns"
    impact: "Users can omit format flag, CLI detects from file content"

tags:
  - format-detection
  - snapgene
  - xdna
  - conversion
  - cli
  - ux-improvement
---

# Phase 12 Plan 01: Format Detection & SnapGene Support Summary

**One-liner:** Magic byte format detection with SnapGene .dna support, filename-based record IDs, and actionable conversion errors

## What Was Built

### Core Features
1. **Format detection module** (`plasmicheck/scripts/format_detection.py`)
   - `detect_format()` — auto-detects SnapGene, xDNA, GenBank, FASTA from file content
   - Magic byte detection for binary formats (SnapGene signature, xDNA header)
   - Text-based fallback for GenBank (LOCUS) and FASTA (>)
   - Actionable error messages listing all supported formats

2. **SnapGene .dna file support**
   - Uses BioPython's native SnapGene parser (1.84+)
   - Converts .dna → FASTA via `SeqIO.parse(file, "snapgene")`
   - Optional GenBank preservation (.dna → .gb → .fasta)

3. **Record ID sanitization**
   - All formats now override `record.id` and `record.name` with `sanitize_filename(Path(input_file).stem)`
   - Fixes xDNA issue where `record.id` was `[author` (from comment field)
   - Ensures valid, predictable IDs across all plasmid formats

4. **Actionable error handling**
   - Detects "version" keyword → suggests GenBank export for unsupported xDNA versions
   - Detects "corrupt"/"truncated" → suggests re-export from cloning software
   - Generic parse failures → lists all supported formats with tool suggestions

5. **CLI improvements**
   - New `--format` flag (primary) with choices: snapgene, xdna, genbank, fasta
   - Auto-detection when `--format` not specified (calls `detect_format()`)
   - Kept `-t/--file_type` as hidden deprecated alias (backward compatibility)

### Configuration Updates
- `config.json` updated with `"snapgene": ".dna"` in `supported_formats`
- `conversion.input_formats` now includes "snapgene"

## Implementation Notes

### Magic Byte Detection
- **SnapGene:** First 13 bytes = `b'\x09\x00\x00\x00\x0ESnapGene'`
- **xDNA:** First byte = 0x00 (version), byte 111 = 0xFF (requires 112 bytes)
- **GenBank/FASTA:** Text-based (first line starts with "LOCUS" or ">")

### Code Quality
- Removed duplicated `sanitize_filename` from `convert_plasmidfile_to_fasta.py`
- Single source in `utils.py`, imported everywhere
- All lint, typecheck, and 297 tests passing

## Deviations from Plan

**None** — plan executed exactly as written.

## Verification Results

All verification criteria met:
- ✓ `format_detection.py` imports successfully
- ✓ `sanitize_filename` imported from utils (no local duplication)
- ✓ `make lint` passes (0 errors)
- ✓ `make typecheck` passes (mypy strict, 0 errors)
- ✓ `make test-fast` passes (297/297 tests)
- ✓ `python -m plasmicheck.cli convert --help` shows `--format` flag with 4 choices

## Commits

| Commit | Type | Description | Files |
|--------|------|-------------|-------|
| dc10c93 | feat | Format detection module and conversion refactor | format_detection.py, convert_plasmidfile_to_fasta.py, config.json |
| dc3bf51 | feat | Add --format flag to convert command with auto-detection | cli.py |

## Success Criteria

All success criteria met:
- [x] New format_detection.py module detects SnapGene, xDNA, GenBank, and FASTA formats
- [x] SnapGene .dna files can be converted via `plasmicheck convert`
- [x] All converted FASTA files use sanitized filename as record ID
- [x] Actionable error messages for unsupported xDNA versions, corrupt files, unrecognized formats
- [x] CLI has --format flag with auto-detection default
- [x] sanitize_filename duplication eliminated
- [x] config.json updated with snapgene format
- [x] Lint, typecheck, and existing tests pass

## Next Phase Readiness

**Ready for 12-02:** Pipeline integration can now use `format_detection.detect_format()` to replace extension-based detection at `run_pipeline.py` lines 512-517.

**Blockers:** None

**Concerns:** None — BioPython's SnapGene parser is stable in 1.84+, no known compatibility issues.

## Usage Examples

```bash
# Auto-detect format (new in 12-01)
plasmicheck convert -i plasmid.dna -o plasmid.fasta

# Explicit format (override detection)
plasmicheck convert -f snapgene -i plasmid.dna -o plasmid.fasta

# Backward compatible (deprecated -t flag still works)
plasmicheck convert -t xdna -i plasmid.xdna -o plasmid.fasta

# Preserve GenBank for gene detection
plasmicheck convert -i plasmid.dna -o plasmid.fasta -pg
```

## Impact

### User-Facing
- **SnapGene users** can now use .dna files directly (no manual GenBank export)
- **xDNA users** get helpful error messages when files fail (instead of cryptic BioPython tracebacks)
- **All users** benefit from auto-detection (no need to remember format names)
- **CLI UX** improved (optional --format flag, sensible defaults)

### Developer-Facing
- **Single source** for filename sanitization (utils.py)
- **Reusable** format detection module for future features
- **Extensible** pattern for adding new binary formats (e.g., SnapGene .dna variants)

### Technical Debt
- **Eliminated:** Duplicated sanitize_filename function
- **Improved:** Error handling with actionable messages

## Metrics

- **Duration:** 6 minutes
- **Files modified:** 4 (created 1, modified 3)
- **Lines added:** ~172
- **Lines removed:** ~18
- **Net change:** +154 lines
- **Tests passing:** 297/297
- **Test coverage:** No regression (existing sanitize_filename tests still pass)

## Related Issues

- Addresses issue #75 (some xDNA files don't work) via:
  - Actionable error messages for unsupported versions
  - Record ID override fixes broken xDNA IDs
  - Format detection reduces user confusion
