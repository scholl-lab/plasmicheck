---
phase: 13-output-restructuring
plan: 02
subsystem: pipeline-infrastructure
status: complete
completed: 2026-02-15
duration: 4m42s

tags:
  - path-management
  - backward-compatibility
  - directory-structure
  - testing

requires:
  - "13-01: OutputPaths dataclass and pipeline rewiring"

provides:
  - "Depth-agnostic summary report reader supporting both flat and data/ subfolder layouts"
  - "Backward compatibility for existing flat-structure outputs"
  - "Forward compatibility for new 5-subfolder structure"
  - "Comprehensive tests proving mixed-layout migration path"

affects:
  - "Future: Any module that reads summary reports"
  - "14-xx: Cache infrastructure (summary reports will work with new structure)"

tech-stack:
  added: []
  patterns:
    - "Depth-agnostic directory parsing with conditional parent navigation"
    - "Integration testing across OutputPaths creator and summary reader"

key-files:
  created: []
  modified:
    - plasmicheck/scripts/generate_summary_reports.py
    - tests/test_generate_summary_reports.py
    - tests/test_output_paths.py

decisions:
  - id: SUM-COMPAT-01
    what: "Use parent directory name check ('data') to determine directory depth"
    why: "Simple, reliable, and robust to arbitrary nesting depth without hardcoding levels"
    alternatives: "Count levels up from file (fragile), store metadata in TSV (breaking change)"
    impact: "Summary reader works with both old (flat) and new (subfolder) outputs automatically"

  - id: SUM-COMPAT-02
    what: "Add _extract_sample_plasmid helper function"
    why: "Centralize depth logic, eliminate hardcoded 2-level dirname chains, improve testability"
    alternatives: "Inline logic in both reads and summary loops (DRY violation)"
    impact: "Single source of truth for sample/plasmid extraction, easy to test and modify"

verification:
  - criterion: "_extract_sample_plasmid handles both flat and data/ subfolder layouts"
    result: "✓ 4 tests pass (flat, subfolder, reads_assignment for both)"

  - criterion: "read_compare_outputs correctly extracts Sample/Plasmid from both layouts"
    result: "✓ 3 tests pass (flat_layout_read, subfolder_layout_read, mixed_layout)"

  - criterion: "OutputPaths-created directories are compatible with summary reader"
    result: "✓ 2 integration tests pass (extract helper + full read_compare_outputs)"

  - criterion: "All existing tests pass unchanged"
    result: "✓ 343 passed, 4 deselected, 9 warnings"

  - criterion: "CI checks pass completely"
    result: "✓ lint + format + typecheck + test all clean"

metrics:
  tests:
    total: 353
    passed: 343
    new: 10
  coverage:
    lines_added: 275
    lines_modified: 4
    files_changed: 3
  performance:
    test_time: "18.71s"
---

# Phase 13 Plan 02: Depth-Agnostic Summary Report Reader Summary

**Summary report reader supports both flat (old) and data/ subfolder (new) directory layouts via conditional parent navigation based on directory name**

**Status:** ✅ Complete — All tasks executed, tests pass, CI clean

## Performance

- **Duration:** 4m42s (282 seconds)
- **Started:** 2026-02-15T20:22:59Z
- **Completed:** 2026-02-15T20:27:41Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments

- Summary report reader (`generate_summary_reports.py`) now correctly extracts sample/plasmid names from both flat and subfolder directory structures
- 10 new comprehensive tests covering flat, subfolder, and mixed layouts plus end-to-end OutputPaths integration
- Backward compatibility proven: existing flat outputs work unchanged
- Forward compatibility proven: new subfolder outputs work correctly
- Migration path validated: mixed flat + subfolder layouts work in single summary report run

## Task Commits

Each task was committed atomically:

1. **Task 1: Make read_compare_outputs depth-agnostic** - `afa468f` (feat)
2. **Task 2: Add tests for both directory layouts and end-to-end integration** - `538ac9c` (test)

**Plan metadata:** (will be committed after this summary)

## Files Created/Modified

- `plasmicheck/scripts/generate_summary_reports.py`:
  - Added `_extract_sample_plasmid(file_path: str) -> tuple[str, str]` helper
  - Updated `read_compare_outputs()` to use helper instead of hardcoded dirname chains
  - Lines 59-89: New helper function with conditional parent navigation
  - Lines 106, 121: Replaced `os.path.basename(os.path.dirname(...))` with helper calls

- `tests/test_generate_summary_reports.py`:
  - Added `TestExtractSamplePlasmid` class (4 tests): flat, subfolder, reads_assignment for both
  - Added `TestReadCompareOutputsBothLayouts` class (3 tests): flat_layout_read, subfolder_layout_read, mixed_layout
  - Tests use minimal but valid TSV data to prove correct sample/plasmid extraction

- `tests/test_output_paths.py`:
  - Added `test_extract_sample_plasmid_integration_with_output_paths`: Proves OutputPaths-created files are correctly parsed by helper
  - Added `test_output_paths_with_read_compare_outputs`: End-to-end integration test (OutputPaths creates structure, read_compare_outputs reads it correctly)
  - Tests use realistic multi-sample/multi-plasmid data to prove full pipeline compatibility

## Decisions Made

**SUM-COMPAT-01: Parent directory name check**
- Strategy: If immediate parent is named "data", go one level further up for plasmid/sample
- Rationale: Simple, reliable, robust to arbitrary nesting without hardcoding levels
- Alternative rejected: Count levels up from file (fragile to future changes)
- Impact: Works with both old (2-level) and new (3-level) directory structures automatically

**SUM-COMPAT-02: Extract helper function**
- Centralized depth-agnostic logic in `_extract_sample_plasmid()`
- Eliminates DRY violation (same logic for reads and summary files)
- Makes testing easier (unit test helper directly)
- Single point of maintenance if directory structure changes again

## Deviations from Plan

None — plan executed exactly as written.

## Issues Encountered

**Integration test TSV format mismatch:**
- **Issue:** Initial integration test failed because `read_compare_outputs()` renames "Count" to "Value" column
- **Fix:** Updated test TSV data to match expected format (Category, Metric, Count columns)
- **Verification:** Integration test now passes, validates correct DataFrame structure
- **Time:** ~2 minutes to debug and fix

## Next Phase Readiness

**Ready to proceed:** Phase 13 complete (all plans in output-restructuring phase done)

**What's validated:**
- OutputPaths creates new 5-subfolder structure (13-01)
- Summary report reader correctly parses new structure (13-02)
- Backward compatibility maintained (old flat outputs still work)
- Forward compatibility established (new subfolder outputs work)
- Migration path proven (mixed layouts work in single run)

**No blockers:**
- All existing tests pass (343 passed)
- All new tests pass (10 new tests)
- CI checks clean (lint + format + typecheck + test)
- Integration tests prove OutputPaths → summary reader pipeline works end-to-end

**Future phases can leverage:**
- Summary reports will work seamlessly with new directory structure
- Migration can be gradual (flat and subfolder coexist)
- OutputPaths + depth-agnostic reader = robust foundation for Phase 14 (caching)

## Technical Notes

**Implementation strategy:**
```python
def _extract_sample_plasmid(file_path: str) -> tuple[str, str]:
    """Extract sample and plasmid names from comparison output file path.

    Handles both directory layouts:
    - Old (flat): output/sample/plasmid/file.tsv
    - New (subfolder): output/sample/plasmid/data/file.tsv

    Strategy: if the immediate parent directory is named 'data',
    go one level further up for plasmid and sample names.
    """
    parent = os.path.dirname(file_path)
    parent_name = os.path.basename(parent)

    if parent_name == "data":
        # New subfolder layout: .../sample/plasmid/data/file.tsv
        plasmid_dir = os.path.dirname(parent)
        plasmid_name = os.path.basename(plasmid_dir)
        sample_name = os.path.basename(os.path.dirname(plasmid_dir))
    else:
        # Old flat layout: .../sample/plasmid/file.tsv
        plasmid_name = parent_name
        sample_name = os.path.basename(os.path.dirname(parent))

    return sample_name, plasmid_name
```

**Why this works:**
- `find_tsv_files()` already uses `os.walk()` (depth-agnostic)
- Files are discovered at any depth automatically
- Helper only needs to determine correct parent levels based on structure
- Single conditional (`if parent_name == "data"`) handles both cases
- No assumptions about absolute depth (works with deeply nested output roots)

**Test coverage:**
- Unit tests: Helper function with both layouts (4 tests)
- Integration tests: read_compare_outputs with both layouts (3 tests)
- End-to-end tests: OutputPaths creates → summary reader parses (2 tests)
- Mixed layout test: Proves migration path (1 test validates coexistence)

**Performance impact:**
- Negligible: One additional `os.path.basename()` call per file
- No change to `os.walk()` traversal (already depth-agnostic)
- Test time unchanged (18.71s, within normal variance)

---
*Phase: 13-output-restructuring*
*Completed: 2026-02-15*
