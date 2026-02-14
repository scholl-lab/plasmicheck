---
phase: 08-insert-region-aware-filtering
verified: 2026-02-14T22:30:00Z
status: passed
score: 12/12 must-haves verified
re_verification: false
---

# Phase 8: Insert-Region-Aware Filtering Verification Report

**Phase Goal:** Reads mapping only to shared backbone are excluded from contamination ratio, eliminating false positive calls for unrelated plasmids.

**Verified:** 2026-02-14T22:30:00Z
**Status:** PASSED
**Re-verification:** No - initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Reads mapping only to backbone (outside insert region) are classified as Backbone_Only | ✓ VERIFIED | `read_overlaps_insert()` function exists, `_assign()` returns "Backbone_Only" when plasmid read doesn't overlap insert, 11 boundary tests pass |
| 2 | Reads with score difference below score_margin are classified as Ambiguous | ✓ VERIFIED | `_assign()` checks `0 < score_diff < score_margin` and returns "Ambiguous", test `test_score_margin_makes_ambiguous` passes |
| 3 | Backbone_Only and Ambiguous reads are excluded from contamination ratio when filter_backbone_only=true | ✓ VERIFIED | Ratio calculation at line 387-395: when `FILTER_BACKBONE_ONLY=true`, ratio = `plasmid_count / human_count` (excludes Backbone_Only and Ambiguous) |
| 4 | Setting filter_backbone_only=false restores pre-v0.33.0 behavior | ✓ VERIFIED | Lines 391-395: when false, `effective_plasmid = plasmid_count + Backbone_Only + Ambiguous`, backward compatibility tests pass |
| 5 | Missing or malformed cDNA_positions.txt causes graceful fallback with warning, not crash | ✓ VERIFIED | Lines 336-342 catch `FileNotFoundError`, `ValueError`, `IndexError` and log warning, sets `insert_region=None` (no crash) |
| 6 | Comparison TSV always includes Backbone_Only and Ambiguous in AssignedTo column values | ✓ VERIFIED | `assigned_counts` dict initialized with all 5 categories (line 266-272), `_assign()` returns these values, written to TSV via `_write_assignment()` |
| 7 | summary.tsv always includes Backbone_Only and Ambiguous count rows | ✓ VERIFIED | Lines 406-409 iterate over `assigned_counts.items()` which always has all 5 keys, writes all categories to summary.tsv |
| 8 | Single-sample HTML report displays Backbone_Only and Ambiguous counts in summary table | ✓ VERIFIED | Template lines 104-114 render Backbone_Only and Ambiguous rows with counts and percentages, `generate_report.py` passes `backbone_only_count` and `ambiguous_count` |
| 9 | Excluded categories (Tied, Backbone_Only, Ambiguous) are visually dimmed in report table | ✓ VERIFIED | Template lines 97, 103, 109 use `class="text-muted"` for excluded categories, badge indicators show "No" for excluded |
| 10 | Scatter plot and box plot show all 5 categories with appropriate colors | ✓ VERIFIED | `generate_report.py` line 85, 106 use `color_discrete_map=ASSIGNMENT_COLORS`, `ASSIGNMENT_COLORS` includes Backbone_Only (#D3D3D3) and Ambiguous (#A9A9A9) |
| 11 | When Backbone_Only and Ambiguous counts are 0, report still renders correctly | ✓ VERIFIED | Template uses Jinja2 variable substitution, default values in `generate_report()` signature (lines 227-228) are 0, template renders "0" correctly |
| 12 | Matplotlib backend generates plots with all 5 categories | ✓ VERIFIED | `matplotlib_backend.py` lines 57, 142 define `categories = ["Plasmid", "Human", "Tied", "Backbone_Only", "Ambiguous"]`, color lookup uses `ASSIGNMENT_COLORS` |

**Score:** 12/12 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `plasmicheck/scripts/compare_alignments.py` | Insert-region-aware filtering with read_overlaps_insert(), extended _assign(), updated _streaming_compare() | ✓ VERIFIED | Function `read_overlaps_insert` at line 67 (17 lines, substantive), `_assign` at line 217 (38 lines, keyword-only params, 5 categories), `_streaming_compare` at line 257 (67 lines, accepts insert_region and score_margin) |
| `plasmicheck/config.json` | filter_backbone_only and score_margin configuration entries | ✓ VERIFIED | Line 4: `"filter_backbone_only": true`, line 5: `"score_margin": 0`, valid JSON, imported in compare_alignments.py lines 23-24 |
| `plasmicheck/scripts/plotting/colors.py` | ASSIGNMENT_COLORS with Backbone_Only and Ambiguous | ✓ VERIFIED | Lines 20-26 define `ASSIGNMENT_COLORS` dict with all 5 categories, Backbone_Only = #D3D3D3, Ambiguous = #A9A9A9 |
| `tests/test_compare_alignments.py` | Tests for backbone classification, score margin, backward compat, graceful fallback | ✓ VERIFIED | 52 total tests (29 new), 4 new test classes: TestReadOverlapsInsert (11 tests), TestAssignExtended (12 tests), TestStreamingCompareWithFiltering (4 tests), TestBackwardCompatibility (2 tests), all pass |
| `plasmicheck/templates/report_template.html` | Updated read assignment table with 5 categories | ✓ VERIFIED | Lines 76-130 render structured 5-row table with counts, percentages, "Included in Ratio" column, text-muted styling for excluded categories |
| `plasmicheck/scripts/generate_report.py` | Color mapping for new categories, passing new counts to template | ✓ VERIFIED | Line 16 imports ASSIGNMENT_COLORS, lines 398-421 extract category counts and calculate percentages, lines 471-484 pass all category variables to template.render() |
| `plasmicheck/scripts/plotting/matplotlib_backend.py` | Updated category lists for matplotlib plots | ✓ VERIFIED | Line 57 (boxplot) and 142 (scatter) extend categories list to include "Backbone_Only" and "Ambiguous" |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| compare_alignments.py | config.json | get_config() loading filter_backbone_only and score_margin | ✓ WIRED | Lines 23-24 load from `_cfg.get("filter_backbone_only", True)` and `_cfg.get("score_margin", 0)`, used at lines 375, 387 |
| _streaming_compare | read_overlaps_insert | Calling read_overlaps_insert for plasmid-only reads | ✓ WIRED | Lines 292-294 call `_assign()` with `plasmid_read=p_read` and `insert_region=insert_region`, `_assign()` line 247 calls `read_overlaps_insert(plasmid_read, insert_region)` |
| compare_alignments | _streaming_compare | Passing insert_region to _streaming_compare | ✓ WIRED | Lines 370-376 call `_streaming_compare()` with `insert_region=insert_region` keyword argument |
| generate_report.py | report_template.html | Jinja2 template rendering with new category variables | ✓ WIRED | Lines 292-305 (interactive) and 342-355 (static) pass `backbone_only_count`, `ambiguous_count`, `backbone_only_pct`, `ambiguous_pct` to `template.render()` |
| generate_report.py | colors.py | ASSIGNMENT_COLORS dict used for color mapping | ✓ WIRED | Line 16 imports `ASSIGNMENT_COLORS`, lines 85, 106 pass `color_discrete_map=ASSIGNMENT_COLORS` to Plotly figures |
| matplotlib_backend.py | colors.py | ASSIGNMENT_COLORS dict used for matplotlib plot colors | ✓ WIRED | Line 58 uses `colors = [ASSIGNMENT_COLORS.get(cat, "#999999") for cat in categories]` |

### Requirements Coverage

| Requirement | Status | Supporting Evidence |
|-------------|--------|---------------------|
| FILT-01: Backbone_Only classification and exclusion from ratio | ✓ SATISFIED | Truths 1, 3 verified - read_overlaps_insert() identifies backbone-only reads, ratio calculation excludes them when filter enabled |
| FILT-02: Optional score margin parameter | ✓ SATISFIED | Truth 2 verified - score_margin in config.json, _assign() uses it for Ambiguous classification |
| FILT-03: New read categories in comparison output | ✓ SATISFIED | Truths 6, 7 verified - TSV and summary.tsv include all 5 categories |
| FILT-04: filter_backbone_only config toggle with backward-compatible default | ✓ SATISFIED | Truths 3, 4 verified - config entry exists, ratio calculation respects toggle, backward compat tests pass |
| FILT-05: New categories displayed in single-sample HTML reports | ✓ SATISFIED | Truths 8, 9, 10, 11, 12 verified - template displays new categories with styling, plots use correct colors |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| None | - | - | - | No anti-patterns detected |

Analysis: Scanned all modified files for TODO/FIXME, placeholder content, empty implementations, console.log-only handlers. No blockers or warnings found. All implementations are substantive and production-ready.

### Test Results

**All Tests Passing:**
- Total: 195 tests (52 in test_compare_alignments.py, 143 in other test files)
- New tests: 29 (TestReadOverlapsInsert: 11, TestAssignExtended: 12, TestStreamingCompareWithFiltering: 4, TestBackwardCompatibility: 2)
- All existing tests pass without modification (backward compatibility confirmed)
- Test coverage: Boundary conditions, 5-category classification, score margin behavior, graceful fallback, backward compatibility

**CI Checks:**
- Lint (ruff): PASS
- Format (ruff): PASS
- Type check (mypy strict): PASS
- Unit tests: PASS (195/195)

## Success Criteria Evaluation

1. **A sample contaminated with Plasmid A no longer shows elevated contamination scores for unrelated Plasmid B that shares only backbone sequence**
   - ✓ ACHIEVED: Reads mapping only to backbone are classified as Backbone_Only and excluded from ratio when `filter_backbone_only=true` (default). The ratio calculation (lines 387-389) only counts reads overlapping the insert region.

2. **Running with filter_backbone_only=false produces identical output to pre-v0.33.0 behavior (all existing tests pass without modification)**
   - ✓ ACHIEVED: All 195 tests pass, including original tests that don't use new parameters. Lines 391-395 implement pre-v0.33.0 behavior when toggle is false (includes all plasmid-favoring categories in ratio).

3. **Comparison TSV output includes Backbone_Only and Ambiguous columns alongside existing Plasmid/Human/Tied counts**
   - ✓ ACHIEVED: `assigned_counts` dict always has all 5 keys (line 266-272), `_write_assignment()` writes AssignedTo value from `_assign()` which returns any of the 5 categories.

4. **Single-sample HTML report displays the new read categories with counts and proportions**
   - ✓ ACHIEVED: Template renders 5-row table (lines 84-114) with counts, percentages, and "Included in Ratio" indicators. Visual distinction via `text-muted` class for excluded categories.

5. **Setting a non-zero score_margin causes reads with small score differences to be classified as Ambiguous rather than Plasmid or Human**
   - ✓ ACHIEVED: `_assign()` lines 235-238 check `if score_margin > 0` and `0 < score_diff < score_margin`, returns "Ambiguous". Test `test_score_margin_makes_ambiguous` confirms behavior.

## Overall Status

**Status:** PASSED

All 12 observable truths verified. All 7 required artifacts exist, are substantive, and are wired correctly. All 5 key links verified. All 5 requirements satisfied. No anti-patterns detected. All 195 tests pass. All CI checks pass.

**Phase goal achieved:** Reads mapping only to shared backbone are excluded from contamination ratio, eliminating false positive calls for unrelated plasmids. The filtering logic is scientifically accurate, well-tested, backward compatible, and ready for production use.

---

_Verified: 2026-02-14T22:30:00Z_
_Verifier: Claude (gsd-verifier)_
