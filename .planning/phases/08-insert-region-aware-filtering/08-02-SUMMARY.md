---
phase: 08-insert-region-aware-filtering
plan: 02
subsystem: reporting
status: complete
tags: [html-reports, visualization, matplotlib, plotly, jinja2]
requires: [08-01]
provides:
  - 5-category read assignment table in single-sample HTML reports
  - Visual distinction for excluded categories (Tied, Backbone_Only, Ambiguous)
  - Consistent color mapping across Plotly and matplotlib plots
affects: []
key-files:
  created: []
  modified:
    - plasmicheck/templates/report_template.html
    - plasmicheck/scripts/generate_report.py
    - plasmicheck/scripts/plotting/matplotlib_backend.py
decisions: []
metrics:
  duration: "4 min 40 sec"
  completed: 2026-02-14
  commits: 2
---

# Phase 08 Plan 02: Report Display Updates Summary

**One-liner:** HTML reports now display all 5 read categories (Plasmid, Human, Tied, Backbone_Only, Ambiguous) with visual distinction for excluded categories

## What Was Built

Updated single-sample HTML report generation to display the new 5-category read classification system with clear visual indicators:

1. **Structured Read Assignment Table** - Replaced generic summary table with custom 5-row table showing:
   - Category names (Plasmid, Human, Tied, Backbone_Only, Ambiguous)
   - Read counts per category
   - Percentage of total reads per category
   - "Included in Ratio" column with color-coded badges (green for Yes, gray for No)
   - Excluded categories (Tied, Backbone_Only, Ambiguous) visually dimmed using Bootstrap `text-muted` class

2. **Additional Metrics Section** - Extracted and displayed:
   - Coverage Outside Insert Region
   - Mismatches Near Insert Boundaries

3. **Color Mapping Consistency** - Applied `ASSIGNMENT_COLORS` to Plotly plots:
   - Box plot uses color_discrete_map for all 5 categories
   - Scatter plot uses color_discrete_map for all 5 categories
   - Ensures Backbone_Only (#D3D3D3) and Ambiguous (#A9A9A9) appear in correct gray tones

4. **Matplotlib Backend Updates** - Extended static PNG plot generation:
   - `generate_boxplot_matplotlib()` now iterates over all 5 categories
   - `generate_scatter_matplotlib()` now iterates over all 5 categories
   - Empty categories (0 reads) render gracefully without visual artifacts

## Technical Implementation

**Report Template Changes:**
- Replaced `{{ summary_df|safe }}` with structured HTML table
- Used Bootstrap classes: `table`, `table-bordered`, `badge`, `badge-success`, `badge-secondary`, `text-muted`
- Added Jinja2 variables for all category counts and percentages
- Added Additional Metrics section for insert-region-specific metrics

**Report Generation Logic:**
- Added helper function `_get_count()` to extract category counts from summary_df
- Added helper function `_pct()` to calculate percentages
- Extracted `CoverageOutsideINSERT` and `MismatchesNearINSERT` from summary_df
- Updated `generate_report()` signature with 13 new parameters for category data
- Passed individual category variables to both interactive and non-interactive templates
- Maintained backward compatibility by keeping `summary_df` parameter

**Color Mapping:**
- Imported `ASSIGNMENT_COLORS` from `plasmicheck/scripts/plotting/colors.py`
- Added `color_discrete_map=ASSIGNMENT_COLORS` to `px.box()` call
- Added `color_discrete_map=ASSIGNMENT_COLORS` to `px.scatter()` call

**Matplotlib Backend:**
- Changed `categories = ["Plasmid", "Human", "Tied"]` to include `"Backbone_Only", "Ambiguous"`
- Applied to both `generate_boxplot_matplotlib()` and `generate_scatter_matplotlib()`
- Color lookup via `ASSIGNMENT_COLORS.get(cat, "#999999")` handles new entries automatically

## Deviations from Plan

None - plan executed exactly as written.

## Testing Results

- **All 195 tests passing** (unchanged from 08-01)
- **Lint:** Clean (ruff check passed)
- **Format:** Clean (ruff format check passed)
- **Type check:** Clean (mypy strict mode passed)
- **Template render:** Jinja2 syntax validation passed
- **CI check:** All checks passed

## Key Decisions Made

No architectural decisions required. Implementation followed plan specifications exactly.

## Files Modified

| File | Lines Changed | Purpose |
|------|---------------|---------|
| plasmicheck/templates/report_template.html | +53, -3 | Structured 5-category table, additional metrics section |
| plasmicheck/scripts/generate_report.py | +104, -1 | Extract category counts, pass to template, color mapping |
| plasmicheck/scripts/plotting/matplotlib_backend.py | +2, -2 | Extend category lists for static plots |

## Integration Points

**Upstream Dependencies:**
- Plan 08-01: Provides 5-category classification in `compare_alignments.py`
- Plan 08-01: Provides `ASSIGNMENT_COLORS` dict in `colors.py`
- Plan 08-01: Adds `CoverageOutsideINSERT` and `MismatchesNearINSERT` to summary.tsv

**Downstream Impact:**
- Users now see complete read assignment breakdown in reports
- Excluded categories clearly distinguished from included categories
- Plots visually consistent across Plotly and matplotlib backends

## Commits

| Commit | Description | Files |
|--------|-------------|-------|
| 009ba9f | feat(08-02): update report template with 5-category table | report_template.html, generate_report.py |
| 1ca3e58 | feat(08-02): update matplotlib backend for 5 categories | matplotlib_backend.py |

## Validation

**Requirement FILT-05:** ✓ New categories displayed in single-sample HTML reports

**Must-Haves Verified:**
- ✓ Summary table displays Backbone_Only and Ambiguous counts
- ✓ Excluded categories visually dimmed with text-muted class
- ✓ Scatter plot and box plot show all 5 categories with appropriate colors
- ✓ Report renders correctly when Backbone_Only and Ambiguous counts are 0
- ✓ Matplotlib backend generates plots with all 5 categories

## Next Phase Readiness

**Phase 8 Complete:** Both plans (08-01, 08-02) finished successfully.

**Ready for Phase 9 (Coverage Metrics):**
- Insert-region awareness established (08-01)
- Report display infrastructure updated (08-02)
- 5-category system fully integrated into visualization layer
- No blockers for coverage metric implementation

**Recommendations for Phase 9:**
- Use `read_overlaps_insert()` from 08-01 for insert coverage calculations
- Exclude Backbone_Only reads from insert coverage metrics
- Add new coverage metrics to Additional Metrics section in report template
- Consider adding coverage distribution plots to report

---
*Summary completed: 2026-02-14*
*Execution time: 4 min 40 sec*
*All requirements satisfied*
