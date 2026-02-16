---
phase: 16
plan: 01
subsystem: reporting
tags: [visualization, plotly, matplotlib, verdict, bar-chart]
requires: [15-02]
provides:
  - VERDICT_COLORS dictionary for consistent verdict color mapping
  - aggregate_verdicts() function with No Result detection for missing sample-plasmid pairs
  - plot_verdict_bar_by_plasmid() for interactive Plotly stacked bar charts
  - plot_verdict_bar_by_sample() for interactive Plotly stacked bar charts
  - generate_verdict_bar_matplotlib() for static PNG fallback
affects: [16-02]
tech-stack:
  added: []
  patterns:
    - Stacked bar chart aggregation with cross-product No Result detection
    - Hover text formatting with name list truncation (max 10 items)
    - Dual-mode plotting (Plotly interactive + matplotlib static fallback)
key-files:
  created: []
  modified:
    - plasmicheck/scripts/plotting/colors.py
    - plasmicheck/scripts/generate_summary_reports.py
    - plasmicheck/scripts/plotting/matplotlib_backend.py
    - tests/test_generate_summary_reports.py
    - tests/test_matplotlib_backend.py
decisions:
  - id: VIS-D01
    context: Verdict aggregation with missing data
    decision: Detect missing sample-plasmid pairs as "No Result" via cross-product left join
    rationale: Some samples may not have verdict data for all plasmids (e.g., failed alignment), showing "No Result" is more accurate than omitting them
  - id: VIS-D02
    context: Hover text readability
    decision: Truncate name lists to 10 items + "...and N more" in hover
    rationale: Long lists make hover tooltips unwieldy, truncation maintains readability
  - id: VIS-D03
    context: Verdict color scheme
    decision: Reuse HEATMAP_COLORS with No Result as lightgray matching Backbone_Only
    rationale: Consistent color scheme across all verdict visualizations, gray for missing/absent data
metrics:
  duration: 440 seconds (~7 minutes)
  commits: 3
  tests-added: 12
  tests-total: 412
  lines-added: 698
completed: 2026-02-16
---

# Phase 16 Plan 01: Verdict Bar Chart Plotting Functions Summary

**One-liner:** Created verdict bar chart aggregation and plotting functions with Plotly interactive and matplotlib static modes, supporting No Result detection for missing sample-plasmid pairs.

## What Was Built

### Core Functionality

**1. VERDICT_COLORS Dictionary (`colors.py`)**
- Maps verdict categories to colors from HEATMAP_COLORS
- Clean → lightblue, Unclear → orange, Contaminated → red
- No Result → #D3D3D3 (lightgray, matches Backbone_Only)

**2. Verdict Aggregation (`aggregate_verdicts()`)**
- Filters summary_df to Category=="Verdict" rows
- Maps full verdict strings to short labels (Contaminated, Unclear, Clean)
- Builds cross-product of all samples × plasmids from full summary_df
- Left-joins verdict data, fills missing pairs with "No Result"
- Returns: (full_verdict_df, by_plasmid_counts, by_sample_counts)
- Ensures column order: Clean, Unclear, Contaminated, No Result

**3. Interactive Plotly Plotting**
- `plot_verdict_bar_by_plasmid()`: Stacked bar chart grouped by plasmid
- `plot_verdict_bar_by_sample()`: Stacked bar chart grouped by sample
- Both use Plotly go.Bar() with barmode='stack'
- Custom hover via customdata with _format_verdict_hover() helper
- Hover shows verdict, count, and list of sample/plasmid names (truncated at 10)
- Outputs HTML to plots/ subdirectory

**4. Static Matplotlib Fallback (`generate_verdict_bar_matplotlib()`)**
- Accepts aggregated DataFrame (group × verdict categories)
- Iterative stacking with np.zeros + cumulative bottom
- Figsize scales with bar count: max(10, len(df) * 0.5)
- Rotates x-labels 45° if >5 bars
- Outputs PNG at 300 DPI

**5. Helper Function (`_format_verdict_hover()`)**
- Builds hover HTML with verdict name, count, and name list
- Truncates at max_display (default 10) with "...and N more"
- Empty list shows verdict + count only

### Implementation Details

**Lazy Imports:**
- Follows existing pattern: plotly.graph_objects imported in function scope
- VERDICT_COLORS imported at top of function or via TYPE_CHECKING guard

**Empty Data Handling:**
- All functions return empty strings or empty DataFrames gracefully
- Matplotlib backend renders "No data available" text when empty

**Static Report Mode:**
- Both Plotly functions check static_report flag
- Conditionally invoke matplotlib backend via lazy import
- Returns tuple (html_path, png_path_or_None)

## Testing

**Unit Tests Added: 12 new tests**

**TestAggregateVerdicts (4 tests):**
- Verdict string mapping to short labels
- No Result detection for missing sample-plasmid pairs
- Column order preservation
- Empty verdict data handling

**TestFormatVerdictHover (3 tests):**
- Small list (≤10) shows all names
- Large list (>10) truncates with "...and N more"
- Empty list shows verdict + count only

**TestPlotVerdictBarByPlasmid (3 tests):**
- HTML file creation with Plotly
- Empty data returns empty strings
- static_report=True with matplotlib backend creates PNG

**TestPlotVerdictBarBySample (3 tests):**
- HTML file creation with Plotly
- Empty data returns empty strings
- static_report=True with matplotlib backend creates PNG

**TestMatplotlibBackend (2 tests in test_matplotlib_backend.py):**
- generate_verdict_bar_matplotlib() creates PNG from aggregated DataFrame
- Empty DataFrame handled gracefully with "No data available" text

**Test Results:**
- 412 unit tests passing (up from 397 at Phase 15)
- All tests marked with @pytest.mark.unit
- make ci-check passes (lint + format + typecheck + tests)

## Deviations from Plan

None - plan executed exactly as written.

## Integration Points

**Exports for Plan 16-02:**
- `aggregate_verdicts(summary_df) -> (full_df, by_plasmid, by_sample)`
- `plot_verdict_bar_by_plasmid(summary_df, output_dir, static_report, plot_backend) -> (html, png)`
- `plot_verdict_bar_by_sample(summary_df, output_dir, static_report, plot_backend) -> (html, png)`
- `generate_verdict_bar_matplotlib(aggregated_df, output_path, title, xlabel) -> None`
- `VERDICT_COLORS` dictionary from plasmicheck.scripts.plotting.colors

**Inputs Expected:**
- summary_df with columns: Sample, Plasmid, Category, Value
- Category=="Verdict" rows with full verdict strings
- All unique sample/plasmid combinations present in summary_df for cross-product

**Outputs Produced:**
- HTML files: plots/verdict_bar_by_plasmid.html, plots/verdict_bar_by_sample.html
- PNG files (if static_report=True): plots/verdict_bar_by_plasmid.png, plots/verdict_bar_by_sample.png

## Next Phase Readiness

**Blockers:** None

**For Plan 16-02 (Template Integration):**
- All plotting functions are public and importable
- Return signatures match expected tuple format (html_path, png_path_or_None)
- Empty data handling prevents crashes
- static_report flag and plot_backend parameter ready for template consumption

**Known Limitations:**
- No current limitation known

## Quality Metrics

- **Type Safety:** mypy strict mode passes, all function signatures typed
- **Linting:** ruff checks pass, no violations
- **Formatting:** ruff format applied consistently
- **Test Coverage:** 12 new unit tests covering aggregation, plotting, and edge cases
- **Documentation:** All functions have complete docstrings with Args/Returns

## Key Learnings

**1. Cross-Product No Result Detection:**
Building the full cross-product of samples × plasmids before left-joining verdict data ensures missing combinations are explicitly marked as "No Result" rather than silently omitted. This is critical for accurate visualization when some samples fail alignment or are skipped.

**2. Hover Data Pre-Aggregation:**
Building the hover_data dict (verdict → group → list of names) from full_verdict_df BEFORE aggregation preserves the raw name lists needed for rich hover tooltips. Aggregation destroys this detail, so pre-compute and store.

**3. Matplotlib Stacking Pattern:**
Using iterative ax.bar() with cumulative bottom works more reliably than stacked bar mode for categorical data. Initialize bottom=np.zeros, then add each layer's counts to bottom for next iteration.

**4. Color Consistency:**
Deriving VERDICT_COLORS from HEATMAP_COLORS ensures consistent color mapping across all verdict visualizations (heatmaps, bar charts, future additions). Single source of truth prevents divergence.

---

**Summary:** Phase 16 Plan 01 complete. Verdict bar chart plotting functions ready for template integration in Plan 16-02.
