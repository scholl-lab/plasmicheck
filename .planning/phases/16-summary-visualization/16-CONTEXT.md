# Phase 16: Summary Visualization - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Multi-sample summary reports gain verdict bar charts — one grouped by plasmid, one grouped by sample — with Plotly interactive and Matplotlib PNG fallback, integrated into the existing summary report template and export packages. No new metrics or data sources; this visualizes existing verdict data in a new way.

</domain>

<decisions>
## Implementation Decisions

### Chart layout & grouping
- Stacked bars: one bar per plasmid (or sample), with colored segments for each verdict
- Vertical bar orientation (names on x-axis, counts on y-axis)
- Y-axis shows absolute counts (not percentages)
- Colors must match existing verdict colors used elsewhere in the report (consistency with heatmaps/tables)

### Interactivity & hover
- Hover shows: verdict name + count + list of sample/plasmid names
- Name list capped at ~10 entries with "...and N more" for larger sets
- No click actions — charts are for overview, details are in tables
- Standard Plotly toolbar (zoom, pan, reset, download PNG) — matching existing report charts

### Template placement
- Charts appear at the **top of the summary report**, before heatmaps and tables
- Section heading: "Verdict Overview"
- Two charts displayed **side by side** (each taking half width)
- Each chart has its own title ("Verdicts by Plasmid", "Verdicts by Sample")

### Edge cases & scaling
- Always show both charts, even with 1 sample or 1 plasmid (single bar is fine)
- For many items (30+): chart grows horizontally with scrollbar if needed (don't compress bars)
- Long names: rotate x-axis labels 45°
- Missing data: add a "No Result" category (4th segment, gray) for sample-plasmid pairs with no verdict

### Claude's Discretion
- Exact chart dimensions and aspect ratio
- Matplotlib fallback styling details
- Scrollbar implementation approach for wide charts
- Gray color shade for "No Result" segment

</decisions>

<specifics>
## Specific Ideas

- Keep the same color scheme as existing report sections for verdict colors (Contaminated/Clean/Unclear)
- "Verdict Overview" section should feel like a dashboard overview at the top of the report

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 16-summary-visualization*
*Context gathered: 2026-02-16*
