# Phase 9: Coverage Metrics - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Compute quantitative depth, breadth, and uniformity metrics per-region (insert vs backbone) for every sample-plasmid pair. Write metrics to summary.tsv and display in per-sample HTML reports. Multi-sample summary visualization is Phase 11.

</domain>

<decisions>
## Implementation Decisions

### Summary.tsv row layout
- Row names use CamelCase to match existing convention (CoverageOutsideINSERT, MismatchesNearINSERT)
- New rows grouped by region: all insert metrics together, then all backbone metrics together, appended after existing rows
- Numeric precision: 2 decimal places for all coverage values
- Breadth values expressed as fractions (0.00–1.00), not percentages

### Single-sample report display
- Coverage metrics appear in the per-sample HTML report (not just summary.tsv)
- Display as a metrics table (no coverage depth plot — that's visualization, Phase 11 territory)
- Placed in a dedicated "Coverage Analysis" section/card, consistent with Phase 8's card-based layout
- Coverage section always visible, even when all values are zero or data is unavailable

### Threshold configuration
- The 5x breadth threshold is configurable via config.json
- Support multiple thresholds as a list (default: [5]). Each threshold produces a breadth row per region
- CV (coefficient of variation) computed for both insert and backbone regions (extends COV-04)
- Metrics computed for insert and backbone regions only — no whole-plasmid aggregate

### Missing data behavior
- When insert region info is unavailable (no cDNA_positions.txt): compute metrics for whole plasmid as one region (fallback), with adjusted labels
- When no reads align: all coverage metrics set to explicit 0.00 (not N/A)
- HTML report shows zeros directly when coverage is zero — no extra note
- HTML report shows a warning when using whole-plasmid fallback: "Insert region not defined — metrics computed for whole plasmid"

### Claude's Discretion
- Exact CamelCase row names (e.g., MeanDepthInsert vs InsertMeanDepth)
- Table styling within the Coverage Analysis card
- Internal architecture for coverage computation functions
- How to structure the config.json entry for breadth thresholds

</decisions>

<specifics>
## Specific Ideas

- Coverage metrics should feel like a natural extension of the existing summary.tsv format — a downstream parser consuming old rows should still work unchanged
- Phase 8 established the insert/backbone region concept and `read_overlaps_insert()` — coverage metrics should build on that same infrastructure
- The HTML coverage card should match the professional redesign from Phase 8 (Bootstrap 5.3 cards, consistent styling)

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 09-coverage-metrics*
*Context gathered: 2026-02-15*
