# Phase 11: Summary Report Integration - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Surface all v0.33.0 metrics (CoverageOutsideINSERT, MismatchesNearINSERT, coverage depth/breadth, resistance genes) in multi-sample summary reports with plots, tables, and exports. Both Plotly (interactive) and matplotlib (static) backends must render all new visualizations. Pre-v0.33.0 data gracefully degrades (missing sections omitted without error).

</domain>

<decisions>
## Implementation Decisions

### Report organization
- Group sections by metric domain, not visualization type
- Section ordering: Contamination (existing) → Coverage → Mismatches → Resistance Genes
- Tables placed inside their parent domain section (verdicts/ratios/p-values under Contamination, resistance gene table under Resistance Genes)
- Collapsible sections — all open by default, users can collapse what they don't need

### Heatmap & chart design
- New coverage heatmaps use the same 3-color scheme as existing ratio heatmap (blue=low/good, red=high/bad) for visual consistency
- Coverage depth heatmap uses linear scale with a sensible cap (values above cap shown as max color)
- CoverageOutsideINSERT gets both heatmap and boxplot (as specified in REPT-04)
- New coverage metrics (depth, breadth) get heatmaps only — no boxplots
- MismatchesNearINSERT visualized as a heatmap (sample × plasmid grid, color = mismatch count)

### Resistance gene presentation
- Flat table: one row per gene per sample-plasmid pair
- Columns: Sample | Plasmid | Gene Name | Mean Depth | Breadth (%) — core metrics only, no coordinates
- When no resistance genes detected for a plasmid, omit silently (no placeholder rows)
- Table only — no additional visual element (no presence/absence heatmap)

### Export structure
- Coverage and mismatch metrics: single combined TSV file (extended_metrics_table.tsv), row-based format matching existing summary.tsv pattern
- Resistance genes: separate TSV file (resistance_genes_table.tsv) due to variable row count
- Single Excel workbook: existing 4 sheets + new sheets for coverage metrics and resistance genes
- Each TSV file corresponds to a sheet in the workbook

### Claude's Discretion
- Linear scale cap value for depth heatmap (pick based on typical PlasmiCheck data ranges)
- Exact collapsible section implementation (Bootstrap accordion or custom)
- Backward compatibility approach for pre-v0.33.0 data (section omission strategy)
- Matplotlib fallback styling details for all new plots

</decisions>

<specifics>
## Specific Ideas

- Follow existing pattern: Plotly figures written as HTML with `include_plotlyjs=False`, optionally converted to PNG for static mode
- Existing summary template uses Bootstrap + DataTables JS — new sections should match
- Current report has ratio heatmap + boxplot + 4 tables; this phase adds ~4 new heatmaps + 1 table

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 11-summary-report-integration*
*Context gathered: 2026-02-15*
