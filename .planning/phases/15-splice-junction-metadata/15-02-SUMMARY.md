---
phase: 15-splice-junction-metadata
plan: 02
subsystem: visualization
tags: [jinja2, html, pandas, report-generation, data-export]

# Dependency graph
requires:
  - phase: 15-01
    provides: Junction metadata in summary.tsv (SpliceJunctionCount, ExonCount, ExonLengths, IntronLengths)
  - phase: 10-resistance-gene-coverage
    provides: Multi-sample summary report structure
  - phase: 09-insert-backbone-coverage
    provides: Coverage metrics display pattern
provides:
  - Splice junction metadata display in per-sample HTML reports
  - Splice junction metadata table in multi-sample summary reports
  - List truncation for long comma-separated values (ExonLengths, IntronLengths)
  - junction_metadata_table.tsv and Junction Metadata Excel sheet
affects: [16-summary-visualization, future-ui-enhancements]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Jinja2 custom filter for list truncation (max 5 items + '...')"
    - "Full values in TSV/Excel, truncated values in HTML with title tooltip"
    - "Bootstrap card collapsible sections for multi-sample reports"
    - "Conditional display (show section only if junction data exists)"

key-files:
  created: []
  modified:
    - plasmicheck/scripts/generate_report.py
    - plasmicheck/scripts/generate_summary_reports.py
    - plasmicheck/templates/report_template.html
    - plasmicheck/templates/summary_template.html

key-decisions:
  - "Hide junction card when SpliceJunctionCount is N/A (unmapped plasmids)"
  - "Truncate comma-separated lists to 5 items in HTML with full value in title attribute"
  - "Add junction categories to extended_categories for automatic inclusion in exports"
  - "Use separate junction_metadata_table.tsv and Junction Metadata Excel sheet"

patterns-established:
  - "List truncation pattern: HTML display truncated, TSV/Excel full values, title attribute for hover"
  - "Conditional section display: {% if field != 'N/A' %}"
  - "Jinja2 filter registration: env.filters['filter_name'] = function"
  - "Pivot table creation for multi-sample metrics: Sample x Plasmid grid"

# Metrics
duration: 7min
completed: 2026-02-16
---

# Phase 15 Plan 02: Splice Junction Metadata Display Summary

**Junction metadata (count, exons, introns) visible in per-sample and multi-sample HTML reports with list truncation and full TSV/Excel export**

## Performance

- **Duration:** 7 min 37 sec
- **Started:** 2026-02-16T10:26:25Z
- **Completed:** 2026-02-16T10:34:03Z
- **Tasks:** 2/2
- **Files modified:** 4

## Accomplishments

- Per-sample HTML reports display Splice Junction Metadata card with SpliceJunctionCount, ExonCount, ExonLengths, IntronLengths
- Multi-sample summary reports include Splice Junction Metadata table (Sample x Plasmid pivot)
- Long comma-separated lists (exon/intron lengths) truncated to 5 items + "..." in HTML with full value in title tooltip
- Junction metadata included in extended_metrics_table.tsv and dedicated junction_metadata_table.tsv
- Junction Metadata sheet added to summary_report_tables.xlsx

## Task Commits

Each task was committed atomically:

1. **Task 1: Per-sample report junction display** - `f158e51` (feat)
   - Extract junction metrics in generate_report.py using _get_metric helper
   - Add truncate_list Jinja2 filter for long comma-separated lists
   - Add Splice Junction Metadata card to report_template.html after Coverage Analysis
   - Register filter with Jinja2 Environment
   - Pass junction_metrics to template rendering

2. **Task 2: Multi-sample summary report junction columns and exports** - `dd92999` (feat)
   - Add junction categories to extended_categories for automatic export
   - Create junction_table_df pivot (Sample x Plasmid) in main()
   - Add _truncate_list_for_html helper for HTML display
   - Update save_tables_as_tsv_and_excel to save junction_metadata_table.tsv and Junction Metadata Excel sheet
   - Convert junction_df to HTML with truncation for ExonLengths/IntronLengths
   - Add Splice Junction Metadata section (Section 6) to summary_template.html
   - Pass junction_table_df to both report generation and export functions

## Files Created/Modified

- **plasmicheck/scripts/generate_report.py** - Extract junction metrics from summary.tsv, add truncate_list Jinja2 filter, pass junction_metrics to template
- **plasmicheck/scripts/generate_summary_reports.py** - Create junction pivot table, add to extended metrics export, generate junction_metadata_table.tsv and Excel sheet
- **plasmicheck/templates/report_template.html** - Add Splice Junction Metadata card with conditional display (hidden when N/A)
- **plasmicheck/templates/summary_template.html** - Add Splice Junction Metadata collapsible section (Section 6) after Resistance Gene Coverage

## Decisions Made

**1. Hide junction card when SpliceJunctionCount is N/A**
- **Context:** Junction metadata is only available for plasmids that align to the human genome. Unmapped plasmids or old data (pre-phase 15) have N/A values.
- **Decision:** Use `{% if junction_metrics and junction_metrics.splice_junction_count != 'N/A' %}` to conditionally display the card.
- **Rationale:** Avoids showing empty or N/A-filled sections to users. Cleaner UI. Aligns with existing pattern for resistance gene coverage.

**2. Truncate to 5 items in HTML, full values in TSV/Excel**
- **Context:** ExonLengths and IntronLengths can be very long comma-separated lists (e.g., 20+ exons for some genes).
- **Decision:** Truncate to 5 items + "..." in HTML display, but preserve full values in title attribute (tooltip) and all TSV/Excel exports.
- **Rationale:** Maintains readability in HTML reports while preserving full data for analysis. Title attribute provides hover access to full values.
- **Implementation:** Jinja2 filter `truncate_list` for templates, Python helper `_truncate_list_for_html` for summary reports.

**3. Add junction categories to extended_categories**
- **Context:** Extended metrics are automatically exported to TSV and Excel.
- **Decision:** Add SpliceJunctionCount, ExonCount, ExonLengths, IntronLengths to extended_categories list.
- **Rationale:** Simplest way to include junction data in extended_metrics_table.tsv and Extended Metrics Excel sheet without duplicating export logic.

**4. Separate junction_metadata_table.tsv file**
- **Context:** Junction metadata is Sample x Plasmid pivot table with 4 columns (SpliceJunctionCount, ExonCount, ExonLengths, IntronLengths).
- **Decision:** Save dedicated junction_metadata_table.tsv and Junction Metadata Excel sheet in addition to extended_metrics_table.tsv.
- **Rationale:** Follows existing pattern for resistance genes (resistance_genes_table.tsv). Easier for users to access junction data directly without filtering extended metrics.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Phase 16 (Summary Visualization) or other visualization work:**
- Junction metadata display complete in both per-sample and multi-sample reports
- Data layer and UI layer connected
- All existing tests pass (397 tests)
- HTML templates follow established patterns (pc-card for per-sample, Bootstrap card for multi-sample)

**Data flow established:**
- summary.tsv (generated in compare_alignments.py) → generate_report.py → report_template.html (per-sample)
- summary.tsv → generate_summary_reports.py → pivot table → summary_template.html (multi-sample)
- Full values preserved in TSV/Excel, truncated values in HTML

**UI patterns:**
- Conditional display for optional data (junction metadata, resistance genes)
- List truncation for long comma-separated values
- Title attribute for full value access on hover
- Consistent section ordering (Coverage Analysis → Junction Metadata → Run Details)

**Known limitations:**
- Junction metadata only shown when available (N/A for unmapped plasmids or old data)
- Truncation threshold (5 items) is hardcoded, not configurable via config.json
- No visualization (charts/plots) for junction distribution across samples (potential future enhancement)

---
*Phase: 15-splice-junction-metadata*
*Completed: 2026-02-16*
