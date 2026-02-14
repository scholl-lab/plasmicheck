# Phase 7: Comparison & Cleanup - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Optimize BAM comparison with faster name grouping (samtools collate), eliminate redundant human reference indexing in batch runs, and add matplotlib as an alternative static plot backend. Interactive HTML reports remain Plotly-based.

</domain>

<decisions>
## Implementation Decisions

### Collate behavior
- Always use samtools collate (not size-dependent threshold)
- Re-sort supplementary alignments within each read group after collate (don't rely on collate ordering)
- Write collated output to temp file (not streaming pipe) for debuggability
- If samtools collate fails or isn't available, fall back to sort -n automatically with a warning logged

### Index deduplication
- Only deduplicate human reference indexing (plasmid indexes are cheap, skip dedup for those)
- Track built indexes by file path matching (not content hash)
- Trust existing indexes — no timestamp checking against reference files
- Log a single INFO line per skipped index: "Skipping index for X (already built)"
- Use upfront planning: PipelinePlan scans all combinations first, builds unique index list, runs all indexing, then processes combinations

### Matplotlib backend
- matplotlib becomes a core dependency (in pyproject.toml requirements, not optional)
- `--plot-backend matplotlib` only affects static PNG generation (`--static-report` required)
- Interactive HTML always uses Plotly.js regardless of `--plot-backend`
- Support all plots: single-sample (score distribution, alignment breakdown) and summary (heatmaps, boxplots)
- Functionally equivalent to Plotly output — same data, same chart types, similar colors; does not need to be pixel-identical

### Batch resilience
- Pipeline continues processing remaining combinations when one fails; failed combinations logged with error details
- Summary report generated for successful combinations only, with clear indication of which were excluded and why
- Log per-combination timing at INFO level (e.g., "Sample A x Plasmid B: 45s")

### Claude's Discretion
- Exact matplotlib style configuration (colors, fonts, figure sizes)
- Temp file cleanup strategy for collate output
- PipelinePlan internal data structure for tracking built indexes
- Error message formatting for failed combinations

</decisions>

<specifics>
## Specific Ideas

- matplotlib should be a regular dependency, not optional — simplifies import handling and user experience
- Upfront index planning phase before combination processing — clear two-phase architecture
- Per-combination timing helps users identify bottlenecks in large batch runs

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 07-comparison-cleanup*
*Context gathered: 2026-02-14*
