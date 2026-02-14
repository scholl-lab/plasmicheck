# Phase 8: Insert-Region-Aware Filtering - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Reads mapping only to shared backbone regions are excluded from the contamination ratio, eliminating false positive calls for unrelated plasmids. New read categories (Backbone_Only, Ambiguous) are tracked and displayed in reports. A config toggle and optional score margin parameter are provided. This phase does NOT add new metrics (Phase 9), resistance gene detection (Phase 10), or summary report integration (Phase 11).

</domain>

<decisions>
## Implementation Decisions

### Default filtering mode
- filter_backbone_only defaults to **true** (ON) -- users get improved accuracy immediately on upgrade
- No runtime upgrade notice -- release notes and documentation are sufficient
- Config toggle lives in config.json only, no CLI flag exposure
- When filter_backbone_only=false: still classify reads as Backbone_Only for visibility, but include them in the ratio (count but don't exclude)

### Report presentation
- Backbone_Only and Ambiguous added as **new rows in the existing** read assignment summary table (same table, not separate section)
- Visual indicator distinguishes included categories (Plasmid, Human) from excluded categories (Tied, Backbone_Only, Ambiguous) -- e.g., dimmed text or "excluded" label
- Existing pie chart / read assignment plot updated to show **all categories** including Backbone_Only and Ambiguous
- New categories use **muted/gray tones** to visually de-emphasize them as excluded reads

### Score margin policy
- Default score_margin = **0 (disabled)** -- backbone filtering is the primary fix; score margin is secondary
- Config.json entry only, minimal documentation -- advanced users find it when they need it
- When margin > 0, Ambiguous is a **single bucket** (no split by leading side)
- Rationale: one behavioral change at a time; no community consensus on margin thresholds; risk of over-filtering legitimate reads

### Missing insert-region handling
- When insert region is unavailable (cDNA_positions.txt missing/empty/malformed): **skip filtering with warning**, fall back to pre-v0.33.0 behavior
- Warning appears in **both terminal log and HTML report** (e.g., "Backbone filtering unavailable for this sample")
- When filtering is skipped, report shows Backbone_Only = 0 and Ambiguous = 0 (consistent table structure)
- Comparison TSV **always includes** Backbone_Only and Ambiguous columns (0 when unavailable) -- consistent schema for downstream parsing

### Claude's Discretion
- Exact implementation of the visual indicator for excluded categories (dimmed text vs footnote vs label styling)
- Specific gray/muted color values for new categories in plots
- Internal code structure for the positional overlap check (insert region vs read alignment position)
- Warning message wording
- How to handle edge cases in cDNA_positions.txt parsing

</decisions>

<specifics>
## Specific Ideas

- Backbone filtering is the scientifically correct approach -- mirrors discrimination strategies from FastQ Screen, ConFindr, ContScout
- Score margin kept at 0 default following bioinformatics convention of shipping with permissive defaults (like samtools MAPQ filter=0)
- Consistent schema matters: always include new columns/rows even when filtering is unavailable, so parsers don't need to handle variable shapes

</specifics>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 08-insert-region-aware-filtering*
*Context gathered: 2026-02-14*
