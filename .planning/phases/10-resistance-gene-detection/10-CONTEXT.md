# Phase 10: Resistance Gene Detection - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Automatically detect resistance genes and reporter genes from plasmid GenBank annotations, compute per-gene coverage metrics (mean depth, breadth, coordinates), and report them in summary.tsv as an independent contamination signal. xDNA files are already converted to GenBank by the existing pipeline. FASTA-only plasmids cannot be analyzed (no annotations). Display in HTML reports is Phase 11 scope.

</domain>

<decisions>
## Implementation Decisions

### Gene identification scope
- Detect **two tiers** of genes from GenBank CDS features:
  - **Tier 1 — Drug resistance:** Bacterial (AmpR, KanR, ChlR, TetR, bla, etc.) and mammalian (PuroR, HygroR, NeoR, BlasticidinR, ZeocinR, etc.)
  - **Tier 2 — Reporter/fluorescent:** GFP, eGFP, mCherry, RFP, dsRed, luciferase, lacZ, etc.
- **Exclude** affinity tags (His-tag, FLAG, Myc, V5) — weak contamination signal, not worth the noise
- **Auto-exclude insert region:** CDS features overlapping the insert region (from cDNA_positions.txt) are skipped — they're the experimental target, not contamination markers

### Detection strategy
- **Curated regex pattern list** matched against GenBank CDS qualifiers (/label, /product, /gene)
- Case-insensitive, boundary-aware matching
- **Skip ambiguous matches** — only report high-confidence hits, no uncertain/low-confidence results
- **Config-extensible:** Default built-in pattern list with optional config.json override for users to add custom patterns for niche plasmids

### Output format (summary.tsv)
- **One row per gene** for detailed metrics: `ResGene_{GeneName}_MeanDepth`, `ResGene_{GeneName}_Breadth`, `ResGene_{GeneName}_Start`, `ResGene_{GeneName}_End`, `ResGene_{GeneName}_Region` (backbone/insert)
- **Summary row:** `ResistanceGenes` with comma-separated list of detected gene names (or `None` if no matches)
- **Same `ResGene_` prefix** for both resistance and reporter genes — no category distinction in TSV keys
- **When no genes detected (but GenBank available):** `ResistanceGenesAvailable: true`, `ResistanceGenes: None`
- **When no GenBank available (FASTA-only):** `ResistanceGenesAvailable: false` + log warning — reports can show "Gene detection not available (no annotations)"
- **Distinct status** between "looked but found nothing" vs "couldn't look" via the ResistanceGenesAvailable flag

### Gene naming
- **GenBank label verbatim** — use whatever /label or /gene qualifier says (e.g., "AmpR", "FRT Hygromycin")
- **Sanitize for TSV row keys:** spaces/special chars replaced with underscores (e.g., "FRT Hygromycin" → `ResGene_FRT_Hygromycin_MeanDepth`)
- Display names in reports/summary row keep original GenBank label

### Claude's Discretion
- Exact regex patterns for the curated list (must achieve >=90% recall on existing GenBank files)
- Which GenBank qualifiers to check beyond /label, /product, /gene
- Config.json structure for custom pattern extension
- How to handle the intermediate .gb file from xDNA conversion (locate it in the pipeline flow)

</decisions>

<specifics>
## Specific Ideas

- xDNA → GenBank conversion already exists in `convert_xdna_to_gb()` in `convert_plasmidfile_to_fasta.py` — the intermediate .gb file should be available for gene detection
- Coverage metrics infrastructure from Phase 9 (`coverage_metrics.py`, `pysam.count_coverage()`) can be reused for per-gene coverage computation
- Existing GenBank files in `reference/plasmid/` contain genes like AmpR, FRT Hygromycin, c-Myc tag, Flag tag — the curated list must detect resistance/reporter genes while ignoring tags
- Per-gene metrics should follow the same pattern as existing coverage metrics (MeanDepth, Breadth naming conventions)

</specifics>

<deferred>
## Deferred Ideas

- BLAST-based gene detection for unannotated plasmids (against CARD/ResFinder databases) — potential future phase for advanced users
- Resistance gene display in HTML single-sample reports — Phase 11 scope
- Resistance gene display in multi-sample summary reports — Phase 11 scope

</deferred>

---

*Phase: 10-resistance-gene-detection*
*Context gathered: 2026-02-15*
