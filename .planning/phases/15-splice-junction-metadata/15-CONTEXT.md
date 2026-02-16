# Phase 15: Splice Junction Metadata - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Reports display splice junction and exon counts extracted from spliced alignments. The spliced alignment maps the plasmid's cDNA insert against the human genome using `minimap2 -ax splice`; N operations in the CIGAR represent introns that were spliced out in the original mRNA. This phase extracts that structural metadata and surfaces it in existing outputs. No changes to the alignment logic itself.

</domain>

<decisions>
## Implementation Decisions

### Data extraction logic
- Count N operations from **primary alignment only** (flag 0 or 16), not supplementary alignments
- SpliceJunctionCount = number of N operations in primary alignment CIGAR
- ExonCount = SpliceJunctionCount + 1
- Store individual exon lengths (M operation lengths) and intron lengths (N operation lengths)
- Use **NA** when no alignment found (unmapped plasmid) — distinguishes "no alignment" from "single-exon gene with 0 junctions"
- Use **0** for single-exon genes (measured zero junctions)
- Lengths reported in **genomic coordinate order** (left-to-right on reference, as they appear in CIGAR) — this is the SAM/GTF/GFF3 standard
- Also store **strand** (+ or -) from primary alignment flag 16

### Output field format
- Keep existing human-readable `key: value` format in cDNA_positions.txt (consistent with current 3 lines)
- New fields added after existing lines:
  ```
  cDNA start position in plasmid: <start>
  cDNA end position in plasmid: <end>
  INSERT_REGION: (<start>, <end>)
  SpliceJunctionCount: <count or NA>
  ExonCount: <count or NA>
  ExonLengths: <comma-separated or NA>
  IntronLengths: <comma-separated or NA>
  Strand: <+ or - or NA>
  Chromosome: <chr or NA>
  GenomicStart: <pos or NA>
  GenomicEnd: <pos or NA>
  ```
- Descriptive key names (SpliceJunctionCount, not Junctions)
- **Always generate** cDNA_positions.txt even for unmapped plasmids (all-NA values)
- Parser: extend `parse_insert_region()` to return a dataclass (`InsertRegionData`) with all fields; new fields are `Optional` with `None` defaults for backward compatibility with old files
- summary.tsv: include SpliceJunctionCount, ExonCount, ExonLengths, and IntronLengths as Category/Count rows

### Report display
- Per-sample HTML report: add SpliceJunctionCount and ExonCount to existing summary statistics table
- Multi-sample summary report: add SpliceJunctionCount, ExonCount, ExonLengths, and IntronLengths as columns
- HTML tables: **truncate long comma-separated lists** (show first 5 values + '...') with tooltip showing full list on hover
- TSV/Excel exports: include full untruncated values

### Multi-transcript handling
- Junction counting: **primary alignment only** — most reliable mapping, represents true gene structure
- Reference extraction: **keep current behavior** using all alignments — different goal (comprehensive coverage for contamination detection)
- These serve different purposes: junction counting needs accuracy, reference extraction needs completeness

### Claude's Discretion
- Exact dataclass field names and implementation details
- How to filter primary alignment in pysam (flag bitmask approach)
- Tooltip CSS/JS implementation in HTML templates
- Error handling for malformed CIGAR strings

</decisions>

<specifics>
## Specific Ideas

- "These are not really splice in the classical sense but large deletions rather" — the insert is a processed mRNA (cDNA), so N operations represent introns from the original gene
- Real data examples: CACNA1H has ~30+ exons (many junctions), KCNJ5 has 1 junction (2 exons), synthetic data has 0 junctions
- Junction metadata is valuable for understanding ratio bias: multi-exon inserts produce more junction-spanning reads assigned to human, potentially lowering contamination ratio
- Store complete genomic provenance: strand, chromosome, start, end — makes cDNA_positions.txt a complete record of the insert's genomic origin

</specifics>

<deferred>
## Deferred Ideas

- Exon-intron structure visualization (diagram in report) — potential future phase
- Using splice junction count to normalize contamination ratio bias — future analysis enhancement
- Transcript-order reporting option (reverse for minus-strand genes) — not needed now, trivial to derive from stored genomic-order + strand

</deferred>

---

*Phase: 15-splice-junction-metadata*
*Context gathered: 2026-02-16*
