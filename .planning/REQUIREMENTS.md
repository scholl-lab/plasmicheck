# Requirements: PlasmiCheck v0.33.0

**Defined:** 2026-02-14
**Core Value:** Accurately detect plasmid contamination through comparative alignment scoring — the contamination verdict must be reliable.

## v0.33.0 Requirements

Requirements for v0.33.0: Scientific & Reporting Enhancements. Each maps to roadmap phases.

### Scoring & Filtering (#82)

- [ ] **FILT-01**: Reads mapping only to backbone (outside insert region) are classified as Backbone_Only and excluded from contamination ratio
- [ ] **FILT-02**: Optional score margin parameter requires minimum score difference for confident Plasmid/Human assignment
- [ ] **FILT-03**: New read categories (Backbone_Only, Ambiguous) tracked in comparison output alongside existing Plasmid/Human/Tied
- [ ] **FILT-04**: `filter_backbone_only` config toggle (default: true) with `false` restoring pre-v0.33.0 behavior
- [ ] **FILT-05**: New categories displayed in single-sample HTML reports

### Coverage Metrics (#65)

- [ ] **COV-01**: Per-region mean and median depth computed for insert and backbone regions
- [ ] **COV-02**: Breadth of coverage (fraction of bases with >=1 read) computed per region
- [ ] **COV-03**: Breadth at 5x threshold computed per region
- [ ] **COV-04**: Coverage uniformity (coefficient of variation) computed for backbone region
- [ ] **COV-05**: All new metrics written as additional rows in summary.tsv (backward compatible)

### Resistance Genes (#64)

- [ ] **RGENE-01**: Resistance genes extracted from GenBank annotations via configurable pattern matching
- [ ] **RGENE-02**: Per-gene coverage metrics (mean depth + breadth) computed using pysam
- [ ] **RGENE-03**: Resistance gene metrics written to summary.tsv per gene
- [ ] **RGENE-04**: Detection works with existing plasmid GenBank files in repository (>=90% recall)

### Report Integration (#58 + visualization)

- [ ] **REPT-04**: CoverageOutsideINSERT parsed and displayed in summary reports with heatmap and boxplot
- [ ] **REPT-05**: MismatchesNearINSERT parsed and displayed in summary reports with chart
- [ ] **REPT-06**: Coverage depth heatmap (Sample x Plasmid, backbone mean depth) in summary reports
- [ ] **REPT-07**: Coverage breadth heatmap (backbone breadth %) in summary reports
- [ ] **REPT-08**: Resistance gene coverage table in summary reports
- [ ] **REPT-09**: Matplotlib fallback for all new summary report plots
- [ ] **REPT-10**: TSV and Excel exports include new metric tables

## Future Requirements

Deferred to v0.34.0+. Tracked but not in current roadmap.

### Scientific Enhancements

- **SCI-01**: Normalize by splice junction count (#68)
- **SCI-02**: Correct cutoffs for exon count/CDS length (#50)
- **SCI-03**: Built-in resistance gene FASTA fallback for unannotated plasmids (#64 extension)

### Visualization

- **VIS-01**: Scatter plot overlay on boxplots (#61)
- **VIS-02**: Additional summary report plots (#51)

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| ML classification (#44) | Exploratory, current scoring works |
| Interactive CLI mode (#25) | Nice-to-have, not blocking |
| Anonymize BAM outputs (#5) | Privacy feature, lower priority |
| IGV session generation (#49) | Power user feature |
| Built-in resistance gene FASTA (#64 extension) | Annotation parsing sufficient for v0.33.0; fallback deferred |
| Benchmarking dataset creation (#72) | Separate effort, not tied to scientific enhancements |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| FILT-01 | — | Pending |
| FILT-02 | — | Pending |
| FILT-03 | — | Pending |
| FILT-04 | — | Pending |
| FILT-05 | — | Pending |
| COV-01 | — | Pending |
| COV-02 | — | Pending |
| COV-03 | — | Pending |
| COV-04 | — | Pending |
| COV-05 | — | Pending |
| RGENE-01 | — | Pending |
| RGENE-02 | — | Pending |
| RGENE-03 | — | Pending |
| RGENE-04 | — | Pending |
| REPT-04 | — | Pending |
| REPT-05 | — | Pending |
| REPT-06 | — | Pending |
| REPT-07 | — | Pending |
| REPT-08 | — | Pending |
| REPT-09 | — | Pending |
| REPT-10 | — | Pending |

**Coverage:**
- v0.33.0 requirements: 21 total
- Mapped to phases: 0
- Unmapped: 21 (pending roadmap)

---
*Requirements defined: 2026-02-14*
*Last updated: 2026-02-14 after initial definition*
