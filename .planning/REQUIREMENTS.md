# Requirements: PlasmiCheck v0.34.0

**Defined:** 2026-02-15
**Core Value:** Accurately detect plasmid contamination through comparative alignment scoring -- the contamination verdict must be reliable.

## v0.34.0 Requirements

Requirements for v0.34.0: Robustness & Infrastructure. Each maps to roadmap phases.

### Format Hardening (#75)

- [x] **FMT-01**: Auto-detect file format via magic bytes to distinguish xDNA vs SnapGene `.dna`
- [x] **FMT-02**: Support SnapGene `.dna` files via BioPython `snapgene` parser
- [x] **FMT-03**: Override xDNA record ID with sanitized filename instead of truncated comment field
- [x] **FMT-04**: Improved error messages for unsupported xDNA versions with conversion suggestions

### Summary Visualization (#51)

- [ ] **VIS-01**: Verdict bar chart grouped by plasmid in summary reports
- [ ] **VIS-02**: Verdict bar chart grouped by sample in summary reports
- [ ] **VIS-03**: Matplotlib fallback for both new verdict bar charts
- [ ] **VIS-04**: New plots integrated in summary HTML template and exports

### Index & Spliced Caching (#39, #43)

- [x] **CACHE-01**: `CacheManager` module with content-addressed cache (`~/.plasmicheck/cache/`)
- [x] **CACHE-02**: Within-run plasmid index deduplication (extend `built_indexes` pattern)
- [x] **CACHE-03**: Cross-run plasmid index caching with manifest (tool versions, input hashes)
- [x] **CACHE-04**: Cross-run spliced alignment caching (BAM, spliced reference, cDNA positions)
- [x] **CACHE-05**: CLI flags: `--cache-dir`, `--no-cache`
- [x] **CACHE-06**: Cache management subcommand: `plasmicheck cache info|clean`

### Output Restructuring (#80)

- [x] **OUT-01**: `OutputPaths` dataclass centralizing all path construction
- [x] **OUT-02**: Subfolder layout: `alignments/`, `references/`, `data/`, `reports/`, `plots/`
- [x] **OUT-03**: Depth-agnostic `read_compare_outputs()` for summary report reader
- [x] **OUT-04**: All pipeline components use `OutputPaths` (no hardcoded path construction)

### Splice Junction Metadata (#68)

- [ ] **SJ-01**: Extract splice junction count and exon count from spliced alignment BAM CIGAR
- [ ] **SJ-02**: Store `SpliceJunctionCount` and `ExonCount` in `cDNA_positions.txt` and `summary.tsv`

## Future Requirements

Deferred to v0.35.0+. Tracked but not in current roadmap.

### Visualization

- **VIS-05**: Violin plot of plasmid scores (per-read data, needs opt-in flag + downsampling)
- **VIS-06**: Facet scatter plot of plasmid vs human scores (per-read data)

### Scientific Enhancements

- **SCI-01**: Normalize contamination ratio by splice junction count (#68 Phase 2)
- **SCI-02**: Calibrate normalized thresholds against known samples (#68 Phase 3)
- **SCI-03**: Correct cutoffs for exon count/CDS length (#50)

### Visualization (continued)

- **VIS-07**: Splice junction metadata display in single-sample and summary reports

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| ML classification (#44) | Exploratory, current scoring works |
| Interactive CLI mode (#25) | Nice-to-have, not blocking |
| Anonymize BAM outputs (#5) | Privacy feature, lower priority |
| IGV session generation (#49) | Power user feature |
| Benchmarking dataset creation (#72) | Separate effort |
| Splice junction normalization (#68 Phases 2-3) | Needs threshold recalibration, deferred |
| Violin + scatter plots (#51 plots 3-4) | Needs memory/opt-in design decision |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| FMT-01 | Phase 12 | Complete |
| FMT-02 | Phase 12 | Complete |
| FMT-03 | Phase 12 | Complete |
| FMT-04 | Phase 12 | Complete |
| OUT-01 | Phase 13 | Complete |
| OUT-02 | Phase 13 | Complete |
| OUT-03 | Phase 13 | Complete |
| OUT-04 | Phase 13 | Complete |
| CACHE-01 | Phase 14 | Complete |
| CACHE-02 | Phase 14 | Complete |
| CACHE-03 | Phase 14 | Complete |
| CACHE-04 | Phase 14 | Complete |
| CACHE-05 | Phase 14 | Complete |
| CACHE-06 | Phase 14 | Complete |
| SJ-01 | Phase 15 | Pending |
| SJ-02 | Phase 15 | Pending |
| VIS-01 | Phase 16 | Pending |
| VIS-02 | Phase 16 | Pending |
| VIS-03 | Phase 16 | Pending |
| VIS-04 | Phase 16 | Pending |

**Coverage:**
- v0.34.0 requirements: 20 total
- Mapped to phases: 20
- Unmapped: 0

**Verification:** 100% coverage (20/20 requirements mapped to phases)

---
*Requirements defined: 2026-02-15*
*Last updated: 2026-02-16 after Phase 14 completion (CACHE-01..06 complete)*
