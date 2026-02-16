# Roadmap: PlasmiCheck

## Milestones

- ✅ **v0.29.0** Bug Fixes & Housekeeping (Phase 1) - shipped 2026-02-13
- ✅ **v0.30.0** Code Quality (Phase 2) - shipped 2026-02-13
- ✅ **v0.31.0** Enhancements & Infrastructure (Phase 3) - shipped 2026-02-13
- ✅ **v0.32.0** Performance Optimization (Phases 4-7) - shipped 2026-02-14
- ✅ **v0.33.0** Scientific & Reporting Enhancements (Phases 8-11) - shipped 2026-02-15
- 🚧 **v0.34.0** Robustness & Infrastructure (Phases 12-16) - in progress

## Phases

<details>
<summary>✅ v0.33.0 Scientific & Reporting Enhancements (Phases 8-11) - SHIPPED 2026-02-15</summary>

### Phase 8: Insert-Region-Aware Filtering
**Goal:** Reads mapping only to shared backbone are excluded from contamination ratio
**Requirements:** FILT-01, FILT-02, FILT-03, FILT-04, FILT-05
**Plans:** 2 plans

Plans:
- [x] 08-01: Core filtering logic, config, pipeline wiring, and unit tests
- [x] 08-02: Report display updates for new read categories

### Phase 9: Coverage Metrics
**Goal:** Quantitative depth, breadth, and uniformity metrics per-region
**Requirements:** COV-01, COV-02, COV-03, COV-04, COV-05
**Plans:** 2 plans

Plans:
- [x] 09-01: Coverage computation engine, config, summary.tsv integration, and unit tests
- [x] 09-02: Coverage Analysis card in per-sample HTML report

### Phase 10: Resistance Gene Detection
**Goal:** Resistance genes automatically detected from GenBank annotations with per-gene coverage
**Requirements:** RGENE-01, RGENE-02, RGENE-03, RGENE-04
**Plans:** 2 plans

Plans:
- [x] 10-01: GenBank preservation, resistance gene detection module, curated patterns, unit tests
- [x] 10-02: Pipeline integration, per-gene coverage, summary.tsv output, integration tests

### Phase 11: Summary Report Integration
**Goal:** All metrics visible in multi-sample summary reports with plots, tables, and exports
**Requirements:** REPT-04, REPT-05, REPT-06, REPT-07, REPT-08, REPT-09, REPT-10
**Plans:** 2 plans

Plans:
- [x] 11-01: New metric plot/table generation functions (Plotly + matplotlib)
- [x] 11-02: Template restructuring, wiring, exports, and tests

</details>

### 🚧 v0.34.0 Robustness & Infrastructure (In Progress)

**Milestone Goal:** Harden file format support, restructure outputs for clarity, add content-addressed caching for batch runs, enrich reports with splice junction metadata, and improve multi-sample visualization.

#### Phase 12: Format Hardening
**Goal:** Users can reliably convert diverse plasmid file formats without manual intervention
**Depends on:** Nothing (first phase of v0.34.0)
**Requirements:** FMT-01, FMT-02, FMT-03, FMT-04
**Success Criteria** (what must be TRUE):
  1. User can convert SnapGene .dna files via `plasmicheck convert`
  2. User receives clear error message for unsupported xDNA versions with conversion instructions
  3. Converted FASTA files use sanitized filename as record ID (not truncated comment field)
  4. Tool automatically distinguishes xDNA vs SnapGene .dna files by magic bytes
**Plans:** 2 plans

Plans:
- [x] 12-01-PLAN.md — Format detection module, SnapGene support, CLI --format flag, sanitized record IDs
- [x] 12-02-PLAN.md — Pipeline integration, quality control update, comprehensive tests

#### Phase 13: Output Restructuring
**Goal:** Pipeline outputs are organized in logical subfolders with predictable paths
**Depends on:** Phase 12
**Requirements:** OUT-01, OUT-02, OUT-03, OUT-04
**Success Criteria** (what must be TRUE):
  1. Pipeline creates outputs in subfolder structure (alignments/, references/, data/, reports/, plots/)
  2. Summary report reader can find comparison files regardless of directory depth
  3. All pipeline components construct paths via OutputPaths dataclass (no hardcoded strings)
  4. Existing summary report TSV/Excel files continue to work with new structure
**Plans:** 2 plans

Plans:
- [x] 13-01-PLAN.md — OutputPaths frozen dataclass, pipeline rewiring, downstream module updates
- [x] 13-02-PLAN.md — Depth-agnostic summary reader and both-layout tests

#### Phase 14: Cache Infrastructure
**Goal:** Batch runs reuse indexes and spliced alignments across samples, dramatically reducing runtime
**Depends on:** Phase 13 (output paths must be stable)
**Requirements:** CACHE-01, CACHE-02, CACHE-03, CACHE-04, CACHE-05, CACHE-06
**Success Criteria** (what must be TRUE):
  1. Identical plasmid files produce same index once per pipeline run (within-run deduplication)
  2. Subsequent pipeline runs reuse cached indexes from ~/.plasmicheck/cache/indexes/
  3. Subsequent pipeline runs reuse cached spliced alignments when human reference and plasmid match
  4. User can disable caching via --no-cache or override cache directory via --cache-dir
  5. User can inspect cache contents and clean old entries via `plasmicheck cache info|clean`
**Plans:** 3 plans

Plans:
- [x] 14-01-PLAN.md — CacheManager module with content-addressed storage, manifest, and unit tests
- [x] 14-02-PLAN.md — Pipeline integration for index and spliced alignment caching
- [x] 14-03-PLAN.md — Cache management CLI subcommand (info, clean)

#### Phase 15: Splice Junction Metadata
**Goal:** Reports display splice junction and exon counts extracted from spliced alignments
**Depends on:** Phase 13 (needs stable output paths)
**Requirements:** SJ-01, SJ-02
**Success Criteria** (what must be TRUE):
  1. cDNA_positions.txt contains SpliceJunctionCount and ExonCount fields
  2. summary.tsv exports include SpliceJunctionCount and ExonCount columns
  3. Splice junction count accurately reflects N operations in spliced alignment CIGAR string
**Plans:** 2 plans

Plans:
- [ ] 15-01-PLAN.md — Junction metadata extraction, cDNA_positions.txt extension, InsertRegionData dataclass, summary.tsv fields
- [ ] 15-02-PLAN.md — Per-sample and multi-sample report display, HTML truncation, TSV/Excel exports

#### Phase 16: Summary Visualization
**Goal:** Multi-sample summary reports include verdict bar charts grouped by plasmid and by sample
**Depends on:** Phase 13 (needs stable output paths)
**Requirements:** VIS-01, VIS-02, VIS-03, VIS-04
**Success Criteria** (what must be TRUE):
  1. Summary report HTML displays verdict bar chart grouped by plasmid
  2. Summary report HTML displays verdict bar chart grouped by sample
  3. Matplotlib fallback PNGs generated for both new verdict bar charts
  4. Exported summary report packages include both new plots
**Plans:** 2 plans (estimated)

Plans:
- [ ] 16-01: Verdict bar chart plotting functions (Plotly + Matplotlib)
- [ ] 16-02: Template integration and export inclusion

## Progress

**Execution Order:**
Phases execute in numeric order: 12 → 13 → 14 → 15 → 16

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 8. Insert-Region-Aware Filtering | v0.33.0 | 2/2 | Complete | 2026-02-14 |
| 9. Coverage Metrics | v0.33.0 | 2/2 | Complete | 2026-02-14 |
| 10. Resistance Gene Detection | v0.33.0 | 2/2 | Complete | 2026-02-15 |
| 11. Summary Report Integration | v0.33.0 | 2/2 | Complete | 2026-02-15 |
| 12. Format Hardening | v0.34.0 | 2/2 | Complete | 2026-02-15 |
| 13. Output Restructuring | v0.34.0 | 2/2 | Complete | 2026-02-15 |
| 14. Cache Infrastructure | v0.34.0 | 3/3 | Complete | 2026-02-16 |
| 15. Splice Junction Metadata | v0.34.0 | 0/2 | Not started | - |
| 16. Summary Visualization | v0.34.0 | 0/2 | Not started | - |

**v0.34.0 Progress:** 7/11 plans complete

---
*Roadmap created: 2026-02-14*
*Last updated: 2026-02-16 after Phase 15 planning*
