# Requirements: PlasmiCheck v0.32.0 Performance Optimization

**Defined:** 2026-02-14
**Core Value:** Accurately detect plasmid contamination through comparative alignment scoring — the contamination verdict must be reliable

## v1 Requirements

Requirements for this milestone. Each maps to roadmap phases.

### Report Optimization

- [ ] **REPT-01**: User can run pipeline without generating static PNG reports (default behavior — no `--static-report` flag)
- [ ] **REPT-02**: User can opt into static PNG report generation with `--static-report` CLI flag
- [ ] **REPT-03**: Interactive HTML reports use shared `plotly.min.js` file (directory mode) instead of embedding 4.6 MB per plot
- [ ] **REPT-04**: User can choose plotly.js inclusion mode via config (`cdn`, `directory`, `embedded`)
- [ ] **REPT-05**: Kaleido uses `start_sync_server()` initialization for faster PNG export when static reports are requested
- [ ] **REPT-06**: Report-related imports (pandas, plotly, jinja2) are lazy-loaded inside functions, not at module level

### Alignment Optimization

- [ ] **ALGN-01**: Plasmid and human alignments run concurrently via ThreadPoolExecutor (~1.8x speedup)
- [ ] **ALGN-02**: Pipeline auto-detects CPU count with cgroup/SLURM awareness for thread tuning
- [ ] **ALGN-03**: User can override thread count with `--threads` CLI flag
- [ ] **ALGN-04**: All `samtools sort` commands use `-m 2G` memory flag for large BAM performance

### Comparison Optimization

- [ ] **COMP-01**: BAM name grouping uses `samtools collate` instead of `samtools sort -n` (30-50% faster)
- [ ] **COMP-02**: Supplementary alignment ordering handled explicitly after collate

### Architectural Cleanup

- [ ] **ARCH-01**: Human reference indexing hoisted out of combination loop in pipeline
- [ ] **ARCH-02**: PipelinePlan tracks which indexes are already built, skipping redundant lock checks
- [ ] **ARCH-03**: User can generate static plots via matplotlib backend (`--plot-backend matplotlib`) without Kaleido dependency

### Testing & Validation

- [ ] **TEST-01**: Regression test suite verifying optimization outputs match pre-optimization baseline (contamination ratios, read assignments)
- [ ] **TEST-02**: Performance benchmark comparing v0.31.0 vs v0.32.0 on synthetic dataset
- [ ] **TEST-03**: Air-gapped environment test (Docker with no network) for directory-mode reports

## v2 Requirements

Deferred to future milestones. Tracked but not in current roadmap.

### Advanced Performance

- **PERF-01**: Batch Kaleido export with `write_images()` for multiple reports
- **PERF-02**: Parallel report generation across multiple sample-plasmid combinations
- **PERF-03**: `lazy_loader` library integration (requires Python 3.11.9+)

### Extended Report Options

- **RPTE-01**: `--profile` flag showing per-step timing breakdown
- **RPTE-02**: Plotly.js version pinning in directory mode for reproducible reports

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Kaleido v0.2.1 downgrade | Incompatible with Plotly 6.x, security vulnerabilities |
| orca for static export | Deprecated by Plotly team |
| Selenium/Playwright for screenshots | Massive overhead for simple PNG task |
| ProcessPoolExecutor | ThreadPoolExecutor is lighter for subprocess-heavy I/O work |
| Global sort -n removal | Coordinate-sorted BAMs still needed for fetch()-based coverage functions |
| Scientific enhancements (#82, #64, #68, #65) | Deferred to v0.33.0+ milestone |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| TEST-01 | Phase 4 | Pending |
| TEST-02 | Phase 4 | Pending |
| TEST-03 | Phase 4 | Pending |
| REPT-01 | Phase 5 | Pending |
| REPT-02 | Phase 5 | Pending |
| REPT-03 | Phase 5 | Pending |
| REPT-04 | Phase 5 | Pending |
| REPT-05 | Phase 5 | Pending |
| REPT-06 | Phase 5 | Pending |
| ALGN-01 | Phase 6 | Pending |
| ALGN-02 | Phase 6 | Pending |
| ALGN-03 | Phase 6 | Pending |
| ALGN-04 | Phase 6 | Pending |
| COMP-01 | Phase 7 | Pending |
| COMP-02 | Phase 7 | Pending |
| ARCH-01 | Phase 7 | Pending |
| ARCH-02 | Phase 7 | Pending |
| ARCH-03 | Phase 7 | Pending |

**Coverage:**
- v1 requirements: 18 total
- Mapped to phases: 18 (100%)
- Unmapped: 0

---
*Requirements defined: 2026-02-14*
*Last updated: 2026-02-14 after roadmap creation*
