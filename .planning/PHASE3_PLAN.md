# Phase 3 Plan: Enhancements & Infrastructure (v0.31.0)

**Date:** 2026-02-13
**Branch:** `feat/phase3-v0.31.0`
**Issues:** #89, #85, #33, #32, #29, #71
**Status:** **ALL IMPLEMENTED** -- verified with real data

---

## Overview

Six issues spanning refactoring, infrastructure, UX, and testing. All implemented and verified.

**Dependency graph:**
```
#71 (test data) ──────────────────────> enables integration tests for all others
#89 (paired-end) ─┐
#85 (memory opt)  ├── independent refactors
#32 (dry-run)    ─┘
#29 (progress)   ──── depends on #32 (dry-run provides the step structure progress wraps)
#33 (Docker)     ──── do last (packages everything above)
```

**Recommended execution order:**
1. #71 — Test data (unblocks integration testing for everything else)
2. #89 — Paired-end refactor
3. #85 — Memory optimization
4. #32 — Dry-run mode
5. #29 — Progress bar
6. #33 — Docker

---

## Issue #71: Synthetic Test Dataset -- DONE

### Goal
Create a small, fast, version-controlled test dataset that enables integration tests in CI without requiring the 134 MB `tests/data/real/` directory or 26 GB reference genome.

### Status: IMPLEMENTED
- Synthetic data in `tests/data/synthetic/` (~50 KB total)
- `generate_test_data.py` script for reproducible data generation
- Integration tests in `tests/test_integration.py`: convert, index, contaminated pipeline, clean pipeline
- `conftest.py` fixture `synthetic_data_dir` shared across tests
- 125 total tests (unit + integration) all passing

### Approach: Tiny Synthetic FASTA + FASTQ

Create minimal synthetic data that exercises the full pipeline in <10 seconds:

**Files to create:**
```
tests/data/synthetic/
├── README.md                          # documents how data was generated
├── human_ref.fasta                    # ~5 KB — 2 fake chromosomes, ~2000 bp each
├── plasmid.gb                         # ~3 KB — fake plasmid with cDNA insert matching one human region
├── contaminated_R1.fastq              # ~10 KB — 200 reads, mix of human + plasmid origin
├── contaminated_R2.fastq              # ~10 KB — paired reads
├── not_contaminated_R1.fastq          # ~10 KB — 200 reads, all human origin
├── not_contaminated_R2.fastq          # ~10 KB — paired reads
└── expected/                          # expected outputs for regression checks
    ├── contaminated.verdict.txt       # "contaminated"
    └── not_contaminated.verdict.txt   # "not contaminated"
```

**Design principles:**
- Human reference: 2 fake chromosomes (`chr_fake1`, `chr_fake2`) of ~2000 bp random sequence
- Plasmid: a GenBank file with a vector backbone + cDNA insert that is a copy of a region from `chr_fake1` (so spliced alignment finds it)
- Contaminated reads: 60% from plasmid, 40% from human reference (with some mutations for realism)
- Clean reads: 100% from human reference
- Read length: 150 bp paired-end, insert size ~300 bp
- Total size: <50 KB (easily version-controlled)

**Script to create:**
- `tests/data/synthetic/generate_test_data.py` — deterministic generation (seeded random) so data is reproducible but the script documents the logic

**Test files to create:**
- `tests/test_integration.py` — pytest markers `@pytest.mark.integration`:
  - `test_pipeline_contaminated` — runs full pipeline on contaminated sample, checks verdict
  - `test_pipeline_not_contaminated` — runs full pipeline on clean sample, checks verdict
  - `test_pipeline_bam_input` — tests BAM input path (align synthetic FASTQs to produce a BAM first)
  - `test_convert_genbank` — convert test plasmid.gb to FASTA, verify output
  - `test_index_creation` — create indexes for synthetic human ref, verify .mmi and .fai exist

**Fixture in conftest.py:**
```python
@pytest.fixture
def synthetic_data_dir() -> Path:
    return Path(__file__).resolve().parent / "data" / "synthetic"
```

**CI update:** Change test command from `-m "not slow"` to `-m "not slow"` (integration tests will be `@pytest.mark.integration` which is separate from `slow`). Or tag them as just `unit` if they're fast enough (<10s).

### Estimated effort: Medium (2-3 hours)

---

## Issue #89: Paired-End File Handling -- DONE

### Goal
Replace the comma-separated paired-end convention with explicit `-sf1`/`-sf2` arguments.

### Status: IMPLEMENTED (legacy removed entirely)
Since PlasmiCheck is in alpha, the legacy `-sf` comma-separated mode was **removed entirely** instead of being deprecated. This simplifies the codebase and avoids maintaining two code paths.

### What was done
- **`cli.py`**: `-sf` argument removed. `-sf1`/`--sequencing_files_r1` is now required. `-sf2`/`--sequencing_files_r2` optional for paired-end.
- **`run_pipeline.py`**: `sequencing_files` parameter removed from `resolve_sequencing_inputs()`, `build_plan()`, and `run_pipeline()`. All comma-parsing logic removed. `import warnings` removed.
- **`run_pipeline.py`**: Added validation — `run_pipeline()` raises `ValueError` if `sequencing_files_r1` is None.
- **`SequencingInput`** frozen dataclass normalizes all input styles (single file, file list, paired-end).
- **`resolve_sequencing_inputs()`** handles: single files, `.txt` file lists, paired R1+R2 with count validation.
- **Tests**: Legacy tests removed. New tests cover `SequencingInput`, `resolve_sequencing_inputs()` (5 cases), `get_file_list`, `read_file_list`, `PipelineStep`, `build_plan`, `print_plan`.
- **Integration tests**: Updated to use `sequencing_files_r1=` / `sequencing_files_r2=` kwargs.

### Bug fix discovered during testing
- **`.fastq.gz` extension handling**: `align_reads.py:34` used `input_file.endswith(".fastq")` which didn't match `.fastq.gz` files. Fixed by using tuple: `input_file.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))`. Same fix applied to QC check in `run_pipeline.py`.

---

## Issue #85: Memory Optimization in compare_alignments.py -- DONE

### Goal
Replace the two in-memory dictionaries (`plasmid_reads`, `human_reads`) with a streaming sorted-merge approach to achieve constant memory usage regardless of BAM size.

### Status: IMPLEMENTED
- Name-sorted BAM merge using `samtools sort -n` + `_iter_reads_by_name()` + two-pointer merge
- Memory: O(1) instead of O(n) — only a few reads in memory at any time
- Name-sorted BAMs cleaned up after comparison (intermediate files)
- Output format (TSV columns, verdict logic) identical to previous implementation
- Verified with real data: contaminated ratio 5.955, clean ratio 0.047 — correct verdicts

---

## Issue #32: Dry-Run Mode -- DONE

### Goal
Add `--dry-run` / `-n` flag to the `pipeline` subcommand that validates inputs, shows the execution plan, and exits without running any processing.

### Status: IMPLEMENTED
- Plan-execute refactor: `build_plan()` → `print_plan()` → execute
- `PipelineStep` and `PipelinePlan` dataclasses with `total_steps`, `skipped_steps`, `combinations` properties
- `--dry-run` / `-n` CLI flag on `pipeline` subcommand
- Verified with real data dry-run — shows correct step counts and file paths
- Unit tests: `TestPipelineStep`, `TestBuildPlan`, `TestPrintPlan`

---

## Issue #29: Progress Bar -- DONE

### Goal
Add visual progress feedback during pipeline execution using Rich.

### Status: IMPLEMENTED
- `rich>=13.0` added to `pyproject.toml`
- Rich progress bar wraps pipeline execution steps
- `--no_progress` CLI flag to disable
- Auto-disabled on non-interactive terminals (`sys.stderr.isatty()`)
- Integrated with the plan-execute refactor from #32

---

## Issue #33: Docker Containerization -- DONE

### Goal
Provide a Dockerfile for reproducible deployment with all dependencies (Python + minimap2 + samtools).

### Status: IMPLEMENTED
- Multi-stage `Dockerfile` (build → runtime → test stages)
- Builds minimap2 and samtools from source with pinned versions
- `.dockerignore` excludes large/unnecessary files
- `.github/workflows/docker.yml` CI workflow for tagged releases
- OCI labels for image metadata

---

## Commit Plan (actual)

Branch: `feat/phase3-v0.31.0`

```
feat: add synthetic test dataset for integration testing (#71)
feat: paired-end refactor, dry-run mode, and progress bar (#89, #32, #29)
feat: streaming name-sorted BAM merge for O(1) memory comparison (#85)
feat: add Docker containerization with CI workflow (#33)
chore: bump version to 0.31.0
chore: update uv.lock
```

Post-testing changes (same branch):
- Removed legacy `-sf` argument entirely (alpha — no backward compatibility needed)
- Fixed `.fastq.gz` extension handling bug in `align_reads.py` and `run_pipeline.py`

---

## Verification Results

All checks passed:

```
make ci-check                    # 125 tests pass, ruff clean, mypy clean

# Real data tests:
# Contaminated sample: ratio 5.955, verdict "contaminated"
# Clean sample:        ratio 0.047, verdict "not contaminated"

# Dry-run smoke test: correct step plan output
# Old -sf flag: properly rejected with argparse error
```

---

## Risk Assessment (post-implementation)

| Issue | Risk | Outcome |
|-------|------|---------|
| #85 memory refactor | Behavior change in edge cases | Verified identical results on real data — contaminated and clean samples produce correct verdicts |
| #89 legacy removal | Users with existing scripts using `-sf` | Alpha stage — no external users. `-sf` properly rejected with clear argparse error |
| #32 plan-execute refactor | Large refactor of run_pipeline.py | 125 tests pass including integration tests. Dry-run verified with real data |
| #33 Docker | minimap2/samtools version mismatch | Versions pinned via ARG. Test stage validates. CI workflow on tagged releases |
| #29 Rich dependency | New runtime dependency | `--no_progress` flag + auto-disable on non-TTY. No impact on CI/testing |
| `.fastq.gz` bug | Blocked all gzipped FASTQ inputs | Fixed in both `align_reads.py` and `run_pipeline.py` QC check |
