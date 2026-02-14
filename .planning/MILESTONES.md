# PlasmiCheck — Milestones

## v0.32.0 Performance Optimization (Shipped: 2026-02-14)

**Delivered:** Pipeline performance optimization delivering 9.5x speedup for small datasets and 1.97x for real-world data through opt-in PNG export, lazy imports, thread auto-detection, and samtools memory tuning.

**Phases completed:** 4-7 (11 plans total)

**Key accomplishments:**

- 9.5x faster pipeline for small datasets (5.5s to 0.577s) by making PNG report generation opt-in
- ~100x faster report step (5.1s to 0.108s) — lazy imports, directory-mode plotly.js, conditional Kaleido
- 1.97x faster alignment on real 3M-read datasets — samtools sort -m 2G delivered 4.8x human alignment speedup
- Auto CPU detection across SLURM/Docker/bare metal with --threads CLI override
- Matplotlib static plot backend — Kaleido-free PNG generation via --plot-backend matplotlib
- Batch resilience with hoisted human indexing, index tracking, and graceful error handling

**Stats:**

- 21 source files created/modified (+42 planning files)
- 3,228 lines added, 237 removed (7,869 total Python LOC)
- 4 phases, 11 plans, 60 commits
- 1 day (2026-02-14), continuing from v0.31.0

**Git range:** `docs(04)` -> `docs(07)`

**Test coverage:** 125 -> 170 tests (+45 new)

**What's next:** Scientific enhancements (v0.33.0+) — improved scoring, additional contamination metrics

---

## v0.31.0 — Enhancements & Infrastructure (Phase 3)

**Completed:** 2026-02-13
**Phases:** 3
**Last phase number:** 3

**What shipped:**
- Synthetic test dataset + integration tests (#71) — 125 tests total
- Paired-end refactor with explicit -sf1/-sf2 (#89) — legacy -sf removed
- O(1) memory streaming BAM comparison (#85)
- Dry-run mode with plan-execute architecture (#32)
- Rich progress bar (#29)
- Docker containerization with CI workflow (#33)
- Bug fix: .fastq.gz extension handling

**Verification:** Contaminated ratio 5.955, clean ratio 0.047 — correct verdicts on real data

## v0.30.0 — Code Quality (Phase 2)

**Completed:** 2026-02-13
**Phases:** 2

**What shipped:**
- Config singleton pattern (#87)
- Centralized CLI arguments (#88)

## v0.29.0 — Bug Fixes & Housekeeping (Phase 1)

**Completed:** 2026-02-13
**Phases:** 1

**What shipped:**
- Closed 6 already-resolved issues (#76, #77, #4, #69, #28, #15)
- Fixed race condition in index creation (#86) with filelock
- Fixed MD5sum file duplicates (#63)
- Fixed template file path resolution (#81) with importlib.resources
- Fixed "wahr" locale bug in Excel output (#79)
