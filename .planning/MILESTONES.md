# PlasmiCheck — Milestones

## v0.29.0 — Bug Fixes & Housekeeping (Phase 1)

**Completed:** 2026-02-13
**Phases:** 1

**What shipped:**
- Closed 6 already-resolved issues (#76, #77, #4, #69, #28, #15)
- Fixed race condition in index creation (#86) with filelock
- Fixed MD5sum file duplicates (#63)
- Fixed template file path resolution (#81) with importlib.resources
- Fixed "wahr" locale bug in Excel output (#79)

## v0.30.0 — Code Quality (Phase 2)

**Completed:** 2026-02-13
**Phases:** 2

**What shipped:**
- Config singleton pattern (#87)
- Centralized CLI arguments (#88)

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
