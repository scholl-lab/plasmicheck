# Phase 6: Alignment Optimization - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Parallelize alignment execution and optimize threading for 1.5-2x speedup on real-world datasets. Covers thread management, CPU detection, and the `--threads` CLI flag. Does not add new alignment algorithms or change alignment logic.

</domain>

<decisions>
## Implementation Decisions

### Thread defaults & limits
- Default thread count: auto-detect all available CPUs (via detection chain)
- Minimum floor: 2 threads (never go below, even on single-core detection)
- Maximum cap: 16 threads (diminishing returns for minimap2 beyond this)
- Thread split between minimap2 and samtools: Claude's discretion informed by best practices — most threads to minimap2, 2-4 to samtools sort (samtools won't use more than needed)

### Parallel execution model
- Plasmid and human alignments run **sequentially**, not concurrently — two concurrent minimap2 instances risk memory pressure
- Each alignment gets full thread allocation (no splitting)
- Fail immediately if first alignment fails — don't proceed to second
- Batch runs (multiple plasmid x sample combinations): sequential, one combination at a time with full threads
- samtools sort `-m` memory flag: configurable via config.json (default 2G)

### CLI flag behavior
- `--threads` controls total system thread budget (split across minimap2 + samtools, not just minimap2)
- CLI `--threads` overrides config.json thread settings; config.json provides defaults
- `--threads` available on `pipeline` subcommand only (not individual subcommands)
- Always log thread count at pipeline start (not just verbose mode)
- Always log detection source (e.g., "Using 8 threads (detected from SLURM_CPUS_PER_TASK)")
- No `--detect-threads` dry-run mode — keep CLI simple

### Environment detection
- Detection chain priority: CLI `--threads` > SLURM (`SLURM_CPUS_PER_TASK`) > cgroup limits > `os.cpu_count()`
- Support Docker (cgroup), SLURM (HPC), and bare metal environments
- Fallback if all detection fails: 4 threads with warning logged
- Detection source always logged for transparency

### Claude's Discretion
- Exact thread split ratio between minimap2 and samtools sort
- cgroup v1 vs v2 detection implementation
- SLURM environment variable selection (SLURM_CPUS_PER_TASK vs alternatives)
- ThreadPoolExecutor vs simple sequential calls for the alignment orchestration
- Config.json schema for sort memory setting

</decisions>

<specifics>
## Specific Ideas

- Memory pressure from concurrent minimap2 is the primary concern — user explicitly chose sequential over concurrent to avoid this
- Thread logging should always show source for HPC debugging (common support issue)
- The 16-thread cap is based on minimap2 scaling benchmarks showing diminishing returns

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 06-alignment-optimization*
*Context gathered: 2026-02-14*
