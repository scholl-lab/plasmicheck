# Phase 6: Alignment Optimization - Research

**Researched:** 2026-02-14
**Domain:** Thread management, CPU detection, subprocess optimization for bioinformatics pipelines
**Confidence:** HIGH

## Summary

This phase optimizes alignment execution through automatic CPU detection and proper thread allocation. Based on user decisions in the CONTEXT document, alignments will run **sequentially** (not concurrently) to avoid memory pressure from multiple minimap2 instances. Performance gains come from: (1) auto-detecting available CPUs across Docker/SLURM/bare-metal environments, (2) properly allocating threads between minimap2 and samtools within each alignment, and (3) adding the `-m 2G` memory flag to samtools sort for large dataset performance.

The standard approach uses a detection chain (CLI --threads > SLURM env vars > cgroup limits > os.cpu_count()) with cgroup v1/v2 awareness for containerized environments. Thread allocation follows the pattern of giving most threads to minimap2 (the bottleneck) with 2-4 to samtools sort. Regression test compatibility requires identical read assignment results, not identical BAM byte sequences.

**Primary recommendation:** Implement sequential alignment with full thread allocation to each step, auto-detect CPUs via detection chain with cgroup awareness, allocate 80-90% of threads to minimap2 and remainder to samtools, enforce 2-thread minimum and 16-thread maximum cap.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| os.cpu_count() | Python 3.10+ stdlib | Bare-metal CPU detection | Built-in, cross-platform, baseline fallback |
| /sys/fs/cgroup/* | Linux kernel | Container CPU limit detection | Standard cgroup v1/v2 interface for Docker/Kubernetes |
| subprocess.run() | Python 3.10+ stdlib | External tool execution | Built-in, synchronous, exception-raising |
| concurrent.futures.ThreadPoolExecutor | Python 3.10+ stdlib | Thread pool management | Built-in, simpler than threading.Thread for this use case |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| psutil | 7.2.2+ (2026) | CPU affinity detection | Optional enhancement for CPU affinity awareness |
| pathlib.Path | Python 3.10+ stdlib | File path operations | Reading cgroup files |
| logging | Python 3.10+ stdlib | Thread detection transparency | Always log detection source and thread count |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| ThreadPoolExecutor | Simple sequential calls | No benefit for sequential execution model; ThreadPoolExecutor adds unnecessary complexity |
| os.cpu_count() | psutil.cpu_count(logical=False) | psutil is external dependency; physical core count not needed for I/O-bound subprocess work |
| Manual cgroup parsing | Python 3.13 -Xcpu_count | Requires Python 3.13+; manual parsing works on 3.10+ |

**Installation:**
```bash
# No additional packages required — uses stdlib only
# Optional psutil for advanced CPU affinity:
pip install psutil>=7.2.2
```

## Architecture Patterns

### Recommended Project Structure
```
plasmicheck/
├── scripts/
│   ├── run_pipeline.py       # Add detect_cpu_count() function here
│   ├── align_reads.py         # Accept thread counts as parameters
│   └── utils.py               # Or add detect_cpu_count() here as shared utility
├── config.json                # Add "alignment": {"samtools_sort_memory": "2G"}
└── config.py                  # Singleton config loader (already exists)
```

### Pattern 1: CPU Detection Chain
**What:** Priority-ordered detection across environments
**When to use:** At pipeline initialization, before any alignment work
**Example:**
```python
# Source: Manual implementation based on cgroup detection patterns
import os
from pathlib import Path

def detect_cpu_count() -> tuple[int, str]:
    """Detect available CPU count with environment awareness.

    Returns:
        (cpu_count, source) where source is detection method for logging
    """
    # Priority 1: SLURM environment (HPC)
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        try:
            count = int(slurm_cpus)
            return (count, f"SLURM_CPUS_PER_TASK={count}")
        except ValueError:
            pass

    # Priority 2: cgroup v2 (modern Docker/Kubernetes)
    cgroup_v2 = Path("/sys/fs/cgroup/cpu.max")
    if cgroup_v2.exists():
        try:
            content = cgroup_v2.read_text().strip()
            parts = content.split()
            if len(parts) == 2 and parts[0] != "max":
                quota = int(parts[0])  # microseconds
                period = int(parts[1])  # microseconds
                count = max(1, quota // period)
                return (count, f"cgroup v2 cpu.max ({quota}/{period})")
        except (ValueError, IOError):
            pass

    # Priority 3: cgroup v1 (older Docker)
    cgroup_v1_quota = Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
    cgroup_v1_period = Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
    if cgroup_v1_quota.exists() and cgroup_v1_period.exists():
        try:
            quota = int(cgroup_v1_quota.read_text().strip())
            period = int(cgroup_v1_period.read_text().strip())
            if quota > 0:
                count = max(1, quota // period)
                return (count, f"cgroup v1 cpu.cfs ({quota}/{period})")
        except (ValueError, IOError):
            pass

    # Priority 4: os.cpu_count() fallback (bare metal)
    count = os.cpu_count()
    if count:
        return (count, f"os.cpu_count()={count}")

    # Final fallback: 4 threads with warning
    return (4, "fallback (detection failed)")
```

### Pattern 2: Thread Allocation with Bounds
**What:** Apply min/max caps and split threads between tools
**When to use:** After detection, before passing to alignment functions
**Example:**
```python
# Source: Based on user CONTEXT decisions and minimap2 scaling research
def allocate_threads(detected_cpus: int) -> tuple[int, int]:
    """Allocate threads to minimap2 and samtools.

    Args:
        detected_cpus: Raw CPU count from detection chain

    Returns:
        (minimap2_threads, samtools_threads)
    """
    # Apply bounds: min 2, max 16
    total = max(2, min(detected_cpus, 16))

    # Allocate 80% to minimap2 (the bottleneck), rest to samtools
    # Ensure samtools gets at least 2, max 4
    minimap2_threads = max(2, int(total * 0.8))
    samtools_threads = max(2, min(4, total - minimap2_threads))

    # Adjust if total would exceed budget
    if minimap2_threads + samtools_threads > total:
        samtools_threads = max(2, total - minimap2_threads)

    return (minimap2_threads, samtools_threads)
```

### Pattern 3: Sequential Alignment with Full Thread Allocation
**What:** Run plasmid then human alignment, each with full thread budget
**When to use:** In run_pipeline.py alignment step (replaces current back-to-back calls)
**Example:**
```python
# Source: User CONTEXT decision for sequential execution
# Current code (run_pipeline.py lines 559-560):
# align_reads(plasmid_index, sequencing_file, plasmid_bam, "plasmid", fastq2)
# align_reads(spliced_index, sequencing_file, spliced_human_bam, "human", fastq2)

# Optimized sequential execution with explicit thread logging:
logging.info(f"Starting alignment phase (minimap2: {mm2_threads}t, samtools: {sam_threads}t)")

# Plasmid alignment with full thread allocation
align_reads(
    plasmid_index, sequencing_file, plasmid_bam, "plasmid", fastq2,
    minimap2_threads=mm2_threads, samtools_threads=sam_threads
)

# Only proceed to human alignment if plasmid succeeded
# (align_reads raises exception on failure via run_command)
align_reads(
    spliced_index, sequencing_file, spliced_human_bam, "human", fastq2,
    minimap2_threads=mm2_threads, samtools_threads=sam_threads
)
```

### Pattern 4: CLI Override of Auto-Detection
**What:** Allow --threads flag to override all detection logic
**When to use:** In cli.py pipeline subcommand argument parsing
**Example:**
```python
# cli.py — add to pipeline subcommand:
parser.add_argument(
    "--threads",
    type=int,
    default=None,
    help="Override auto-detected CPU count (default: auto-detect via SLURM/cgroup/os.cpu_count())"
)

# In run_pipeline.py:
def run_pipeline(..., threads: int | None = None, ...):
    if threads is not None:
        total_threads = threads
        source = f"CLI --threads={threads}"
    else:
        total_threads, source = detect_cpu_count()

    logging.info(f"Using {total_threads} threads (detected from {source})")
    mm2_threads, sam_threads = allocate_threads(total_threads)
```

### Anti-Patterns to Avoid
- **Concurrent minimap2 execution:** Two simultaneous minimap2 processes cause memory pressure and potential thrashing (60-70% of minimap2 runtime is sequential chaining work that doesn't parallelize well)
- **Ignoring cgroup limits:** Using os.cpu_count() in containers leads to oversubscription (Docker --cpus=4 on 16-core host reports 16, not 4)
- **Thread splitting between concurrent alignments:** Giving each of two concurrent alignments half the threads is slower than sequential with full threads (due to minimap2's sequential bottleneck)
- **Not logging detection source:** HPC users need to debug thread allocation; always log where the thread count came from

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| CPU detection in containers | Custom Docker API queries | cgroup file parsing | Cgroup is kernel interface, works across Docker/Kubernetes/Podman; Docker API requires socket access |
| Thread pool management | Manual threading.Thread list | concurrent.futures.ThreadPoolExecutor | Handles exceptions, join semantics, and cleanup automatically |
| SLURM CPU detection | Parse scontrol output | SLURM_CPUS_PER_TASK env var | Env var is official interface; scontrol requires parsing unstable output |
| Subprocess execution | os.system() or shell=True everywhere | subprocess.run() with list args | Prevents shell injection, better error handling |

**Key insight:** Container CPU limits are managed by the kernel via cgroups, not the container runtime. Reading cgroup files directly is more portable than querying Docker/Podman/containerd APIs.

## Common Pitfalls

### Pitfall 1: os.cpu_count() Reports Host CPUs in Containers
**What goes wrong:** Pipeline detects 16 CPUs on a Docker container limited to 4, spawns 16 threads, causes resource contention
**Why it happens:** Python's os.cpu_count() reads /sys/devices/system/cpu which reflects the host system, not cgroup limits
**How to avoid:** Always check cgroup files before falling back to os.cpu_count()
**Warning signs:** Users report "pipeline very slow in Docker" or "OOM killed in Kubernetes" despite setting resource limits

### Pitfall 2: cgroup v1 vs v2 File Path Differences
**What goes wrong:** Detection code only checks cgroup v1 paths, fails silently on modern systems using cgroup v2
**Why it happens:** Docker/Kubernetes migrated to cgroup v2 in 2021-2023; different file structure (cpu.max vs cpu.cfs_quota_us/cpu.cfs_period_us)
**How to avoid:** Check both v2 (/sys/fs/cgroup/cpu.max) and v1 (/sys/fs/cgroup/cpu/cpu.cfs_*) paths in detection chain
**Warning signs:** Auto-detection works on developer's local Docker but fails on HPC cluster or CI system

### Pitfall 3: SLURM Environment Variable Confusion
**What goes wrong:** Using SLURM_CPUS_ON_NODE instead of SLURM_CPUS_PER_TASK, over-allocating threads
**Why it happens:** Multiple SLURM variables exist; SLURM_CPUS_ON_NODE is total CPUs on the node, not CPUs allocated to this task
**How to avoid:** Use SLURM_CPUS_PER_TASK which reflects the --cpus-per-task allocation for the job
**Warning signs:** HPC jobs fail with "thread spawn error" or job scheduler kills task for exceeding CPU allocation

### Pitfall 4: Not Enforcing Thread Minimums
**What goes wrong:** Detection returns 1 CPU in single-core container, samtools sort fails with "thread count must be >= 1" (since -@ is additional threads)
**Why it happens:** Some tools interpret thread count differently (-t is total, -@ is additional)
**How to avoid:** Enforce minimum of 2 threads total, ensure at least 2 go to minimap2 and 2 to samtools
**Warning signs:** Intermittent failures with "invalid thread count" errors on small cloud instances

### Pitfall 5: Assuming ThreadPoolExecutor Improves Sequential Code
**What goes wrong:** Developer wraps sequential alignment calls in ThreadPoolExecutor expecting speedup, adds overhead instead
**Why it happens:** Misunderstanding that ThreadPoolExecutor only helps with concurrent tasks
**How to avoid:** For sequential execution model (user's decision), use simple function calls; ThreadPoolExecutor adds complexity without benefit
**Warning signs:** Code complexity increases, no performance gain measured

### Pitfall 6: samtools sort -m Flag Confusion
**What goes wrong:** Setting -m 2G with 8 threads allocates 16GB total (2G per thread), causes OOM
**Why it happens:** The -m flag is per-thread, not total; documentation ambiguous
**How to avoid:** Calculate total memory as (threads × per-thread memory); use config.json to make it tunable
**Warning signs:** Large BAM sorts trigger OOM killer despite "plenty of RAM available"

## Code Examples

Verified patterns from official sources and performance analysis:

### Thread Detection with Logging
```python
# Source: Combined from cgroup detection research and user CONTEXT requirements
import logging
import os
from pathlib import Path

def detect_cpu_count() -> tuple[int, str]:
    """Detect CPU count across environments with transparency."""
    # Check SLURM first (HPC common)
    if "SLURM_CPUS_PER_TASK" in os.environ:
        try:
            count = int(os.environ["SLURM_CPUS_PER_TASK"])
            return (count, "SLURM_CPUS_PER_TASK")
        except ValueError:
            pass

    # Check cgroup v2 (modern containers)
    cpu_max = Path("/sys/fs/cgroup/cpu.max")
    if cpu_max.exists():
        try:
            quota_str, period_str = cpu_max.read_text().split()
            if quota_str != "max":
                count = max(1, int(quota_str) // int(period_str))
                return (count, "cgroup v2")
        except (ValueError, IOError, OSError):
            pass

    # Check cgroup v1 (older Docker)
    quota_file = Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
    period_file = Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
    if quota_file.exists() and period_file.exists():
        try:
            quota = int(quota_file.read_text().strip())
            period = int(period_file.read_text().strip())
            if quota > 0:
                count = max(1, quota // period)
                return (count, "cgroup v1")
        except (ValueError, IOError, OSError):
            pass

    # Fallback to os.cpu_count()
    count = os.cpu_count()
    if count:
        return (count, "os.cpu_count()")

    # Final fallback
    logging.warning("CPU detection failed, using fallback of 4 threads")
    return (4, "fallback")

def get_thread_allocation(cli_threads: int | None = None) -> tuple[int, int, str]:
    """Get minimap2 and samtools thread counts.

    Returns:
        (minimap2_threads, samtools_threads, detection_source)
    """
    if cli_threads is not None:
        total = cli_threads
        source = f"CLI --threads={cli_threads}"
    else:
        total, detection = detect_cpu_count()
        source = detection

    # Apply bounds: min 2, max 16
    bounded = max(2, min(total, 16))

    # Allocate ~80% to minimap2, rest to samtools (max 4)
    mm2 = max(2, int(bounded * 0.8))
    sam = min(4, max(2, bounded - mm2))

    return (mm2, sam, source)
```

### Modified align_reads.py Signature
```python
# Source: Based on current align_reads.py with thread parameterization
def align_reads(
    reference_index: str,
    input_file: str,
    output_bam: str,
    alignment_type: str,
    fastq2: str | None = None,
    minimap2_threads: int = 8,  # NEW: accept as parameter
    samtools_threads: int = 4,  # NEW: accept as parameter
    samtools_sort_memory: str = "2G",  # NEW: configurable memory
) -> None:
    logging.info(
        f"Aligning {input_file} to {reference_index} as {alignment_type} "
        f"(minimap2: {minimap2_threads}t, samtools: {samtools_threads}t, sort mem: {samtools_sort_memory})"
    )

    # ... input validation ...

    if input_file.endswith(".bam"):
        command = (
            f"samtools fasta {input_file} | "
            f"minimap2 -t {minimap2_threads} -ax sr {reference_index} - | "
            f"samtools view -@ {samtools_threads} -h -F 4 - | "
            f"samtools sort -@ {samtools_threads} -m {samtools_sort_memory} -o {output_bam}"
        )
    elif fastq2:
        command = (
            f"minimap2 -t {minimap2_threads} -ax sr {reference_index} {input_file} {fastq2} | "
            f"samtools view -@ {samtools_threads} -h -F 4 - | "
            f"samtools sort -@ {samtools_threads} -m {samtools_sort_memory} -o {output_bam}"
        )
    # ... rest of function
```

### Pipeline Integration with Logging
```python
# Source: User CONTEXT requirement for transparent logging
# In run_pipeline.py:

def run_pipeline(
    human_fasta: str,
    plasmid_files: str,
    output_folder: str = "output",
    # ... existing args ...
    threads: int | None = None,  # NEW: CLI override
    # ... rest of args ...
) -> None:
    # Get thread allocation at pipeline start
    mm2_threads, sam_threads, source = get_thread_allocation(threads)

    # Always log (not just verbose mode)
    logging.info(f"Thread allocation: {mm2_threads + sam_threads} total ({source})")
    logging.info(f"  minimap2: {mm2_threads} threads")
    logging.info(f"  samtools: {sam_threads} threads")

    # Load sort memory from config
    cfg = get_config()
    sort_memory = cfg.get("alignment", {}).get("samtools_sort_memory", "2G")

    # ... quality control, indexing ...

    # Step 5: Align reads (sequential, full thread allocation)
    logging.info("Aligning reads to plasmid and human references (sequential)...")

    align_reads(
        plasmid_index, sequencing_file, plasmid_bam, "plasmid", fastq2,
        minimap2_threads=mm2_threads,
        samtools_threads=sam_threads,
        samtools_sort_memory=sort_memory,
    )

    align_reads(
        spliced_index, sequencing_file, spliced_human_bam, "human", fastq2,
        minimap2_threads=mm2_threads,
        samtools_threads=sam_threads,
        samtools_sort_memory=sort_memory,
    )
```

### Config.json Update
```json
{
  "alignment": {
    "minimap2_threads": 8,
    "samtools_threads": 4,
    "samtools_sort_memory": "2G",
    "fasta_extensions": [".fasta", ".fa", ".fna", ".fsa", ".ffn"]
  }
}
```

**Note:** The config.json values become defaults when --threads is not provided AND detection fails. In normal operation, auto-detection overrides these static values.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Static thread config (minimap2_threads: 8) | Auto-detection with cgroup awareness | 2024-2026 (containerization rise) | Prevents oversubscription in containers |
| Concurrent alignments with split threads | Sequential alignments with full threads | Phase 6 (user decision 2026-02-14) | Avoids memory pressure, simpler error handling |
| samtools sort default memory (768MB) | -m 2G for large datasets | 2020+ (large dataset optimization) | 10-19% speedup on large BAMs |
| os.cpu_count() only | Detection chain (SLURM > cgroup > os) | 2023+ (HPC/cloud adoption) | Respects scheduler/container limits |
| Hardcoded thread counts | CLI --threads override | Standard practice | User control for debugging/optimization |

**Deprecated/outdated:**
- **ProcessPoolExecutor for subprocess orchestration:** ThreadPoolExecutor is sufficient and lighter for I/O-bound subprocess work (GIL doesn't affect subprocesses)
- **shell=True with string commands:** Security risk; prefer list-form args with shell=False (current code already uses shell pipes only where necessary)
- **cgroup v1-only detection:** Modern systems use cgroup v2; must support both

## Open Questions

Things that couldn't be fully resolved:

1. **psutil vs manual cgroup parsing**
   - What we know: psutil 7.2.2 has cpu_affinity() but cpu_count() still doesn't respect cgroup limits automatically
   - What's unclear: Whether psutil.Process().cpu_affinity() returns container-aware CPU list on all platforms
   - Recommendation: Start with manual cgroup parsing (stdlib only, verified approach); psutil can be added later if affinity detection proves valuable

2. **Exact thread split ratio (80/20 vs 75/25 vs 90/10)**
   - What we know: minimap2 is the bottleneck; samtools sort benefits from 2-4 threads; diminishing returns above 4
   - What's unclear: Whether 80/20 is optimal across all dataset sizes and read lengths
   - Recommendation: Implement 80/20 split with samtools cap at 4 threads; make ratio tunable in config.json for future optimization

3. **Thread count impact on memory usage**
   - What we know: minimap2 memory scales with -K batch size more than thread count; samtools sort memory is per-thread (-m flag)
   - What's unclear: Exact memory consumption formula for {threads} × {-m value} with specific dataset sizes
   - Recommendation: Document memory calculation in user docs; log warning if detected threads × 2G > 50% system RAM

4. **Optimal maximum thread cap**
   - What we know: User set cap at 16 based on "diminishing returns" research; some benchmarks show gains up to 24 threads
   - What's unclear: Whether cap should be architecture-dependent (more threads on high-memory systems)
   - Recommendation: Keep 16 as default cap per user decision; make tunable via config.json for advanced users

## Sources

### Primary (HIGH confidence)
- [samtools sort manual](http://www.htslib.org/doc/samtools-sort.html) — official samtools documentation for -@ and -m flags
- [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html) — official minimap2 documentation for -t threading option
- [SLURM environment variables](https://docs.hpc.shef.ac.uk/en/latest/referenceinfo/scheduler/SLURM/SLURM-environment-variables.html) — SLURM_CPUS_PER_TASK and related variables
- [Python subprocess documentation](https://docs.python.org/3/library/subprocess.html) — subprocess.run() best practices
- [Python concurrent.futures documentation](https://docs.python.org/3/library/concurrent.futures.html) — ThreadPoolExecutor usage

### Secondary (MEDIUM confidence)
- [How to get the number of CPU cores inside a container](https://donghao.org/2022/01/20/how-to-get-the-number-of-cpu-cores-inside-a-container/) — cgroup v1/v2 detection patterns (verified approach)
- [Python CPU detection in containers - Issue #36054](https://bugs.python.org/issue36054) — os.cpu_count() cgroup awareness discussion
- [Kubernetes cgroup v2 CPU conversion](https://kubernetes.io/blog/2026/01/30/new-cgroup-v1-to-v2-cpu-conversion-formula/) — cgroup v2 migration context
- [samtools threading guidance - Issue #1931](https://github.com/samtools/samtools/issues/1931) — community discussion on thread allocation
- [minimap2 threading performance](https://www.biostars.org/p/9543993/) — community benchmarks

### Tertiary (LOW confidence)
- [ThreadPoolExecutor vs ProcessPoolExecutor comparison](https://superfastpython.com/threadpoolexecutor-vs-processpoolexecutor/) — general Python concurrency patterns (not bioinformatics-specific)
- [Python threading and subprocesses explained](https://www.infoworld.com/article/2257425/python-threading-and-subprocesses-explained.html) — GIL behavior with subprocess

### Project-Specific
- `/mnt/c/development/scholl-lab/plasmicheck/.planning/PERFORMANCE_ANALYSIS.md` — measured performance data showing alignment takes 0.34s for 200 reads
- `/mnt/c/development/scholl-lab/plasmicheck/plasmicheck/scripts/align_reads.py` — current implementation with hardcoded thread config
- `/mnt/c/development/scholl-lab/plasmicheck/.planning/phases/06-alignment-optimization/06-CONTEXT.md` — user decisions constraining this phase

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all stdlib components, verified cgroup file paths
- Architecture patterns: HIGH — detection chain verified across sources, thread allocation based on tool documentation
- Thread split ratio: MEDIUM — 80/20 split is informed estimate, not benchmarked for this specific pipeline
- Pitfalls: HIGH — based on real-world container/HPC issues documented in Python tracker and SLURM docs
- Memory usage formulas: MEDIUM — samtools docs confirm per-thread memory, but interaction with large datasets not fully profiled

**Research date:** 2026-02-14
**Valid until:** 60 days (2026-04-15) — thread management patterns are stable; cgroup v2 adoption may continue evolving
