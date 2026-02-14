"""Thread detection and allocation for alignment operations.

Detects available CPUs from SLURM, cgroups (v1/v2), or OS, then allocates
threads between minimap2 and samtools with configurable bounds.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path


def detect_cpu_count() -> tuple[int, str]:
    """Detect available CPU count with fallback chain.

    Detection priority:
    1. SLURM_CPUS_PER_TASK environment variable
    2. cgroup v2 cpu.max file
    3. cgroup v1 cpu.cfs_quota_us + cpu.cfs_period_us
    4. os.cpu_count()
    5. Fallback to 4

    Returns:
        (cpu_count, source_description): CPU count and human-readable source name
    """
    # SLURM environment
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        try:
            cpu_count = int(slurm_cpus)
            if cpu_count > 0:
                return (cpu_count, "SLURM_CPUS_PER_TASK")
        except (ValueError, TypeError):
            pass

    # cgroup v2
    cgroup_v2_path = Path("/sys/fs/cgroup/cpu.max")
    try:
        if cgroup_v2_path.exists():
            content = cgroup_v2_path.read_text().strip()
            parts = content.split()
            if len(parts) == 2 and parts[0] != "max":
                quota = int(parts[0])
                period = int(parts[1])
                if quota > 0 and period > 0:
                    cpu_count = quota // period
                    if cpu_count > 0:
                        return (cpu_count, "cgroup v2")
    except (ValueError, OSError):
        pass

    # cgroup v1
    cgroup_v1_quota_path = Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
    cgroup_v1_period_path = Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
    try:
        if cgroup_v1_quota_path.exists() and cgroup_v1_period_path.exists():
            quota = int(cgroup_v1_quota_path.read_text().strip())
            period = int(cgroup_v1_period_path.read_text().strip())
            if quota > 0 and quota != -1 and period > 0:
                cpu_count = quota // period
                if cpu_count > 0:
                    return (cpu_count, "cgroup v1")
    except (ValueError, OSError):
        pass

    # os.cpu_count()
    try:
        detected = os.cpu_count()
        if detected is not None and detected > 0:
            return (detected, "os.cpu_count()")
    except Exception:
        pass

    # Fallback
    logging.warning(
        "Could not detect CPU count from SLURM, cgroups, or os.cpu_count(). "
        "Falling back to 4 CPUs."
    )
    return (4, "fallback")


def allocate_threads(total_cpus: int) -> tuple[int, int]:
    """Allocate threads between minimap2 and samtools.

    Applies bounds (min 2, max 16), then allocates ~80% to minimap2,
    remainder to samtools (max 4).

    Args:
        total_cpus: Total CPUs available

    Returns:
        (minimap2_threads, samtools_threads): Thread counts for each tool
    """
    # Clamp to bounds
    total_cpus = max(2, min(16, total_cpus))

    # Allocate ~80% to minimap2, at least 2
    minimap2_threads = max(2, int(total_cpus * 0.8))

    # Remainder to samtools, at least 2, max 4
    samtools_threads = max(2, min(4, total_cpus - minimap2_threads))

    # If samtools min requirement pushes us over total, reduce samtools
    if minimap2_threads + samtools_threads > total_cpus:
        samtools_threads = max(2, total_cpus - minimap2_threads)

    return (minimap2_threads, samtools_threads)
