"""Unit tests for thread_config module."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from plasmicheck.thread_config import allocate_threads, detect_cpu_count


@pytest.mark.unit
def test_detect_cpu_count_slurm(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test CPU detection from SLURM_CPUS_PER_TASK env var."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    cpu_count, source = detect_cpu_count()
    assert cpu_count == 16
    assert source == "SLURM_CPUS_PER_TASK"


@pytest.mark.unit
def test_detect_cpu_count_cgroup_v2(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Test CPU detection from cgroup v2 cpu.max file."""
    # Remove SLURM env var
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)

    # Create fake cgroup v2 file
    fake_cgroup = tmp_path / "cpu.max"
    fake_cgroup.write_text("200000 100000\n")  # 200000/100000 = 2 CPUs

    # Monkeypatch the Path used in detect_cpu_count
    original_path = Path

    class MockPath(type(Path())):  # type: ignore[misc]
        def __new__(cls, *args, **kwargs):  # type: ignore[no-untyped-def]
            if args and args[0] == "/sys/fs/cgroup/cpu.max":
                return original_path(fake_cgroup)
            return original_path(*args, **kwargs)

    monkeypatch.setattr("plasmicheck.thread_config.Path", MockPath)

    cpu_count, source = detect_cpu_count()
    assert cpu_count == 2
    assert source == "cgroup v2"


@pytest.mark.unit
def test_detect_cpu_count_cgroup_v2_unlimited(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test cgroup v2 with unlimited quota (max) falls through."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)

    # Create fake cgroup v2 file with "max" quota
    fake_cgroup = tmp_path / "cpu.max"
    fake_cgroup.write_text("max 100000\n")

    original_path = Path

    class MockPath(type(Path())):  # type: ignore[misc]
        def __new__(cls, *args, **kwargs):  # type: ignore[no-untyped-def]
            if args and args[0] == "/sys/fs/cgroup/cpu.max":
                return original_path(fake_cgroup)
            return original_path(*args, **kwargs)

    monkeypatch.setattr("plasmicheck.thread_config.Path", MockPath)

    # Should fall through to os.cpu_count or fallback
    _, source = detect_cpu_count()
    assert source in ("os.cpu_count()", "fallback")


@pytest.mark.unit
def test_detect_cpu_count_cgroup_v1(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Test CPU detection from cgroup v1 quota and period files."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)

    # Create fake cgroup v1 files
    fake_quota = tmp_path / "cpu.cfs_quota_us"
    fake_period = tmp_path / "cpu.cfs_period_us"
    fake_quota.write_text("400000\n")  # 400000/100000 = 4 CPUs
    fake_period.write_text("100000\n")

    original_path = Path

    class MockPath(type(Path())):  # type: ignore[misc]
        def __new__(cls, *args, **kwargs):  # type: ignore[no-untyped-def]
            path_str = str(args[0]) if args else ""
            if "/sys/fs/cgroup/cpu.max" in path_str:
                # v2 doesn't exist
                p = original_path(tmp_path / "nonexistent_v2")
                return p
            elif "/cpu.cfs_quota_us" in path_str:
                return original_path(fake_quota)
            elif "/cpu.cfs_period_us" in path_str:
                return original_path(fake_period)
            return original_path(*args, **kwargs)

    monkeypatch.setattr("plasmicheck.thread_config.Path", MockPath)

    cpu_count, source = detect_cpu_count()
    assert cpu_count == 4
    assert source == "cgroup v1"


@pytest.mark.unit
def test_detect_cpu_count_cgroup_v1_unlimited(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Test cgroup v1 with unlimited quota (-1) falls through."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)

    # Create fake cgroup v1 files with -1 quota (unlimited)
    fake_quota = tmp_path / "cpu.cfs_quota_us"
    fake_period = tmp_path / "cpu.cfs_period_us"
    fake_quota.write_text("-1\n")
    fake_period.write_text("100000\n")

    original_path = Path

    class MockPath(type(Path())):  # type: ignore[misc]
        def __new__(cls, *args, **kwargs):  # type: ignore[no-untyped-def]
            path_str = str(args[0]) if args else ""
            if "/sys/fs/cgroup/cpu.max" in path_str:
                p = original_path(tmp_path / "nonexistent_v2")
                return p
            elif "/cpu.cfs_quota_us" in path_str:
                return original_path(fake_quota)
            elif "/cpu.cfs_period_us" in path_str:
                return original_path(fake_period)
            return original_path(*args, **kwargs)

    monkeypatch.setattr("plasmicheck.thread_config.Path", MockPath)

    # Should fall through to os.cpu_count or fallback
    _, source = detect_cpu_count()
    assert source in ("os.cpu_count()", "fallback")


@pytest.mark.unit
def test_detect_cpu_count_os(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Test CPU detection falls back to os.cpu_count()."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)

    # Mock Path to make cgroup files not exist
    original_path = Path

    class MockPath(type(Path())):  # type: ignore[misc]
        def __new__(cls, *args, **kwargs):  # type: ignore[no-untyped-def]
            path_str = str(args[0]) if args else ""
            if "cgroup" in path_str or "cpu.max" in path_str or "cpu.cfs" in path_str:
                p = original_path(tmp_path / "nonexistent")
                return p
            return original_path(*args, **kwargs)

    monkeypatch.setattr("plasmicheck.thread_config.Path", MockPath)

    cpu_count, source = detect_cpu_count()
    # Should be from os.cpu_count()
    assert source == "os.cpu_count()"
    assert cpu_count > 0


@pytest.mark.unit
def test_detect_cpu_count_fallback(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Test CPU detection falls back to 4 when all detection fails."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)

    # Mock Path to make cgroup files not exist
    original_path = Path

    class MockPath(type(Path())):  # type: ignore[misc]
        def __new__(cls, *args, **kwargs):  # type: ignore[no-untyped-def]
            path_str = str(args[0]) if args else ""
            if "cgroup" in path_str or "cpu.max" in path_str or "cpu.cfs" in path_str:
                p = original_path(tmp_path / "nonexistent")
                return p
            return original_path(*args, **kwargs)

    monkeypatch.setattr("plasmicheck.thread_config.Path", MockPath)

    # Mock os.cpu_count to return None
    monkeypatch.setattr(os, "cpu_count", lambda: None)

    cpu_count, source = detect_cpu_count()
    assert cpu_count == 4
    assert source == "fallback"


@pytest.mark.unit
def test_allocate_threads_minimum() -> None:
    """Test thread allocation with minimum CPUs (2)."""
    minimap2, samtools = allocate_threads(2)
    assert minimap2 == 2
    assert samtools == 1


@pytest.mark.unit
def test_allocate_threads_below_minimum() -> None:
    """Test thread allocation with below-minimum CPUs (1) clamps to 2."""
    minimap2, samtools = allocate_threads(1)
    assert minimap2 == 2
    assert samtools == 1


@pytest.mark.unit
def test_allocate_threads_4_cpus() -> None:
    """Test thread allocation with 4 CPUs."""
    minimap2, samtools = allocate_threads(4)
    # 80% of 4 = 3.2 -> 3, remainder: 4-3=1
    assert minimap2 == 3
    assert samtools == 1
    assert minimap2 + samtools <= 4


@pytest.mark.unit
def test_allocate_threads_8_cpus() -> None:
    """Test thread allocation with 8 CPUs."""
    minimap2, samtools = allocate_threads(8)
    # 80% of 8 = 6.4 -> 6
    # Remainder: 8-6=2
    assert minimap2 == 6
    assert samtools == 2


@pytest.mark.unit
def test_allocate_threads_16_cpus() -> None:
    """Test thread allocation with 16 CPUs (at max bound)."""
    minimap2, samtools = allocate_threads(16)
    # 80% of 16 = 12.8 -> 12
    # Remainder: 16-12=4
    assert minimap2 == 12
    assert samtools == 4


@pytest.mark.unit
def test_allocate_threads_32_cpus() -> None:
    """Test thread allocation with 32 CPUs (capped to 16)."""
    minimap2, samtools = allocate_threads(32)
    # Clamped to 16
    # 80% of 16 = 12.8 -> 12
    # Remainder: 16-12=4
    assert minimap2 == 12
    assert samtools == 4


@pytest.mark.unit
def test_allocate_threads_bounds() -> None:
    """Test thread allocation respects bounds across various inputs."""
    for total in range(1, 33):
        minimap2, samtools = allocate_threads(total)
        # minimap2 at least 2, samtools at least 1
        assert minimap2 >= 2, f"minimap2={minimap2} < 2 for total={total}"
        assert samtools >= 1, f"samtools={samtools} < 1 for total={total}"
        # samtools capped at 4
        assert samtools <= 4, f"samtools={samtools} > 4 for total={total}"
        # minimap2 should be in reasonable range (capped at 16 total)
        assert minimap2 <= 16, f"minimap2={minimap2} > 16 for total={total}"
        # Sum should not exceed clamped total (except total=2 where mm2 min forces 3)
        clamped = max(2, min(16, total))
        assert minimap2 + samtools <= clamped + 1, (
            f"sum={minimap2 + samtools} > {clamped}+1 for total={total}"
        )
