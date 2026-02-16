"""CLI smoke tests for plasmicheck."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from plasmicheck.version import __version__

SUBCOMMANDS = [
    "convert",
    "index",
    "align",
    "compare",
    "spliced",
    "pipeline",
    "report",
    "summary_reports",
    "cache",
]

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)


class TestCLI:
    @pytest.mark.unit
    def test_help_exits_zero(self) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", "--help"],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        assert "plasmicheck" in result.stdout.lower()

    @pytest.mark.unit
    @pytest.mark.parametrize("subcommand", SUBCOMMANDS)
    def test_subcommand_help(self, subcommand: str) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", subcommand, "--help"],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0

    @pytest.mark.unit
    def test_version_output(self) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", "--version"],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        assert __version__ in result.stdout

    @pytest.mark.unit
    @pytest.mark.parametrize("subcommand", ["pipeline", "report", "summary_reports"])
    def test_report_flags_in_help(self, subcommand: str) -> None:
        """Verify --static-report and --plotly-mode appear in subcommand help."""
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", subcommand, "--help"],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        assert "--static-report" in result.stdout
        assert "--plotly-mode" in result.stdout
        assert "directory" in result.stdout  # Default value should appear

    @pytest.mark.unit
    def test_generate_report_lazy_imports(self) -> None:
        """Verify that importing generate_report doesn't trigger heavy imports."""
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "import sys; "
                "before = set(sys.modules.keys()); "
                "import plasmicheck.scripts.generate_report; "
                "after = set(sys.modules.keys()); "
                "new_modules = after - before; "
                "heavy = {'pandas', 'plotly', 'jinja2', 'kaleido'}; "
                "loaded_heavy = heavy & {m.split('.')[0] for m in new_modules}; "
                "print(','.join(sorted(loaded_heavy)) if loaded_heavy else 'NONE'); ",
            ],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        # pandas, plotly, jinja2, kaleido should NOT be loaded on import
        assert result.stdout.strip() == "NONE", (
            f"Heavy modules loaded on import: {result.stdout.strip()}"
        )

    @pytest.mark.unit
    def test_generate_summary_reports_lazy_imports(self) -> None:
        """Verify that importing generate_summary_reports doesn't trigger heavy imports."""
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "import sys; "
                "before = set(sys.modules.keys()); "
                "import plasmicheck.scripts.generate_summary_reports; "
                "after = set(sys.modules.keys()); "
                "new_modules = after - before; "
                "heavy = {'pandas', 'plotly', 'jinja2', 'kaleido', 'numpy', 'scipy', 'statsmodels'}; "
                "loaded_heavy = heavy & {m.split('.')[0] for m in new_modules}; "
                "print(','.join(sorted(loaded_heavy)) if loaded_heavy else 'NONE'); ",
            ],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        # pandas, plotly, jinja2, kaleido, numpy, scipy, statsmodels should NOT be loaded on import
        assert result.stdout.strip() == "NONE", (
            f"Heavy modules loaded on import: {result.stdout.strip()}"
        )

    @pytest.mark.unit
    def test_cache_info_empty(self, tmp_path: Path) -> None:
        """Test cache info with empty cache."""
        cache_dir = tmp_path / "cache"
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", "cache", "info", "--cache-dir", str(cache_dir)],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        assert "Total entries:   0" in result.stdout
        assert "Total size:      0 B" in result.stdout
        assert f"Cache directory: {cache_dir}" in result.stdout

    @pytest.mark.unit
    def test_cache_info_with_entries(self, tmp_path: Path) -> None:
        """Test cache info with cached entries."""
        from plasmicheck.cache_manager import CacheManager

        cache_dir = tmp_path / "cache"
        cache = CacheManager(cache_dir)

        # Create dummy index files
        dummy_fasta = tmp_path / "test.fa"
        dummy_fasta.write_text(">seq1\nACGT\n")
        dummy_mmi = tmp_path / "test.mmi"
        dummy_mmi.write_bytes(b"dummy_mmi_content")
        dummy_fai = tmp_path / "test.fai"
        dummy_fai.write_bytes(b"seq1\t4\t6\t4\t5\n")

        # Put entry in cache
        cache.put_index(str(dummy_fasta), str(dummy_mmi), str(dummy_fai))

        # Run cache info
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", "cache", "info", "--cache-dir", str(cache_dir)],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        assert "Total entries:   1" in result.stdout
        assert "Total size:" in result.stdout
        # Size should be non-zero (check it's not exactly "0 B")
        assert "Total size:      0 B" not in result.stdout

    @pytest.mark.unit
    def test_cache_clean_all(self, tmp_path: Path) -> None:
        """Test cache clean --all removes all entries."""
        from plasmicheck.cache_manager import CacheManager

        cache_dir = tmp_path / "cache"
        cache = CacheManager(cache_dir)

        # Create dummy index files
        dummy_fasta = tmp_path / "test.fa"
        dummy_fasta.write_text(">seq1\nACGT\n")
        dummy_mmi = tmp_path / "test.mmi"
        dummy_mmi.write_bytes(b"dummy_mmi_content")
        dummy_fai = tmp_path / "test.fai"
        dummy_fai.write_bytes(b"seq1\t4\t6\t4\t5\n")

        # Put entry in cache
        cache.put_index(str(dummy_fasta), str(dummy_mmi), str(dummy_fai))

        # Run cache clean
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", "cache", "clean", "--all", "--cache-dir", str(cache_dir)],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        assert "Removed 1 cache entries" in result.stdout

        # Verify cache is now empty
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", "cache", "info", "--cache-dir", str(cache_dir)],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        assert "Total entries:   0" in result.stdout

    @pytest.mark.unit
    def test_cache_help(self) -> None:
        """Test cache subcommand help."""
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", "cache", "--help"],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        assert "info" in result.stdout
        assert "clean" in result.stdout

    @pytest.mark.unit
    def test_cache_no_minimap2_required(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test that cache commands work without minimap2 on PATH."""
        import os

        cache_dir = tmp_path / "cache"

        # Remove all directories from PATH to ensure minimap2/samtools can't be found
        monkeypatch.setenv("PATH", "")

        # cache info should still work
        result = subprocess.run(
            [sys.executable, "-m", "plasmicheck.cli", "cache", "info", "--cache-dir", str(cache_dir)],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
            env={**os.environ, "PATH": ""},
        )
        assert result.returncode == 0
        assert "Total entries:   0" in result.stdout
        # Should NOT contain error about missing tools
        assert "minimap2" not in result.stderr.lower()
        assert "samtools" not in result.stderr.lower()
