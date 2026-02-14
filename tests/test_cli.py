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
