"""CLI smoke tests for plasmicheck."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

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
        assert "0.31.0" in result.stdout
