"""CLI parsing, dispatch, and entry-point smoke tests."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from plasmicheck import cli
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


@pytest.mark.unit
def test_version_without_requirements(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(sys, "argv", ["plasmicheck", "--version"])
    with (
        patch.object(cli, "check_requirements") as requirements,
        pytest.raises(SystemExit) as exc,
    ):
        cli.main()
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out
    requirements.assert_not_called()


@pytest.mark.unit
@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["pipeline"], "--human_fasta"),
        (["compare", "-p", "p.bam", "-m", "h.bam", "-o", "out", "-t", "bad"], "invalid float"),
    ],
)
def test_invalid_arguments_exit_before_requirements(
    arguments: list[str],
    message: str,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(sys, "argv", ["plasmicheck", *arguments])
    with (
        patch.object(cli, "check_requirements") as requirements,
        pytest.raises(SystemExit) as exc,
    ):
        cli.main()
    assert exc.value.code == 2
    assert message in capsys.readouterr().err
    requirements.assert_not_called()


@pytest.mark.unit
@pytest.mark.parametrize(("threshold_args", "expected"), [([], 0.42), (["-t", "0.25"], 0.25)])
def test_compare_forwards_threshold(
    threshold_args: list[str], expected: float, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        ["plasmicheck", "compare", "-p", "p.bam", "-m", "h.bam", "-o", "out", *threshold_args],
    )
    with (
        patch.object(cli, "setup_logging"),
        patch.object(cli, "check_requirements"),
        patch("plasmicheck.scripts.compare_alignments.compare_alignments") as compare,
    ):
        cli.main(default_threshold=0.42)
    compare.assert_called_once_with("p.bam", "h.bam", "out", expected)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("interactive", "progress_args", "expected_progress"),
    [(True, [], True), (True, ["--no_progress"], False), (False, [], False)],
)
def test_pipeline_forwards_run_options(
    interactive: bool,
    progress_args: list[str],
    expected_progress: bool,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plasmicheck",
            "pipeline",
            "-hf",
            "human.fa",
            "-pf",
            "plasmid.gb",
            "-sf1",
            "r1.fastq",
            "-sf2",
            "r2.fastq",
            "-o",
            "output",
            "--no-keep_intermediate",
            "--threshold",
            "0.25",
            "--dry_run",
            "--threads",
            "3",
            "--cache-dir",
            "custom-cache",
            "--no-cache",
            "--static-report",
            "--plotly-mode",
            "embedded",
            "--plot-backend",
            "matplotlib",
            *progress_args,
        ],
    )
    monkeypatch.setattr(sys.stderr, "isatty", lambda: interactive)
    with (
        patch.object(cli, "setup_logging"),
        patch.object(cli, "check_requirements"),
        patch("plasmicheck.scripts.run_pipeline.run_pipeline") as pipeline,
    ):
        cli.main()
    pipeline.assert_called_once_with(
        "human.fa",
        "plasmid.gb",
        "output",
        False,
        500,
        False,
        False,
        1000,
        0.25,
        "intermediate",
        None,
        False,
        sequencing_files_r1="r1.fastq",
        sequencing_files_r2="r2.fastq",
        dry_run=True,
        progress=expected_progress,
        static_report=True,
        plotly_mode="embedded",
        plot_backend="matplotlib",
        threads=3,
        cache_dir="custom-cache",
        no_cache=True,
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    ("tools_missing", "package_missing"), [(True, False), (False, True), (True, True)]
)
def test_missing_requirements_stop_dispatch(
    tools_missing: bool,
    package_missing: bool,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    monkeypatch.setattr(
        sys, "argv", ["plasmicheck", "compare", "-p", "p.bam", "-m", "h.bam", "-o", "out"]
    )
    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(cli, "REQUIRED_TOOLS", ["minimap2"] if tools_missing else [])
    monkeypatch.setattr(
        cli,
        "REQUIRED_PYTHON_PACKAGES",
        ["sys", "_plasmicheck_missing_package"] if package_missing else ["sys"],
    )
    with (
        patch.object(cli, "setup_logging"),
        patch("plasmicheck.scripts.compare_alignments.compare_alignments") as compare,
        pytest.raises(SystemExit) as exc,
    ):
        cli.main()
    assert exc.value.code == 1
    assert ("Missing tools: minimap2" in caplog.text) == tools_missing
    assert (
        "Missing Python packages: _plasmicheck_missing_package" in caplog.text
    ) == package_missing
    compare.assert_not_called()


@pytest.mark.unit
@pytest.mark.parametrize("cdna_output", [None, "positions.txt"])
def test_spliced_optional_cdna_output(
    cdna_output: str | None, monkeypatch: pytest.MonkeyPatch
) -> None:
    args = [
        "plasmicheck",
        "spliced",
        "-i",
        "human.mmi",
        "-p",
        "plasmid.fa",
        "-o",
        "regions.fa",
        "-b",
        "spliced.bam",
        "-hf",
        "human.fa",
    ]
    if cdna_output:
        args.extend(["--cDNA_output", cdna_output])
    monkeypatch.setattr(sys, "argv", args)
    with (
        patch.object(cli, "setup_logging"),
        patch.object(cli, "check_requirements"),
        patch("plasmicheck.scripts.spliced_alignment.spliced_alignment", return_value=[]) as align,
        patch("plasmicheck.scripts.spliced_alignment.extract_human_reference") as extract,
        patch("plasmicheck.scripts.spliced_alignment.extract_plasmid_cdna_positions") as cdna,
    ):
        cli.main()

    align.assert_called_once_with("human.mmi", "plasmid.fa", "spliced.bam", 1000)
    extract.assert_called_once_with("human.fa", [], "regions.fa")
    if cdna_output:
        cdna.assert_called_once_with("plasmid.fa", "spliced.bam", "positions.txt")
    else:
        cdna.assert_not_called()


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
    def test_subcommand_help(
        self, subcommand: str, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
    ) -> None:
        monkeypatch.setattr(sys, "argv", ["plasmicheck", subcommand, "--help"])
        with (
            patch.object(cli, "check_requirements") as requirements,
            pytest.raises(SystemExit) as exc,
        ):
            cli.main()
        assert exc.value.code == 0
        assert subcommand in capsys.readouterr().out
        requirements.assert_not_called()

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
    def test_report_flags_in_help(
        self, subcommand: str, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Verify --static-report and --plotly-mode appear in subcommand help."""
        monkeypatch.setattr(sys, "argv", ["plasmicheck", subcommand, "--help"])
        with pytest.raises(SystemExit) as exc:
            cli.main()
        assert exc.value.code == 0
        output = capsys.readouterr().out
        assert "--static-report" in output
        assert "--plotly-mode" in output
        assert "directory" in output

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
            [
                sys.executable,
                "-m",
                "plasmicheck.cli",
                "cache",
                "info",
                "--cache-dir",
                str(cache_dir),
            ],
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
            [
                sys.executable,
                "-m",
                "plasmicheck.cli",
                "cache",
                "info",
                "--cache-dir",
                str(cache_dir),
            ],
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
            [
                sys.executable,
                "-m",
                "plasmicheck.cli",
                "cache",
                "clean",
                "--all",
                "--cache-dir",
                str(cache_dir),
            ],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0
        assert "Removed 1 cache entries" in result.stdout

        # Verify cache is now empty
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "plasmicheck.cli",
                "cache",
                "info",
                "--cache-dir",
                str(cache_dir),
            ],
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
    def test_cache_no_minimap2_required(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test that cache commands work without minimap2 on PATH."""
        import os

        cache_dir = tmp_path / "cache"

        # Remove all directories from PATH to ensure minimap2/samtools can't be found
        monkeypatch.setenv("PATH", "")

        # cache info should still work
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "plasmicheck.cli",
                "cache",
                "info",
                "--cache-dir",
                str(cache_dir),
            ],
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
