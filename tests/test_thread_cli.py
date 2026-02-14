"""Tests for CLI --threads flag and pipeline thread integration."""

from __future__ import annotations

import logging
import subprocess
from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch

import pytest

from plasmicheck.thread_config import allocate_threads

if TYPE_CHECKING:
    from _pytest.logging import LogCaptureFixture


@pytest.mark.unit
def test_pipeline_help_shows_threads_flag() -> None:
    """Test that plasmicheck pipeline --help includes --threads flag."""
    result = subprocess.run(
        ["plasmicheck", "pipeline", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "--threads" in result.stdout
    assert "auto-detect" in result.stdout.lower()


@pytest.mark.unit
def test_pipeline_threads_override(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test run_pipeline with --threads override (threads=8)."""
    # Mock all pipeline dependencies
    mock_align_reads = MagicMock()
    mock_quality_control = MagicMock()
    mock_convert = MagicMock()
    mock_create_indexes = MagicMock()
    mock_spliced_alignment = MagicMock(return_value=[])
    mock_extract_human_reference = MagicMock()
    mock_extract_cdna = MagicMock()
    mock_compare_alignments = MagicMock()
    mock_generate_report = MagicMock()
    mock_write_md5sum = MagicMock()

    # Patch all the imports
    with (
        patch("plasmicheck.scripts.run_pipeline.align_reads", mock_align_reads),
        patch("plasmicheck.scripts.run_pipeline.quality_control", mock_quality_control),
        patch(
            "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
            mock_convert,
        ),
        patch("plasmicheck.scripts.run_pipeline.create_indexes", mock_create_indexes),
        patch(
            "plasmicheck.scripts.run_pipeline.spliced_alignment",
            mock_spliced_alignment,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.extract_human_reference",
            mock_extract_human_reference,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions",
            mock_extract_cdna,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.compare_alignments",
            mock_compare_alignments,
        ),
        patch("plasmicheck.scripts.run_pipeline.generate_report", mock_generate_report),
        patch("plasmicheck.scripts.run_pipeline.write_md5sum", mock_write_md5sum),
        patch("plasmicheck.scripts.run_pipeline.os.makedirs"),
        patch("plasmicheck.scripts.run_pipeline.os.path.exists", return_value=False),
    ):
        from plasmicheck.scripts.run_pipeline import run_pipeline

        # Call run_pipeline with threads=8
        run_pipeline(
            human_fasta="test.fasta",
            plasmid_files="test_plasmid.gb",
            output_folder="test_output",
            sequencing_files_r1="test_reads.fastq",
            threads=8,
        )

        # Verify align_reads was called with correct thread allocation
        mm2_threads, sam_threads = allocate_threads(8)
        assert mock_align_reads.call_count == 2  # plasmid + human alignments

        # Check that both align_reads calls got the right thread parameters
        for call_args in mock_align_reads.call_args_list:
            kwargs = call_args.kwargs
            assert kwargs["minimap2_threads"] == mm2_threads
            assert kwargs["samtools_threads"] == sam_threads
            assert kwargs["samtools_sort_memory"] == "2G"  # from config


@pytest.mark.unit
def test_pipeline_threads_autodetect(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test run_pipeline without --threads (auto-detection)."""
    # Mock detect_cpu_count to return predictable value
    mock_detect = MagicMock(return_value=(4, "test-mock"))

    # Mock all pipeline dependencies
    mock_align_reads = MagicMock()
    mock_quality_control = MagicMock()
    mock_convert = MagicMock()
    mock_create_indexes = MagicMock()
    mock_spliced_alignment = MagicMock(return_value=[])
    mock_extract_human_reference = MagicMock()
    mock_extract_cdna = MagicMock()
    mock_compare_alignments = MagicMock()
    mock_generate_report = MagicMock()
    mock_write_md5sum = MagicMock()

    with (
        patch("plasmicheck.scripts.run_pipeline.detect_cpu_count", mock_detect),
        patch("plasmicheck.scripts.run_pipeline.align_reads", mock_align_reads),
        patch("plasmicheck.scripts.run_pipeline.quality_control", mock_quality_control),
        patch(
            "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
            mock_convert,
        ),
        patch("plasmicheck.scripts.run_pipeline.create_indexes", mock_create_indexes),
        patch(
            "plasmicheck.scripts.run_pipeline.spliced_alignment",
            mock_spliced_alignment,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.extract_human_reference",
            mock_extract_human_reference,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions",
            mock_extract_cdna,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.compare_alignments",
            mock_compare_alignments,
        ),
        patch("plasmicheck.scripts.run_pipeline.generate_report", mock_generate_report),
        patch("plasmicheck.scripts.run_pipeline.write_md5sum", mock_write_md5sum),
        patch("plasmicheck.scripts.run_pipeline.os.makedirs"),
        patch("plasmicheck.scripts.run_pipeline.os.path.exists", return_value=False),
    ):
        from plasmicheck.scripts.run_pipeline import run_pipeline

        # Call run_pipeline with threads=None (auto-detect)
        run_pipeline(
            human_fasta="test.fasta",
            plasmid_files="test_plasmid.gb",
            output_folder="test_output",
            sequencing_files_r1="test_reads.fastq",
            threads=None,
        )

        # Verify detect_cpu_count was called
        mock_detect.assert_called_once()

        # Verify align_reads received threads matching allocate_threads(4)
        mm2_threads, sam_threads = allocate_threads(4)
        assert mock_align_reads.call_count == 2

        for call_args in mock_align_reads.call_args_list:
            kwargs = call_args.kwargs
            assert kwargs["minimap2_threads"] == mm2_threads
            assert kwargs["samtools_threads"] == sam_threads
            assert kwargs["samtools_sort_memory"] == "2G"


@pytest.mark.unit
def test_pipeline_thread_logging(caplog: LogCaptureFixture) -> None:
    """Test that thread allocation is logged correctly."""
    # Mock all pipeline dependencies
    mock_align_reads = MagicMock()
    mock_quality_control = MagicMock()
    mock_convert = MagicMock()
    mock_create_indexes = MagicMock()
    mock_spliced_alignment = MagicMock(return_value=[])
    mock_extract_human_reference = MagicMock()
    mock_extract_cdna = MagicMock()
    mock_compare_alignments = MagicMock()
    mock_generate_report = MagicMock()
    mock_write_md5sum = MagicMock()

    with (
        patch("plasmicheck.scripts.run_pipeline.align_reads", mock_align_reads),
        patch("plasmicheck.scripts.run_pipeline.quality_control", mock_quality_control),
        patch(
            "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
            mock_convert,
        ),
        patch("plasmicheck.scripts.run_pipeline.create_indexes", mock_create_indexes),
        patch(
            "plasmicheck.scripts.run_pipeline.spliced_alignment",
            mock_spliced_alignment,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.extract_human_reference",
            mock_extract_human_reference,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions",
            mock_extract_cdna,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.compare_alignments",
            mock_compare_alignments,
        ),
        patch("plasmicheck.scripts.run_pipeline.generate_report", mock_generate_report),
        patch("plasmicheck.scripts.run_pipeline.write_md5sum", mock_write_md5sum),
        patch("plasmicheck.scripts.run_pipeline.os.makedirs"),
        patch("plasmicheck.scripts.run_pipeline.os.path.exists", return_value=False),
        caplog.at_level(logging.INFO),
    ):
        from plasmicheck.scripts.run_pipeline import run_pipeline

        # Call with threads=8 to get predictable logging
        run_pipeline(
            human_fasta="test.fasta",
            plasmid_files="test_plasmid.gb",
            output_folder="test_output",
            sequencing_files_r1="test_reads.fastq",
            threads=8,
        )

        # Check that logging includes thread info
        log_text = caplog.text
        assert "Using 8 threads (CLI --threads=8)" in log_text
        assert "minimap2:" in log_text
        assert "samtools:" in log_text
        assert "samtools sort memory:" in log_text


@pytest.mark.unit
def test_pipeline_thread_logging_autodetect(caplog: LogCaptureFixture) -> None:
    """Test that auto-detected thread count is logged with source."""
    # Mock detect_cpu_count
    mock_detect = MagicMock(return_value=(6, "SLURM_CPUS_PER_TASK"))

    # Mock all pipeline dependencies
    mock_align_reads = MagicMock()
    mock_quality_control = MagicMock()
    mock_convert = MagicMock()
    mock_create_indexes = MagicMock()
    mock_spliced_alignment = MagicMock(return_value=[])
    mock_extract_human_reference = MagicMock()
    mock_extract_cdna = MagicMock()
    mock_compare_alignments = MagicMock()
    mock_generate_report = MagicMock()
    mock_write_md5sum = MagicMock()

    with (
        patch("plasmicheck.scripts.run_pipeline.detect_cpu_count", mock_detect),
        patch("plasmicheck.scripts.run_pipeline.align_reads", mock_align_reads),
        patch("plasmicheck.scripts.run_pipeline.quality_control", mock_quality_control),
        patch(
            "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
            mock_convert,
        ),
        patch("plasmicheck.scripts.run_pipeline.create_indexes", mock_create_indexes),
        patch(
            "plasmicheck.scripts.run_pipeline.spliced_alignment",
            mock_spliced_alignment,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.extract_human_reference",
            mock_extract_human_reference,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions",
            mock_extract_cdna,
        ),
        patch(
            "plasmicheck.scripts.run_pipeline.compare_alignments",
            mock_compare_alignments,
        ),
        patch("plasmicheck.scripts.run_pipeline.generate_report", mock_generate_report),
        patch("plasmicheck.scripts.run_pipeline.write_md5sum", mock_write_md5sum),
        patch("plasmicheck.scripts.run_pipeline.os.makedirs"),
        patch("plasmicheck.scripts.run_pipeline.os.path.exists", return_value=False),
        caplog.at_level(logging.INFO),
    ):
        from plasmicheck.scripts.run_pipeline import run_pipeline

        # Call with threads=None for auto-detection
        run_pipeline(
            human_fasta="test.fasta",
            plasmid_files="test_plasmid.gb",
            output_folder="test_output",
            sequencing_files_r1="test_reads.fastq",
            threads=None,
        )

        # Check that logging includes auto-detected source
        log_text = caplog.text
        assert "Using 6 threads (SLURM_CPUS_PER_TASK)" in log_text
        assert "minimap2:" in log_text
        assert "samtools:" in log_text
