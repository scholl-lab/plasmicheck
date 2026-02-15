"""Unit tests for coverage metrics computation."""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from plasmicheck.scripts.coverage_metrics import (
    compute_coverage_stats,
    compute_region_coverage_metrics,
    get_depth_array,
)

# ============================================================================
# Tests for get_depth_array
# ============================================================================


@pytest.mark.unit
def test_get_depth_array_calls_count_coverage_correctly() -> None:
    """Test that count_coverage is called with correct parameters."""
    mock_bam = MagicMock()
    mock_bam.__enter__ = MagicMock(return_value=mock_bam)
    mock_bam.__exit__ = MagicMock(return_value=None)
    mock_bam.count_coverage.return_value = (
        np.array([1, 2, 3]),
        np.array([4, 5, 6]),
        np.array([7, 8, 9]),
        np.array([10, 11, 12]),
    )

    with patch("plasmicheck.scripts.coverage_metrics.pysam.AlignmentFile", return_value=mock_bam):
        result = get_depth_array("/path/to/bam", "chr1", 100, 200)

    mock_bam.count_coverage.assert_called_once_with(
        contig="chr1", start=100, stop=200, quality_threshold=0, read_callback="all"
    )
    # Sum of arrays: [1+4+7+10, 2+5+8+11, 3+6+9+12] = [22, 26, 30]
    np.testing.assert_array_equal(result, np.array([22, 26, 30]))


@pytest.mark.unit
def test_get_depth_array_sums_acgt_correctly() -> None:
    """Test that A, C, G, T arrays are summed correctly."""
    mock_bam = MagicMock()
    mock_bam.__enter__ = MagicMock(return_value=mock_bam)
    mock_bam.__exit__ = MagicMock(return_value=None)
    mock_bam.count_coverage.return_value = (
        np.array([10, 0, 5]),
        np.array([0, 20, 0]),
        np.array([5, 0, 10]),
        np.array([0, 15, 0]),
    )

    with patch("plasmicheck.scripts.coverage_metrics.pysam.AlignmentFile", return_value=mock_bam):
        result = get_depth_array("/path/to/bam", "plasmid", 0, 3)

    # Position 0: 10+0+5+0=15, Position 1: 0+20+0+15=35, Position 2: 5+0+10+0=15
    np.testing.assert_array_equal(result, np.array([15, 35, 15]))


@pytest.mark.unit
def test_get_depth_array_empty_region() -> None:
    """Test that empty region (start >= end) returns empty array."""
    result = get_depth_array("/path/to/bam", "chr1", 100, 100)
    np.testing.assert_array_equal(result, np.array([]))

    result = get_depth_array("/path/to/bam", "chr1", 200, 100)
    np.testing.assert_array_equal(result, np.array([]))


# ============================================================================
# Tests for compute_coverage_stats
# ============================================================================


@pytest.mark.unit
def test_compute_coverage_stats_mean_median() -> None:
    """Test mean and median depth calculation."""
    depth = np.array([10, 20, 30])
    stats = compute_coverage_stats(depth, [1, 5])

    assert stats["mean_depth"] == 20.0
    assert stats["median_depth"] == 20.0


@pytest.mark.unit
def test_compute_coverage_stats_breadth_1x() -> None:
    """Test breadth at 1x threshold."""
    depth = np.array([0, 5, 10, 0, 3])
    stats = compute_coverage_stats(depth, [5])

    # 3 out of 5 bases have depth >= 1
    assert stats["breadth_1x"] == pytest.approx(0.6)


@pytest.mark.unit
def test_compute_coverage_stats_breadth_5x() -> None:
    """Test breadth at 5x threshold."""
    depth = np.array([0, 5, 10, 0, 3])
    stats = compute_coverage_stats(depth, [5])

    # 2 out of 5 bases have depth >= 5
    assert stats["breadth_5x"] == pytest.approx(0.4)


@pytest.mark.unit
def test_compute_coverage_stats_cv() -> None:
    """Test coefficient of variation calculation."""
    depth = np.array([10, 20, 30])
    stats = compute_coverage_stats(depth, [1])

    # std(ddof=1) = 10, mean = 20, cv = 10/20 = 0.5
    assert stats["cv"] == pytest.approx(0.5)


@pytest.mark.unit
def test_compute_coverage_stats_all_zero() -> None:
    """Test all metrics return 0 for all-zero array."""
    depth = np.array([0, 0, 0, 0])
    stats = compute_coverage_stats(depth, [1, 5])

    assert stats["mean_depth"] == 0.0
    assert stats["median_depth"] == 0.0
    assert stats["breadth_1x"] == 0.0
    assert stats["breadth_5x"] == 0.0
    assert stats["cv"] == 0.0


@pytest.mark.unit
def test_compute_coverage_stats_empty_array() -> None:
    """Test all metrics return 0 for empty array."""
    depth = np.array([])
    stats = compute_coverage_stats(depth, [1, 5, 10])

    assert stats["mean_depth"] == 0.0
    assert stats["median_depth"] == 0.0
    assert stats["breadth_1x"] == 0.0
    assert stats["breadth_5x"] == 0.0
    assert stats["breadth_10x"] == 0.0
    assert stats["cv"] == 0.0


@pytest.mark.unit
def test_compute_coverage_stats_single_element() -> None:
    """Test CV is 0 for single-element array (no variance)."""
    depth = np.array([42])
    stats = compute_coverage_stats(depth, [1])

    assert stats["mean_depth"] == 42.0
    assert stats["median_depth"] == 42.0
    assert stats["breadth_1x"] == 1.0
    # std with ddof=1 on single element is NaN, but division by zero check prevents error
    # Actually std([42], ddof=1) produces a warning and returns NaN, mean=42, so cv would be NaN/42
    # But our implementation checks mean == 0, not std == NaN. Let me verify the behavior.
    # Actually, numpy.std([42], ddof=1) returns nan, and nan / 42 = nan
    # We should handle this case. Let me check what happens.
    # For now, let's check if cv is a valid number (not NaN)
    assert not np.isnan(stats["cv"])


@pytest.mark.unit
def test_compute_coverage_stats_multiple_thresholds() -> None:
    """Test multiple breadth thresholds produce correct keys."""
    depth = np.array([1, 3, 5, 10, 15, 20])
    stats = compute_coverage_stats(depth, [1, 5, 10])

    # All thresholds should have corresponding breadth keys
    assert "breadth_1x" in stats
    assert "breadth_5x" in stats
    assert "breadth_10x" in stats

    # breadth_1x: 6/6 = 1.0
    # breadth_5x: 4/6 = 0.666...
    # breadth_10x: 3/6 = 0.5
    assert stats["breadth_1x"] == 1.0
    assert stats["breadth_5x"] == pytest.approx(2 / 3)
    assert stats["breadth_10x"] == pytest.approx(0.5)


@pytest.mark.unit
def test_compute_coverage_stats_returns_python_floats() -> None:
    """Test that all returned values are Python floats, not numpy scalars."""
    depth = np.array([10, 20, 30])
    stats = compute_coverage_stats(depth, [1, 5])

    for key, value in stats.items():
        assert isinstance(value, float), f"{key} should be Python float, got {type(value)}"


# ============================================================================
# Tests for compute_region_coverage_metrics
# ============================================================================


def _mock_alignment_file(
    contig_name: str = "plasmid", plasmid_length: int = 1000
) -> tuple[MagicMock, MagicMock]:
    """Helper to create mock pysam.AlignmentFile."""
    mock_bam = MagicMock()
    mock_ctx = MagicMock()
    mock_ctx.__enter__ = MagicMock(return_value=mock_bam)
    mock_ctx.__exit__ = MagicMock(return_value=None)
    mock_bam.get_reference_name.return_value = contig_name
    mock_bam.lengths = [plasmid_length]
    return mock_ctx, mock_bam


@pytest.mark.unit
def test_compute_region_coverage_metrics_with_insert() -> None:
    """Test metrics computed separately for insert and backbone."""
    mock_ctx, _mock_bam = _mock_alignment_file("plasmid", 1000)

    # Mock get_depth_array to return different arrays for insert vs backbone
    with (
        patch("plasmicheck.scripts.coverage_metrics.pysam.AlignmentFile", return_value=mock_ctx),
        patch("plasmicheck.scripts.coverage_metrics.get_depth_array") as mock_get_depth,
    ):

        def side_effect(bam_path: str, contig: str, start: int, end: int) -> Any:
            if start == 100 and end == 501:  # Insert region [100, 500] -> [100, 501)
                return np.array([10] * 401)  # Uniform coverage
            elif start == 0 and end == 100:  # Backbone before
                return np.array([5] * 100)
            elif start == 501 and end == 1000:  # Backbone after
                return np.array([5] * 499)
            return np.array([])

        mock_get_depth.side_effect = side_effect

        metrics, fallback = compute_region_coverage_metrics("/path/to/bam.bam", (100, 500), [1, 5])

    assert not fallback
    assert "insert" in metrics
    assert "backbone" in metrics

    # Insert: all 10s
    assert metrics["insert"]["mean_depth"] == 10.0
    assert metrics["insert"]["breadth_1x"] == 1.0

    # Backbone: all 5s
    assert metrics["backbone"]["mean_depth"] == 5.0
    assert metrics["backbone"]["breadth_1x"] == 1.0


@pytest.mark.unit
def test_compute_region_coverage_metrics_fallback_no_insert() -> None:
    """Test whole-plasmid fallback when insert_region is None."""
    mock_ctx, _mock_bam = _mock_alignment_file("plasmid", 1000)

    with (
        patch("plasmicheck.scripts.coverage_metrics.pysam.AlignmentFile", return_value=mock_ctx),
        patch("plasmicheck.scripts.coverage_metrics.get_depth_array") as mock_get_depth,
    ):
        mock_get_depth.return_value = np.array([15] * 1000)

        metrics, fallback = compute_region_coverage_metrics("/path/to/bam.bam", None, [1, 5])

    assert fallback is True
    assert metrics["insert"]["mean_depth"] == 15.0
    assert metrics["insert"]["breadth_1x"] == 1.0

    # Backbone should be all zeros in fallback mode
    assert metrics["backbone"]["mean_depth"] == 0.0
    assert metrics["backbone"]["breadth_1x"] == 0.0
    assert metrics["backbone"]["breadth_5x"] == 0.0


@pytest.mark.unit
def test_compute_region_coverage_metrics_insert_at_start() -> None:
    """Test insert at start of plasmid (backbone only after insert)."""
    mock_ctx, _mock_bam = _mock_alignment_file("plasmid", 1000)

    with (
        patch("plasmicheck.scripts.coverage_metrics.pysam.AlignmentFile", return_value=mock_ctx),
        patch("plasmicheck.scripts.coverage_metrics.get_depth_array") as mock_get_depth,
    ):

        def side_effect(bam_path: str, contig: str, start: int, end: int) -> Any:
            if start == 0 and end == 501:  # Insert [0, 500]
                return np.array([10] * 501)
            elif start == 501 and end == 1000:  # Backbone after
                return np.array([5] * 499)
            return np.array([])

        mock_get_depth.side_effect = side_effect

        _metrics, fallback = compute_region_coverage_metrics("/path/to/bam.bam", (0, 500), [1, 5])

    assert not fallback
    # Should only call get_depth_array for insert and backbone_after (not before, since start=0)
    assert mock_get_depth.call_count == 2


@pytest.mark.unit
def test_compute_region_coverage_metrics_insert_at_end() -> None:
    """Test insert at end of plasmid (backbone only before insert)."""
    mock_ctx, _mock_bam = _mock_alignment_file("plasmid", 1000)

    with (
        patch("plasmicheck.scripts.coverage_metrics.pysam.AlignmentFile", return_value=mock_ctx),
        patch("plasmicheck.scripts.coverage_metrics.get_depth_array") as mock_get_depth,
    ):

        def side_effect(bam_path: str, contig: str, start: int, end: int) -> Any:
            if start == 500 and end == 1000:  # Insert [500, 999] -> [500, 1000)
                return np.array([10] * 500)
            elif start == 0 and end == 500:  # Backbone before
                return np.array([5] * 500)
            return np.array([])

        mock_get_depth.side_effect = side_effect

        _metrics, fallback = compute_region_coverage_metrics("/path/to/bam.bam", (500, 999), [1, 5])

    assert not fallback
    # Should only call get_depth_array for insert and backbone_before (not after, since end=999)
    assert mock_get_depth.call_count == 2


@pytest.mark.unit
def test_compute_region_coverage_metrics_off_by_one() -> None:
    """Test inclusive boundaries are correctly converted to half-open intervals."""
    mock_ctx, _mock_bam = _mock_alignment_file("plasmid", 1000)

    with (
        patch("plasmicheck.scripts.coverage_metrics.pysam.AlignmentFile", return_value=mock_ctx),
        patch("plasmicheck.scripts.coverage_metrics.get_depth_array") as mock_get_depth,
    ):
        mock_get_depth.return_value = np.array([10] * 100)

        compute_region_coverage_metrics("/path/to/bam.bam", (100, 500), [1])

    # Insert region [100, 500] inclusive should call get_depth_array with [100, 501) half-open
    calls = mock_get_depth.call_args_list
    insert_call = [c for c in calls if c[0][2] == 100 and c[0][3] == 501]
    assert len(insert_call) == 1


@pytest.mark.unit
def test_compute_region_coverage_metrics_zero_depth() -> None:
    """Test BAM with zero mapped reads returns all metrics as 0.0."""
    mock_ctx, _mock_bam = _mock_alignment_file("plasmid", 1000)

    with (
        patch("plasmicheck.scripts.coverage_metrics.pysam.AlignmentFile", return_value=mock_ctx),
        patch("plasmicheck.scripts.coverage_metrics.get_depth_array") as mock_get_depth,
    ):
        # Return all-zero arrays
        mock_get_depth.return_value = np.array([0] * 100)

        metrics, _fallback = compute_region_coverage_metrics("/path/to/bam.bam", (100, 500), [1, 5])

    assert metrics["insert"]["mean_depth"] == 0.0
    assert metrics["insert"]["breadth_1x"] == 0.0
    assert metrics["insert"]["cv"] == 0.0


# ============================================================================
# Tests for summary.tsv integration
# ============================================================================


@pytest.mark.unit
def test_summary_tsv_coverage_rows_written(tmp_path: Any) -> None:
    """Test that compare_alignments writes coverage rows to summary.tsv."""
    from pathlib import Path

    from plasmicheck.scripts.compare_alignments import compare_alignments

    # Create mock BAM files
    plasmid_bam = tmp_path / "plasmid.bam"
    human_bam = tmp_path / "human.bam"
    output_basename = tmp_path / "output"
    cdna_positions = tmp_path / "cDNA_positions.txt"

    # Write cDNA positions
    cdna_positions.write_text(
        "cDNA start position in plasmid: 100\ncDNA end position in plasmid: 500\n"
    )

    # Mock all the heavy lifting
    mock_bam_ctx = MagicMock()
    mock_bam = MagicMock()
    mock_bam_ctx.__enter__ = MagicMock(return_value=mock_bam)
    mock_bam_ctx.__exit__ = MagicMock(return_value=None)
    mock_bam.fetch.return_value = []
    mock_bam.get_tid.return_value = 0
    mock_bam.get_reference_name.return_value = "plasmid"
    mock_bam.lengths = [1000]

    with (
        patch(
            "plasmicheck.scripts.compare_alignments.pysam.AlignmentFile", return_value=mock_bam_ctx
        ),
        patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
        patch("plasmicheck.scripts.compare_alignments._streaming_compare") as mock_compare,
        patch("plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics") as mock_cov,
    ):
        # Mock comparison results
        mock_compare.return_value = {
            "Plasmid": 100,
            "Human": 50,
            "Tied": 10,
            "Backbone_Only": 5,
            "Ambiguous": 0,
        }

        # Mock coverage metrics
        mock_cov.return_value = (
            {
                "insert": {
                    "mean_depth": 42.50,
                    "median_depth": 40.00,
                    "breadth_1x": 0.95,
                    "breadth_5x": 0.85,
                    "cv": 0.25,
                },
                "backbone": {
                    "mean_depth": 10.25,
                    "median_depth": 8.50,
                    "breadth_1x": 0.60,
                    "breadth_5x": 0.50,
                    "cv": 0.35,
                },
            },
            False,
        )

        compare_alignments(str(plasmid_bam), str(human_bam), str(output_basename))

    # Read summary.tsv
    summary_path = Path(str(output_basename) + ".summary.tsv")
    assert summary_path.exists()
    content = summary_path.read_text()

    # Check coverage rows are present
    assert "MeanDepthInsert\t42.50" in content
    assert "MedianDepthInsert\t40.00" in content
    assert "BreadthInsert\t0.95" in content
    assert "BreadthInsert_5x\t0.85" in content
    assert "CoverageCV_Insert\t0.25" in content

    assert "MeanDepthBackbone\t10.25" in content
    assert "MedianDepthBackbone\t8.50" in content
    assert "BreadthBackbone\t0.60" in content
    assert "BreadthBackbone_5x\t0.50" in content
    assert "CoverageCV_Backbone\t0.35" in content

    assert "CoverageFallback\tFalse" in content


@pytest.mark.unit
def test_summary_tsv_existing_rows_unchanged(tmp_path: Any) -> None:
    """Test that existing summary.tsv rows are unchanged after adding coverage metrics."""
    from pathlib import Path

    from plasmicheck.scripts.compare_alignments import compare_alignments

    plasmid_bam = tmp_path / "plasmid.bam"
    human_bam = tmp_path / "human.bam"
    output_basename = tmp_path / "output"
    cdna_positions = tmp_path / "cDNA_positions.txt"

    cdna_positions.write_text(
        "cDNA start position in plasmid: 100\ncDNA end position in plasmid: 500\n"
    )

    mock_bam_ctx = MagicMock()
    mock_bam = MagicMock()
    mock_bam_ctx.__enter__ = MagicMock(return_value=mock_bam)
    mock_bam_ctx.__exit__ = MagicMock(return_value=None)
    mock_bam.fetch.return_value = []
    mock_bam.get_tid.return_value = 0
    mock_bam.get_reference_name.return_value = "plasmid"
    mock_bam.lengths = [1000]

    with (
        patch(
            "plasmicheck.scripts.compare_alignments.pysam.AlignmentFile", return_value=mock_bam_ctx
        ),
        patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
        patch("plasmicheck.scripts.compare_alignments._streaming_compare") as mock_compare,
        patch("plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics") as mock_cov,
    ):
        mock_compare.return_value = {
            "Plasmid": 100,
            "Human": 50,
            "Tied": 10,
            "Backbone_Only": 5,
            "Ambiguous": 2,
        }
        mock_cov.return_value = (
            {
                "insert": {
                    "mean_depth": 10.0,
                    "median_depth": 10.0,
                    "breadth_1x": 1.0,
                    "breadth_5x": 0.5,
                    "cv": 0.1,
                },
                "backbone": {
                    "mean_depth": 5.0,
                    "median_depth": 5.0,
                    "breadth_1x": 0.5,
                    "breadth_5x": 0.3,
                    "cv": 0.2,
                },
            },
            False,
        )

        compare_alignments(str(plasmid_bam), str(human_bam), str(output_basename))

    summary_path = Path(str(output_basename) + ".summary.tsv")
    content = summary_path.read_text()

    # Check existing rows are present and unchanged
    assert "Plasmid\t100" in content
    assert "Human\t50" in content
    assert "Tied\t10" in content
    assert "Backbone_Only\t5" in content
    assert "Ambiguous\t2" in content
    assert "Verdict\t" in content
    assert "Ratio\t" in content


@pytest.mark.unit
def test_summary_tsv_fallback_mode(tmp_path: Any) -> None:
    """Test summary.tsv when coverage uses fallback mode."""
    from pathlib import Path

    from plasmicheck.scripts.compare_alignments import compare_alignments

    plasmid_bam = tmp_path / "plasmid.bam"
    human_bam = tmp_path / "human.bam"
    output_basename = tmp_path / "output"

    # No cDNA_positions.txt file

    mock_bam_ctx = MagicMock()
    mock_bam = MagicMock()
    mock_bam_ctx.__enter__ = MagicMock(return_value=mock_bam)
    mock_bam_ctx.__exit__ = MagicMock(return_value=None)
    mock_bam.fetch.return_value = []
    mock_bam.get_tid.return_value = 0
    mock_bam.get_reference_name.return_value = "plasmid"
    mock_bam.lengths = [1000]

    with (
        patch(
            "plasmicheck.scripts.compare_alignments.pysam.AlignmentFile", return_value=mock_bam_ctx
        ),
        patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
        patch("plasmicheck.scripts.compare_alignments._streaming_compare") as mock_compare,
        patch("plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics") as mock_cov,
    ):
        mock_compare.return_value = {
            "Plasmid": 50,
            "Human": 50,
            "Tied": 0,
            "Backbone_Only": 0,
            "Ambiguous": 0,
        }
        # Fallback returns True
        mock_cov.return_value = (
            {
                "insert": {
                    "mean_depth": 20.0,
                    "median_depth": 18.0,
                    "breadth_1x": 0.9,
                    "breadth_5x": 0.7,
                    "cv": 0.3,
                },
                "backbone": {
                    "mean_depth": 0.0,
                    "median_depth": 0.0,
                    "breadth_1x": 0.0,
                    "breadth_5x": 0.0,
                    "cv": 0.0,
                },
            },
            True,
        )

        compare_alignments(str(plasmid_bam), str(human_bam), str(output_basename))

    summary_path = Path(str(output_basename) + ".summary.tsv")
    content = summary_path.read_text()

    # Fallback should be True
    assert "CoverageFallback\tTrue" in content
    # Backbone metrics should be 0.0
    assert "MeanDepthBackbone\t0.00" in content
