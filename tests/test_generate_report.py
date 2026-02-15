"""Unit tests for plasmicheck.scripts.generate_report."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from plasmicheck.scripts.generate_report import (
    downsample_data,
    extract_verdict_from_summary,
    load_data,
    main,
)


class TestDownsampleData:
    @pytest.mark.unit
    def test_under_limit_unchanged(self, sample_reads_assignment_df: pd.DataFrame) -> None:
        df = sample_reads_assignment_df  # 5 rows
        result, downsampled = downsample_data(df, 100)
        assert not downsampled
        assert len(result) == len(df)

    @pytest.mark.unit
    def test_over_limit_sampled(self) -> None:
        df = pd.DataFrame({"A": range(100)})
        result, downsampled = downsample_data(df, 10)
        assert downsampled
        assert len(result) == 10

    @pytest.mark.unit
    def test_exactly_at_limit(self) -> None:
        df = pd.DataFrame({"A": range(50)})
        result, downsampled = downsample_data(df, 50)
        assert not downsampled
        assert len(result) == 50

    @pytest.mark.unit
    def test_reproducible(self) -> None:
        df = pd.DataFrame({"A": range(100)})
        r1, _ = downsample_data(df, 10)
        r2, _ = downsample_data(df, 10)
        pd.testing.assert_frame_equal(r1.reset_index(drop=True), r2.reset_index(drop=True))


class TestExtractVerdictFromSummary:
    @pytest.mark.unit
    def test_has_verdict_row(self, sample_summary_df: pd.DataFrame) -> None:
        verdict = extract_verdict_from_summary(sample_summary_df)
        assert "contaminated" in verdict.lower()

    @pytest.mark.unit
    def test_no_verdict_row(self) -> None:
        df = pd.DataFrame(
            {
                "Category": ["Plasmid", "Human"],
                "Count": [100, 50],
            }
        )
        verdict = extract_verdict_from_summary(df)
        assert verdict == "Verdict not found in summary file"

    @pytest.mark.unit
    def test_multiple_verdict_rows(self) -> None:
        df = pd.DataFrame(
            {
                "Category": ["Verdict", "Verdict"],
                "Count": ["first verdict", "second verdict"],
            }
        )
        verdict = extract_verdict_from_summary(df)
        assert verdict == "first verdict"


class TestLoadData:
    @pytest.mark.unit
    def test_valid_tsv_files(self, tmp_path: Path) -> None:
        reads_file = tmp_path / "reads.tsv"
        summary_file = tmp_path / "summary.tsv"

        reads_file.write_text(
            "ReadID\tAssignedTo\tPlasmidScore\tHumanScore\tPlasmidCIGAR\tHumanCIGAR\tPlasmidMapQ\tHumanMapQ\n"
            "read1\tPlasmid\t70\t30\t100M\t50M50S\t60\t20\n"
        )
        summary_file.write_text(
            "Category\tCount\nPlasmid\t100\nHuman\t50\nVerdict\tSample is contaminated\n"
        )

        reads_df, summary_df = load_data(str(reads_file), str(summary_file))
        assert "ReadID" in reads_df.columns
        assert "Category" in summary_df.columns
        assert len(reads_df) == 1
        assert len(summary_df) == 3


class TestCoverageMetricsInReport:
    """Tests for coverage metrics display in HTML reports (Phase 9 Plan 02)."""

    @pytest.mark.unit
    def test_coverage_metrics_in_report(
        self,
        tmp_path: Path,
        sample_reads_assignment_df: pd.DataFrame,
        sample_summary_df_with_coverage: pd.DataFrame,
    ) -> None:
        """Verify coverage metrics are passed through to the template."""
        reads_file = tmp_path / "reads_assignment.tsv"
        summary_file = tmp_path / "summary.tsv"
        output_folder = tmp_path / "output"
        output_folder.mkdir()

        sample_reads_assignment_df.to_csv(reads_file, sep="\t", index=False)
        sample_summary_df_with_coverage.to_csv(summary_file, sep="\t", index=False)

        # Mock plot generation to avoid plotly/kaleido dependencies
        with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
            mock_plots.return_value = (
                str(output_folder / "plots" / "box.html"),
                None,
                str(output_folder / "plots" / "scatter.html"),
                None,
            )
            # Create dummy plot files so template can read them
            plots_dir = output_folder / "plots"
            plots_dir.mkdir()
            (plots_dir / "box.html").write_text("<div>box</div>")
            (plots_dir / "scatter.html").write_text("<div>scatter</div>")

            main(
                str(reads_file),
                str(summary_file),
                str(output_folder),
                command_line="test command",
            )

        # Verify interactive report was created
        report_path = output_folder / "report_interactive.html"
        assert report_path.exists()
        html_content = report_path.read_text()

        # Verify Coverage Analysis card is present
        assert "Coverage Analysis" in html_content

        # Verify at least one coverage metric value appears in the report
        assert "45.25" in html_content  # MeanDepthInsert from fixture

    @pytest.mark.unit
    def test_coverage_metrics_missing_graceful(
        self,
        tmp_path: Path,
        sample_reads_assignment_df: pd.DataFrame,
        sample_summary_df: pd.DataFrame,
    ) -> None:
        """Verify report renders without error when summary.tsv has no coverage rows."""
        reads_file = tmp_path / "reads_assignment.tsv"
        summary_file = tmp_path / "summary.tsv"
        output_folder = tmp_path / "output"
        output_folder.mkdir()

        sample_reads_assignment_df.to_csv(reads_file, sep="\t", index=False)
        sample_summary_df.to_csv(summary_file, sep="\t", index=False)

        # Mock plot generation
        with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
            mock_plots.return_value = (
                str(output_folder / "plots" / "box.html"),
                None,
                str(output_folder / "plots" / "scatter.html"),
                None,
            )
            plots_dir = output_folder / "plots"
            plots_dir.mkdir()
            (plots_dir / "box.html").write_text("<div>box</div>")
            (plots_dir / "scatter.html").write_text("<div>scatter</div>")

            main(
                str(reads_file),
                str(summary_file),
                str(output_folder),
                command_line="test command",
            )

        # Verify report was created and contains Coverage Analysis card
        report_path = output_folder / "report_interactive.html"
        assert report_path.exists()
        html_content = report_path.read_text()

        assert "Coverage Analysis" in html_content
        # Default values should be shown (0.00)
        assert "0.00" in html_content

    @pytest.mark.unit
    def test_coverage_fallback_warning(
        self,
        tmp_path: Path,
        sample_reads_assignment_df: pd.DataFrame,
    ) -> None:
        """Verify fallback warning appears when CoverageFallback is True."""
        reads_file = tmp_path / "reads_assignment.tsv"
        summary_file = tmp_path / "summary.tsv"
        output_folder = tmp_path / "output"
        output_folder.mkdir()

        sample_reads_assignment_df.to_csv(reads_file, sep="\t", index=False)

        # Create summary with CoverageFallback=True
        summary_df_fallback = pd.DataFrame(
            {
                "Category": [
                    "Plasmid",
                    "Human",
                    "Tied",
                    "Backbone_Only",
                    "Ambiguous",
                    "Verdict",
                    "Ratio",
                    "CoverageOutsideINSERT",
                    "MismatchesNearINSERT",
                    "MeanDepthInsert",
                    "MedianDepthInsert",
                    "BreadthInsert",
                    "BreadthInsert_5x",
                    "CoverageCV_Insert",
                    "MeanDepthBackbone",
                    "MedianDepthBackbone",
                    "BreadthBackbone",
                    "BreadthBackbone_5x",
                    "CoverageCV_Backbone",
                    "CoverageFallback",
                ],
                "Count": [
                    100,
                    50,
                    10,
                    15,
                    5,
                    "Sample is contaminated with plasmid DNA",
                    2.0,
                    0.1234,
                    "{'with_mismatches_or_clipping': 5, 'without_mismatches_or_clipping': 10}",
                    45.25,
                    42.00,
                    0.95,
                    0.85,
                    0.22,
                    0.00,
                    0.00,
                    0.00,
                    0.00,
                    0.00,
                    "True",  # CoverageFallback is True
                ],
            }
        )
        summary_df_fallback.to_csv(summary_file, sep="\t", index=False)

        # Mock plot generation
        with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
            mock_plots.return_value = (
                str(output_folder / "plots" / "box.html"),
                None,
                str(output_folder / "plots" / "scatter.html"),
                None,
            )
            plots_dir = output_folder / "plots"
            plots_dir.mkdir()
            (plots_dir / "box.html").write_text("<div>box</div>")
            (plots_dir / "scatter.html").write_text("<div>scatter</div>")

            main(
                str(reads_file),
                str(summary_file),
                str(output_folder),
                command_line="test command",
            )

        # Verify fallback warning appears
        report_path = output_folder / "report_interactive.html"
        assert report_path.exists()
        html_content = report_path.read_text()

        assert "Insert region not defined" in html_content
