"""Unit tests for plasmicheck.scripts.generate_report."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from plasmicheck.scripts.generate_report import (
    downsample_data,
    extract_verdict_from_summary,
    load_data,
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
