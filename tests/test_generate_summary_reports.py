"""Unit tests for plasmicheck.scripts.generate_summary_reports."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from plasmicheck.scripts.generate_summary_reports import (
    apply_sorting,
    calculate_variations,
    find_tsv_files,
)


class TestApplySorting:
    @pytest.mark.unit
    def test_none_config_unchanged(self) -> None:
        df = pd.DataFrame({"A": [3, 1, 2], "B": ["c", "a", "b"]})
        result = apply_sorting(df, None)
        pd.testing.assert_frame_equal(result, df)

    @pytest.mark.unit
    def test_single_column_ascending(self) -> None:
        df = pd.DataFrame({"A": [3, 1, 2]})
        result = apply_sorting(df, {"columns": ["A"], "ascending": [True]})
        assert list(result["A"]) == [1, 2, 3]

    @pytest.mark.unit
    def test_multi_column_mixed(self) -> None:
        df = pd.DataFrame(
            {
                "A": ["x", "x", "y", "y"],
                "B": [2, 1, 4, 3],
            }
        )
        result = apply_sorting(df, {"columns": ["A", "B"], "ascending": [True, False]})
        assert list(result["A"]) == ["x", "x", "y", "y"]
        assert list(result["B"]) == [2, 1, 4, 3]

    @pytest.mark.unit
    def test_empty_dataframe(self) -> None:
        df = pd.DataFrame({"A": [], "B": []})
        result = apply_sorting(df, {"columns": ["A"], "ascending": [True]})
        assert len(result) == 0


class TestFindTsvFiles:
    @pytest.mark.unit
    def test_finds_matching_files(self, tmp_path: Path) -> None:
        sub = tmp_path / "sample" / "plasmid"
        sub.mkdir(parents=True)
        (sub / "test.reads_assignment.tsv").write_text("header\n")
        (sub / "test.summary.tsv").write_text("header\n")
        (sub / "other.txt").write_text("ignore")

        result = find_tsv_files(str(tmp_path), ".reads_assignment.tsv")
        assert len(result) == 1
        assert result[0].endswith(".reads_assignment.tsv")

    @pytest.mark.unit
    def test_ignores_non_matching(self, tmp_path: Path) -> None:
        (tmp_path / "file.txt").write_text("nope")
        result = find_tsv_files(str(tmp_path), ".reads_assignment.tsv")
        assert result == []

    @pytest.mark.unit
    def test_empty_dir(self, tmp_path: Path) -> None:
        result = find_tsv_files(str(tmp_path), ".summary.tsv")
        assert result == []

    @pytest.mark.unit
    def test_recursive_discovery(self, tmp_path: Path) -> None:
        deep = tmp_path / "a" / "b" / "c"
        deep.mkdir(parents=True)
        (deep / "data.summary.tsv").write_text("header\n")
        result = find_tsv_files(str(tmp_path), ".summary.tsv")
        assert len(result) == 1


class TestCalculateVariations:
    @pytest.mark.unit
    def test_multiple_samples(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s2", "s3", "s1", "s2", "s3"],
                "Plasmid": ["p1", "p1", "p1", "p2", "p2", "p2"],
                "Value": [0.5, 0.6, 0.7, 0.1, 0.2, 0.3],
            }
        )
        boxplot_data, p_values_df = calculate_variations(df)
        assert "p1" in boxplot_data.columns
        assert "p2" in boxplot_data.columns
        assert not p_values_df.empty
        assert "p_value" in p_values_df.columns
        assert "p_value_corrected" in p_values_df.columns

    @pytest.mark.unit
    def test_single_sample_skipped(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Value": [0.5],
            }
        )
        _, p_values_df = calculate_variations(df)
        assert p_values_df.empty

    @pytest.mark.unit
    def test_fdr_correction_applied(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s2", "s3", "s4"],
                "Plasmid": ["p1", "p1", "p1", "p1"],
                "Value": [0.1, 0.2, 0.3, 10.0],
            }
        )
        _, p_values_df = calculate_variations(df)
        if not p_values_df.empty:
            assert "p_value_corrected" in p_values_df.columns
            # Corrected p-values should be >= raw p-values (FDR only inflates)
            for _, row in p_values_df.iterrows():
                assert row["p_value_corrected"] >= row["p_value"] - 1e-10
