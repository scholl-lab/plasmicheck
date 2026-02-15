"""Unit tests for plasmicheck.scripts.generate_summary_reports."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from plasmicheck.scripts.generate_summary_reports import (
    _sanitize_for_excel,
    apply_sorting,
    calculate_variations,
    create_metric_heatmap_data,
    create_resistance_gene_table,
    find_tsv_files,
    plot_coverage_outside_insert_boxplot,
    plot_metric_heatmap,
    save_tables_as_tsv_and_excel,
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


class TestSanitizeForExcel:
    @pytest.mark.unit
    def test_converts_bool_columns_to_strings(self) -> None:
        df = pd.DataFrame({"name": ["a", "b"], "flag": [True, False]})
        result = _sanitize_for_excel(df)
        # pandas 2.x may return StringDtype instead of object; check values
        assert result["flag"].dtype != bool
        assert list(result["flag"]) == ["True", "False"]

    @pytest.mark.unit
    def test_preserves_non_bool_columns(self) -> None:
        df = pd.DataFrame({"val": [1, 2], "name": ["a", "b"]})
        result = _sanitize_for_excel(df)
        pd.testing.assert_frame_equal(result, df)

    @pytest.mark.unit
    def test_does_not_mutate_original(self) -> None:
        df = pd.DataFrame({"flag": [True, False]})
        _sanitize_for_excel(df)
        assert df["flag"].dtype == bool


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


class TestCreateMetricHeatmapData:
    @pytest.mark.unit
    def test_extracts_metric_and_pivots(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s1", "s2", "s2"],
                "Plasmid": ["p1", "p2", "p1", "p2"],
                "Category": [
                    "MeanDepthBackbone",
                    "MeanDepthBackbone",
                    "MeanDepthBackbone",
                    "MeanDepthBackbone",
                ],
                "Value": [10.5, 20.3, 15.2, 25.8],
            }
        )
        result = create_metric_heatmap_data(df, "MeanDepthBackbone")
        assert result.shape == (2, 2)  # 2 samples x 2 plasmids
        assert result.loc["s1", "p1"] == 10.5
        assert result.loc["s2", "p2"] == 25.8

    @pytest.mark.unit
    def test_empty_for_missing_category(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["MeanDepthBackbone"],
                "Value": [10.0],
            }
        )
        result = create_metric_heatmap_data(df, "NonExistent")
        assert result.empty

    @pytest.mark.unit
    def test_coerces_non_numeric(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s2"],
                "Plasmid": ["p1", "p1"],
                "Category": ["MeanDepthBackbone", "MeanDepthBackbone"],
                "Value": [10.0, "bad"],
            }
        )
        result = create_metric_heatmap_data(df, "MeanDepthBackbone")
        assert pd.isna(result.loc["s2", "p1"])

    @pytest.mark.unit
    def test_handles_single_sample(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s1"],
                "Plasmid": ["p1", "p2"],
                "Category": ["BreadthBackbone", "BreadthBackbone"],
                "Value": [95.0, 98.5],
            }
        )
        result = create_metric_heatmap_data(df, "BreadthBackbone")
        assert result.shape == (1, 2)  # 1 sample, 2 plasmids


class TestCreateResistanceGeneTable:
    @pytest.mark.unit
    def test_extracts_genes_correctly(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": [
                    "s1",
                    "s1",
                    "s1",
                    "s1",
                    "s1",
                    "s1",
                ],
                "Plasmid": ["p1", "p1", "p1", "p1", "p1", "p1"],
                "Category": [
                    "ResistanceGenesAvailable",
                    "ResistanceGenes",
                    "ResGene_AmpR_MeanDepth",
                    "ResGene_AmpR_Breadth",
                    "ResGene_KanR_MeanDepth",
                    "ResGene_KanR_Breadth",
                ],
                "Value": ["true", "AmpR,KanR", "50.5", "95.2", "30.1", "88.0"],
            }
        )
        result = create_resistance_gene_table(df)
        assert len(result) == 2
        assert list(result.columns) == [
            "Sample",
            "Plasmid",
            "Gene Name",
            "Mean Depth",
            "Breadth (%)",
        ]
        assert result.iloc[0]["Gene Name"] == "AmpR"
        assert result.iloc[0]["Mean Depth"] == 50.5
        assert result.iloc[1]["Gene Name"] == "KanR"
        assert result.iloc[1]["Breadth (%)"] == 88.0

    @pytest.mark.unit
    def test_empty_when_no_genbank(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["ResistanceGenesAvailable"],
                "Value": ["false"],
            }
        )
        result = create_resistance_gene_table(df)
        assert result.empty
        assert list(result.columns) == [
            "Sample",
            "Plasmid",
            "Gene Name",
            "Mean Depth",
            "Breadth (%)",
        ]

    @pytest.mark.unit
    def test_empty_when_no_resistance_genes_available_row(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["SomeOtherCategory"],
                "Value": ["123"],
            }
        )
        result = create_resistance_gene_table(df)
        assert result.empty

    @pytest.mark.unit
    def test_skips_samples_with_no_genes(self) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s1", "s1", "s1", "s2", "s2"],
                "Plasmid": ["p1", "p1", "p1", "p1", "p1", "p1"],
                "Category": [
                    "ResistanceGenesAvailable",
                    "ResistanceGenes",
                    "ResGene_AmpR_MeanDepth",
                    "ResGene_AmpR_Breadth",
                    "ResistanceGenesAvailable",
                    "ResistanceGenes",
                ],
                "Value": ["true", "AmpR", "40.0", "90.0", "true", ""],
            }
        )
        result = create_resistance_gene_table(df)
        assert len(result) == 1
        assert result.iloc[0]["Sample"] == "s1"


class TestPlotMetricHeatmap:
    @pytest.mark.unit
    def test_creates_html_file(self, tmp_path: Path) -> None:
        df = pd.DataFrame({"p1": [10.0, 15.0], "p2": [20.0, 25.0]}, index=["s1", "s2"])
        html_file, png_file = plot_metric_heatmap(
            df, str(tmp_path), "Test Metric", static_report=False
        )
        assert Path(html_file).exists()
        assert png_file is None

    @pytest.mark.unit
    def test_empty_data_returns_empty_strings(self, tmp_path: Path) -> None:
        df = pd.DataFrame()
        html_file, png_file = plot_metric_heatmap(
            df, str(tmp_path), "Test Metric", static_report=False
        )
        assert html_file == ""
        assert png_file is None

    @pytest.mark.unit
    def test_vmax_clips_data(self, tmp_path: Path) -> None:
        df = pd.DataFrame({"p1": [100.0, 200.0]}, index=["s1", "s2"])
        html_file, png_file = plot_metric_heatmap(
            df, str(tmp_path), "Test Metric", vmax=50.0, static_report=False
        )
        assert Path(html_file).exists()
        assert png_file is None


class TestPlotCoverageOutsideInsertBoxplot:
    @pytest.mark.unit
    def test_creates_html_file(self, tmp_path: Path) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s2"],
                "Plasmid": ["p1", "p1"],
                "Category": ["CoverageOutsideINSERT", "CoverageOutsideINSERT"],
                "Value": [5.0, 10.0],
            }
        )
        html_file, png_file = plot_coverage_outside_insert_boxplot(
            df, str(tmp_path), static_report=False
        )
        assert Path(html_file).exists()
        assert png_file is None

    @pytest.mark.unit
    def test_empty_returns_empty(self, tmp_path: Path) -> None:
        df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["SomeOtherCategory"],
                "Value": [10.0],
            }
        )
        html_file, png_file = plot_coverage_outside_insert_boxplot(
            df, str(tmp_path), static_report=False
        )
        assert html_file == ""
        assert png_file is None


class TestSaveTablesExtended:
    @pytest.mark.unit
    def test_extended_metrics_tsv_created(self, tmp_path: Path) -> None:
        combined_df = pd.DataFrame(
            {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Plasmid"], "Value": [100]}
        )
        verdict_df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["Verdict"],
                "Value": ["Contaminated"],
            }
        )
        ratio_df = pd.DataFrame(
            {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Ratio"], "Value": [1.5]}
        )
        p_values_df = pd.DataFrame(
            {"Plasmid": ["p1"], "p_value": [0.05], "p_value_corrected": [0.05]}
        )
        extended_df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["MeanDepthBackbone"],
                "Value": [50.0],
            }
        )

        save_tables_as_tsv_and_excel(
            combined_df,
            verdict_df,
            ratio_df,
            p_values_df,
            str(tmp_path),
            extended_metrics_df=extended_df,
        )

        tsv_path = tmp_path / "data" / "extended_metrics_table.tsv"
        assert tsv_path.exists()

    @pytest.mark.unit
    def test_resistance_genes_tsv_created(self, tmp_path: Path) -> None:
        combined_df = pd.DataFrame(
            {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Plasmid"], "Value": [100]}
        )
        verdict_df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["Verdict"],
                "Value": ["Contaminated"],
            }
        )
        ratio_df = pd.DataFrame(
            {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Ratio"], "Value": [1.5]}
        )
        p_values_df = pd.DataFrame(
            {"Plasmid": ["p1"], "p_value": [0.05], "p_value_corrected": [0.05]}
        )
        resistance_df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Gene Name": ["AmpR"],
                "Mean Depth": [40.0],
                "Breadth (%)": [90.0],
            }
        )

        save_tables_as_tsv_and_excel(
            combined_df,
            verdict_df,
            ratio_df,
            p_values_df,
            str(tmp_path),
            resistance_gene_df=resistance_df,
        )

        tsv_path = tmp_path / "data" / "resistance_genes_table.tsv"
        assert tsv_path.exists()

    @pytest.mark.unit
    def test_excel_has_new_sheets(self, tmp_path: Path) -> None:
        combined_df = pd.DataFrame(
            {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Plasmid"], "Value": [100]}
        )
        verdict_df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["Verdict"],
                "Value": ["Contaminated"],
            }
        )
        ratio_df = pd.DataFrame(
            {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Ratio"], "Value": [1.5]}
        )
        p_values_df = pd.DataFrame(
            {"Plasmid": ["p1"], "p_value": [0.05], "p_value_corrected": [0.05]}
        )
        extended_df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["MeanDepthBackbone"],
                "Value": [50.0],
            }
        )
        resistance_df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Gene Name": ["AmpR"],
                "Mean Depth": [40.0],
                "Breadth (%)": [90.0],
            }
        )

        save_tables_as_tsv_and_excel(
            combined_df,
            verdict_df,
            ratio_df,
            p_values_df,
            str(tmp_path),
            extended_metrics_df=extended_df,
            resistance_gene_df=resistance_df,
        )

        excel_path = tmp_path / "data" / "summary_report_tables.xlsx"
        assert excel_path.exists()

        # Read Excel and check sheets exist
        import openpyxl

        wb = openpyxl.load_workbook(excel_path)
        sheet_names = wb.sheetnames
        assert "Extended Metrics" in sheet_names
        assert "Resistance Genes" in sheet_names

    @pytest.mark.unit
    def test_backward_compatible_without_new_args(self, tmp_path: Path) -> None:
        combined_df = pd.DataFrame(
            {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Plasmid"], "Value": [100]}
        )
        verdict_df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["Verdict"],
                "Value": ["Contaminated"],
            }
        )
        ratio_df = pd.DataFrame(
            {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Ratio"], "Value": [1.5]}
        )
        p_values_df = pd.DataFrame(
            {"Plasmid": ["p1"], "p_value": [0.05], "p_value_corrected": [0.05]}
        )

        # Call without new optional args (like pre-v0.33.0 code)
        save_tables_as_tsv_and_excel(combined_df, verdict_df, ratio_df, p_values_df, str(tmp_path))

        data_dir = tmp_path / "data"
        assert (data_dir / "combined_table.tsv").exists()
        assert (data_dir / "verdict_table.tsv").exists()
        assert (data_dir / "ratio_table.tsv").exists()
        assert (data_dir / "p_value_table.tsv").exists()
        assert (data_dir / "summary_report_tables.xlsx").exists()
        # New files should NOT exist
        assert not (data_dir / "extended_metrics_table.tsv").exists()
        assert not (data_dir / "resistance_genes_table.tsv").exists()
