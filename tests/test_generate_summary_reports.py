"""Unit tests for plasmicheck.scripts.generate_summary_reports."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from plasmicheck.scripts.generate_summary_reports import (
    PLOT_CONFIG,
    UNCLEAR_RANGE,
    _add_plasmid_flag,
    _extract_sample_plasmid,
    _format_verdict_hover,
    _reestimate_cohort_depth_fractions,
    _sanitize_for_excel,
    aggregate_verdicts,
    apply_sorting,
    calculate_Metrics,
    calculate_variations,
    categorize_value,
    compute_plasmid_flag,
    create_metric_heatmap_data,
    create_resistance_gene_table,
    find_tsv_files,
    generate_report,
    plot_boxplot,
    plot_coverage_outside_insert_boxplot,
    plot_heatmap,
    plot_metric_heatmap,
    plot_verdict_bar_by_plasmid,
    plot_verdict_bar_by_sample,
    read_compare_outputs,
    save_tables_as_tsv_and_excel,
)


@pytest.mark.unit
def test_ratio_plots_disable_cloud_sharing(tmp_path: Path) -> None:
    ratios = pd.DataFrame({"Sample": ["s1"], "Plasmid": ["p1"], "Value": [0.1]})
    pivot = ratios.pivot(index="Sample", columns="Plasmid", values="Value")
    outputs = [
        plot_boxplot(pivot, str(tmp_path)),
        plot_heatmap(ratios, str(tmp_path), 0.8, UNCLEAR_RANGE, PLOT_CONFIG),
    ]
    for html_file, _ in outputs:
        assert '"showSendToCloud": false' in Path(html_file).read_text()


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
        assert '"showSendToCloud": false' in Path(html_file).read_text()
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
        assert '"showSendToCloud": false' in Path(html_file).read_text()
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
        assert '"showSendToCloud": false' in Path(html_file).read_text()
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


class TestExtractSamplePlasmid:
    @pytest.mark.unit
    def test_flat_layout(self) -> None:
        """Test extraction from flat directory structure (old)."""
        file_path = "/output/sampleA/plasmidX/comparison_result.summary.tsv"
        sample, plasmid = _extract_sample_plasmid(file_path)
        assert sample == "sampleA"
        assert plasmid == "plasmidX"

    @pytest.mark.unit
    def test_subfolder_layout(self) -> None:
        """Test extraction from data/ subfolder structure (new)."""
        file_path = "/output/sampleA/plasmidX/data/comparison_result.summary.tsv"
        sample, plasmid = _extract_sample_plasmid(file_path)
        assert sample == "sampleA"
        assert plasmid == "plasmidX"

    @pytest.mark.unit
    def test_reads_assignment_flat(self) -> None:
        """Test extraction from flat layout with reads_assignment file."""
        file_path = "/output/sample1/plasmid1/comparison_result.reads_assignment.tsv"
        sample, plasmid = _extract_sample_plasmid(file_path)
        assert sample == "sample1"
        assert plasmid == "plasmid1"

    @pytest.mark.unit
    def test_reads_assignment_subfolder(self) -> None:
        """Test extraction from subfolder layout with reads_assignment file."""
        file_path = "/output/sample1/plasmid1/data/comparison_result.reads_assignment.tsv"
        sample, plasmid = _extract_sample_plasmid(file_path)
        assert sample == "sample1"
        assert plasmid == "plasmid1"


class TestReadCompareOutputsBothLayouts:
    @pytest.mark.unit
    def test_flat_layout_read(self, tmp_path: Path) -> None:
        """Test read_compare_outputs with flat directory structure."""
        # Create flat structure: output/sample/plasmid/file.tsv
        sample_dir = tmp_path / "sampleA" / "plasmidX"
        sample_dir.mkdir(parents=True)

        # Write minimal reads_assignment TSV
        reads_file = sample_dir / "comparison_result.reads_assignment.tsv"
        reads_file.write_text("ReadName\tAssignedTo\nread1\tPlasmid\n")

        # Write minimal summary TSV
        summary_file = sample_dir / "comparison_result.summary.tsv"
        summary_file.write_text("Category\tMetric\tCount\nVerdict\tContaminated\t1\n")

        # Read outputs
        summary_df = read_compare_outputs(str(tmp_path))

        assert len(summary_df) == 1
        assert summary_df.iloc[0]["Sample"] == "sampleA"
        assert summary_df.iloc[0]["Plasmid"] == "plasmidX"

    @pytest.mark.unit
    def test_subfolder_layout_read(self, tmp_path: Path) -> None:
        """Test read_compare_outputs with data/ subfolder structure."""
        # Create subfolder structure: output/sample/plasmid/data/file.tsv
        data_dir = tmp_path / "sampleB" / "plasmidY" / "data"
        data_dir.mkdir(parents=True)

        # Write minimal reads_assignment TSV
        reads_file = data_dir / "comparison_result.reads_assignment.tsv"
        reads_file.write_text("ReadName\tAssignedTo\nread2\tHuman\n")

        # Write minimal summary TSV
        summary_file = data_dir / "comparison_result.summary.tsv"
        summary_file.write_text("Category\tMetric\tCount\nVerdict\tNotContaminated\t1\n")

        # Read outputs
        summary_df = read_compare_outputs(str(tmp_path))

        assert len(summary_df) == 1
        assert summary_df.iloc[0]["Sample"] == "sampleB"
        assert summary_df.iloc[0]["Plasmid"] == "plasmidY"

    @pytest.mark.unit
    def test_mixed_layout(self, tmp_path: Path) -> None:
        """Test read_compare_outputs with mixed flat and subfolder layouts."""
        # Create flat layout for sample1
        flat_dir = tmp_path / "sample1" / "plasmid1"
        flat_dir.mkdir(parents=True)
        (flat_dir / "comparison_result.reads_assignment.tsv").write_text(
            "ReadName\tAssignedTo\nread1\tPlasmid\n"
        )
        (flat_dir / "comparison_result.summary.tsv").write_text(
            "Category\tMetric\tCount\nVerdict\tContaminated\t1\n"
        )

        # Create subfolder layout for sample2
        subfolder_dir = tmp_path / "sample2" / "plasmid2" / "data"
        subfolder_dir.mkdir(parents=True)
        (subfolder_dir / "comparison_result.reads_assignment.tsv").write_text(
            "ReadName\tAssignedTo\nread2\tHuman\n"
        )
        (subfolder_dir / "comparison_result.summary.tsv").write_text(
            "Category\tMetric\tCount\nVerdict\tNotContaminated\t1\n"
        )

        # Read outputs
        summary_df = read_compare_outputs(str(tmp_path))

        assert len(summary_df) == 2
        assert set(summary_df["Sample"]) == {"sample1", "sample2"}
        assert set(summary_df["Plasmid"]) == {"plasmid1", "plasmid2"}


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
        openpyxl = pytest.importorskip("openpyxl")

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


class TestAggregateVerdicts:
    @pytest.mark.unit
    def test_maps_verdicts_correctly(self) -> None:
        """Test verdict string mapping to short labels."""
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s2", "s3"],
                "Plasmid": ["p1", "p1", "p1"],
                "Category": ["Verdict", "Verdict", "Verdict"],
                "Value": [
                    "Sample is contaminated with plasmid DNA",
                    "Sample contamination status is unclear",
                    "Sample is not contaminated with plasmid DNA",
                ],
            }
        )

        full_df, by_plasmid, _by_sample = aggregate_verdicts(df)

        # Check mappings
        assert set(full_df["VerdictShort"]) == {"Contaminated", "Unclear", "Clean"}
        assert len(full_df) == 3

        # Check by_plasmid aggregation
        assert by_plasmid.loc["p1", "Contaminated"] == 1
        assert by_plasmid.loc["p1", "Unclear"] == 1
        assert by_plasmid.loc["p1", "Clean"] == 1
        assert by_plasmid.loc["p1", "No Result"] == 0

    @pytest.mark.unit
    def test_detects_missing_pairs_as_no_result(self) -> None:
        """Test No Result detection for missing sample-plasmid combinations."""
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s1", "s2"],
                "Plasmid": ["p1", "p2", "p1"],
                "Category": ["Verdict", "Verdict", "Verdict"],
                "Value": [
                    "Sample is contaminated with plasmid DNA",
                    "Sample is not contaminated with plasmid DNA",
                    "Sample is not contaminated with plasmid DNA",
                ],
            }
        )

        full_df, _by_plasmid, _by_sample = aggregate_verdicts(df)

        # s2 x p2 should be missing -> No Result
        assert len(full_df) == 4  # 2 samples x 2 plasmids
        no_result_rows = full_df[full_df["VerdictShort"] == "No Result"]
        assert len(no_result_rows) == 1
        assert no_result_rows.iloc[0]["Sample"] == "s2"
        assert no_result_rows.iloc[0]["Plasmid"] == "p2"

    @pytest.mark.unit
    def test_column_order_preserved(self) -> None:
        """Test that verdict columns are in correct order."""
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s2"],
                "Plasmid": ["p1", "p1"],
                "Category": ["Verdict", "Verdict"],
                "Value": [
                    "Sample is contaminated with plasmid DNA",
                    "Sample is not contaminated with plasmid DNA",
                ],
            }
        )

        _, by_plasmid, by_sample = aggregate_verdicts(df)

        expected_order = ["Clean", "Unclear", "Contaminated", "No Result"]
        assert list(by_plasmid.columns) == expected_order
        assert list(by_sample.columns) == expected_order

    @pytest.mark.unit
    def test_empty_verdict_data_returns_empty(self) -> None:
        """Test that empty verdict data returns empty DataFrames."""
        df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["SomeOtherCategory"],
                "Value": ["123"],
            }
        )

        full_df, by_plasmid, by_sample = aggregate_verdicts(df)

        assert full_df.empty
        assert by_plasmid.empty
        assert by_sample.empty


class TestFormatVerdictHover:
    @pytest.mark.unit
    def test_small_list_shows_all(self) -> None:
        """Test hover with <= 10 names shows all names."""
        names = ["sample1", "sample2", "sample3"]
        result = _format_verdict_hover("Clean", 3, names)
        assert "Clean" in result
        assert "Count: 3" in result
        assert "sample1" in result
        assert "sample2" in result
        assert "sample3" in result
        assert "...and" not in result

    @pytest.mark.unit
    def test_large_list_truncated(self) -> None:
        """Test hover with > 10 names shows first 10 + truncation."""
        names = [f"sample{i}" for i in range(15)]
        result = _format_verdict_hover("Contaminated", 15, names, max_display=10)
        assert "Contaminated" in result
        assert "Count: 15" in result
        assert "sample0" in result
        assert "sample9" in result
        assert "...and 5 more" in result
        assert "sample14" not in result

    @pytest.mark.unit
    def test_empty_list_shows_only_verdict_and_count(self) -> None:
        """Test hover with empty name list shows verdict and count only."""
        result = _format_verdict_hover("Unclear", 0, [])
        assert "Unclear" in result
        assert "Count: 0" in result
        # Should not have extra lines
        assert result.count("<br>") == 1  # Only one <br> between verdict and count


class TestPlotVerdictBarByPlasmid:
    @pytest.mark.unit
    def test_creates_html_file(self, tmp_path: Path) -> None:
        """Test plot_verdict_bar_by_plasmid creates HTML output."""
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s2", "s3"],
                "Plasmid": ["p1", "p1", "p1"],
                "Category": ["Verdict", "Verdict", "Verdict"],
                "Value": [
                    "Sample is contaminated with plasmid DNA",
                    "Sample contamination status is unclear",
                    "Sample is not contaminated with plasmid DNA",
                ],
            }
        )

        html_file, png_file = plot_verdict_bar_by_plasmid(df, str(tmp_path), static_report=False)

        assert Path(html_file).exists()
        assert '"showSendToCloud": false' in Path(html_file).read_text()
        assert png_file is None

    @pytest.mark.unit
    def test_empty_verdict_data_returns_empty_string(self, tmp_path: Path) -> None:
        """Test empty verdict data returns empty paths."""
        df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["SomeOther"],
                "Value": ["123"],
            }
        )

        html_file, png_file = plot_verdict_bar_by_plasmid(df, str(tmp_path), static_report=False)

        assert html_file == ""
        assert png_file is None

    @pytest.mark.unit
    def test_static_report_triggers_matplotlib(self, tmp_path: Path) -> None:
        """Test static_report=True with matplotlib backend creates PNG."""
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s2"],
                "Plasmid": ["p1", "p1"],
                "Category": ["Verdict", "Verdict"],
                "Value": [
                    "Sample is contaminated with plasmid DNA",
                    "Sample is not contaminated with plasmid DNA",
                ],
            }
        )

        html_file, png_file = plot_verdict_bar_by_plasmid(
            df, str(tmp_path), static_report=True, plot_backend="matplotlib"
        )

        assert Path(html_file).exists()
        assert '"showSendToCloud": false' in Path(html_file).read_text()
        assert png_file is not None
        assert Path(png_file).exists()


class TestPlotVerdictBarBySample:
    @pytest.mark.unit
    def test_creates_html_file(self, tmp_path: Path) -> None:
        """Test plot_verdict_bar_by_sample creates HTML output."""
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s1", "s1"],
                "Plasmid": ["p1", "p2", "p3"],
                "Category": ["Verdict", "Verdict", "Verdict"],
                "Value": [
                    "Sample is contaminated with plasmid DNA",
                    "Sample contamination status is unclear",
                    "Sample is not contaminated with plasmid DNA",
                ],
            }
        )

        html_file, png_file = plot_verdict_bar_by_sample(df, str(tmp_path), static_report=False)

        assert Path(html_file).exists()
        assert '"showSendToCloud": false' in Path(html_file).read_text()
        assert png_file is None

    @pytest.mark.unit
    def test_empty_verdict_data_returns_empty_string(self, tmp_path: Path) -> None:
        """Test empty verdict data returns empty paths."""
        df = pd.DataFrame(
            {
                "Sample": ["s1"],
                "Plasmid": ["p1"],
                "Category": ["SomeOther"],
                "Value": ["123"],
            }
        )

        html_file, png_file = plot_verdict_bar_by_sample(df, str(tmp_path), static_report=False)

        assert html_file == ""
        assert png_file is None

    @pytest.mark.unit
    def test_static_report_triggers_matplotlib(self, tmp_path: Path) -> None:
        """Test static_report=True with matplotlib backend creates PNG."""
        df = pd.DataFrame(
            {
                "Sample": ["s1", "s2"],
                "Plasmid": ["p1", "p1"],
                "Category": ["Verdict", "Verdict"],
                "Value": [
                    "Sample is contaminated with plasmid DNA",
                    "Sample is not contaminated with plasmid DNA",
                ],
            }
        )

        html_file, png_file = plot_verdict_bar_by_sample(
            df, str(tmp_path), static_report=True, plot_backend="matplotlib"
        )

        assert Path(html_file).exists()
        assert '"showSendToCloud": false' in Path(html_file).read_text()
        assert png_file is not None
        assert Path(png_file).exists()


@pytest.mark.unit
def test_calculate_metrics_one_sample_and_zero_variance() -> None:
    """Test calculate_Metrics handles single-sample and zero variance gracefully."""
    # 1. Single sample should return empty p-values DataFrame with expected columns
    df_single = pd.DataFrame({"Sample": ["s1"], "Plasmid": ["p1"], "Value": [0.5]})
    _boxplot_data, p_values_df = calculate_Metrics(df_single)
    assert p_values_df.empty
    assert list(p_values_df.columns) == ["Sample", "Plasmid", "p_value", "p_value_corrected"]

    # 2. Multi-sample where p1 has zero variance (identical values) and p2 has varying values
    df_multi = pd.DataFrame(
        {
            "Sample": ["s1", "s2", "s3", "s1", "s2", "s3"],
            "Plasmid": ["p1", "p1", "p1", "p2", "p2", "p2"],
            "Value": [1.0, 1.0, 1.0, 0.1, 0.4, 0.8],
        }
    )
    _boxplot_data, p_values_df = calculate_Metrics(df_multi)
    assert not p_values_df.empty
    p1_rows = p_values_df[p_values_df["Plasmid"] == "p1"]
    p2_rows = p_values_df[p_values_df["Plasmid"] == "p2"]

    # For p1, zero variance leaves p_value as NaN
    assert p1_rows["p_value"].isna().all()
    assert p1_rows["p_value_corrected"].isna().all()

    # For p2, p_value and p_value_corrected are computed and valid
    assert p2_rows["p_value"].notna().all()
    assert p2_rows["p_value_corrected"].notna().all()


@pytest.mark.unit
def test_mismatches_near_insert_dict_string_parsed() -> None:
    """Test that MismatchesNearINSERT dict string is parsed to numeric mismatch count."""
    df = pd.DataFrame(
        {
            "Sample": ["s1", "s2"],
            "Plasmid": ["p1", "p1"],
            "Category": ["MismatchesNearINSERT", "MismatchesNearINSERT"],
            "Value": [
                "{'with_mismatches_or_clipping': 8, 'without_mismatches_or_clipping': 12}",
                "{'with_mismatches_or_clipping': 3, 'without_mismatches_or_clipping': 17}",
            ],
        }
    )
    result = create_metric_heatmap_data(df, "MismatchesNearINSERT")
    assert result.loc["s1", "p1"] == 8.0
    assert result.loc["s2", "p1"] == 3.0


@pytest.mark.unit
def test_create_resistance_gene_table_sanitizes_gene_names() -> None:
    """Test that resistance gene table sanitizes gene names when looking up metrics."""
    df = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s1", "s1"],
            "Plasmid": ["p1", "p1", "p1", "p1"],
            "Category": [
                "ResistanceGenesAvailable",
                "ResistanceGenes",
                "ResGene_bla_TEM-1__MeanDepth",
                "ResGene_bla_TEM-1__Breadth",
            ],
            "Value": ["true", "bla(TEM-1)", "42.0", "99.5"],
        }
    )
    result = create_resistance_gene_table(df)
    assert len(result) == 1
    assert result.iloc[0]["Gene Name"] == "bla(TEM-1)"
    assert result.iloc[0]["Mean Depth"] == 42.0
    assert result.iloc[0]["Breadth (%)"] == 99.5


@pytest.mark.unit
def test_categorize_value_boundary_consistency() -> None:
    """Test categorize_value matches verdict boundaries with inclusive unclear range."""
    unclear = {"lower_bound": 0.05, "upper_bound": 0.10}
    threshold = 0.10

    assert categorize_value(0.04, threshold, unclear) == "Clean"
    assert categorize_value(0.05, threshold, unclear) == "Unclear"
    assert categorize_value(0.08, threshold, unclear) == "Unclear"
    assert categorize_value(0.10, threshold, unclear) == "Unclear"
    assert categorize_value(0.11, threshold, unclear) == "Contaminated"

    # Default config: threshold=0.8 is midpoint of unclear_range=[0.6, 1.0]
    default_unclear = {"lower_bound": 0.6, "upper_bound": 1.0}
    assert categorize_value(0.59, 0.8, default_unclear) == "Clean"
    assert categorize_value(0.60, 0.8, default_unclear) == "Unclear"
    assert categorize_value(0.80, 0.8, default_unclear) == "Unclear"
    assert categorize_value(0.95, 0.8, default_unclear) == "Unclear"
    assert categorize_value(1.00, 0.8, default_unclear) == "Unclear"
    assert categorize_value(1.01, 0.8, default_unclear) == "Contaminated"


@pytest.mark.unit
def test_summary_report_autoescape(tmp_path: Path) -> None:
    """Test that summary_template.html autoescapes potentially dangerous content."""
    heatmap_html = tmp_path / "heatmap.html"
    boxplot_html = tmp_path / "boxplot.html"
    heatmap_html.write_text("<div>Heatmap</div>")
    boxplot_html.write_text("<div>Boxplot</div>")

    combined_df = pd.DataFrame(
        {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Verdict"], "Value": ["Clean"]}
    )
    verdict_df = pd.DataFrame(
        {"Sample": ["s1"], "Plasmid": ["p1"], "Category": ["Verdict"], "Value": ["Clean"]}
    )
    ratio_df = pd.DataFrame({"Sample": ["s1"], "Plasmid": ["p1"], "Value": [0.01]})
    p_value_tsv = tmp_path / "p_values.tsv"
    p_value_tsv.write_text("Sample\tPlasmid\tp_value\tp_value_corrected\n")

    generate_report(
        combined_df=combined_df,
        verdict_df=verdict_df,
        ratio_df=ratio_df,
        heatmap_filename_interactive=str(heatmap_html),
        heatmap_filename_png=None,
        boxplot_filename_interactive=str(boxplot_html),
        boxplot_filename_png=None,
        p_value_table_filename=str(p_value_tsv),
        output_folder=str(tmp_path),
        threshold=0.8,
        unclear_range={"lower_bound": 0.5, "upper_bound": 0.8},
        command_line="<script>alert('cmd')</script>",
    )

    report_file = tmp_path / "summary_report_interactive.html"
    assert report_file.exists()
    content = report_file.read_text()
    assert "<script>alert('cmd')</script>" not in content
    assert "&lt;script&gt;alert(&#39;cmd&#39;)&lt;/script&gt;" in content


@pytest.mark.unit
def test_truncate_list_for_html_semantics() -> None:
    """Verify _truncate_list_for_html normalizes missing values to N/A."""
    from plasmicheck.scripts.generate_summary_reports import _truncate_list_for_html

    assert _truncate_list_for_html("") == "N/A"
    assert _truncate_list_for_html("NA") == "N/A"
    assert _truncate_list_for_html("nan") == "N/A"
    assert _truncate_list_for_html("1,2,3,4,5") == "1,2,3,4,5"
    assert _truncate_list_for_html("1,2,3,4,5,6") == "1,2,3,4,5,..."


@pytest.mark.unit
def test_read_compare_outputs_returns_single_dataframe(tmp_path: Path) -> None:
    """Verify read_compare_outputs returns a single DataFrame rather than a tuple."""
    sample_dir = tmp_path / "sampleA" / "plasmidX"
    sample_dir.mkdir(parents=True)
    summary_file = sample_dir / "comparison_result.summary.tsv"
    summary_file.write_text("Category\tMetric\tCount\nVerdict\tClean\t1\n")

    result: object = read_compare_outputs(str(tmp_path))
    assert not isinstance(result, tuple)
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 1


@pytest.mark.unit
def test_plot_boxplot_custom_threshold(tmp_path: Path) -> None:
    """Verify plot_boxplot accepts custom threshold and unclear range."""
    boxplot_data = pd.DataFrame({"Plasmid1": [0.1, 0.5, 0.9]}, index=["s1", "s2", "s3"])
    unclear = {"lower_bound": 0.5, "upper_bound": 1.2}

    interactive_html, _png_path = plot_boxplot(
        boxplot_data,
        str(tmp_path),
        threshold=1.2,
        unclear_range=unclear,
    )
    assert Path(interactive_html).exists()
    content = Path(interactive_html).read_text()
    assert "Unclear (0.5-1.2)" in content
    assert any(
        s in content
        for s in ("Contaminated (> 1.2)", "Contaminated (&gt; 1.2)", r"Contaminated (\u003e 1.2)")
    )


@pytest.mark.unit
def test_summary_reports_main_effective_threshold_reclassification(tmp_path: Path) -> None:
    """Verify main() reclassifies stored verdicts using effective threshold from TSV."""
    from plasmicheck.scripts.generate_summary_reports import main

    sample_dir = tmp_path / "sample1" / "plasmid1"
    sample_dir.mkdir(parents=True)
    summary_file = sample_dir / "comparison_result.summary.tsv"
    # Stored verdict is Contaminated (based on default 0.8), Ratio is 0.9, but Threshold in TSV is 1.2
    summary_file.write_text(
        "Category\tMetric\tCount\n"
        "Ratio\tRatio\t0.9\n"
        "Threshold\tThreshold\t1.2\n"
        "Verdict\tVerdict\tSample is contaminated with plasmid DNA\n"
    )

    out_dir = tmp_path / "output"
    out_dir.mkdir()

    main(str(tmp_path), str(out_dir))

    report_file = out_dir / "summary_report_interactive.html"
    assert report_file.exists()
    content = report_file.read_text()
    # With threshold 1.2 and unclear 0.6..1.2, ratio 0.9 is Unclear
    assert "Unclear" in content


@pytest.mark.unit
def test_summary_report_cohort_matrix_and_links(tmp_path: Path) -> None:
    """Verify summary report renders unified cohort matrix, sample report links, and readable junction headers."""
    from plasmicheck.scripts.generate_summary_reports import main

    sample_dir = tmp_path / "sample1_merged" / "plasmid1"
    reports_dir = sample_dir / "reports"
    reports_dir.mkdir(parents=True)
    (reports_dir / "report_interactive.html").write_text("<html>sample report</html>")
    summary_file = sample_dir / "comparison_result.summary.tsv"
    summary_file.write_text(
        "Category\tCount\n"
        "Plasmid\t1200\n"
        "Human\t300\n"
        "Tied\t50\n"
        "Ratio\t4.0\n"
        "Threshold\t0.8\n"
        "Verdict\tSample is contaminated with plasmid DNA\n"
        "MeanDepthBackbone\t150.5\n"
        "BreadthBackbone\t0.98\n"
        "SpliceJunctionCount\t5\n"
        "ExonCount\t6\n"
        "SpliceJunctionReads\t1100\n"
        "SpliceJunctionEvents\t2200\n"
        "SpliceJunctionsCovered\t5\n"
        "ExonIntronBoundaryReads\t40\n"
        "EstimatedContaminationFraction\t0.9650\n"
    )

    out_dir = tmp_path / "summary"
    out_dir.mkdir()

    main(str(tmp_path), str(out_dir), substring_to_remove="_merged")

    report_file = out_dir / "summary_report_interactive.html"
    assert report_file.exists()
    content = report_file.read_text()
    assert "cohortOverviewTable" in content
    assert "96.5%" in content
    assert "150.5x (98%)" in content
    assert "2,200 / 40" in content
    assert "../sample1_merged/plasmid1/reports/report_interactive.html" in content
    assert "Plasmid Span Events (EE)" in content
    assert "Human Intron Clips (EI)" in content


class TestPlasmidFlagAndCohortDepth:
    @pytest.mark.unit
    def test_contaminated_verdict_is_positive_with_zero_f_plasmid(self) -> None:
        """Contaminated verdict -> 'positive' even when f_plasmid == 0.0 (single-exon construct without target BED)."""
        assert (
            compute_plasmid_flag(
                verdict="Sample is contaminated with plasmid DNA",
                f_plasmid=0.0,
            )
            == "positive"
        )
        assert compute_plasmid_flag(verdict="Contaminated", f_plasmid=0.0) == "positive"

    @pytest.mark.unit
    def test_unclear_verdict_is_weak_at_low_ratio(self) -> None:
        """Unclear verdict -> 'weak' even at ratio = 0.62 where 0.5 * ratio / (1 + ratio) < 0.20."""
        summary_df = pd.DataFrame(
            {
                "Sample": ["APA25-N", "APA25-N"],
                "Plasmid": ["CACNA1H", "CACNA1H"],
                "Category": ["Ratio", "Verdict"],
                "Value": ["0.62", "Sample contamination status is unclear"],
            }
        )
        verdict_df = summary_df[summary_df["Category"] == "Verdict"]
        flagged = _add_plasmid_flag(verdict_df, summary_df)
        assert list(flagged["plasmid_flag"]) == ["weak"]

    @pytest.mark.unit
    def test_f_plasmid_thresholds(self) -> None:
        """f_plasmid >= 0.50 -> 'positive'; 0.20 <= f_plasmid < 0.50 -> 'weak'."""
        assert (
            compute_plasmid_flag("Sample is not contaminated with plasmid DNA", f_plasmid=0.50)
            == "positive"
        )
        assert (
            compute_plasmid_flag("Sample is not contaminated with plasmid DNA", f_plasmid=0.75)
            == "positive"
        )
        assert (
            compute_plasmid_flag("Sample is not contaminated with plasmid DNA", f_plasmid=0.20)
            == "weak"
        )
        assert (
            compute_plasmid_flag("Sample is not contaminated with plasmid DNA", f_plasmid=0.49)
            == "weak"
        )
        assert (
            compute_plasmid_flag("Sample is not contaminated with plasmid DNA", f_plasmid=0.19)
            == "negative"
        )

    @pytest.mark.unit
    def test_structural_junction_and_bridge_reads_trigger_weak(self) -> None:
        """SpliceJunctionReads > 0, VectorInsertJunctionReads > 0, VectorInsertBridgePairs > 0 -> 'weak' on NotContaminated verdict."""
        summary_df = pd.DataFrame(
            {
                "Sample": ["s1", "s1", "s2", "s2", "s3", "s3"],
                "Plasmid": ["p1", "p1", "p1", "p1", "p1", "p1"],
                "Category": [
                    "Verdict",
                    "SpliceJunctionReads",
                    "Verdict",
                    "VectorInsertJunctionReads",
                    "Verdict",
                    "VectorInsertBridgePairs",
                ],
                "Value": ["NotContaminated", "2", "NotContaminated", "1", "NotContaminated", "3"],
            }
        )
        verdict_df = summary_df[summary_df["Category"] == "Verdict"]
        flagged = _add_plasmid_flag(verdict_df, summary_df)
        assert list(flagged["plasmid_flag"]) == ["weak", "weak", "weak"]

    @pytest.mark.unit
    def test_backbone_only_negative_control(self) -> None:
        """Backbone_Only > 0 and MeanDepthBackbone > 0 with NotContaminated verdict and 0 junction/bridge reads -> 'negative'."""
        summary_df = pd.DataFrame(
            {
                "Sample": ["s1"] * 7,
                "Plasmid": ["pWPI-KCNJ5-WT"] * 7,
                "Category": [
                    "Verdict",
                    "Backbone_Only",
                    "MeanDepthBackbone",
                    "BreadthBackbone",
                    "SpliceJunctionReads",
                    "VectorInsertJunctionReads",
                    "VectorInsertBridgePairs",
                ],
                "Value": ["NotContaminated", "420", "85.5", "0.95", "0", "0", "0"],
            }
        )
        verdict_df = summary_df[summary_df["Category"] == "Verdict"]
        flagged = _add_plasmid_flag(verdict_df, summary_df)
        assert list(flagged["plasmid_flag"]) == ["negative"]

    @pytest.mark.unit
    def test_variant_cooccurrence_and_cdna_artefacts(self) -> None:
        """cooccurrence=True -> 'positive'; cdna_artefacts > 0 -> 'weak'."""
        assert compute_plasmid_flag("NotContaminated", cooccurrence=True) == "positive"
        assert compute_plasmid_flag("NotContaminated", cdna_artefacts=1.0) == "weak"

        verdict_df = pd.DataFrame(
            {
                "Sample": ["s_cooc", "s_cdna", "s_clean"],
                "Plasmid": ["p1", "p1", "p1"],
                "Category": ["Verdict", "Verdict", "Verdict"],
                "Value": ["NotContaminated", "NotContaminated", "NotContaminated"],
            }
        )
        flagged = _add_plasmid_flag(
            verdict_df,
            verdict_df,
            cooccurrence_pairs={("s_cooc", "p1")},
            cdna_artefact_counts={("s_cdna", "p1"): 2},
        )
        assert list(flagged["plasmid_flag"]) == ["positive", "weak", "negative"]

    @pytest.mark.unit
    def test_verdict_table_tsv_includes_plasmid_flag_column(self, tmp_path: Path) -> None:
        """Verify verdict_table.tsv written by save_tables_as_tsv_and_excel includes plasmid_flag and preserves columns."""
        combined_df = pd.DataFrame(
            {
                "Sample": ["s1", "s2"],
                "Plasmid": ["p1", "p1"],
                "Category": ["Plasmid", "Plasmid"],
                "Value": [100, 0],
            }
        )
        verdict_df = pd.DataFrame(
            {
                "Sample": ["s1", "s2"],
                "Plasmid": ["p1", "p1"],
                "Category": ["Verdict", "Verdict"],
                "Value": [
                    "Sample is contaminated with plasmid DNA",
                    "Sample is not contaminated with plasmid DNA",
                ],
            }
        )
        ratio_df = pd.DataFrame(
            {
                "Sample": ["s1", "s2"],
                "Plasmid": ["p1", "p1"],
                "Category": ["Ratio", "Ratio"],
                "Value": [5.0, 0.0],
            }
        )
        p_values_df = pd.DataFrame(
            {
                "Sample": ["s1", "s2"],
                "Plasmid": ["p1", "p1"],
                "p_value": [0.01, 0.9],
                "p_value_corrected": [0.02, 0.9],
            }
        )

        save_tables_as_tsv_and_excel(combined_df, verdict_df, ratio_df, p_values_df, str(tmp_path))

        verdict_tsv = pd.read_csv(tmp_path / "data" / "verdict_table.tsv", sep="\t")
        assert list(verdict_tsv.columns) == [
            "Sample",
            "Plasmid",
            "Category",
            "Value",
            "plasmid_flag",
        ]
        by_sample = dict(zip(verdict_tsv["Sample"], verdict_tsv["plasmid_flag"], strict=True))
        assert by_sample == {"s1": "positive", "s2": "negative"}

    @pytest.mark.unit
    def test_reestimate_cohort_depth_fractions_multi_sample(self) -> None:
        """Verify _reestimate_cohort_depth_fractions updates r_cohort, f_excess, and f_plasmid across a cohort."""
        # 3 clean samples (Ratio < 1.5) with dp_locus / dp_sample_median = [0.70, 0.80, 0.90] -> median r_cohort = 0.80
        # 1 contaminated single-exon sample (s4, Ratio = 3.0) with dp_sample_median = 100, dp_locus = 160 (2x expected 80)
        # -> f_excess = 1 - (0.80 * 100) / 160 = 0.50, f_plasmid = 0.50
        rows = []
        sample_specs = [
            ("s1", 100.0, 70.0, 0.01, 0.60, 0.0),
            ("s2", 100.0, 80.0, 0.02, 0.60, 0.0),
            ("s3", 100.0, 90.0, 0.05, 0.60, 0.0),
            ("s4", 100.0, 160.0, 3.00, 0.60, 0.0),
        ]
        for s, dp_med, dp_loc, ratio, r_init, f_junc in sample_specs:
            rows.extend(
                [
                    {
                        "Sample": s,
                        "Plasmid": "KCNJ5_construct",
                        "Category": "dp_sample_median",
                        "Value": dp_med,
                    },
                    {
                        "Sample": s,
                        "Plasmid": "KCNJ5_construct",
                        "Category": "dp_locus",
                        "Value": dp_loc,
                    },
                    {
                        "Sample": s,
                        "Plasmid": "KCNJ5_construct",
                        "Category": "Ratio",
                        "Value": ratio,
                    },
                    {
                        "Sample": s,
                        "Plasmid": "KCNJ5_construct",
                        "Category": "r_cohort",
                        "Value": r_init,
                    },
                    {
                        "Sample": s,
                        "Plasmid": "KCNJ5_construct",
                        "Category": "f_excess",
                        "Value": 0.0,
                    },
                    {
                        "Sample": s,
                        "Plasmid": "KCNJ5_construct",
                        "Category": "f_junction",
                        "Value": f_junc,
                    },
                    {
                        "Sample": s,
                        "Plasmid": "KCNJ5_construct",
                        "Category": "f_plasmid",
                        "Value": 0.0,
                    },
                ]
            )
        summary_df = pd.DataFrame(rows)
        updated = _reestimate_cohort_depth_fractions(summary_df)

        def _val(sample: str, cat: str) -> float:
            sub = updated[(updated["Sample"] == sample) & (updated["Category"] == cat)]
            return float(sub["Value"].iloc[0])

        for s in ("s1", "s2", "s3", "s4"):
            assert _val(s, "r_cohort") == pytest.approx(0.80)

        # Clean samples s1 (70 <= 80) and s2 (80 <= 80) have f_excess == 0.0
        assert _val("s1", "f_excess") == pytest.approx(0.0)
        assert _val("s2", "f_excess") == pytest.approx(0.0)
        # Spiked sample s4 has 2x expected locus depth -> f_excess == 0.50 and f_plasmid == 0.50
        assert _val("s4", "f_excess") == pytest.approx(0.50)
        assert _val("s4", "f_plasmid") == pytest.approx(0.50)

    @pytest.mark.unit
    def test_main_with_vcf_file_writes_flagged_variants(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Verify main(..., vcf_file=...) writes plasmid_flagged_variants.tsv and feeds variant evidence into plasmid_flag."""
        from plasmicheck.scripts.generate_summary_reports import main

        for s, p, cvar in (
            ("s1-T", "p_g151r", "128911724:G>A"),
            ("s1-T", "p_l168r", "128911776:T>G"),
            ("s2-T", "p_cacna1h", "NA"),
        ):
            sdir = tmp_path / "input" / s / p / "data"
            sdir.mkdir(parents=True)
            (sdir / "comparison_result.summary.tsv").write_text(
                "Category\tCount\n"
                "Plasmid\t0\n"
                "Human\t500\n"
                "Tied\t10\n"
                "Ratio\t0.0\n"
                "Threshold\t0.8\n"
                "Verdict\tSample is not contaminated with plasmid DNA\n"
                "SpliceJunctionReads\t1\n"
                f"Chromosome\t{'chr11' if 's1' in s else 'chr16'}\n"
                f"GenomicStart\t{128910000 if 's1' in s else 1200000}\n"
                f"GenomicEnd\t{128920000 if 's1' in s else 1220000}\n"
                f"GenomicIntrons\t{'NA' if 's1' in s else '1207118-1207210'}\n"
                f"ConstructVariants\t{cvar}\n"
                "f_plasmid\t0.10\n"
            )

        vcf_tsv_path = tmp_path / "cohort_variants.tsv"
        vcf_tsv_path.write_text(
            "sample\tchr\tpos\tref\talt\tad\tdp\tvaf\n"
            "s1_To\tchr11\t128911724\tG\tA\t15\t100\t0.15\n"
            "s1_To\tchr11\t128911776\tT\tG\t49\t100\t0.49\n"
            "s2_To\tchr16\t1207114\tACGTACGTACGTACG\tA\t8\t100\t0.08\n"
        )
        out_dir = tmp_path / "summary_out"
        out_dir.mkdir()

        main(str(tmp_path / "input"), str(out_dir), vcf_file=str(vcf_tsv_path))

        flagged_variants_tsv = out_dir / "data" / "plasmid_flagged_variants.tsv"
        assert flagged_variants_tsv.exists()
        verdict_tsv = pd.read_csv(out_dir / "data" / "verdict_table.tsv", sep="\t")
        by_pair = {
            (str(r["Sample"]), str(r["Plasmid"])): str(r["plasmid_flag"])
            for _, r in verdict_tsv.iterrows()
        }
        assert by_pair[("s1-T", "p_g151r")] == "positive"
        assert by_pair[("s1-T", "p_l168r")] == "positive"
        assert by_pair[("s2-T", "p_cacna1h")] == "weak"
