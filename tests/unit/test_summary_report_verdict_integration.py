"""Integration-style unit tests for verdict bar chart wiring in summary reports."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from jinja2 import Environment, FileSystemLoader

from plasmicheck.resources import get_resource_path


@pytest.mark.unit
def test_template_renders_verdict_overview_interactive(tmp_path: Path) -> None:
    """Test that summary template renders Verdict Overview section with interactive charts."""
    from plasmicheck.config import get_config

    cfg = get_config()
    template_dir = cfg["paths"]["template_dir"]

    # Build Jinja2 environment
    env = Environment(loader=FileSystemLoader(str(get_resource_path(template_dir))))
    template = env.get_template("summary_template.html")

    # Minimal template variables
    template_vars = {
        "interactive": True,
        "verdict_by_plasmid_content": "<div>mock plotly chart plasmid</div>",
        "verdict_by_sample_content": "<div>mock plotly chart sample</div>",
        "combined_df": "<table></table>",
        "verdict_df": "<table></table>",
        "ratio_df": "<table></table>",
        "heatmap_content": "<div></div>",
        "boxplot_content": "<div></div>",
        "p_value_table": "<table></table>",
        "logo_base64": "data:image/png;base64,abc",
        "threshold": 0.8,
        "unclear_range": {"lower_bound": 0.6, "upper_bound": 1.0},
        "version": "0.34.0",
        "run_date": "2026-02-16",
        "command_line": "test",
        "plotly_mode": "embedded",
        "plotly_version": "5.0.0",
        "plotly_js_path": "",
        "plotly_js_inline": "<script></script>",
    }

    # Render template
    html = template.render(**template_vars)

    # Assertions
    assert "Verdict Overview" in html
    assert "Verdicts by Plasmid" in html
    assert "Verdicts by Sample" in html
    assert "mock plotly chart plasmid" in html
    assert "mock plotly chart sample" in html


@pytest.mark.unit
def test_template_renders_verdict_overview_static(tmp_path: Path) -> None:
    """Test that summary template renders Verdict Overview section with static PNG charts."""
    from plasmicheck.config import get_config

    cfg = get_config()
    template_dir = cfg["paths"]["template_dir"]

    # Build Jinja2 environment
    env = Environment(loader=FileSystemLoader(str(get_resource_path(template_dir))))
    template = env.get_template("summary_template.html")

    # Minimal template variables for static report
    template_vars = {
        "interactive": False,
        "verdict_by_plasmid_content": "data:image/png;base64,plasmid123",
        "verdict_by_sample_content": "data:image/png;base64,sample456",
        "combined_df": "<table></table>",
        "verdict_df": "<table></table>",
        "ratio_df": "<table></table>",
        "heatmap_content": "data:image/png;base64,heatmap",
        "boxplot_content": "data:image/png;base64,boxplot",
        "p_value_table": "<table></table>",
        "logo_base64": "data:image/png;base64,logo",
        "threshold": 0.8,
        "unclear_range": {"lower_bound": 0.6, "upper_bound": 1.0},
        "version": "0.34.0",
        "run_date": "2026-02-16",
        "command_line": "test",
        "plotly_mode": "embedded",
        "plotly_version": "5.0.0",
        "plotly_js_path": "",
        "plotly_js_inline": "<script></script>",
    }

    # Render template
    html = template.render(**template_vars)

    # Assertions
    assert "Verdict Overview" in html
    assert 'src="data:image/png;base64,plasmid123"' in html
    assert 'src="data:image/png;base64,sample456"' in html


@pytest.mark.unit
def test_template_renders_without_verdict_data() -> None:
    """Test that template renders without error when verdict data is missing (backward compat)."""
    from plasmicheck.config import get_config

    cfg = get_config()
    template_dir = cfg["paths"]["template_dir"]

    # Build Jinja2 environment
    env = Environment(loader=FileSystemLoader(str(get_resource_path(template_dir))))
    template = env.get_template("summary_template.html")

    # Minimal template variables WITHOUT verdict fields
    template_vars = {
        "interactive": True,
        "combined_df": "<table></table>",
        "verdict_df": "<table></table>",
        "ratio_df": "<table></table>",
        "heatmap_content": "<div></div>",
        "boxplot_content": "<div></div>",
        "p_value_table": "<table></table>",
        "logo_base64": "data:image/png;base64,abc",
        "threshold": 0.8,
        "unclear_range": {"lower_bound": 0.6, "upper_bound": 1.0},
        "version": "0.34.0",
        "run_date": "2026-02-16",
        "command_line": "test",
        "plotly_mode": "embedded",
        "plotly_version": "5.0.0",
        "plotly_js_path": "",
        "plotly_js_inline": "<script></script>",
    }

    # Render template - should not raise
    html = template.render(**template_vars)

    # Verdict Overview section should still be present
    assert "Verdict Overview" in html
    # But the plot divs will be empty (Jinja2 undefined defaults to empty string)


@pytest.mark.unit
def test_generate_report_passes_verdict_data_to_template(tmp_path: Path) -> None:
    """Test that generate_report() passes verdict chart data to template.render()."""
    from plasmicheck.scripts.generate_summary_reports import generate_report

    # Create minimal DataFrames matching the expected structure
    combined_df = pd.DataFrame(
        {"Sample": ["S1"], "Plasmid": ["P1"], "Category": ["Plasmid"], "Value": [100]}
    )
    verdict_df = pd.DataFrame(
        {"Sample": ["S1"], "Plasmid": ["P1"], "Category": ["Verdict"], "Value": ["Clean"]}
    )
    ratio_df = pd.DataFrame(
        {"Sample": ["S1"], "Plasmid": ["P1"], "Category": ["Ratio"], "Value": [0.1]}
    )

    # Create temporary plot files
    verdict_plasmid_html = tmp_path / "verdict_plasmid.html"
    verdict_plasmid_html.write_text("<div>verdict by plasmid chart</div>")

    verdict_sample_html = tmp_path / "verdict_sample.html"
    verdict_sample_html.write_text("<div>verdict by sample chart</div>")

    heatmap_html = tmp_path / "heatmap.html"
    heatmap_html.write_text("<div>heatmap</div>")

    boxplot_html = tmp_path / "boxplot.html"
    boxplot_html.write_text("<div>boxplot</div>")

    p_value_file = tmp_path / "p_value.tsv"
    p_value_file.write_text("Sample\tPlasmid\tp_value\tp_value_corrected\nS1\tP1\t0.05\t0.05\n")

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Mock the jinja2.Environment class
    with patch("jinja2.Environment") as mock_env_class:
        mock_env = MagicMock()
        mock_template = MagicMock()
        mock_env.get_template.return_value = mock_template
        mock_env_class.return_value = mock_env

        # Mock template.render to return empty string
        mock_template.render.return_value = "<html></html>"

        # Call generate_report with verdict data
        generate_report(
            combined_df=combined_df,
            verdict_df=verdict_df,
            ratio_df=ratio_df,
            heatmap_filename_interactive=str(heatmap_html),
            heatmap_filename_png=None,
            boxplot_filename_interactive=str(boxplot_html),
            boxplot_filename_png=None,
            p_value_table_filename=str(p_value_file),
            output_folder=str(output_dir),
            threshold=0.8,
            unclear_range={"lower_bound": 0.6, "upper_bound": 1.0},
            command_line="test",
            static_report=False,
            verdict_by_plasmid_interactive=str(verdict_plasmid_html),
            verdict_by_sample_interactive=str(verdict_sample_html),
        )

        # Assert template.render was called with verdict data
        assert mock_template.render.called
        call_kwargs = mock_template.render.call_args[1]
        assert "verdict_by_plasmid_content" in call_kwargs
        assert "verdict by plasmid chart" in call_kwargs["verdict_by_plasmid_content"]
        assert "verdict_by_sample_content" in call_kwargs
        assert "verdict by sample chart" in call_kwargs["verdict_by_sample_content"]


@pytest.mark.unit
def test_main_calls_verdict_plot_functions(tmp_path: Path) -> None:
    """Test that main() calls plot_verdict_bar_by_plasmid and plot_verdict_bar_by_sample."""
    from plasmicheck.scripts.generate_summary_reports import main

    # Create minimal input directory structure
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    sample_dir = input_dir / "S1" / "P1" / "data"
    sample_dir.mkdir(parents=True)

    # Create minimal TSV files
    reads_file = sample_dir / "comparison.reads_assignment.tsv"
    reads_file.write_text("ReadID\tAssignment\tPlasmidScore\tHumanScore\nR1\tPlasmid\t100\t50\n")

    summary_file = sample_dir / "comparison.summary.tsv"
    summary_file.write_text(
        "Category\tCount\nPlasmid\t100\nHuman\t50\nTied\t10\nVerdict\tSample is not contaminated with plasmid DNA\nRatio\t0.1\n"
    )

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Mock all plot functions and generate_report
    with (
        patch("plasmicheck.scripts.generate_summary_reports.plot_boxplot") as mock_boxplot,
        patch("plasmicheck.scripts.generate_summary_reports.plot_heatmap") as mock_heatmap,
        patch(
            "plasmicheck.scripts.generate_summary_reports.plot_metric_heatmap"
        ) as mock_metric_heatmap,
        patch(
            "plasmicheck.scripts.generate_summary_reports.plot_coverage_outside_insert_boxplot"
        ) as mock_cov_boxplot,
        patch(
            "plasmicheck.scripts.generate_summary_reports.create_resistance_gene_table"
        ) as mock_res_genes,
        patch(
            "plasmicheck.scripts.generate_summary_reports.plot_verdict_bar_by_plasmid"
        ) as mock_verdict_plasmid,
        patch(
            "plasmicheck.scripts.generate_summary_reports.plot_verdict_bar_by_sample"
        ) as mock_verdict_sample,
        patch("plasmicheck.scripts.generate_summary_reports.generate_report") as mock_gen_report,
        patch("plasmicheck.scripts.generate_summary_reports.save_tables_as_tsv_and_excel"),
    ):
        # Configure mocks to return empty paths
        mock_boxplot.return_value = ("", None)
        mock_heatmap.return_value = ("", None)
        mock_metric_heatmap.return_value = ("", None)
        mock_cov_boxplot.return_value = ("", None)
        mock_res_genes.return_value = pd.DataFrame()
        mock_verdict_plasmid.return_value = ("/path/to/verdict_plasmid.html", None)
        mock_verdict_sample.return_value = ("/path/to/verdict_sample.html", None)

        # Call main
        main(
            input_dir=str(input_dir),
            output_dir=str(output_dir),
            threshold=0.8,
            unclear_range={"lower_bound": 0.6, "upper_bound": 1.0},
            static_report=False,
            plot_backend="plotly",
        )

        # Assert verdict plot functions were called
        assert mock_verdict_plasmid.called
        assert mock_verdict_sample.called

        # Assert generate_report was called with verdict data
        assert mock_gen_report.called
        call_kwargs = mock_gen_report.call_args[1]
        assert "verdict_by_plasmid_interactive" in call_kwargs
        assert call_kwargs["verdict_by_plasmid_interactive"] == "/path/to/verdict_plasmid.html"
        assert "verdict_by_sample_interactive" in call_kwargs
        assert call_kwargs["verdict_by_sample_interactive"] == "/path/to/verdict_sample.html"
