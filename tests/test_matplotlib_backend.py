"""Tests for matplotlib plotting backend."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd
import pytest


@pytest.mark.unit
def test_generate_boxplot_matplotlib_creates_png(tmp_path: Path) -> None:
    """Test that generate_boxplot_matplotlib creates a PNG file."""
    from plasmicheck.scripts.plotting.matplotlib_backend import generate_boxplot_matplotlib

    # Create test DataFrame
    reads_df = pd.DataFrame(
        {
            "ReadID": ["read1", "read2", "read3", "read4", "read5", "read6"],
            "AssignedTo": ["Plasmid", "Plasmid", "Human", "Human", "Tied", "Tied"],
            "PlasmidScore": [50, 60, 20, 25, 40, 45],
            "HumanScore": [10, 15, 55, 60, 38, 42],
        }
    )

    output_path = tmp_path / "box.png"
    title = "Test Boxplot (Total Reads: 6)"

    generate_boxplot_matplotlib(reads_df, output_path, title, width=900, height=600)

    assert output_path.exists()
    assert output_path.stat().st_size > 0


@pytest.mark.unit
def test_generate_scatter_matplotlib_creates_png(tmp_path: Path) -> None:
    """Test that generate_scatter_matplotlib creates a PNG file."""
    from plasmicheck.scripts.plotting.matplotlib_backend import generate_scatter_matplotlib

    # Create test DataFrame
    reads_df = pd.DataFrame(
        {
            "ReadID": ["read1", "read2", "read3", "read4", "read5", "read6"],
            "AssignedTo": ["Plasmid", "Plasmid", "Human", "Human", "Tied", "Tied"],
            "PlasmidScore": [50, 60, 20, 25, 40, 45],
            "HumanScore": [10, 15, 55, 60, 38, 42],
        }
    )

    output_path = tmp_path / "scatter.png"
    title = "Test Scatter (Total Reads: 6)"

    generate_scatter_matplotlib(reads_df, output_path, title, width=900, height=600)

    assert output_path.exists()
    assert output_path.stat().st_size > 0


@pytest.mark.unit
def test_generate_heatmap_matplotlib_creates_png(tmp_path: Path) -> None:
    """Test that generate_heatmap_matplotlib creates a PNG file."""
    from plasmicheck.scripts.plotting.matplotlib_backend import generate_heatmap_matplotlib

    # Create pivoted DataFrame (Sample x Plasmid)
    ratio_data = pd.DataFrame(
        {
            "Plasmid1": [0.5, 0.8, 1.2],
            "Plasmid2": [0.3, 1.5, 0.9],
            "Plasmid3": [0.1, 0.4, 2.0],
        },
        index=["Sample1", "Sample2", "Sample3"],
    )

    output_path = tmp_path / "heatmap.png"
    threshold = 0.8
    plot_config = {}

    generate_heatmap_matplotlib(ratio_data, output_path, threshold, plot_config)

    assert output_path.exists()
    assert output_path.stat().st_size > 0


@pytest.mark.unit
def test_generate_summary_boxplot_matplotlib_creates_png(tmp_path: Path) -> None:
    """Test that generate_summary_boxplot_matplotlib creates a PNG file."""
    from plasmicheck.scripts.plotting.matplotlib_backend import (
        generate_summary_boxplot_matplotlib,
    )

    # Create pivoted DataFrame (Sample x Plasmid)
    boxplot_data = pd.DataFrame(
        {
            "Plasmid1": [0.5, 0.8, 1.2],
            "Plasmid2": [0.3, 1.5, 0.9],
            "Plasmid3": [0.1, 0.4, 2.0],
        },
        index=["Sample1", "Sample2", "Sample3"],
    )

    output_path = tmp_path / "summary_boxplot.png"

    generate_summary_boxplot_matplotlib(boxplot_data, output_path, log_offset=0.001)

    assert output_path.exists()
    assert output_path.stat().st_size > 0


@pytest.mark.unit
def test_matplotlib_colors_match_plotly() -> None:
    """Test that ASSIGNMENT_COLORS has the expected keys and hex values."""
    from plasmicheck.scripts.plotting.colors import ASSIGNMENT_COLORS

    assert "Plasmid" in ASSIGNMENT_COLORS
    assert "Human" in ASSIGNMENT_COLORS
    assert "Tied" in ASSIGNMENT_COLORS

    # Verify they are hex strings
    assert ASSIGNMENT_COLORS["Plasmid"].startswith("#")
    assert ASSIGNMENT_COLORS["Human"].startswith("#")
    assert ASSIGNMENT_COLORS["Tied"].startswith("#")

    # Verify they match expected Plotly colors
    assert ASSIGNMENT_COLORS["Plasmid"] == "#636EFA"
    assert ASSIGNMENT_COLORS["Human"] == "#EF553B"
    assert ASSIGNMENT_COLORS["Tied"] == "#00CC96"


@pytest.mark.unit
def test_empty_dataframe_handled_gracefully(tmp_path: Path) -> None:
    """Test that matplotlib functions handle empty DataFrames without crashing."""
    from plasmicheck.scripts.plotting.matplotlib_backend import (
        generate_boxplot_matplotlib,
        generate_heatmap_matplotlib,
        generate_scatter_matplotlib,
        generate_summary_boxplot_matplotlib,
    )

    empty_df = pd.DataFrame()

    # Test boxplot
    output_box = tmp_path / "empty_box.png"
    generate_boxplot_matplotlib(empty_df, output_box, "Empty Boxplot")
    assert output_box.exists()

    # Test scatter
    output_scatter = tmp_path / "empty_scatter.png"
    generate_scatter_matplotlib(empty_df, output_scatter, "Empty Scatter")
    assert output_scatter.exists()

    # Test heatmap
    output_heatmap = tmp_path / "empty_heatmap.png"
    generate_heatmap_matplotlib(empty_df, output_heatmap, 0.8, {})
    assert output_heatmap.exists()

    # Test summary boxplot
    output_summary = tmp_path / "empty_summary.png"
    generate_summary_boxplot_matplotlib(empty_df, output_summary)
    assert output_summary.exists()


@pytest.mark.unit
def test_plot_backend_cli_flag_exists() -> None:
    """Test that --plot-backend flag exists in CLI help output."""
    # Test pipeline command
    result_pipeline = subprocess.run(
        ["plasmicheck", "pipeline", "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--plot-backend" in result_pipeline.stdout

    # Test report command
    result_report = subprocess.run(
        ["plasmicheck", "report", "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--plot-backend" in result_report.stdout

    # Test summary_reports command
    result_summary = subprocess.run(
        ["plasmicheck", "summary_reports", "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--plot-backend" in result_summary.stdout


@pytest.mark.unit
def test_plot_backend_default_is_plotly() -> None:
    """Test that the default value for --plot-backend is 'plotly'."""
    import inspect

    from plasmicheck.scripts.generate_report import main

    sig = inspect.signature(main)
    plot_backend_param = sig.parameters.get("plot_backend")

    assert plot_backend_param is not None
    assert plot_backend_param.default == "plotly"


@pytest.mark.unit
def test_matplotlib_backend_no_kaleido_import(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Test that matplotlib backend doesn't import kaleido."""
    # Block kaleido import by making it unavailable
    import sys

    original_kaleido = sys.modules.get("kaleido")

    # Remove kaleido from sys.modules if it exists
    if "kaleido" in sys.modules:
        del sys.modules["kaleido"]

    # Create a mock module that raises ImportError
    class MockKaleido:
        def __getattr__(self, name: str) -> None:
            raise ImportError("kaleido is not available")

    monkeypatch.setitem(sys.modules, "kaleido", MockKaleido())  # type: ignore[arg-type]

    try:
        from plasmicheck.scripts.generate_report import generate_plots

        # Create test DataFrame
        reads_df = pd.DataFrame(
            {
                "ReadID": ["read1", "read2"],
                "AssignedTo": ["Plasmid", "Human"],
                "PlasmidScore": [50, 20],
                "HumanScore": [10, 55],
            }
        )

        # Call generate_plots with matplotlib backend
        # Should succeed without importing kaleido
        output_folder = str(tmp_path)
        os.makedirs(output_folder, exist_ok=True)

        (
            _boxplot_interactive,
            boxplot_png,
            _scatter_interactive,
            scatter_png,
        ) = generate_plots(reads_df, output_folder, static_report=True, plot_backend="matplotlib")

        # Verify PNG files were created
        assert boxplot_png is not None
        assert os.path.exists(boxplot_png)
        assert scatter_png is not None
        assert os.path.exists(scatter_png)

    finally:
        # Restore original kaleido module if it existed
        if original_kaleido is not None:
            sys.modules["kaleido"] = original_kaleido
        elif "kaleido" in sys.modules:
            del sys.modules["kaleido"]
