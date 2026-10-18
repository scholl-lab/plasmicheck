from __future__ import annotations

import logging
import os
import sys
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any

from plasmicheck.config import get_config
from plasmicheck.resources import get_resource_path
from plasmicheck.version import __version__ as VERSION

from .utils import (
    categorize_value as categorize_value,
)
from .utils import (
    encode_image_to_base64,
    ensure_plotly_assets,
    truncate_comma_separated,
)
from .utils import (
    verdict_bounds as verdict_bounds,
)

if TYPE_CHECKING:
    import pandas as pd

_cfg = get_config()
DEFAULT_THRESHOLD: float = _cfg["default_threshold"]
UNCLEAR_RANGE: dict[str, float] = _cfg["unclear_range"]
PLOT_CONFIG: dict[str, Any] = _cfg["plot_summary"]
TEMPLATE_DIR: str = _cfg["paths"]["template_dir"]
LOGO_PATH: str = _cfg["paths"]["logo_path"]
TABLE_SORTING: dict[str, Any] = _cfg["table_sorting"]
ROUND_DECIMALS: int = PLOT_CONFIG.get("round_decimals", 3)
LOG_OFFSET: float = PLOT_CONFIG.get("log_offset", 1e-9)
MARKER_STYLE: dict[str, Any] = PLOT_CONFIG.get("marker_style", {"size": 7, "opacity": 0.75})
DEPTH_HEATMAP_CAP: float = PLOT_CONFIG.get("depth_heatmap_cap", 100.0)
_PLOTLY_CONFIG: dict[str, Any] = {
    "showSendToCloud": False,
    "displayModeBar": "hover",
    "displaylogo": False,
    "responsive": True,
}
_SUMMARY_FONT = '-apple-system, BlinkMacSystemFont, "Inter", "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif'
_SUMMARY_PLOTLY_LAYOUT: dict[str, Any] = {
    "template": "plotly_white",
    "autosize": True,
    "font": {"family": _SUMMARY_FONT, "size": 12, "color": "#0f172a"},
    "paper_bgcolor": "white",
    "plot_bgcolor": "white",
    "margin": {"l": 64, "r": 24, "t": 28, "b": 96, "pad": 4},
    "hoverlabel": {
        "bgcolor": "#0f172a",
        "bordercolor": "#334155",
        "font": {"size": 12, "family": _SUMMARY_FONT, "color": "#f8fafc"},
    },
}

# Setup logging for the libraries used in this script
logging.getLogger("matplotlib").setLevel(logging.ERROR)
logging.getLogger("seaborn").setLevel(logging.ERROR)
logging.getLogger("jinja2").setLevel(logging.ERROR)
logging.getLogger("fontTools").setLevel(logging.ERROR)


# Helper Functions
def find_tsv_files(input_dir: str, pattern: str) -> list[str]:
    """Recursively find all tsv files matching the pattern in the input directory."""
    tsv_files: list[str] = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(pattern):
                tsv_files.append(os.path.join(root, file))
    return tsv_files


def _extract_sample_plasmid(file_path: str) -> tuple[str, str]:
    """Extract sample and plasmid names from a comparison output file path.

    Handles both directory layouts:
    - Old (flat): output/sample/plasmid/file.tsv
    - New (subfolder): output/sample/plasmid/data/file.tsv

    Strategy: if the immediate parent directory is named 'data',
    go one level further up for plasmid and sample names.

    Args:
        file_path: Path to the comparison file

    Returns:
        Tuple of (sample_name, plasmid_name)
    """
    parent = os.path.dirname(file_path)
    parent_name = os.path.basename(parent)

    if parent_name == "data":
        # New subfolder layout: .../sample/plasmid/data/file.tsv
        plasmid_dir = os.path.dirname(parent)
        plasmid_name = os.path.basename(plasmid_dir)
        sample_name = os.path.basename(os.path.dirname(plasmid_dir))
    else:
        # Old flat layout: .../sample/plasmid/file.tsv
        plasmid_name = parent_name
        sample_name = os.path.basename(os.path.dirname(parent))

    return sample_name, plasmid_name


def read_compare_outputs(input_dir: str, substring_to_remove: str | None = None) -> pd.DataFrame:
    import pandas as pd

    logging.info(f"Reading compare outputs from {input_dir}")
    summary_files: list[str] = find_tsv_files(input_dir, ".summary.tsv")

    if not summary_files:
        logging.error("No valid TSV files found in the specified directory structure.")
        raise ValueError("No valid TSV files found in the specified directory structure.")

    summary_df_list: list[pd.DataFrame] = []
    for file in summary_files:
        sample_name, plasmid_name = _extract_sample_plasmid(file)

        # Remove the substring from the sample name if provided
        if substring_to_remove:
            sample_name = sample_name.replace(substring_to_remove, "")

        df = pd.read_csv(
            file, sep="\t", low_memory=False
        )  # Disable low_memory mode to avoid DtypeWarning
        df["Sample"] = sample_name
        df["Plasmid"] = plasmid_name
        df = df.rename(columns={"Count": "Value"})  # Rename Count to Value
        summary_df_list.append(df)

    summary_df = pd.concat(summary_df_list, ignore_index=True)

    return summary_df


def apply_sorting(df: pd.DataFrame, sorting_config: dict[str, Any] | None) -> pd.DataFrame:
    """Apply sorting to a DataFrame based on the given sorting configuration."""
    if sorting_config:
        columns: list[str] = sorting_config.get("columns", [])
        ascending: list[bool] = sorting_config.get("ascending", [True] * len(columns))
        df = df.sort_values(by=columns, ascending=ascending)
    return df


def _sanitize_for_excel(df: pd.DataFrame) -> pd.DataFrame:
    """Convert boolean columns to strings to avoid locale-dependent rendering."""
    result = df.copy()
    bool_cols = result.select_dtypes(include=["bool"]).columns
    for col in bool_cols:
        result[col] = result[col].astype(str)
    return result


def _truncate_list_for_html(value: str, max_items: int = 5) -> str:
    """Truncate comma-separated list for HTML table display."""
    if not value or value in ("N/A", "NA", "nan"):
        return "N/A"
    return truncate_comma_separated(value, max_items)


def save_tables_as_tsv_and_excel(
    combined_df: pd.DataFrame,
    verdict_df: pd.DataFrame,
    ratio_df: pd.DataFrame,
    p_values_df: pd.DataFrame,
    output_dir: str,
    extended_metrics_df: pd.DataFrame | None = None,
    resistance_gene_df: pd.DataFrame | None = None,
    junction_df: pd.DataFrame | None = None,
) -> None:
    """Save the combined, verdict, ratio, and p-value dataframes as TSV and Excel files."""
    import pandas as pd

    data_dir: str = os.path.join(output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)

    # Apply sorting to each DataFrame
    combined_df = apply_sorting(combined_df, TABLE_SORTING.get("combined"))
    verdict_df = apply_sorting(verdict_df, TABLE_SORTING.get("verdict"))
    ratio_df = apply_sorting(ratio_df, TABLE_SORTING.get("ratio"))
    p_values_df = apply_sorting(p_values_df, TABLE_SORTING.get("p_value"))

    # Save TSV files
    combined_df.to_csv(os.path.join(data_dir, "combined_table.tsv"), sep="\t", index=False)
    verdict_df.to_csv(os.path.join(data_dir, "verdict_table.tsv"), sep="\t", index=False)
    ratio_df.to_csv(os.path.join(data_dir, "ratio_table.tsv"), sep="\t", index=False)
    p_values_df.to_csv(os.path.join(data_dir, "p_value_table.tsv"), sep="\t", index=False)

    # Save extended metrics TSV if provided
    if extended_metrics_df is not None and not extended_metrics_df.empty:
        extended_metrics_df.to_csv(
            os.path.join(data_dir, "extended_metrics_table.tsv"), sep="\t", index=False
        )

    # Save resistance genes TSV if provided
    if resistance_gene_df is not None and not resistance_gene_df.empty:
        resistance_gene_df.to_csv(
            os.path.join(data_dir, "resistance_genes_table.tsv"), sep="\t", index=False
        )

    # Save junction metadata TSV if provided
    if junction_df is not None and not junction_df.empty:
        junction_df.to_csv(
            os.path.join(data_dir, "junction_metadata_table.tsv"), sep="\t", index=False
        )

    # Save Excel file with multiple sheets
    excel_file: str = os.path.join(data_dir, "summary_report_tables.xlsx")
    with pd.ExcelWriter(excel_file, engine="xlsxwriter") as writer:
        _sanitize_for_excel(combined_df).to_excel(writer, sheet_name="Combined", index=False)
        _sanitize_for_excel(verdict_df).to_excel(writer, sheet_name="Verdicts", index=False)
        _sanitize_for_excel(ratio_df).to_excel(writer, sheet_name="Ratios", index=False)
        _sanitize_for_excel(p_values_df).to_excel(writer, sheet_name="P-Values", index=False)

        # Add extended metrics sheet if provided
        if extended_metrics_df is not None and not extended_metrics_df.empty:
            _sanitize_for_excel(extended_metrics_df).to_excel(
                writer, sheet_name="Extended Metrics", index=False
            )

        # Add resistance genes sheet if provided
        if resistance_gene_df is not None and not resistance_gene_df.empty:
            _sanitize_for_excel(resistance_gene_df).to_excel(
                writer, sheet_name="Resistance Genes", index=False
            )

        # Add junction metadata sheet if provided
        if junction_df is not None and not junction_df.empty:
            _sanitize_for_excel(junction_df).to_excel(
                writer, sheet_name="Junction Metadata", index=False
            )

    logging.info("TSV and Excel files have been saved.")


# Data Manipulation Functions
def calculate_variations(ratio_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate p-values and apply FDR correction for contamination ratios by plasmid."""
    import numpy as np
    import pandas as pd
    from scipy.stats import false_discovery_control
    from scipy.stats import t as student_t

    logging.info("Calculating variations for contamination ratios.")

    # Convert 'Value' to numeric, ensuring all entries are valid numbers.
    # Use .assign() to avoid in-place dtype conflict with StringDtype (pandas 2.2+).
    ratio_df = ratio_df.assign(Value=pd.to_numeric(ratio_df["Value"], errors="coerce"))
    boxplot_data: pd.DataFrame = ratio_df.pivot(index="Sample", columns="Plasmid", values="Value")

    # Drop NaNs
    boxplot_data = boxplot_data.dropna(how="all", axis=1)

    p_values_list: list[dict[str, Any]] = []
    for plasmid in boxplot_data.columns:
        # Ensure the data is numeric and drop NaN values
        data = pd.to_numeric(boxplot_data[plasmid], errors="coerce").dropna()
        n = len(data)
        if n > 1:
            x = data.to_numpy(dtype=np.float64)
            samples = data.index.tolist()
            if n < 3:
                for s in samples:
                    p_values_list.append({"Sample": s, "Plasmid": plasmid, "p_value": np.nan})
            else:
                s_tot = float(np.sum(x))
                q_tot = float(np.dot(x, x))
                m = n - 1
                mean_others = (s_tot - x) / m
                ss_others = np.maximum((q_tot - x * x) - m * (mean_others * mean_others), 0.0)
                var_others = ss_others / (m - 1)
                se = np.sqrt(var_others / m)
                with np.errstate(divide="ignore", invalid="ignore"):
                    t_stat = (mean_others - x) / se
                    p_vals = np.where(
                        var_others > 0.0,
                        student_t.sf(np.abs(t_stat), df=m - 1) * 2.0,
                        np.nan,
                    )
                for s, pv in zip(samples, p_vals, strict=True):
                    p_values_list.append(
                        {
                            "Sample": s,
                            "Plasmid": plasmid,
                            "p_value": float(pv) if not np.isnan(pv) else np.nan,
                        }
                    )

    p_values_df = pd.DataFrame(p_values_list, columns=["Sample", "Plasmid", "p_value"])
    p_values_df["p_value_corrected"] = np.nan
    valid = p_values_df["p_value"].notna()
    if valid.any():
        p_values_df.loc[valid, "p_value_corrected"] = false_discovery_control(
            p_values_df.loc[valid, "p_value"].to_numpy(dtype=float), method="bh"
        )

    return boxplot_data, p_values_df


calculate_Metrics = calculate_variations  # noqa: N816


# Plotting Functions
def plot_boxplot(
    boxplot_data: pd.DataFrame,
    output_dir: str,
    static_report: bool = False,
    plot_backend: str = "matplotlib",
    threshold: float = DEFAULT_THRESHOLD,
    unclear_range: dict[str, float] | None = None,
) -> tuple[str, str | None]:
    """Generate and save boxplots for contamination ratios by plasmid."""
    import numpy as np
    import pandas as pd
    import plotly.express as px

    logging.info("Generating boxplots for contamination ratios.")

    df_reset = boxplot_data.reset_index()
    if "Sample" not in df_reset.columns:
        df_reset = df_reset.rename(columns={df_reset.columns[0]: "Sample"})
    boxplot_df = pd.melt(df_reset, id_vars=["Sample"], var_name="Plasmid", value_name="Value")
    boxplot_df["log_value"] = boxplot_df["Value"].apply(lambda x: np.log10(x + LOG_OFFSET))

    fig = px.box(
        boxplot_df,
        x="Plasmid",
        y="log_value",
        points="all",
        hover_data={"Sample": True, "Value": ":.4f", "log_value": False},
    )

    if unclear_range is None:
        unclear_range = UNCLEAR_RANGE
    low, high = verdict_bounds(threshold, unclear_range)
    thresh_log = float(np.log10(high + LOG_OFFSET))
    low_log = float(np.log10(low + LOG_OFFSET))
    high_log = float(np.log10(high + LOG_OFFSET))

    if high > low:
        fig.add_hrect(
            y0=low_log,
            y1=high_log,
            fillcolor="rgba(245, 158, 11, 0.12)",
            line_width=0,
            annotation_text=f"Unclear ({low}-{high})",
            annotation_position="bottom right",
            annotation_font={"size": 10, "color": "#b45309"},
        )
    fig.add_hline(
        y=thresh_log,
        line_dash="dash",
        line_color="#dc2626",
        line_width=1.5,
        annotation_text=f"Contaminated (> {high})",
        annotation_position="top right",
        annotation_font={"size": 10, "color": "#dc2626"},
    )

    fig.update_layout(
        **_SUMMARY_PLOTLY_LAYOUT,
        xaxis_title="Plasmid",
        yaxis_title="log₁₀(Contamination Ratio)",
        height=420,
        xaxis={
            "tickangle": -30,
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
            "showgrid": False,
            "linecolor": "#cbd5e1",
        },
        yaxis={
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
            "gridcolor": "#f1f5f9",
            "linecolor": "#cbd5e1",
        },
    )

    fig.update_traces(
        marker={**MARKER_STYLE, "color": "#2563eb"},
        line={"color": "#1d4ed8", "width": 1.6},
        fillcolor="rgba(37, 99, 235, 0.18)",
    )

    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    boxplot_filename: str = os.path.join(plots_dir, "boxplot_contamination_ratios.html")
    fig.write_html(boxplot_filename, include_plotlyjs=False, full_html=False, config=_PLOTLY_CONFIG)

    boxplot_filename_png: str | None = None
    if static_report:
        boxplot_filename_png = os.path.join(plots_dir, "boxplot_contamination_ratios.png")
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_summary_boxplot_matplotlib

            generate_summary_boxplot_matplotlib(
                boxplot_data,
                boxplot_filename_png,
                LOG_OFFSET,
                "outliers",
                threshold=threshold,
                unclear_range=unclear_range,
            )
        else:
            fig.write_image(boxplot_filename_png, width=1000, height=420)

    logging.info("Boxplots generated.")

    return boxplot_filename, boxplot_filename_png


def plot_heatmap(
    ratio_df: pd.DataFrame,
    output_dir: str,
    threshold: float,
    unclear_range: dict[str, float],
    plot_config: dict[str, Any],
    static_report: bool = False,
    plot_backend: str = "matplotlib",
) -> tuple[str, str | None]:
    """Generate and save heatmap for contamination ratios by plasmid."""
    import pandas as pd
    import plotly.express as px

    logging.info("Generating heatmap for contamination ratios.")

    # Use .assign() to avoid dtype conflict with StringDtype (pandas 2.2+)
    ratio_df = ratio_df.assign(Value=pd.to_numeric(ratio_df["Value"], errors="coerce"))

    # Pivot the DataFrame for the heatmap
    ratio_data: pd.DataFrame = ratio_df.pivot(
        index="Sample", columns="Plasmid", values="Value"
    ).round(ROUND_DECIMALS)

    # Define a custom color scale using Plotly's built-in scaling
    max_val = float(ratio_data.max().max()) if not ratio_data.empty else 0.0
    if max_val > 0:
        lower_norm = min(unclear_range["lower_bound"] / max_val, 0.99)
        upper_norm = min(unclear_range["upper_bound"] / max_val, 0.995)
    else:
        lower_norm = 0.5
        upper_norm = 0.75
    color_scale: list[list[float | str]] = [
        [0.0, "#e0f2fe"],
        [lower_norm, "#fde68a"],
        [upper_norm, "#fca5a5"],
        [1.0, "#dc2626"],
    ]

    n_samples, n_plasmids = ratio_data.shape
    chart_height = max(420, min(640, 180 + n_samples * 28))
    show_annot = n_samples <= 20 and n_plasmids <= 15

    fig = px.imshow(
        ratio_data,
        color_continuous_scale=color_scale,
        height=chart_height,
        aspect="auto",
        text_auto=".2f" if show_annot else False,
    )
    fig.update_traces(xgap=2, ygap=2, textfont={"size": 11})

    fig.update_layout(
        **_SUMMARY_PLOTLY_LAYOUT,
        xaxis_title="Plasmid",
        yaxis_title=None,
        xaxis={
            "tickangle": -30,
            "tickfont": {"size": 11, "color": "#334155"},
            "tickmode": "linear",
            "automargin": True,
        },
        yaxis={
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
        },
    )

    fig.update_coloraxes(
        showscale=True,
        colorbar={
            "title": {"text": "Ratio", "font": {"size": 11}},
            "thickness": 10,
            "len": 0.82,
            "tickfont": {"size": 10},
            "outlinewidth": 0,
        },
    )

    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    heatmap_filename_interactive: str = os.path.join(
        plots_dir, plot_config["output_filename"].replace(".png", ".html")
    )
    fig.write_html(
        heatmap_filename_interactive, include_plotlyjs=False, full_html=False, config=_PLOTLY_CONFIG
    )

    heatmap_filename_png: str | None = None
    if static_report:
        heatmap_filename_png = os.path.join(plots_dir, plot_config["output_filename"])
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_heatmap_matplotlib

            generate_heatmap_matplotlib(ratio_data, heatmap_filename_png, threshold, plot_config)
        else:
            fig.write_image(heatmap_filename_png, width=1000, height=chart_height)

    logging.info("Heatmap generated.")

    return heatmap_filename_interactive, heatmap_filename_png


def create_metric_heatmap_data(summary_df: pd.DataFrame, metric_category: str) -> pd.DataFrame:
    """Extract and pivot a single metric category from summary_df to Sample x Plasmid grid.

    Args:
        summary_df: DataFrame with columns Sample, Plasmid, Category, Value
        metric_category: Category name to filter (e.g., "MeanDepthInsert", "CoverageOutsideINSERT")

    Returns:
        Pivoted DataFrame with samples as rows, plasmids as columns, values as cells.
        Returns empty DataFrame if metric_category not found (graceful degradation).
    """
    import pandas as pd

    # Filter for the specific metric category
    filtered = summary_df[summary_df["Category"] == metric_category].copy()

    if filtered.empty:
        logging.warning(f"Metric category '{metric_category}' not found in summary data")
        return pd.DataFrame()

    if metric_category == "MismatchesNearINSERT":
        from plasmicheck.scripts.generate_report import _parse_mismatches

        def _extract_mismatches(val: Any) -> float:
            if isinstance(val, (int, float)):
                return float(val)
            w, _ = _parse_mismatches(str(val))
            return float(w)

        filtered["Value"] = filtered["Value"].apply(_extract_mismatches)
    else:
        # Convert Value to numeric
        filtered["Value"] = pd.to_numeric(filtered["Value"], errors="coerce")

    # Pivot to Sample x Plasmid grid
    metric_data = filtered.pivot(index="Sample", columns="Plasmid", values="Value").round(
        ROUND_DECIMALS
    )

    return metric_data


def plot_metric_heatmap(
    metric_data: pd.DataFrame,
    output_dir: str,
    metric_name: str,
    vmax: float | None = None,
    static_report: bool = False,
    plot_backend: str = "matplotlib",
) -> tuple[str, str | None]:
    """Generate Plotly heatmap for any Sample x Plasmid metric grid.

    Args:
        metric_data: Pivoted DataFrame (Sample x Plasmid)
        output_dir: Output directory for plots
        metric_name: Display name for the metric (e.g., "Mean Depth Insert")
        vmax: Optional maximum value for color scale (data will be clipped)
        static_report: If True, also generate PNG
        plot_backend: "plotly" or "matplotlib" for PNG generation

    Returns:
        Tuple of (html_filename, png_filename_or_None)
    """
    import plotly.express as px

    if metric_data.empty:
        logging.warning(f"No data for {metric_name} heatmap")
        return ("", None)

    logging.info(f"Generating heatmap for {metric_name}.")

    # Clip data if vmax provided
    plot_data = metric_data.clip(upper=vmax) if vmax is not None else metric_data

    color_scale: list[list[float | str]] = [
        [0.0, "#f0f9ff"],
        [0.35, "#7dd3fc"],
        [0.7, "#0284c7"],
        [1.0, "#be123c"],
    ]

    n_samples, n_plasmids = plot_data.shape
    chart_height = max(420, min(640, 180 + n_samples * 28))
    show_annot = n_samples <= 20 and n_plasmids <= 15
    fmt = ".0f" if "Mismatches" in metric_name else ".2f"

    fig = px.imshow(
        plot_data,
        color_continuous_scale=color_scale,
        height=chart_height,
        aspect="auto",
        text_auto=fmt if show_annot else False,
    )
    fig.update_traces(xgap=2, ygap=2, textfont={"size": 11})

    # Update layout
    fig.update_layout(
        **_SUMMARY_PLOTLY_LAYOUT,
        xaxis_title="Plasmid",
        yaxis_title=None,
        xaxis={
            "tickangle": -30,
            "tickfont": {"size": 11, "color": "#334155"},
            "tickmode": "linear",
            "automargin": True,
        },
        yaxis={
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
        },
    )

    fig.update_coloraxes(
        showscale=True,
        colorbar={
            "thickness": 10,
            "len": 0.82,
            "tickfont": {"size": 10},
            "outlinewidth": 0,
        },
    )

    # Save interactive HTML
    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    slug = metric_name.lower().replace(" ", "_")
    html_filename: str = os.path.join(plots_dir, f"heatmap_{slug}.html")
    fig.write_html(html_filename, include_plotlyjs=False, full_html=False, config=_PLOTLY_CONFIG)

    # Save static PNG if requested
    png_filename: str | None = None
    if static_report:
        png_filename = os.path.join(plots_dir, f"heatmap_{slug}.png")
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_metric_heatmap_matplotlib

            generate_metric_heatmap_matplotlib(metric_data, png_filename, metric_name, vmax)
        else:
            fig.write_image(png_filename, width=1000, height=chart_height)

    logging.info(f"{metric_name} heatmap generated.")

    return html_filename, png_filename


def plot_coverage_outside_insert_boxplot(
    summary_df: pd.DataFrame,
    output_dir: str,
    static_report: bool = False,
    plot_backend: str = "matplotlib",
) -> tuple[str, str | None]:
    """Generate boxplot for CoverageOutsideINSERT distribution across plasmids.

    Args:
        summary_df: DataFrame with columns Sample, Plasmid, Category, Value
        output_dir: Output directory for plots
        static_report: If True, also generate PNG
        plot_backend: "plotly" or "matplotlib" for PNG generation

    Returns:
        Tuple of (html_filename, png_filename_or_None)
    """
    import pandas as pd
    import plotly.express as px

    # Extract CoverageOutsideINSERT rows
    coverage_df = summary_df[summary_df["Category"] == "CoverageOutsideINSERT"].copy()

    if coverage_df.empty:
        logging.warning("No CoverageOutsideINSERT data found")
        return ("", None)

    # Convert Value to numeric
    coverage_df["Value"] = pd.to_numeric(coverage_df["Value"], errors="coerce")

    logging.info("Generating CoverageOutsideINSERT boxplot.")

    # Pivot to Sample x Plasmid, then melt for boxplot
    coverage_pivot = coverage_df.pivot(index="Sample", columns="Plasmid", values="Value")
    boxplot_data = coverage_pivot.melt(
        var_name="Plasmid", value_name="Value", ignore_index=False
    ).reset_index()

    # Create boxplot
    fig = px.box(
        boxplot_data,
        x="Plasmid",
        y="Value",
        points="all",
        hover_data={"Sample": True, "Value": ":.4f"},
    )

    fig.update_layout(
        **_SUMMARY_PLOTLY_LAYOUT,
        xaxis_title="Plasmid",
        yaxis_title="Coverage Outside Insert",
        height=420,
        xaxis={
            "tickangle": -30,
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
            "showgrid": False,
            "linecolor": "#cbd5e1",
        },
        yaxis={
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
            "gridcolor": "#f1f5f9",
            "linecolor": "#cbd5e1",
        },
    )
    fig.update_traces(
        marker={**MARKER_STYLE, "color": "#0284c7"},
        line={"color": "#0369a1", "width": 1.6},
        fillcolor="rgba(2, 132, 199, 0.18)",
    )

    # Save interactive HTML
    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    html_filename: str = os.path.join(plots_dir, "boxplot_coverage_outside_insert.html")
    fig.write_html(html_filename, include_plotlyjs=False, full_html=False, config=_PLOTLY_CONFIG)

    # Save static PNG if requested
    png_filename: str | None = None
    if static_report:
        png_filename = os.path.join(plots_dir, "boxplot_coverage_outside_insert.png")
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_coverage_boxplot_matplotlib

            generate_coverage_boxplot_matplotlib(coverage_pivot, png_filename)
        else:
            fig.write_image(png_filename, width=1000, height=420)

    logging.info("CoverageOutsideINSERT boxplot generated.")

    return html_filename, png_filename


def create_resistance_gene_table(summary_df: pd.DataFrame) -> pd.DataFrame:
    """Extract resistance gene data from summary_df into a flat table.

    Args:
        summary_df: DataFrame with columns Sample, Plasmid, Category, Value

    Returns:
        DataFrame with columns: Sample, Plasmid, Gene Name, Mean Depth, Breadth (%)
        Returns empty DataFrame with correct columns if no genes found.
    """
    import pandas as pd

    from plasmicheck.scripts.resistance_gene_detection import sanitize_gene_name

    # Check if resistance genes are available
    available_rows = summary_df[summary_df["Category"] == "ResistanceGenesAvailable"]
    if available_rows.empty or (available_rows["Value"] == "false").all():
        # No genes available - return empty DataFrame with correct structure
        return pd.DataFrame(columns=["Sample", "Plasmid", "Gene Name", "Mean Depth", "Breadth (%)"])

    # Find all ResistanceGenes rows
    gene_rows = summary_df[summary_df["Category"] == "ResistanceGenes"]
    if gene_rows.empty:
        return pd.DataFrame(columns=["Sample", "Plasmid", "Gene Name", "Breadth (%)"])

    # Pre-index ResGene_ metrics for O(1) lookup (768x faster on large cohorts)
    res_rows = summary_df[summary_df["Category"].str.startswith("ResGene_", na=False)]
    res_lookup: dict[tuple[Any, Any, Any], Any] = {
        (s, p, c): v
        for s, p, c, v in zip(
            res_rows["Sample"],
            res_rows["Plasmid"],
            res_rows["Category"],
            res_rows["Value"],
            strict=True,
        )
    }

    # Build flat table
    table_rows: list[dict[str, Any]] = []

    for sample, plasmid, gene_names_str in zip(
        gene_rows["Sample"], gene_rows["Plasmid"], gene_rows["Value"], strict=True
    ):
        if not isinstance(gene_names_str, str) or not gene_names_str or gene_names_str == "None":
            continue

        # Split gene names by comma
        gene_names = [g.strip() for g in gene_names_str.split(",") if g.strip()]

        for gene_name in gene_names:
            # Sanitize gene name for lookup
            sanitized_name = sanitize_gene_name(gene_name)

            # Look up depth and breadth in O(1)
            depth_val = res_lookup.get((sample, plasmid, f"ResGene_{sanitized_name}_MeanDepth"))
            breadth_val = res_lookup.get((sample, plasmid, f"ResGene_{sanitized_name}_Breadth"))

            mean_depth = float(depth_val) if depth_val is not None else 0.0
            breadth = float(breadth_val) if breadth_val is not None else 0.0

            table_rows.append(
                {
                    "Sample": sample,
                    "Plasmid": plasmid,
                    "Gene Name": gene_name,
                    "Mean Depth": mean_depth,
                    "Breadth (%)": breadth,
                }
            )

    return pd.DataFrame(table_rows)


def _format_verdict_hover(
    verdict: str, count: int, name_list: list[str], max_display: int = 10
) -> str:
    """Build hover HTML string with verdict, count, and truncated name list.

    Args:
        verdict: Verdict category (e.g., "Clean", "Contaminated")
        count: Number of items in this category
        name_list: List of sample or plasmid names
        max_display: Maximum names to display before truncating

    Returns:
        HTML-formatted string for hover display
    """
    lines = [f"<b>{verdict}</b>", f"Count: {count}"]

    if name_list:
        if len(name_list) <= max_display:
            lines.extend(name_list)
        else:
            lines.extend(name_list[:max_display])
            lines.append(f"...and {len(name_list) - max_display} more")

    return "<br>".join(lines)


def aggregate_verdicts(
    summary_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Aggregate verdict data by plasmid and by sample with No Result detection.

    Args:
        summary_df: DataFrame with columns Sample, Plasmid, Category, Value

    Returns:
        Tuple of (full_verdict_df, by_plasmid_counts, by_sample_counts):
        - full_verdict_df: DataFrame with [Sample, Plasmid, VerdictShort] including No Result
        - by_plasmid_counts: Plasmids as index, verdict categories as columns
        - by_sample_counts: Samples as index, verdict categories as columns
    """
    import pandas as pd

    # Map full verdict strings to short labels
    verdict_map = {
        "Sample is contaminated with plasmid DNA": "Contaminated",
        "Sample contamination status is unclear": "Unclear",
        "Sample is not contaminated with plasmid DNA": "Clean",
    }

    # Extract verdict rows
    verdict_rows = summary_df[summary_df["Category"] == "Verdict"].copy()

    if verdict_rows.empty:
        # Return empty DataFrames with correct structure
        empty_df = pd.DataFrame(columns=["Sample", "Plasmid", "VerdictShort"])
        empty_counts = pd.DataFrame(columns=["Clean", "Unclear", "Contaminated", "No Result"])
        return empty_df, empty_counts, empty_counts

    # Map to short labels
    verdict_rows["VerdictShort"] = verdict_rows["Value"].map(verdict_map)

    # Build cross-product of all unique samples x plasmids from FULL summary_df
    all_samples = summary_df["Sample"].unique()
    all_plasmids = summary_df["Plasmid"].unique()

    # Create full cross-product
    full_cross = pd.DataFrame(
        [(s, p) for s in all_samples for p in all_plasmids], columns=["Sample", "Plasmid"]
    )

    # Left join verdict data
    full_verdict_df = full_cross.merge(
        verdict_rows[["Sample", "Plasmid", "VerdictShort"]], on=["Sample", "Plasmid"], how="left"
    )

    # Fill missing with "No Result"
    full_verdict_df["VerdictShort"] = full_verdict_df["VerdictShort"].fillna("No Result")

    # Aggregate by plasmid
    by_plasmid = full_verdict_df.groupby(["Plasmid", "VerdictShort"]).size().unstack(fill_value=0)

    # Aggregate by sample
    by_sample = full_verdict_df.groupby(["Sample", "VerdictShort"]).size().unstack(fill_value=0)

    # Ensure column order
    verdict_order = ["Clean", "Unclear", "Contaminated", "No Result"]
    by_plasmid_counts = by_plasmid.reindex(columns=verdict_order, fill_value=0)
    by_sample_counts = by_sample.reindex(columns=verdict_order, fill_value=0)

    return full_verdict_df, by_plasmid_counts, by_sample_counts


def plot_verdict_bar_by_plasmid(
    summary_df: pd.DataFrame,
    output_dir: str,
    static_report: bool = False,
    plot_backend: str = "matplotlib",
) -> tuple[str, str | None]:
    """Generate stacked bar chart of verdicts grouped by plasmid.

    Args:
        summary_df: DataFrame with columns Sample, Plasmid, Category, Value
        output_dir: Output directory for plots
        static_report: If True, also generate PNG via matplotlib
        plot_backend: Backend for static plots ("plotly" or "matplotlib")

    Returns:
        Tuple of (html_path, png_path_or_None)
    """
    import plotly.graph_objects as go

    from .plotting.colors import VERDICT_COLORS

    logging.info("Generating verdict bar chart by plasmid.")

    # Get aggregated data
    full_verdict_df, by_plasmid_counts, _ = aggregate_verdicts(summary_df)

    if full_verdict_df.empty:
        logging.warning("No verdict data available for plotting.")
        return "", None

    # Build hover data: verdict -> plasmid -> list of sample names
    hover_data: dict[str, dict[str, list[str]]] = {
        verdict: {} for verdict in ["Clean", "Unclear", "Contaminated", "No Result"]
    }

    for _, row in full_verdict_df.iterrows():
        verdict = row["VerdictShort"]
        plasmid = row["Plasmid"]
        sample = row["Sample"]

        if plasmid not in hover_data[verdict]:
            hover_data[verdict][plasmid] = []
        hover_data[verdict][plasmid].append(sample)

    # Create traces
    fig = go.Figure()
    verdict_order = ["Clean", "Unclear", "Contaminated", "No Result"]

    for verdict in verdict_order:
        if verdict not in by_plasmid_counts.columns:
            continue

        counts = by_plasmid_counts[verdict].values
        plasmids = by_plasmid_counts.index.tolist()

        # Build customdata for hover
        customdata = [
            _format_verdict_hover(verdict, int(counts[i]), hover_data[verdict].get(plasmids[i], []))
            for i in range(len(plasmids))
        ]

        fig.add_trace(
            go.Bar(
                x=plasmids,
                y=counts,
                name=verdict,
                marker_color=VERDICT_COLORS[verdict],
                customdata=customdata,
                hovertemplate="%{customdata}<extra></extra>",
            )
        )

    max_count = int(by_plasmid_counts.sum(axis=1).max()) if not by_plasmid_counts.empty else 1
    fig.update_layout(
        **_SUMMARY_PLOTLY_LAYOUT,
        barmode="stack",
        xaxis_title="Plasmid",
        yaxis_title="Sample Count",
        xaxis={
            "tickangle": -30,
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
            "showgrid": False,
            "linecolor": "#cbd5e1",
        },
        yaxis={
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
            "gridcolor": "#f1f5f9",
            "linecolor": "#cbd5e1",
            "dtick": 1 if max_count <= 10 else None,
        },
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "right",
            "x": 1.0,
            "font": {"size": 11, "color": "#334155"},
        },
        height=420,
    )

    # Write HTML
    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    html_path: str = os.path.join(plots_dir, "verdict_bar_by_plasmid.html")
    fig.write_html(html_path, include_plotlyjs=False, full_html=False, config=_PLOTLY_CONFIG)

    # Generate PNG if requested
    png_path: str | None = None
    if static_report:
        png_path = os.path.join(plots_dir, "verdict_bar_by_plasmid.png")
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_verdict_bar_matplotlib

            generate_verdict_bar_matplotlib(
                by_plasmid_counts, png_path, "Verdicts by Plasmid", "Plasmid"
            )
        else:
            # Fallback to plotly kaleido
            fig.write_image(png_path, width=1000, height=420)

    return html_path, png_path


def plot_verdict_bar_by_sample(
    summary_df: pd.DataFrame,
    output_dir: str,
    static_report: bool = False,
    plot_backend: str = "matplotlib",
) -> tuple[str, str | None]:
    """Generate stacked bar chart of verdicts grouped by sample.

    Args:
        summary_df: DataFrame with columns Sample, Plasmid, Category, Value
        output_dir: Output directory for plots
        static_report: If True, also generate PNG via matplotlib
        plot_backend: Backend for static plots ("plotly" or "matplotlib")

    Returns:
        Tuple of (html_path, png_path_or_None)
    """
    import plotly.graph_objects as go

    from .plotting.colors import VERDICT_COLORS

    logging.info("Generating verdict bar chart by sample.")

    # Get aggregated data
    full_verdict_df, _, by_sample_counts = aggregate_verdicts(summary_df)

    if full_verdict_df.empty:
        logging.warning("No verdict data available for plotting.")
        return "", None

    # Build hover data: verdict -> sample -> list of plasmid names
    hover_data: dict[str, dict[str, list[str]]] = {
        verdict: {} for verdict in ["Clean", "Unclear", "Contaminated", "No Result"]
    }

    for _, row in full_verdict_df.iterrows():
        verdict = row["VerdictShort"]
        sample = row["Sample"]
        plasmid = row["Plasmid"]

        if sample not in hover_data[verdict]:
            hover_data[verdict][sample] = []
        hover_data[verdict][sample].append(plasmid)

    # Create traces
    fig = go.Figure()
    verdict_order = ["Clean", "Unclear", "Contaminated", "No Result"]

    for verdict in verdict_order:
        if verdict not in by_sample_counts.columns:
            continue

        counts = by_sample_counts[verdict].values
        samples = by_sample_counts.index.tolist()

        # Build customdata for hover
        customdata = [
            _format_verdict_hover(verdict, int(counts[i]), hover_data[verdict].get(samples[i], []))
            for i in range(len(samples))
        ]

        fig.add_trace(
            go.Bar(
                x=samples,
                y=counts,
                name=verdict,
                marker_color=VERDICT_COLORS[verdict],
                customdata=customdata,
                hovertemplate="%{customdata}<extra></extra>",
            )
        )

    max_count = int(by_sample_counts.sum(axis=1).max()) if not by_sample_counts.empty else 1
    fig.update_layout(
        **_SUMMARY_PLOTLY_LAYOUT,
        barmode="stack",
        xaxis_title="Sample",
        yaxis_title="Plasmid Count",
        xaxis={
            "tickangle": -30,
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
            "showgrid": False,
            "linecolor": "#cbd5e1",
        },
        yaxis={
            "tickfont": {"size": 11, "color": "#334155"},
            "automargin": True,
            "gridcolor": "#f1f5f9",
            "linecolor": "#cbd5e1",
            "dtick": 1 if max_count <= 10 else None,
        },
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "right",
            "x": 1.0,
            "font": {"size": 11, "color": "#334155"},
        },
        height=420,
    )

    # Write HTML
    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    html_path: str = os.path.join(plots_dir, "verdict_bar_by_sample.html")
    fig.write_html(html_path, include_plotlyjs=False, full_html=False, config=_PLOTLY_CONFIG)

    # Generate PNG if requested
    png_path: str | None = None
    if static_report:
        png_path = os.path.join(plots_dir, "verdict_bar_by_sample.png")
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_verdict_bar_matplotlib

            generate_verdict_bar_matplotlib(
                by_sample_counts, png_path, "Verdicts by Sample", "Sample"
            )
        else:
            # Fallback to plotly kaleido
            fig.write_image(png_path, width=1000, height=420)

    return html_path, png_path


def _format_table_float(val: Any, decimals: int = 4, sci_threshold: float = 0.0) -> Any:
    """Format floating-point numbers cleanly for HTML table display."""
    import math

    if val is None or (isinstance(val, float) and math.isnan(val)):
        return "N/A"
    if isinstance(val, (int, bool)):
        return val
    try:
        fval = float(val)
        if math.isnan(fval):
            return "N/A"
        if fval == int(fval) and isinstance(val, str) and "." not in val:
            return val
        if sci_threshold > 0 and 0 < abs(fval) < sci_threshold:
            return f"{fval:.2e}"
        return f"{fval:.{decimals}f}"
    except (ValueError, TypeError):
        return val


# Report Generation Functions
def generate_report(
    combined_df: pd.DataFrame,
    verdict_df: pd.DataFrame,
    ratio_df: pd.DataFrame,
    heatmap_filename_interactive: str,
    heatmap_filename_png: str | None,
    boxplot_filename_interactive: str,
    boxplot_filename_png: str | None,
    p_value_table_filename: str,
    output_folder: str,
    threshold: float,
    unclear_range: dict[str, float],
    command_line: str,
    static_report: bool = False,
    plotly_mode: str = "directory",
    plotly_js_path: str = "",
    plotly_version: str = "",
    plotly_js_inline: str = "",
    coverage_outside_heatmap_interactive: str = "",
    coverage_outside_heatmap_png: str | None = None,
    coverage_outside_boxplot_interactive: str = "",
    coverage_outside_boxplot_png: str | None = None,
    mismatches_heatmap_interactive: str = "",
    mismatches_heatmap_png: str | None = None,
    depth_heatmap_interactive: str = "",
    depth_heatmap_png: str | None = None,
    breadth_heatmap_interactive: str = "",
    breadth_heatmap_png: str | None = None,
    resistance_gene_df: pd.DataFrame | None = None,
    junction_df: pd.DataFrame | None = None,
    verdict_by_plasmid_interactive: str = "",
    verdict_by_plasmid_png: str | None = None,
    verdict_by_sample_interactive: str = "",
    verdict_by_sample_png: str | None = None,
) -> None:
    import pandas as pd
    from jinja2 import Environment, FileSystemLoader

    logging.info("Generating summary reports")

    # Apply sorting to each DataFrame
    combined_df = apply_sorting(combined_df, TABLE_SORTING.get("combined"))
    verdict_df = apply_sorting(verdict_df, TABLE_SORTING.get("verdict"))
    if os.path.exists(p_value_table_filename) and os.path.getsize(p_value_table_filename) > 0:
        try:
            p_values_df = pd.read_csv(p_value_table_filename, sep="\t")
        except pd.errors.EmptyDataError:
            p_values_df = pd.DataFrame(
                columns=["Sample", "Plasmid", "p_value", "p_value_corrected"]
            )
    else:
        p_values_df = pd.DataFrame(columns=["Sample", "Plasmid", "p_value", "p_value_corrected"])
    p_values_df = apply_sorting(p_values_df, TABLE_SORTING.get("p_value"))

    env = Environment(
        loader=FileSystemLoader(str(get_resource_path(TEMPLATE_DIR))), autoescape=True
    )
    template = env.get_template("summary_template.html")

    # Read the Plotly HTML files content for interactive report
    with open(heatmap_filename_interactive) as f:
        heatmap_html_content = f.read()

    with open(boxplot_filename_interactive) as f:
        boxplot_html_content = f.read()

    # Read new plot HTML files (if they exist)
    coverage_outside_heatmap_html = ""
    if coverage_outside_heatmap_interactive:
        with open(coverage_outside_heatmap_interactive) as f:
            coverage_outside_heatmap_html = f.read()

    coverage_outside_boxplot_html = ""
    if coverage_outside_boxplot_interactive:
        with open(coverage_outside_boxplot_interactive) as f:
            coverage_outside_boxplot_html = f.read()

    mismatches_heatmap_html = ""
    if mismatches_heatmap_interactive:
        with open(mismatches_heatmap_interactive) as f:
            mismatches_heatmap_html = f.read()

    depth_heatmap_html = ""
    if depth_heatmap_interactive:
        with open(depth_heatmap_interactive) as f:
            depth_heatmap_html = f.read()

    breadth_heatmap_html = ""
    if breadth_heatmap_interactive:
        with open(breadth_heatmap_interactive) as f:
            breadth_heatmap_html = f.read()

    verdict_by_plasmid_html = ""
    if verdict_by_plasmid_interactive:
        with open(verdict_by_plasmid_interactive) as f:
            verdict_by_plasmid_html = f.read()

    verdict_by_sample_html = ""
    if verdict_by_sample_interactive:
        with open(verdict_by_sample_interactive) as f:
            verdict_by_sample_html = f.read()

    # Read the PNG files as base64 for non-interactive report (only if static_report=True)
    heatmap_png_base64 = ""
    boxplot_png_base64 = ""
    if static_report and heatmap_filename_png and boxplot_filename_png:
        heatmap_png_base64 = encode_image_to_base64(heatmap_filename_png)
        boxplot_png_base64 = encode_image_to_base64(boxplot_filename_png)

    # Read new plot PNG files (if they exist and static_report=True)
    coverage_outside_heatmap_png_base64 = ""
    coverage_outside_boxplot_png_base64 = ""
    mismatches_heatmap_png_base64 = ""
    depth_heatmap_png_base64 = ""
    breadth_heatmap_png_base64 = ""
    verdict_by_plasmid_png_base64 = ""
    verdict_by_sample_png_base64 = ""
    if static_report:
        if coverage_outside_heatmap_png:
            coverage_outside_heatmap_png_base64 = encode_image_to_base64(
                coverage_outside_heatmap_png
            )
        if coverage_outside_boxplot_png:
            coverage_outside_boxplot_png_base64 = encode_image_to_base64(
                coverage_outside_boxplot_png
            )
        if mismatches_heatmap_png:
            mismatches_heatmap_png_base64 = encode_image_to_base64(mismatches_heatmap_png)
        if depth_heatmap_png:
            depth_heatmap_png_base64 = encode_image_to_base64(depth_heatmap_png)
        if breadth_heatmap_png:
            breadth_heatmap_png_base64 = encode_image_to_base64(breadth_heatmap_png)
        if verdict_by_plasmid_png:
            verdict_by_plasmid_png_base64 = encode_image_to_base64(verdict_by_plasmid_png)
        if verdict_by_sample_png:
            verdict_by_sample_png_base64 = encode_image_to_base64(verdict_by_sample_png)

    # Format display copies of DataFrames for clean HTML rendering
    combined_display_df = combined_df.copy()
    if "Value" in combined_display_df.columns:
        combined_display_df["Value"] = combined_display_df["Value"].apply(
            lambda v: _format_table_float(v, decimals=4)
        )

    ratio_display_df = ratio_df.copy()
    if "Value" in ratio_display_df.columns:
        ratio_display_df["Value"] = ratio_display_df["Value"].apply(
            lambda v: _format_table_float(v, decimals=4)
        )

    p_values_display_df = p_values_df.copy()
    for col in ("p_value", "p_value_corrected"):
        if col in p_values_display_df.columns:
            p_values_display_df[col] = p_values_display_df[col].apply(
                lambda v: _format_table_float(v, decimals=4, sci_threshold=1e-3)
            )

    # Convert sorted DataFrames to HTML
    combined_html = combined_display_df.to_html(
        classes="table table-striped table-bordered", index=False
    )
    verdict_html = verdict_df.to_html(classes="table table-striped table-bordered", index=False)
    ratio_html = ratio_display_df.to_html(classes="table table-striped table-bordered", index=False)
    p_value_html = p_values_display_df.to_html(
        classes="table table-striped table-bordered", index=False
    )

    # Convert resistance gene DataFrame to HTML if provided
    resistance_gene_html = ""
    if resistance_gene_df is not None and not resistance_gene_df.empty:
        res_display_df = resistance_gene_df.copy()
        for col in ("Mean Depth", "Breadth (%)"):
            if col in res_display_df.columns:
                res_display_df[col] = res_display_df[col].apply(
                    lambda v: _format_table_float(v, decimals=2)
                )
        resistance_gene_html = res_display_df.to_html(
            classes="table table-striped table-bordered", index=False
        )

    # Convert junction metadata DataFrame to HTML with truncation
    junction_table_html = ""
    if junction_df is not None and not junction_df.empty:
        # Truncate long comma-separated values for HTML display
        display_df = junction_df.copy()
        display_df.columns.name = None
        for col in ["ExonLengths", "IntronLengths"]:
            if col in display_df.columns:
                display_df[col] = display_df[col].apply(
                    lambda v: _truncate_list_for_html(str(v)) if pd.notna(v) else "N/A"
                )
        junction_table_html = display_df.to_html(
            classes="table table-striped table-bordered", index=False
        )

    # Encode the logo
    logo_base64 = encode_image_to_base64(str(get_resource_path(LOGO_PATH)))

    # Render interactive HTML report
    html_content_interactive = template.render(
        combined_df=combined_html,
        verdict_df=verdict_html,
        ratio_df=ratio_html,
        heatmap_content=heatmap_html_content,
        boxplot_content=boxplot_html_content,
        p_value_table=p_value_html,
        coverage_outside_heatmap_content=coverage_outside_heatmap_html,
        coverage_outside_boxplot_content=coverage_outside_boxplot_html,
        mismatches_heatmap_content=mismatches_heatmap_html,
        depth_heatmap_content=depth_heatmap_html,
        breadth_heatmap_content=breadth_heatmap_html,
        resistance_gene_table=resistance_gene_html,
        junction_table=junction_table_html,
        verdict_by_plasmid_content=verdict_by_plasmid_html,
        verdict_by_sample_content=verdict_by_sample_html,
        logo_base64=logo_base64,
        threshold=threshold,
        unclear_range=unclear_range,
        unclear_lower=unclear_range.get("lower_bound", 0.6),
        unclear_upper=unclear_range.get("upper_bound", 1.0),
        contam_cutoff=unclear_range.get("upper_bound", threshold),
        version=VERSION,
        run_date=datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line,
        interactive=True,
        plotly_mode=plotly_mode,
        plotly_version=plotly_version,
        plotly_js_path=plotly_js_path,
        plotly_js_inline=plotly_js_inline,
    )

    # Save interactive HTML report
    html_report_interactive = os.path.join(output_folder, "summary_report_interactive.html")
    with open(html_report_interactive, "w") as f:
        f.write(html_content_interactive)

    # Save non-interactive HTML report (only if static_report=True)
    if static_report and heatmap_png_base64 and boxplot_png_base64:
        html_content_non_interactive = template.render(
            combined_df=combined_html,
            verdict_df=verdict_html,
            ratio_df=ratio_html,
            heatmap_content=heatmap_png_base64,
            boxplot_content=boxplot_png_base64,
            p_value_table=p_value_html,
            coverage_outside_heatmap_content=coverage_outside_heatmap_png_base64,
            coverage_outside_boxplot_content=coverage_outside_boxplot_png_base64,
            mismatches_heatmap_content=mismatches_heatmap_png_base64,
            depth_heatmap_content=depth_heatmap_png_base64,
            breadth_heatmap_content=breadth_heatmap_png_base64,
            resistance_gene_table=resistance_gene_html,
            junction_table=junction_table_html,
            verdict_by_plasmid_content=verdict_by_plasmid_png_base64,
            verdict_by_sample_content=verdict_by_sample_png_base64,
            logo_base64=logo_base64,
            threshold=threshold,
            unclear_range=unclear_range,
            unclear_lower=unclear_range.get("lower_bound", 0.6),
            unclear_upper=unclear_range.get("upper_bound", 1.0),
            contam_cutoff=unclear_range.get("upper_bound", threshold),
            version=VERSION,
            run_date=datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
            command_line=command_line,
            interactive=False,
            plotly_mode=plotly_mode,
            plotly_version=plotly_version,
            plotly_js_path=plotly_js_path,
            plotly_js_inline=plotly_js_inline,
        )

        html_report_non_interactive = os.path.join(
            output_folder, "summary_report_non_interactive.html"
        )
        with open(html_report_non_interactive, "w") as f:
            f.write(html_content_non_interactive)


def main(
    input_dir: str,
    output_dir: str,
    threshold: float = DEFAULT_THRESHOLD,
    unclear_range: dict[str, float] = UNCLEAR_RANGE,
    substring_to_remove: str | None = None,
    static_report: bool = False,
    plotly_mode: str = "directory",
    plot_backend: str = "matplotlib",
) -> None:
    # Initialize kaleido once before plotting if static_report is True and using plotly backend
    if static_report and plot_backend == "plotly":
        import kaleido  # type: ignore[import-untyped]

        kaleido.start_sync_server()

    # Setup plotly.js asset management
    plotly_js_path = ""
    plotly_version = ""
    plotly_js_inline = ""

    if plotly_mode in ("directory", "cdn"):
        import plotly

        plotly_version = plotly.__version__

    if plotly_mode == "directory":
        ensure_plotly_assets(output_dir)
        plotly_js_path = "assets/plotly.min.js"  # Relative from output_dir

    if plotly_mode == "embedded":
        import plotly as _plotly

        plotly_js_source = os.path.join(
            os.path.dirname(_plotly.__file__), "package_data", "plotly.min.js"
        )
        with open(plotly_js_source) as f:
            plotly_js_inline = f"<script>{f.read()}</script>"

    import pandas as pd

    summary_df = read_compare_outputs(input_dir, substring_to_remove=substring_to_remove)

    # Derive effective threshold and bounds
    effective_threshold = threshold
    thresh_rows = summary_df[summary_df["Category"] == "Threshold"]
    if not thresh_rows.empty:
        try:
            val_col = "Value" if "Value" in thresh_rows.columns else "Count"
            effective_threshold = float(thresh_rows[val_col].iloc[0])
        except (ValueError, TypeError):
            effective_threshold = threshold

    unclear_lower, unclear_upper = verdict_bounds(effective_threshold, unclear_range)
    effective_unclear = {"lower_bound": unclear_lower, "upper_bound": unclear_upper}

    # Recompute display verdicts from valid stored Ratio values under effective threshold/bounds
    category_to_full = {
        "Contaminated": "Sample is contaminated with plasmid DNA",
        "Unclear": "Sample contamination status is unclear",
        "Clean": "Sample is not contaminated with plasmid DNA",
    }
    ratio_rows = summary_df[summary_df["Category"] == "Ratio"]
    if not ratio_rows.empty:
        ratio_map: dict[tuple[str, str], str] = {}
        for _, r_row in ratio_rows.iterrows():
            val = pd.to_numeric(r_row["Value"], errors="coerce")
            if pd.notna(val):
                cat = categorize_value(float(val), effective_threshold, effective_unclear)
                ratio_map[(str(r_row["Sample"]), str(r_row["Plasmid"]))] = category_to_full.get(
                    cat, cat
                )

        def _recalc_verdict(v_row: pd.Series) -> Any:
            key = (str(v_row["Sample"]), str(v_row["Plasmid"]))
            return ratio_map.get(key, v_row["Value"])

        verdict_mask = summary_df["Category"] == "Verdict"
        if verdict_mask.any():
            summary_df.loc[verdict_mask, "Value"] = summary_df.loc[verdict_mask].apply(
                _recalc_verdict, axis=1
            )

    # Calculate variations
    boxplot_data, p_values_df = calculate_variations(summary_df[summary_df["Category"] == "Ratio"])

    # Generate plots
    boxplot_filename_interactive, boxplot_filename_png = plot_boxplot(
        boxplot_data,
        output_dir,
        static_report=static_report,
        plot_backend=plot_backend,
        threshold=effective_threshold,
        unclear_range=effective_unclear,
    )
    heatmap_filename_interactive, heatmap_filename_png = plot_heatmap(
        summary_df[summary_df["Category"] == "Ratio"],
        output_dir,
        effective_threshold,
        effective_unclear,
        PLOT_CONFIG,
        static_report=static_report,
        plot_backend=plot_backend,
    )

    # --- New: Coverage Outside Insert (REPT-04) ---
    cov_outside_data = create_metric_heatmap_data(summary_df, "CoverageOutsideINSERT")
    cov_outside_heatmap_interactive, cov_outside_heatmap_png = plot_metric_heatmap(
        cov_outside_data,
        output_dir,
        "Coverage Outside Insert",
        static_report=static_report,
        plot_backend=plot_backend,
    )
    cov_outside_boxplot_interactive, cov_outside_boxplot_png = plot_coverage_outside_insert_boxplot(
        summary_df, output_dir, static_report=static_report, plot_backend=plot_backend
    )

    # --- New: MismatchesNearINSERT (REPT-05) ---
    mismatches_data = create_metric_heatmap_data(summary_df, "MismatchesNearINSERT")
    mismatches_heatmap_interactive, mismatches_heatmap_png = plot_metric_heatmap(
        mismatches_data,
        output_dir,
        "Mismatches Near Insert",
        static_report=static_report,
        plot_backend=plot_backend,
    )

    # --- New: Coverage Depth (REPT-06) ---
    depth_data = create_metric_heatmap_data(summary_df, "MeanDepthBackbone")
    depth_heatmap_interactive, depth_heatmap_png = plot_metric_heatmap(
        depth_data,
        output_dir,
        "Backbone Mean Depth",
        vmax=DEPTH_HEATMAP_CAP,
        static_report=static_report,
        plot_backend=plot_backend,
    )

    # --- New: Coverage Breadth (REPT-07) ---
    breadth_data = create_metric_heatmap_data(summary_df, "BreadthBackbone")
    breadth_vmax = (
        1.0 if (not breadth_data.empty and float(breadth_data.max().max()) <= 1.0) else 100.0
    )
    breadth_heatmap_interactive, breadth_heatmap_png = plot_metric_heatmap(
        breadth_data,
        output_dir,
        "Backbone Breadth",
        vmax=breadth_vmax,
        static_report=static_report,
        plot_backend=plot_backend,
    )

    # --- New: Resistance Gene Table (REPT-08) ---
    resistance_gene_df = create_resistance_gene_table(summary_df)

    # --- New: Verdict Bar Charts (VIS-01, VIS-02) ---
    verdict_by_plasmid_interactive, verdict_by_plasmid_png = plot_verdict_bar_by_plasmid(
        summary_df, output_dir, static_report=static_report, plot_backend=plot_backend
    )
    verdict_by_sample_interactive, verdict_by_sample_png = plot_verdict_bar_by_sample(
        summary_df, output_dir, static_report=static_report, plot_backend=plot_backend
    )

    # --- New: Splice Junction Metadata Table (SJ-02) ---
    junction_categories = [
        "SpliceJunctionCount",
        "ExonCount",
        "SpliceJunctionReads",
        "SpliceJunctionEvents",
        "SpliceJunctionsCovered",
        "ExonIntronBoundaryReads",
        "EstimatedContaminationFraction",
        "VectorInsertJunctionReads",
        "VectorInsertBridgePairs",
        "ExonLengths",
        "IntronLengths",
    ]
    junction_df = summary_df[summary_df["Category"].isin(junction_categories)]

    junction_table_df = None
    if not junction_df.empty:
        junction_pivot = junction_df.pivot_table(
            index=["Sample", "Plasmid"],
            columns="Category",
            values="Value",
            aggfunc="first",
        ).reset_index()
        junction_pivot.columns.name = None
        # Reorder columns
        col_order = ["Sample", "Plasmid"] + [
            c for c in junction_categories if c in junction_pivot.columns
        ]
        junction_table_df = junction_pivot[col_order]

    # --- Build extended metrics DataFrame for export (REPT-10) ---
    extended_categories = [
        "CoverageOutsideINSERT",
        "MismatchesNearINSERT",
        "VectorInsertJunctionReads",
        "VectorInsertBridgePairs",
        "MeanDepthBackbone",
        "BreadthBackbone",
        "SpliceJunctionCount",
        "ExonCount",
        "SpliceJunctionReads",
        "SpliceJunctionEvents",
        "SpliceJunctionsCovered",
        "ExonIntronBoundaryReads",
        "EstimatedContaminationFraction",
        "ExonLengths",
        "IntronLengths",
    ]
    extended_metrics_df = summary_df[summary_df["Category"].isin(extended_categories)]

    # Save tables
    save_tables_as_tsv_and_excel(
        summary_df[summary_df["Category"].isin(["Plasmid", "Human", "Tied"])],
        summary_df[summary_df["Category"] == "Verdict"],
        summary_df[summary_df["Category"] == "Ratio"],
        p_values_df,
        output_dir,
        extended_metrics_df=extended_metrics_df,
        resistance_gene_df=resistance_gene_df,
        junction_df=junction_table_df,
    )

    # Generate report
    generate_report(
        summary_df[summary_df["Category"].isin(["Plasmid", "Human", "Tied"])],
        summary_df[summary_df["Category"] == "Verdict"],
        summary_df[summary_df["Category"] == "Ratio"],
        heatmap_filename_interactive,
        heatmap_filename_png,
        boxplot_filename_interactive,
        boxplot_filename_png,
        os.path.join(output_dir, "data", "p_value_table.tsv"),
        output_dir,
        effective_threshold,
        effective_unclear,
        " ".join(sys.argv),
        static_report=static_report,
        plotly_mode=plotly_mode,
        plotly_js_path=plotly_js_path,
        plotly_version=plotly_version,
        plotly_js_inline=plotly_js_inline,
        coverage_outside_heatmap_interactive=cov_outside_heatmap_interactive,
        coverage_outside_heatmap_png=cov_outside_heatmap_png,
        coverage_outside_boxplot_interactive=cov_outside_boxplot_interactive,
        coverage_outside_boxplot_png=cov_outside_boxplot_png,
        mismatches_heatmap_interactive=mismatches_heatmap_interactive,
        mismatches_heatmap_png=mismatches_heatmap_png,
        depth_heatmap_interactive=depth_heatmap_interactive,
        depth_heatmap_png=depth_heatmap_png,
        breadth_heatmap_interactive=breadth_heatmap_interactive,
        breadth_heatmap_png=breadth_heatmap_png,
        resistance_gene_df=resistance_gene_df,
        junction_df=junction_table_df,
        verdict_by_plasmid_interactive=verdict_by_plasmid_interactive,
        verdict_by_plasmid_png=verdict_by_plasmid_png,
        verdict_by_sample_interactive=verdict_by_sample_interactive,
        verdict_by_sample_png=verdict_by_sample_png,
    )
