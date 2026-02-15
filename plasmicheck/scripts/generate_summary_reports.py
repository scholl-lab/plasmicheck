from __future__ import annotations

import base64
import logging
import os
import sys
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any

from plasmicheck.config import get_config
from plasmicheck.resources import get_resource_path
from plasmicheck.version import __version__ as VERSION

from .utils import add_logging_args, configure_logging_from_args

if TYPE_CHECKING:
    import pandas as pd

_cfg = get_config()
DEFAULT_THRESHOLD: float = _cfg["default_threshold"]
UNCLEAR_RANGE: dict[str, float] = _cfg["unclear_range"]
PLOT_CONFIG: dict[str, Any] = _cfg["plot_summary"]
TEMPLATE_DIR: str = _cfg["paths"]["template_dir"]
LOGO_PATH: str = _cfg["paths"]["logo_path"]
TABLE_SORTING: dict[str, Any] = _cfg["table_sorting"]
PLOT_DIMENSIONS: dict[str, int] = PLOT_CONFIG.get(
    "plot_dimensions", {"width": 1200, "height": 1200}
)
ROUND_DECIMALS: int = PLOT_CONFIG.get("round_decimals", 3)
LOG_OFFSET: float = PLOT_CONFIG.get("log_offset", 1e-9)
MARKER_STYLE: dict[str, Any] = PLOT_CONFIG.get("marker_style", {"size": 8, "opacity": 0.7})
DEPTH_HEATMAP_CAP: float = PLOT_CONFIG.get("depth_heatmap_cap", 100.0)

# Setup logging for the libraries used in this script
logging.getLogger("matplotlib").setLevel(logging.ERROR)
logging.getLogger("seaborn").setLevel(logging.ERROR)
logging.getLogger("jinja2").setLevel(logging.ERROR)
logging.getLogger("fontTools").setLevel(logging.ERROR)


# Helper Functions
def encode_image_to_base64(image_path: str) -> str:
    logging.info(f"Encoding image {image_path} to base64")
    with open(image_path, "rb") as image_file:
        encoded_string: str = base64.b64encode(image_file.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded_string}"


def find_tsv_files(input_dir: str, pattern: str) -> list[str]:
    """Recursively find all tsv files matching the pattern in the input directory."""
    tsv_files: list[str] = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(pattern):
                tsv_files.append(os.path.join(root, file))
    return tsv_files


def read_compare_outputs(
    input_dir: str, substring_to_remove: str | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    import pandas as pd

    logging.info(f"Reading compare outputs from {input_dir}")
    reads_assignment_files: list[str] = find_tsv_files(input_dir, ".reads_assignment.tsv")
    summary_files: list[str] = find_tsv_files(input_dir, ".summary.tsv")

    if not reads_assignment_files or not summary_files:
        logging.error("No valid TSV files found in the specified directory structure.")
        raise ValueError("No valid TSV files found in the specified directory structure.")

    reads_df_list: list[pd.DataFrame] = []
    for file in reads_assignment_files:
        sample_name: str = os.path.basename(os.path.dirname(os.path.dirname(file)))
        plasmid_name: str = os.path.basename(os.path.dirname(file))

        # Remove the substring from the sample name if provided
        if substring_to_remove:
            sample_name = sample_name.replace(substring_to_remove, "")

        df = pd.read_csv(
            file, sep="\t", low_memory=False
        )  # Disable low_memory mode to avoid DtypeWarning
        df["Sample"] = sample_name
        df["Plasmid"] = plasmid_name
        reads_df_list.append(df)

    summary_df_list: list[pd.DataFrame] = []
    for file in summary_files:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(file)))
        plasmid_name = os.path.basename(os.path.dirname(file))

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

    reads_df = pd.concat(reads_df_list, ignore_index=True)
    summary_df = pd.concat(summary_df_list, ignore_index=True)

    return reads_df, summary_df


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


def save_tables_as_tsv_and_excel(
    combined_df: pd.DataFrame,
    verdict_df: pd.DataFrame,
    ratio_df: pd.DataFrame,
    p_values_df: pd.DataFrame,
    output_dir: str,
    extended_metrics_df: pd.DataFrame | None = None,
    resistance_gene_df: pd.DataFrame | None = None,
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

    logging.info("TSV and Excel files have been saved.")


# Data Manipulation Functions
def calculate_variations(ratio_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate p-values and apply FDR correction for contamination ratios by plasmid."""
    import pandas as pd
    import scipy.stats as stats
    from statsmodels.stats.multitest import multipletests

    logging.info("Calculating variations for contamination ratios.")

    # Convert 'Value' to numeric, ensuring all entries are valid numbers
    ratio_df.loc[:, "Value"] = pd.to_numeric(ratio_df["Value"], errors="coerce")
    boxplot_data: pd.DataFrame = ratio_df.pivot(index="Sample", columns="Plasmid", values="Value")

    # Drop NaNs
    boxplot_data = boxplot_data.dropna(how="all", axis=1)

    p_values_list: list[dict[str, Any]] = []
    for plasmid in boxplot_data.columns:
        # Ensure the data is numeric and drop NaN values
        data = pd.to_numeric(boxplot_data[plasmid], errors="coerce").dropna()

        # Only run the t-test if there are enough data points
        if len(data) > 1:
            for i in range(len(data)):
                # Exclude the current sample and run the one-sample t-test
                try:
                    _, p_value = stats.ttest_1samp(data.drop(data.index[i]), data.iloc[i])
                    p_values_list.append(
                        {
                            "Sample": data.index[i],
                            "Plasmid": plasmid,
                            "p_value": p_value,
                        }
                    )
                except Exception as e:
                    logging.error(f"Error during t-test for {plasmid}, sample {data.index[i]}: {e}")

    p_values_df: pd.DataFrame = pd.DataFrame(p_values_list)
    if not p_values_df.empty:
        p_values_df["p_value_corrected"] = multipletests(p_values_df["p_value"], method="fdr_bh")[1]

    return boxplot_data, p_values_df


# Plotting Functions
def plot_boxplot(
    boxplot_data: pd.DataFrame,
    output_dir: str,
    static_report: bool = False,
    plot_backend: str = "plotly",
) -> tuple[str, str | None]:
    """Generate and save boxplots for contamination ratios by plasmid."""
    import numpy as np
    import pandas as pd
    import plotly.express as px

    logging.info("Generating boxplots for contamination ratios.")

    boxplot_df = pd.melt(
        boxplot_data.reset_index(), id_vars=["Sample"], var_name="Plasmid", value_name="Value"
    )
    boxplot_df["log_value"] = boxplot_df["Value"].apply(lambda x: np.log10(x + LOG_OFFSET))

    fig = px.box(
        boxplot_df,
        x="Plasmid",
        y="log_value",
        points="all",
        hover_data={"Sample": True, "Value": True, "log_value": False},
        title="Contamination Ratios by Plasmid (Log Scale)",
    )

    # Increase the height of the boxplot
    fig.update_layout(
        yaxis_title="Log of Value",
        height=800,  # Adjust this value to make the boxplot taller
    )

    fig.update_traces(marker=MARKER_STYLE)

    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    boxplot_filename: str = os.path.join(plots_dir, "boxplot_contamination_ratios.html")
    fig.write_html(boxplot_filename, include_plotlyjs=False, full_html=False)

    boxplot_filename_png: str | None = None
    if static_report:
        boxplot_filename_png = os.path.join(plots_dir, "boxplot_contamination_ratios.png")
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_summary_boxplot_matplotlib

            generate_summary_boxplot_matplotlib(
                boxplot_data, boxplot_filename_png, LOG_OFFSET, "outliers"
            )
        else:
            fig.write_image(boxplot_filename_png)

    logging.info("Boxplots generated.")

    return boxplot_filename, boxplot_filename_png


def plot_heatmap(
    ratio_df: pd.DataFrame,
    output_dir: str,
    threshold: float,
    unclear_range: dict[str, float],
    plot_config: dict[str, Any],
    static_report: bool = False,
    plot_backend: str = "plotly",
) -> tuple[str, str | None]:
    """Generate and save heatmap for contamination ratios by plasmid without the color bar."""
    import pandas as pd
    import plotly.express as px

    logging.info("Generating heatmap for contamination ratios.")

    # Use .loc[] to avoid SettingWithCopyWarning
    ratio_df.loc[:, "Value"] = pd.to_numeric(ratio_df["Value"], errors="coerce")

    # Pivot the DataFrame for the heatmap
    ratio_data: pd.DataFrame = ratio_df.pivot(
        index="Sample", columns="Plasmid", values="Value"
    ).round(ROUND_DECIMALS)

    # Define a custom color scale using Plotly's built-in scaling
    max_val = ratio_data.max().max()
    if max_val > 0:
        lower_norm = min(unclear_range["lower_bound"] / max_val, 0.99)
        upper_norm = min(unclear_range["upper_bound"] / max_val, 0.995)
    else:
        lower_norm = 0.5
        upper_norm = 0.75
    color_scale: list[list[float | str]] = [
        [0.0, plot_config["colors"]["not_contaminated"]],  # Start
        [lower_norm, plot_config["colors"]["unclear"]],  # Lower bound
        [upper_norm, plot_config["colors"]["unclear"]],  # Upper bound
        [1.0, plot_config["colors"]["contaminated"]],  # End
    ]

    fig = px.imshow(
        ratio_data,
        color_continuous_scale=color_scale,
        height=PLOT_DIMENSIONS["height"],
        aspect="auto",  # Adjust the aspect ratio to fit the plot dimensions
    )

    # Update layout for axis labels and title
    fig.update_layout(
        title=plot_config["title"],
        xaxis_title="Plasmid",
        yaxis_title="Sample",
        xaxis={
            "tickangle": plot_config["xticks_rotation"],
            "tickfont": {"size": 9},
            "tickmode": "linear",
        },
        yaxis={"tickfont": {"size": 9}},
    )

    # Disable the color bar
    fig.update_coloraxes(showscale=False)

    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    heatmap_filename_interactive: str = os.path.join(
        plots_dir, plot_config["output_filename"].replace(".png", ".html")
    )
    fig.write_html(heatmap_filename_interactive, include_plotlyjs=False, full_html=False)

    heatmap_filename_png: str | None = None
    if static_report:
        heatmap_filename_png = os.path.join(plots_dir, plot_config["output_filename"])
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_heatmap_matplotlib

            generate_heatmap_matplotlib(ratio_data, heatmap_filename_png, threshold, plot_config)
        else:
            fig.write_image(heatmap_filename_png)

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
    plot_backend: str = "plotly",
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

    # Blue-to-red color scale matching existing ratio heatmap
    color_scale: list[list[float | str]] = [
        [0.0, PLOT_CONFIG["colors"]["not_contaminated"]],
        [1.0, PLOT_CONFIG["colors"]["contaminated"]],
    ]

    fig = px.imshow(
        plot_data,
        color_continuous_scale=color_scale,
        height=PLOT_DIMENSIONS["height"],
        aspect="auto",
    )

    # Update layout
    fig.update_layout(
        title=f"Heatmap of {metric_name} by Sample and Plasmid",
        xaxis_title="Plasmid",
        yaxis_title="Sample",
        xaxis={
            "tickangle": PLOT_CONFIG["xticks_rotation"],
            "tickfont": {"size": 9},
            "tickmode": "linear",
        },
        yaxis={"tickfont": {"size": 9}},
    )

    # Disable color bar
    fig.update_coloraxes(showscale=False)

    # Save interactive HTML
    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    slug = metric_name.lower().replace(" ", "_")
    html_filename: str = os.path.join(plots_dir, f"heatmap_{slug}.html")
    fig.write_html(html_filename, include_plotlyjs=False, full_html=False)

    # Save static PNG if requested
    png_filename: str | None = None
    if static_report:
        png_filename = os.path.join(plots_dir, f"heatmap_{slug}.png")
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_metric_heatmap_matplotlib

            generate_metric_heatmap_matplotlib(metric_data, png_filename, metric_name, vmax)
        else:
            fig.write_image(png_filename)

    logging.info(f"{metric_name} heatmap generated.")

    return html_filename, png_filename


def plot_coverage_outside_insert_boxplot(
    summary_df: pd.DataFrame,
    output_dir: str,
    static_report: bool = False,
    plot_backend: str = "plotly",
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
        hover_data={"Sample": True},
        title="CoverageOutsideINSERT Distribution by Plasmid",
    )

    fig.update_layout(
        yaxis_title="Coverage Outside Insert",
        height=800,
    )

    # Save interactive HTML
    plots_dir: str = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    html_filename: str = os.path.join(plots_dir, "boxplot_coverage_outside_insert.html")
    fig.write_html(html_filename, include_plotlyjs=False, full_html=False)

    # Save static PNG if requested
    png_filename: str | None = None
    if static_report:
        png_filename = os.path.join(plots_dir, "boxplot_coverage_outside_insert.png")
        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import generate_coverage_boxplot_matplotlib

            generate_coverage_boxplot_matplotlib(coverage_pivot, png_filename)
        else:
            fig.write_image(png_filename)

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

    # Check if resistance genes are available
    available_rows = summary_df[summary_df["Category"] == "ResistanceGenesAvailable"]
    if available_rows.empty or (available_rows["Value"] == "false").all():
        # No genes available - return empty DataFrame with correct structure
        return pd.DataFrame(columns=["Sample", "Plasmid", "Gene Name", "Mean Depth", "Breadth (%)"])

    # Find all ResistanceGenes rows
    gene_rows = summary_df[summary_df["Category"] == "ResistanceGenes"]
    if gene_rows.empty:
        return pd.DataFrame(columns=["Sample", "Plasmid", "Gene Name", "Breadth (%)"])

    # Build flat table
    table_rows: list[dict[str, Any]] = []

    for _, row in gene_rows.iterrows():
        sample = row["Sample"]
        plasmid = row["Plasmid"]
        gene_names_str = row["Value"]

        if not isinstance(gene_names_str, str) or not gene_names_str:
            continue

        # Split gene names by comma
        gene_names = [g.strip() for g in gene_names_str.split(",")]

        for gene_name in gene_names:
            # Sanitize gene name for lookup (replace "/" and " " with "_")
            sanitized_name = gene_name.replace("/", "_").replace(" ", "_")

            # Look up depth and breadth
            depth_category = f"ResGene_{sanitized_name}_MeanDepth"
            breadth_category = f"ResGene_{sanitized_name}_Breadth"

            depth_row = summary_df[
                (summary_df["Sample"] == sample)
                & (summary_df["Plasmid"] == plasmid)
                & (summary_df["Category"] == depth_category)
            ]
            breadth_row = summary_df[
                (summary_df["Sample"] == sample)
                & (summary_df["Plasmid"] == plasmid)
                & (summary_df["Category"] == breadth_category)
            ]

            mean_depth = float(depth_row.iloc[0]["Value"]) if not depth_row.empty else 0.0
            breadth = float(breadth_row.iloc[0]["Value"]) if not breadth_row.empty else 0.0

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
) -> None:
    import pandas as pd
    from jinja2 import Environment, FileSystemLoader

    logging.info("Generating summary reports")

    # Apply sorting to each DataFrame
    combined_df = apply_sorting(combined_df, TABLE_SORTING.get("combined"))
    verdict_df = apply_sorting(verdict_df, TABLE_SORTING.get("verdict"))
    ratio_df = apply_sorting(ratio_df, TABLE_SORTING.get("ratio"))
    p_values_df = pd.read_csv(p_value_table_filename, sep="\t")
    p_values_df = apply_sorting(p_values_df, TABLE_SORTING.get("p_value"))

    env = Environment(loader=FileSystemLoader(str(get_resource_path(TEMPLATE_DIR))))
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

    # Convert sorted DataFrames to HTML
    combined_html = combined_df.to_html(classes="table table-striped table-bordered", index=False)
    verdict_html = verdict_df.to_html(classes="table table-striped table-bordered", index=False)
    ratio_html = ratio_df.to_html(classes="table table-striped table-bordered", index=False)
    p_value_html = p_values_df.to_html(classes="table table-striped table-bordered", index=False)

    # Convert resistance gene DataFrame to HTML if provided
    resistance_gene_html = ""
    if resistance_gene_df is not None and not resistance_gene_df.empty:
        resistance_gene_html = resistance_gene_df.to_html(
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
        logo_base64=logo_base64,
        threshold=threshold,
        unclear_range=unclear_range,
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
            logo_base64=logo_base64,
            threshold=threshold,
            unclear_range=unclear_range,
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
    plot_backend: str = "plotly",
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
        # Copy plotly.min.js to output_dir/assets/
        assets_dir = os.path.join(output_dir, "assets")
        os.makedirs(assets_dir, exist_ok=True)
        dest = os.path.join(assets_dir, "plotly.min.js")
        if not os.path.exists(dest):
            import shutil

            import plotly as _plotly

            plotly_dir = os.path.dirname(_plotly.__file__)
            source = os.path.join(plotly_dir, "package_data", "plotly.min.js")
            shutil.copy2(source, dest)
        plotly_js_path = "assets/plotly.min.js"  # Relative from output_dir

    if plotly_mode == "embedded":
        import plotly as _plotly

        plotly_js_source = os.path.join(
            os.path.dirname(_plotly.__file__), "package_data", "plotly.min.js"
        )
        with open(plotly_js_source) as f:
            plotly_js_inline = f"<script>{f.read()}</script>"

    _, summary_df = read_compare_outputs(input_dir, substring_to_remove=substring_to_remove)

    # Calculate variations
    boxplot_data, p_values_df = calculate_variations(summary_df[summary_df["Category"] == "Ratio"])

    # Generate plots
    boxplot_filename_interactive, boxplot_filename_png = plot_boxplot(
        boxplot_data, output_dir, static_report=static_report, plot_backend=plot_backend
    )
    heatmap_filename_interactive, heatmap_filename_png = plot_heatmap(
        summary_df[summary_df["Category"] == "Ratio"],
        output_dir,
        threshold,
        unclear_range,
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
    breadth_heatmap_interactive, breadth_heatmap_png = plot_metric_heatmap(
        breadth_data,
        output_dir,
        "Backbone Breadth",
        vmax=100.0,
        static_report=static_report,
        plot_backend=plot_backend,
    )

    # --- New: Resistance Gene Table (REPT-08) ---
    resistance_gene_df = create_resistance_gene_table(summary_df)

    # --- Build extended metrics DataFrame for export (REPT-10) ---
    extended_categories = [
        "CoverageOutsideINSERT",
        "MismatchesNearINSERT",
        "MeanDepthBackbone",
        "BreadthBackbone",
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
        threshold,
        unclear_range,
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
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Merge and generate summary reports for multiple samples and plasmids."
    )
    parser.add_argument(
        "-i", "--input_dir", help="Directory containing compare outputs", required=True
    )
    parser.add_argument(
        "-o", "--output_dir", help="Directory to save the plots and reports", required=True
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})",
    )
    parser.add_argument(
        "--substring_to_remove", help="Substring to remove from sample names", default=None
    )
    parser.add_argument(
        "--static-report",
        action="store_true",
        help="Generate PNG files and non-interactive HTML report (slower, requires kaleido)",
    )
    parser.add_argument(
        "--plotly-mode",
        choices=["directory", "cdn", "embedded"],
        default="directory",
        help="Plotly.js inclusion mode: directory (shared file), cdn (online), embedded (inline)",
    )
    add_logging_args(parser)

    args = parser.parse_args()

    configure_logging_from_args(args)

    main(
        args.input_dir,
        args.output_dir,
        args.threshold,
        unclear_range=UNCLEAR_RANGE,
        substring_to_remove=args.substring_to_remove,
        static_report=args.static_report,
        plotly_mode=args.plotly_mode,
    )
