from __future__ import annotations

import base64
import logging
import os
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd

from plasmicheck.config import get_config
from plasmicheck.resources import get_resource_path
from plasmicheck.version import __version__ as VERSION

from .plotting.colors import ASSIGNMENT_COLORS
from .utils import add_logging_args, configure_logging_from_args

_cfg = get_config()
DEFAULT_THRESHOLD: float = _cfg["default_threshold"]
UNCLEAR_RANGE: dict[str, float] = _cfg["unclear_range"]
PLOT_SAMPLE_REPORT: dict[str, Any] = _cfg["plot_sample_report"]
TEMPLATE_DIR: str = _cfg["paths"]["template_dir"]
LOGO_PATH: str = _cfg["paths"]["logo_path"]
PLOT_DIMENSIONS: dict[str, int] = PLOT_SAMPLE_REPORT.get("figsize", {"width": 1000, "height": 400})
DOWNSAMPLE_LIMIT: int = _cfg.get("downsample_limit", 5000)

# Setup logging for the libraries used in this script
logging.getLogger("jinja2").setLevel(logging.ERROR)
logging.getLogger("fontTools").setLevel(logging.ERROR)


def load_data(reads_assignment_file: str, summary_file: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    import pandas as pd

    logging.info(f"Loading data from {reads_assignment_file} and {summary_file}")
    reads_df = pd.read_csv(reads_assignment_file, sep="\t")
    summary_df = pd.read_csv(summary_file, sep="\t")
    return reads_df, summary_df


def downsample_data(df: pd.DataFrame, limit: int) -> tuple[pd.DataFrame, bool]:
    if len(df) > limit:
        logging.info(f"Downsampling data from {len(df)} to {limit} points.")
        df = df.sample(n=limit, random_state=1)  # Random downsample to the limit
        downsampled = True
    else:
        downsampled = False
    return df, downsampled


_CATEGORY_ORDER: list[str] = ["Plasmid", "Human", "Tied", "Backbone_Only", "Ambiguous"]

_PLOTLY_FONT = (
    '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif'
)

_PLOTLY_LAYOUT: dict[str, Any] = {
    "font": {"family": _PLOTLY_FONT, "size": 13, "color": "#2d3436"},
    "paper_bgcolor": "white",
    "plot_bgcolor": "white",
    "margin": {"l": 60, "r": 30, "t": 10, "b": 50, "pad": 4},
    "xaxis": {
        "showgrid": False,
        "showline": True,
        "linewidth": 1,
        "linecolor": "#e9ecef",
        "zeroline": False,
        "tickfont": {"size": 12, "color": "#636e72"},
        "title_font": {"size": 13, "color": "#636e72"},
    },
    "yaxis": {
        "showgrid": True,
        "gridwidth": 0.5,
        "gridcolor": "#f0f0f0",
        "showline": True,
        "linewidth": 1,
        "linecolor": "#e9ecef",
        "zeroline": False,
        "tickfont": {"size": 12, "color": "#636e72"},
        "title_font": {"size": 13, "color": "#636e72"},
    },
    "legend": {
        "bgcolor": "rgba(255,255,255,0.9)",
        "bordercolor": "#e9ecef",
        "borderwidth": 1,
        "font": {"size": 12},
    },
    "hoverlabel": {
        "bgcolor": "white",
        "bordercolor": "#dee2e6",
        "font": {"size": 12, "family": _PLOTLY_FONT},
    },
}

_PLOTLY_CONFIG: dict[str, Any] = {
    "displayModeBar": "hover",
    "displaylogo": False,
    "modeBarButtonsToRemove": ["pan2d", "lasso2d", "select2d", "autoScale2d"],
    "responsive": True,
    "toImageButtonOptions": {
        "format": "png",
        "filename": "plasmicheck_chart",
        "height": 800,
        "width": 1200,
        "scale": 2,
    },
}


def generate_plots(
    reads_df: pd.DataFrame,
    output_folder: str,
    static_report: bool = False,
    plot_backend: str = "plotly",
) -> tuple[str, str | None, str, str | None]:
    import plotly.express as px

    logging.info("Generating plots")
    width = PLOT_SAMPLE_REPORT.get("figsize", {}).get("width", 1000)
    height = PLOT_SAMPLE_REPORT.get("figsize", {}).get("height", 400)
    if not isinstance(width, int):
        width = 1000
    if not isinstance(height, int):
        height = 400

    plots_dir = os.path.join(output_folder, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # Consistent category ordering across all plots
    cat_order = {"AssignedTo": _CATEGORY_ORDER}

    # ── Box plot ──────────────────────────────────────────
    fig_box = px.box(
        reads_df,
        x="AssignedTo",
        y="PlasmidScore",
        points="all",
        color="AssignedTo",
        color_discrete_map=ASSIGNMENT_COLORS,
        category_orders=cat_order,
        labels={
            "PlasmidScore": PLOT_SAMPLE_REPORT["box_plot_y_label"],
            "AssignedTo": "",
        },
        width=width,
        height=height,
    )
    fig_box.update_layout(**_PLOTLY_LAYOUT, showlegend=False, hovermode="closest")
    fig_box.update_traces(
        marker={"size": 3, "opacity": 0.5, "line_width": 0},
        line={"width": 1.5},
    )

    boxplot_filename_interactive = os.path.join(
        plots_dir, PLOT_SAMPLE_REPORT["output_box_plot_filename"].replace(".png", ".html")
    )
    fig_box.write_html(
        boxplot_filename_interactive,
        include_plotlyjs=False,
        full_html=False,
        config=_PLOTLY_CONFIG,
    )

    # ── Scatter plot ──────────────────────────────────────
    fig_scatter = px.scatter(
        reads_df,
        x="PlasmidScore",
        y="HumanScore",
        color="AssignedTo",
        color_discrete_map=ASSIGNMENT_COLORS,
        category_orders=cat_order,
        labels={
            "PlasmidScore": PLOT_SAMPLE_REPORT["scatter_plot_x_label"],
            "HumanScore": PLOT_SAMPLE_REPORT["scatter_plot_y_label"],
        },
        width=width,
        height=height,
        render_mode="webgl",
    )
    fig_scatter.update_layout(**_PLOTLY_LAYOUT, showlegend=True, hovermode="closest")
    fig_scatter.update_traces(marker={"size": 4, "opacity": 0.5, "line_width": 0})

    # Diagonal reference line (equal scores)
    import plotly.graph_objects as go

    score_max = max(
        reads_df["PlasmidScore"].max(),
        reads_df["HumanScore"].max(),
        1,
    )
    fig_scatter.add_trace(
        go.Scattergl(
            x=[0, score_max],
            y=[0, score_max],
            mode="lines",
            line={"color": "#dee2e6", "width": 1, "dash": "dash"},
            showlegend=False,
            hoverinfo="skip",
        )
    )

    scatter_filename_interactive = os.path.join(
        plots_dir, PLOT_SAMPLE_REPORT["output_scatter_plot_filename"].replace(".png", ".html")
    )
    fig_scatter.write_html(
        scatter_filename_interactive,
        include_plotlyjs=False,
        full_html=False,
        config=_PLOTLY_CONFIG,
    )

    # Only generate PNGs when static_report=True
    boxplot_filename_png = None
    scatter_filename_png = None
    if static_report:
        boxplot_filename_png = os.path.join(
            plots_dir, PLOT_SAMPLE_REPORT["output_box_plot_filename"]
        )
        scatter_filename_png = os.path.join(
            plots_dir, PLOT_SAMPLE_REPORT["output_scatter_plot_filename"]
        )

        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import (
                generate_boxplot_matplotlib,
                generate_scatter_matplotlib,
            )

            # Generate plots using matplotlib
            boxplot_title = f"{PLOT_SAMPLE_REPORT['title_box_plot']} (Total Reads: {len(reads_df)})"
            scatter_title = (
                f"{PLOT_SAMPLE_REPORT['title_scatter_plot']} (Total Reads: {len(reads_df)})"
            )

            generate_boxplot_matplotlib(
                reads_df, boxplot_filename_png, boxplot_title, width, height
            )
            generate_scatter_matplotlib(
                reads_df, scatter_filename_png, scatter_title, width, height
            )
        else:
            # Default plotly/kaleido backend
            import kaleido  # type: ignore[import-untyped]

            kaleido.start_sync_server()

            fig_box.write_image(boxplot_filename_png)
            fig_scatter.write_image(scatter_filename_png)

    logging.info("Plots generated.")

    return (
        boxplot_filename_interactive,
        boxplot_filename_png,
        scatter_filename_interactive,
        scatter_filename_png,
    )


def _ensure_plotly_assets(output_root: str) -> None:
    """Copy plotly.min.js to {output_root}/assets/ if not already present."""
    import plotly

    plotly_dir = os.path.dirname(plotly.__file__)
    plotly_js_source = os.path.join(plotly_dir, "package_data", "plotly.min.js")

    assets_dir = os.path.join(output_root, "assets")
    os.makedirs(assets_dir, exist_ok=True)

    dest = os.path.join(assets_dir, "plotly.min.js")
    if not os.path.exists(dest):
        import shutil

        shutil.copy2(plotly_js_source, dest)
        logging.info(f"Copied plotly.min.js to {dest}")


def encode_image_to_base64(image_path: str) -> str:
    logging.info(f"Encoding image {image_path} to base64")
    with open(image_path, "rb") as image_file:
        encoded_string: str = base64.b64encode(image_file.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded_string}"


def extract_verdict_from_summary(summary_df: pd.DataFrame) -> str:
    verdict_row: pd.DataFrame = summary_df[summary_df["Category"] == "Verdict"]
    if not verdict_row.empty:
        return str(verdict_row["Count"].values[0])
    else:
        return "Verdict not found in summary file"


def _derive_verdict_info(verdict: str) -> tuple[str, str]:
    """Return (verdict_short, verdict_class) from the full verdict string."""
    if "not contaminated" in verdict.lower():
        return "Clean", "clean"
    elif "unclear" in verdict.lower():
        return "Unclear", "unclear"
    else:
        return "Contaminated", "contaminated"


def _compute_gauge_position(
    ratio: float,
    unclear_lower: float,
    unclear_upper: float,
) -> float:
    """Map ratio to 0-100 gauge position using piecewise linear scale.

    Zones: [0, unclear_lower] -> 0-40%,  [unclear_lower, unclear_upper] -> 40-60%,
    [unclear_upper, 10*unclear_upper] -> 60-95%, clamped at 97.
    """
    if ratio <= 0:
        return 0.0
    if ratio <= unclear_lower:
        return (ratio / unclear_lower) * 40.0
    if ratio <= unclear_upper:
        return 40.0 + ((ratio - unclear_lower) / (unclear_upper - unclear_lower)) * 20.0
    upper_max = 10.0 * unclear_upper
    if ratio <= upper_max:
        return 60.0 + ((ratio - unclear_upper) / (upper_max - unclear_upper)) * 35.0
    return 97.0


def _parse_mismatches(raw: str) -> tuple[int, int]:
    """Parse mismatches_near_insert from its raw string representation.

    Returns (with_mismatches, without_mismatches).
    """
    import ast

    if raw in ("N/A", ""):
        return 0, 0
    try:
        d = ast.literal_eval(raw) if isinstance(raw, str) else raw
        return (
            int(d.get("with_mismatches_or_clipping", 0)),
            int(d.get("without_mismatches_or_clipping", 0)),
        )
    except (ValueError, SyntaxError, AttributeError):
        return 0, 0


def generate_report(
    summary_df: pd.DataFrame,
    output_folder: str,
    verdict: str,
    ratio: float,
    threshold: float,
    unclear: dict[str, float],
    command_line: str,
    human_fasta: str = "None",
    plasmid_gb: str = "None",
    sequencing_file: str = "None",
    boxplot_filename_interactive: str = "",
    boxplot_filename_png: str | None = None,
    scatter_filename_interactive: str = "",
    scatter_filename_png: str | None = None,
    downsampled: bool = False,
    static_report: bool = False,
    plotly_mode: str = "directory",
    plotly_version: str = "",
    plotly_js_path: str = "",
    plotly_js_inline: str = "",
    plasmid_count: int = 0,
    human_count: int = 0,
    tied_count: int = 0,
    backbone_only_count: int = 0,
    ambiguous_count: int = 0,
    total_reads: int = 0,
    plasmid_pct: str = "0.0",
    human_pct: str = "0.0",
    tied_pct: str = "0.0",
    backbone_only_pct: str = "0.0",
    ambiguous_pct: str = "0.0",
    coverage_outside_insert: str = "N/A",
    mismatches_near_insert: str = "N/A",
    backbone_warning: str = "",
    verdict_short: str = "",
    verdict_class: str = "",
    gauge_position: float = 50.0,
    gauge_clean_end: float = 40.0,
    gauge_unclear_end: float = 60.0,
    gauge_threshold_pos: float = 50.0,
    mismatches_with: int = 0,
    mismatches_without: int = 0,
    sample_name: str = "",
    plasmid_name: str = "",
    coverage_metrics: dict[str, str] | None = None,
    coverage_fallback: bool = False,
) -> None:
    from jinja2 import Environment, FileSystemLoader

    logging.info("Generating report")
    env = Environment(loader=FileSystemLoader(str(get_resource_path(TEMPLATE_DIR))))
    template = env.get_template("report_template.html")

    logo_base64 = encode_image_to_base64(str(get_resource_path(LOGO_PATH)))

    # Prepare downsample message if data was downsampled
    downsample_message = (
        f"Plots show a random subset of {DOWNSAMPLE_LIMIT:,} reads for rendering performance."
        " The verdict and all statistics are computed from the complete dataset."
        if downsampled
        else ""
    )

    # Read the interactive HTML plot content
    with open(boxplot_filename_interactive) as f:
        boxplot_html_content = f.read()

    with open(scatter_filename_interactive) as f:
        scatter_html_content = f.read()

    # Render the interactive HTML report
    html_content_interactive = template.render(
        box_plot=boxplot_html_content,
        scatter_plot=scatter_html_content,
        verdict=verdict,
        verdict_short=verdict_short,
        verdict_class=verdict_class,
        ratio=f"{ratio:.3f}",
        threshold=threshold,
        unclear_range=unclear,
        version=VERSION,
        logo_base64=logo_base64,
        human_fasta=human_fasta,
        plasmid_gb=plasmid_gb,
        sequencing_file=sequencing_file,
        run_date=datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line,
        interactive=True,
        downsample_message=downsample_message,
        plotly_mode=plotly_mode,
        plotly_version=plotly_version,
        plotly_js_path=plotly_js_path,
        plotly_js_inline=plotly_js_inline,
        plasmid_count=plasmid_count,
        human_count=human_count,
        tied_count=tied_count,
        backbone_only_count=backbone_only_count,
        ambiguous_count=ambiguous_count,
        total_reads=total_reads,
        plasmid_pct=plasmid_pct,
        human_pct=human_pct,
        tied_pct=tied_pct,
        backbone_only_pct=backbone_only_pct,
        ambiguous_pct=ambiguous_pct,
        assignment_colors=ASSIGNMENT_COLORS,
        coverage_outside_insert=coverage_outside_insert,
        mismatches_with=mismatches_with,
        mismatches_without=mismatches_without,
        backbone_warning=backbone_warning,
        gauge_position=gauge_position,
        gauge_clean_end=gauge_clean_end,
        gauge_unclear_end=gauge_unclear_end,
        gauge_threshold_pos=gauge_threshold_pos,
        sample_name=sample_name,
        plasmid_name=plasmid_name,
        coverage_metrics=coverage_metrics or {},
        coverage_fallback=coverage_fallback,
    )

    # Save interactive HTML report
    html_report_interactive = os.path.join(output_folder, "report_interactive.html")
    with open(html_report_interactive, "w") as f:
        f.write(html_content_interactive)

    # Only generate non-interactive HTML report when static_report=True
    if static_report and boxplot_filename_png and scatter_filename_png:
        # Convert the static PNG files to Base64
        boxplot_png_base64 = encode_image_to_base64(boxplot_filename_png)
        scatter_png_base64 = encode_image_to_base64(scatter_filename_png)

        # Render the non-interactive HTML report
        html_content_non_interactive = template.render(
            box_plot=boxplot_png_base64,
            scatter_plot=scatter_png_base64,
            verdict=verdict,
            verdict_short=verdict_short,
            verdict_class=verdict_class,
            ratio=f"{ratio:.3f}",
            threshold=threshold,
            unclear_range=unclear,
            version=VERSION,
            logo_base64=logo_base64,
            human_fasta=human_fasta,
            plasmid_gb=plasmid_gb,
            sequencing_file=sequencing_file,
            run_date=datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
            command_line=command_line,
            interactive=False,
            downsample_message=downsample_message,
            plotly_mode=plotly_mode,
            plotly_version=plotly_version,
            plotly_js_path=plotly_js_path,
            plotly_js_inline=plotly_js_inline,
            plasmid_count=plasmid_count,
            human_count=human_count,
            tied_count=tied_count,
            backbone_only_count=backbone_only_count,
            ambiguous_count=ambiguous_count,
            total_reads=total_reads,
            plasmid_pct=plasmid_pct,
            human_pct=human_pct,
            tied_pct=tied_pct,
            backbone_only_pct=backbone_only_pct,
            ambiguous_pct=ambiguous_pct,
            assignment_colors=ASSIGNMENT_COLORS,
            coverage_outside_insert=coverage_outside_insert,
            mismatches_with=mismatches_with,
            mismatches_without=mismatches_without,
            backbone_warning=backbone_warning,
            gauge_position=gauge_position,
            gauge_clean_end=gauge_clean_end,
            gauge_unclear_end=gauge_unclear_end,
            gauge_threshold_pos=gauge_threshold_pos,
            sample_name=sample_name,
            plasmid_name=plasmid_name,
            coverage_metrics=coverage_metrics or {},
            coverage_fallback=coverage_fallback,
        )

        # Save non-interactive HTML report
        html_report_non_interactive = os.path.join(output_folder, "report_non_interactive.html")
        with open(html_report_non_interactive, "w") as f:
            f.write(html_content_non_interactive)


def main(
    reads_assignment_file: str,
    summary_file: str,
    output_folder: str,
    threshold: float = DEFAULT_THRESHOLD,
    unclear: dict[str, float] = UNCLEAR_RANGE,
    human_fasta: str = "None",
    plasmid_gb: str = "None",
    sequencing_file: str = "None",
    command_line: str = "",
    static_report: bool = False,
    plotly_mode: str = "directory",
    plot_backend: str = "plotly",
    output_root: str | None = None,
) -> None:
    reads_df, summary_df = load_data(reads_assignment_file, summary_file)

    # Downsample data if necessary
    reads_df, downsampled = downsample_data(reads_df, DOWNSAMPLE_LIMIT)

    # Generate plots
    (
        boxplot_filename_interactive,
        boxplot_filename_png,
        scatter_filename_interactive,
        scatter_filename_png,
    ) = generate_plots(
        reads_df, output_folder, static_report=static_report, plot_backend=plot_backend
    )

    # Extract the verdict from the summary file
    verdict = extract_verdict_from_summary(summary_df)

    # Extract category counts from summary_df
    def _get_count(df: pd.DataFrame, category: str) -> int:
        rows = df[df["Category"] == category]
        return int(rows["Count"].values[0]) if not rows.empty else 0

    plasmid_count = _get_count(summary_df, "Plasmid")
    human_count = _get_count(summary_df, "Human")
    tied_count = _get_count(summary_df, "Tied")
    backbone_only_count = _get_count(summary_df, "Backbone_Only")
    ambiguous_count = _get_count(summary_df, "Ambiguous")
    total_reads = plasmid_count + human_count + tied_count + backbone_only_count + ambiguous_count

    # Calculate percentages
    def _pct(count: int, total: int) -> str:
        return f"{count / total * 100:.1f}" if total > 0 else "0.0"

    # Extract additional metrics
    cov_rows = summary_df[summary_df["Category"] == "CoverageOutsideINSERT"]
    coverage_outside_insert = str(cov_rows["Count"].values[0]) if not cov_rows.empty else "N/A"
    mis_rows = summary_df[summary_df["Category"] == "MismatchesNearINSERT"]
    mismatches_near_insert = str(mis_rows["Count"].values[0]) if not mis_rows.empty else "N/A"

    # Extract coverage metrics (Phase 9)
    def _get_metric(df: pd.DataFrame, category: str, default: str = "0.00") -> str:
        rows = df[df["Category"] == category]
        return str(rows["Count"].values[0]) if not rows.empty else default

    coverage_metrics = {
        "mean_depth_insert": _get_metric(summary_df, "MeanDepthInsert"),
        "median_depth_insert": _get_metric(summary_df, "MedianDepthInsert"),
        "breadth_insert": _get_metric(summary_df, "BreadthInsert"),
        "breadth_insert_5x": _get_metric(summary_df, "BreadthInsert_5x"),
        "cv_insert": _get_metric(summary_df, "CoverageCV_Insert"),
        "mean_depth_backbone": _get_metric(summary_df, "MeanDepthBackbone"),
        "median_depth_backbone": _get_metric(summary_df, "MedianDepthBackbone"),
        "breadth_backbone": _get_metric(summary_df, "BreadthBackbone"),
        "breadth_backbone_5x": _get_metric(summary_df, "BreadthBackbone_5x"),
        "cv_backbone": _get_metric(summary_df, "CoverageCV_Backbone"),
    }

    # Detect coverage fallback mode
    coverage_fallback_rows = summary_df[summary_df["Category"] == "CoverageFallback"]
    coverage_fallback = (
        not coverage_fallback_rows.empty
        and str(coverage_fallback_rows["Count"].values[0]).lower() == "true"
    )

    # Determine if backbone filtering was unavailable
    backbone_warning = ""

    ratio = plasmid_count / human_count if human_count != 0 else float("inf")

    # Derive new template variables
    verdict_short, verdict_class = _derive_verdict_info(verdict)
    unclear_lower = float(UNCLEAR_RANGE.get("lower_bound", 0.6))
    unclear_upper = float(UNCLEAR_RANGE.get("upper_bound", 1.0))
    gauge_position = _compute_gauge_position(ratio, unclear_lower, unclear_upper)
    gauge_clean_end = _compute_gauge_position(unclear_lower, unclear_lower, unclear_upper)
    gauge_unclear_end = _compute_gauge_position(unclear_upper, unclear_lower, unclear_upper)
    gauge_threshold_pos = _compute_gauge_position(threshold, unclear_lower, unclear_upper)
    mismatches_with, mismatches_without = _parse_mismatches(mismatches_near_insert)

    # Derive sample and plasmid display names from file paths
    sample_name = (
        os.path.splitext(os.path.basename(sequencing_file))[0]
        if sequencing_file not in ("None", "")
        else ""
    )
    plasmid_name = (
        os.path.splitext(os.path.basename(plasmid_gb))[0] if plasmid_gb not in ("None", "") else ""
    )

    # Determine output_root for shared assets
    effective_root = output_root if output_root else output_folder

    plotly_js_path = ""
    plotly_version = ""
    plotly_js_inline = ""

    if plotly_mode in ("directory", "cdn"):
        import plotly

        plotly_version = plotly.__version__
    if plotly_mode == "directory":
        _ensure_plotly_assets(effective_root)
        # Compute relative path from output_folder to assets/plotly.min.js
        assets_js = os.path.join(effective_root, "assets", "plotly.min.js")
        plotly_js_path = os.path.relpath(assets_js, output_folder)
    elif plotly_mode == "embedded":
        import plotly as _plotly

        plotly_js_source = os.path.join(
            os.path.dirname(_plotly.__file__), "package_data", "plotly.min.js"
        )
        with open(plotly_js_source) as f:
            plotly_js_inline = f"<script>{f.read()}</script>"

    # Generate report
    generate_report(
        summary_df,
        output_folder,
        verdict,
        ratio,
        threshold,
        UNCLEAR_RANGE,
        command_line,
        human_fasta,
        plasmid_gb,
        sequencing_file,
        boxplot_filename_interactive=boxplot_filename_interactive,
        boxplot_filename_png=boxplot_filename_png,
        scatter_filename_interactive=scatter_filename_interactive,
        scatter_filename_png=scatter_filename_png,
        downsampled=downsampled,
        static_report=static_report,
        plotly_mode=plotly_mode,
        plotly_version=plotly_version,
        plotly_js_path=plotly_js_path,
        plotly_js_inline=plotly_js_inline,
        plasmid_count=plasmid_count,
        human_count=human_count,
        tied_count=tied_count,
        backbone_only_count=backbone_only_count,
        ambiguous_count=ambiguous_count,
        total_reads=total_reads,
        plasmid_pct=_pct(plasmid_count, total_reads),
        human_pct=_pct(human_count, total_reads),
        tied_pct=_pct(tied_count, total_reads),
        backbone_only_pct=_pct(backbone_only_count, total_reads),
        ambiguous_pct=_pct(ambiguous_count, total_reads),
        coverage_outside_insert=coverage_outside_insert,
        mismatches_near_insert=mismatches_near_insert,
        backbone_warning=backbone_warning,
        verdict_short=verdict_short,
        verdict_class=verdict_class,
        gauge_position=gauge_position,
        gauge_clean_end=gauge_clean_end,
        gauge_unclear_end=gauge_unclear_end,
        gauge_threshold_pos=gauge_threshold_pos,
        mismatches_with=mismatches_with,
        mismatches_without=mismatches_without,
        sample_name=sample_name,
        plasmid_name=plasmid_name,
        coverage_metrics=coverage_metrics,
        coverage_fallback=coverage_fallback,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate a visualized HTML report from alignment comparison results"
    )
    parser.add_argument(
        "-r",
        "--reads_assignment_file",
        help="Reads assignment file (reads_assignment.tsv)",
        required=True,
    )
    parser.add_argument("-s", "--summary_file", help="Summary file (summary.tsv)", required=True)
    parser.add_argument(
        "-o", "--output_folder", help="Folder to write the report and plots", required=True
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})",
    )
    parser.add_argument("-hf", "--human_fasta", default="None", help="Human reference FASTA file")
    parser.add_argument("-pg", "--plasmid_gb", default="None", help="GenBank plasmid file")
    parser.add_argument(
        "-sf",
        "--sequencing_file",
        default="None",
        help="Sequencing file (BAM, interleaved FASTQ, or first FASTQ file for paired FASTQ)",
    )
    parser.add_argument(
        "--static-report",
        action="store_true",
        help="Generate static PNG plots and non-interactive HTML report (slow)",
    )
    parser.add_argument(
        "--plotly-mode",
        choices=["directory", "cdn", "embedded"],
        default="directory",
        help="Plotly.js inclusion mode: directory (local file), cdn (CDN), or embedded (inline)",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help="Root directory for shared assets (plotly.min.js). Defaults to output_folder.",
    )
    add_logging_args(parser)

    import sys

    args = parser.parse_args()
    command_line = " ".join(sys.argv)

    configure_logging_from_args(args)

    main(
        args.reads_assignment_file,
        args.summary_file,
        args.output_folder,
        args.threshold,
        UNCLEAR_RANGE,
        args.human_fasta,
        args.plasmid_gb,
        args.sequencing_file,
        command_line,
        args.static_report,
        args.plotly_mode,
        args.output_root,
    )
