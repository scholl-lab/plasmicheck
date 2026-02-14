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


def generate_plots(
    reads_df: pd.DataFrame,
    output_folder: str,
    static_report: bool = False,
) -> tuple[str, str | None, str, str | None]:
    import plotly.express as px

    logging.info("Generating plots")
    # Ensure default values for width and height
    width = PLOT_SAMPLE_REPORT.get("figsize", {}).get("width", 1000)  # Convert inches to pixels
    height = PLOT_SAMPLE_REPORT.get("figsize", {}).get("height", 400)  # Convert inches to pixels

    # Validate that width and height are integers
    if not isinstance(width, int):
        logging.warning(f"Invalid width value: {width}. Defaulting to 1000.")
        width = 1000
    if not isinstance(height, int):
        logging.warning(f"Invalid height value: {height}. Defaulting to 400.")
        height = 400

    # Create directory for plots if not exist
    plots_dir = os.path.join(output_folder, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # Box plot using Plotly
    boxplot_df = reads_df.copy()
    fig_box = px.box(
        boxplot_df,
        x="AssignedTo",
        y="PlasmidScore",
        points="all",
        color="AssignedTo",
        title=f"{PLOT_SAMPLE_REPORT['title_box_plot']}<br>(Total Reads: {len(reads_df)})",
        labels={
            "PlasmidScore": PLOT_SAMPLE_REPORT["box_plot_y_label"],
            "AssignedTo": PLOT_SAMPLE_REPORT["box_plot_x_label"],
        },
        width=width,
        height=height,
    )

    boxplot_filename_interactive = os.path.join(
        plots_dir, PLOT_SAMPLE_REPORT["output_box_plot_filename"].replace(".png", ".html")
    )
    fig_box.write_html(boxplot_filename_interactive, include_plotlyjs=False, full_html=False)

    # Scatter plot using Plotly
    fig_scatter = px.scatter(
        reads_df,
        x="PlasmidScore",
        y="HumanScore",
        color="AssignedTo",
        title=f"{PLOT_SAMPLE_REPORT['title_scatter_plot']}<br>(Total Reads: {len(reads_df)})",
        labels={
            "PlasmidScore": PLOT_SAMPLE_REPORT["scatter_plot_x_label"],
            "HumanScore": PLOT_SAMPLE_REPORT["scatter_plot_y_label"],
        },
        width=width,
        height=height,
        render_mode="webgl",  # Use WebGL for better performance with large datasets
    )

    scatter_filename_interactive = os.path.join(
        plots_dir, PLOT_SAMPLE_REPORT["output_scatter_plot_filename"].replace(".png", ".html")
    )
    fig_scatter.write_html(scatter_filename_interactive, include_plotlyjs=False, full_html=False)

    # Only generate PNGs when static_report=True
    boxplot_filename_png = None
    scatter_filename_png = None
    if static_report:
        import kaleido

        kaleido.start_sync_server()

        boxplot_filename_png = os.path.join(plots_dir, PLOT_SAMPLE_REPORT["output_box_plot_filename"])
        fig_box.write_image(boxplot_filename_png)

        scatter_filename_png = os.path.join(
            plots_dir, PLOT_SAMPLE_REPORT["output_scatter_plot_filename"]
        )
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
) -> None:
    from jinja2 import Environment, FileSystemLoader

    logging.info("Generating report")
    env = Environment(loader=FileSystemLoader(str(get_resource_path(TEMPLATE_DIR))))
    template = env.get_template("report_template.html")

    logo_base64 = encode_image_to_base64(str(get_resource_path(LOGO_PATH)))

    verdict_color = (
        "green"
        if "Sample is not contaminated with plasmid DNA" in verdict
        else "orange"
        if "Sample contamination status is unclear" in verdict
        else "red"
    )

    # Prepare downsample message if data was downsampled
    downsample_message = (
        "Data was downsampled to improve performance. Only a subset of reads is displayed."
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
        summary_df=summary_df.to_html(classes="table table-striped"),
        box_plot=boxplot_html_content,
        scatter_plot=scatter_html_content,
        verdict=verdict,
        ratio=f"{ratio:.3f}",
        threshold=threshold,
        unclear_range=unclear,
        verdict_color=verdict_color,
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
            summary_df=summary_df.to_html(classes="table table-striped"),
            box_plot=boxplot_png_base64,
            scatter_plot=scatter_png_base64,
            verdict=verdict,
            ratio=f"{ratio:.3f}",
            threshold=threshold,
            unclear_range=unclear,
            verdict_color=verdict_color,
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
    ) = generate_plots(reads_df, output_folder, static_report=static_report)

    # Extract the verdict from the summary file
    verdict = extract_verdict_from_summary(summary_df)

    plasmid_count = int(summary_df[summary_df["Category"] == "Plasmid"]["Count"].values[0])
    human_count = int(summary_df[summary_df["Category"] == "Human"]["Count"].values[0])
    ratio = plasmid_count / human_count if human_count != 0 else float("inf")

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
