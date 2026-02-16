"""Matplotlib backend for static PNG plot generation."""

from __future__ import annotations

import gc
import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # Headless backend — must be set before importing pyplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

from .colors import ASSIGNMENT_COLORS

logger = logging.getLogger(__name__)

# Set seaborn theme once at module level (repeated calls accumulate state)
sns.set_theme(style="whitegrid")


def _close_and_collect(fig: matplotlib.figure.Figure) -> None:
    """Close a matplotlib figure and run garbage collection to free memory."""
    plt.close(fig)
    gc.collect()


def generate_boxplot_matplotlib(
    reads_df: pd.DataFrame,
    output_path: Path | str,
    title: str,
    width: int = 900,
    height: int = 600,
) -> None:
    """Generate single-sample score distribution boxplot using matplotlib.

    Args:
        reads_df: DataFrame with columns [ReadID, AssignedTo, PlasmidScore, HumanScore]
        output_path: Output PNG file path
        title: Plot title (typically includes total read count)
        width: Plot width in pixels (default: 900)
        height: Plot height in pixels (default: 600)
    """
    if reads_df.empty:
        logger.warning("Empty DataFrame provided to generate_boxplot_matplotlib")
        fig, ax = plt.subplots(figsize=(width / 100, height / 100))
        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="center",
            va="center",
            fontsize=12,
        )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        _close_and_collect(fig)
        return

    fig, ax = plt.subplots(figsize=(width / 100, height / 100))

    # Create boxplot with custom colors
    categories = ["Plasmid", "Human", "Tied", "Backbone_Only", "Ambiguous"]
    colors = [ASSIGNMENT_COLORS.get(cat, "#999999") for cat in categories]

    # Prepare data for boxplot
    data_list = [
        reads_df[reads_df["AssignedTo"] == cat]["PlasmidScore"].tolist() for cat in categories
    ]

    bp = ax.boxplot(
        data_list,
        patch_artist=True,
        widths=0.6,
    )

    # Set x-axis labels
    ax.set_xticks(range(1, len(categories) + 1))
    ax.set_xticklabels(categories)

    # Apply colors to boxes
    for patch, color in zip(bp["boxes"], colors, strict=True):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Add jittered scatter points matching Plotly's points="all"
    for i, cat in enumerate(categories):
        cat_data = reads_df[reads_df["AssignedTo"] == cat]["PlasmidScore"].tolist()
        if len(cat_data) > 0:
            # Add jitter to x-coordinates
            x = np.random.normal(i + 1, 0.04, size=len(cat_data))
            ax.scatter(
                x,
                cat_data,
                alpha=0.3,
                s=2,
                color=colors[i],
            )

    ax.set_xlabel("AssignedTo")
    ax.set_ylabel("PlasmidScore")
    ax.set_title(title)

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved boxplot to {output_path}")


def generate_scatter_matplotlib(
    reads_df: pd.DataFrame,
    output_path: Path | str,
    title: str,
    width: int = 900,
    height: int = 600,
) -> None:
    """Generate single-sample score scatter plot using matplotlib.

    Args:
        reads_df: DataFrame with columns [ReadID, AssignedTo, PlasmidScore, HumanScore]
        output_path: Output PNG file path
        title: Plot title (typically includes total read count)
        width: Plot width in pixels (default: 900)
        height: Plot height in pixels (default: 600)
    """
    if reads_df.empty:
        logger.warning("Empty DataFrame provided to generate_scatter_matplotlib")
        fig, ax = plt.subplots(figsize=(width / 100, height / 100))
        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="center",
            va="center",
            fontsize=12,
        )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        _close_and_collect(fig)
        return

    fig, ax = plt.subplots(figsize=(width / 100, height / 100))

    # Plot each category with its color
    categories = ["Plasmid", "Human", "Tied", "Backbone_Only", "Ambiguous"]
    for cat in categories:
        cat_data = reads_df[reads_df["AssignedTo"] == cat]
        if not cat_data.empty:
            ax.scatter(
                cat_data["PlasmidScore"],
                cat_data["HumanScore"],
                alpha=0.3,
                s=2,
                color=ASSIGNMENT_COLORS.get(cat, "#999999"),
                label=cat,
            )

    ax.set_xlabel("PlasmidScore")
    ax.set_ylabel("HumanScore")
    ax.set_title(title)
    ax.legend()

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved scatter plot to {output_path}")


def generate_heatmap_matplotlib(
    ratio_data: pd.DataFrame,
    output_path: Path | str,
    threshold: float,
    plot_config: dict[str, Any],
) -> None:
    """Generate summary contamination ratio heatmap using matplotlib.

    Args:
        ratio_data: Pivoted DataFrame (Sample x Plasmid) with contamination ratios
        output_path: Output PNG file path
        threshold: Contamination threshold (not directly used, kept for API consistency)
        plot_config: Plot configuration dict (kept for API consistency)
    """
    if ratio_data.empty:
        logger.warning("Empty DataFrame provided to generate_heatmap_matplotlib")
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="center",
            va="center",
            fontsize=12,
        )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        _close_and_collect(fig)
        return

    # Ensure numeric dtype for heatmap rendering
    ratio_data = ratio_data.apply(pd.to_numeric, errors="coerce")

    # Determine figure size based on data dimensions
    n_samples, n_plasmids = ratio_data.shape
    fig_width = max(10, n_plasmids * 0.8)
    fig_height = max(6, n_samples * 0.4)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Create heatmap
    sns.heatmap(
        ratio_data,
        cmap="viridis",
        annot=False,
        cbar=False,
        ax=ax,
    )

    # Rotate x-axis labels to match Plotly layout
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved heatmap to {output_path}")


def generate_summary_boxplot_matplotlib(
    boxplot_data: pd.DataFrame,
    output_path: Path | str,
    log_offset: float = 0.001,
    marker_style: str = "outliers",
) -> None:
    """Generate summary contamination ratio boxplot (log scale) using matplotlib.

    Args:
        boxplot_data: Pivoted DataFrame (Sample x Plasmid) with contamination ratios
        output_path: Output PNG file path
        log_offset: Offset for log10 transformation to handle zeros (default: 0.001)
        marker_style: Marker style for outliers (kept for API consistency)
    """
    if boxplot_data.empty:
        logger.warning("Empty DataFrame provided to generate_summary_boxplot_matplotlib")
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="center",
            va="center",
            fontsize=12,
        )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        _close_and_collect(fig)
        return

    # Melt the data for boxplot
    melted = boxplot_data.melt(var_name="Plasmid", value_name="Ratio")
    melted["Ratio"] = pd.to_numeric(melted["Ratio"], errors="coerce")

    # Apply log10 transformation
    melted["LogRatio"] = np.log10(melted["Ratio"] + log_offset)

    # Height=800px equivalent (~8 inches at 100dpi)
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create boxplot by Plasmid
    plasmids = melted["Plasmid"].unique().tolist()
    data_list = [melted[melted["Plasmid"] == p]["LogRatio"].tolist() for p in plasmids]

    ax.boxplot(
        data_list,
        patch_artist=True,
    )

    # Set x-axis labels
    ax.set_xticks(range(1, len(plasmids) + 1))
    ax.set_xticklabels(plasmids)

    ax.set_xlabel("Plasmid")
    ax.set_ylabel(f"log10(Contamination Ratio + {log_offset})")
    ax.set_title("Contamination Ratio Distribution by Plasmid")

    # Rotate x-axis labels if many plasmids
    if len(plasmids) > 5:
        ax.set_xticklabels(plasmids, rotation=45, ha="right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved summary boxplot to {output_path}")


def generate_metric_heatmap_matplotlib(
    metric_data: pd.DataFrame,
    output_path: Path | str,
    metric_name: str,
    vmax: float | None = None,
) -> None:
    """Generate heatmap for any metric using matplotlib.

    Args:
        metric_data: Pivoted DataFrame (Sample x Plasmid) with metric values
        output_path: Output PNG file path
        metric_name: Display name for the metric (e.g., "Mean Depth Insert")
        vmax: Optional maximum value for color scale (data will be clipped)
    """
    if metric_data.empty:
        logger.warning(
            f"Empty DataFrame provided to generate_metric_heatmap_matplotlib for {metric_name}"
        )
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="center",
            va="center",
            fontsize=12,
        )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        _close_and_collect(fig)
        return

    # Ensure numeric dtype for heatmap rendering
    metric_data = metric_data.apply(pd.to_numeric, errors="coerce")

    # Clip data if vmax provided
    plot_data = metric_data.clip(upper=vmax) if vmax is not None else metric_data

    # Create custom colormap from lightblue to red
    cmap = LinearSegmentedColormap.from_list("coverage", ["lightblue", "red"], N=100)

    # Determine figure size based on data dimensions
    n_samples, n_plasmids = metric_data.shape
    figwidth = max(10, n_plasmids * 0.8)
    figheight = max(6, n_samples * 0.4)

    fig, ax = plt.subplots(figsize=(figwidth, figheight))

    # Generate heatmap
    sns.heatmap(
        plot_data,
        annot=False,
        cbar=False,
        cmap=cmap,
        vmin=0,
        vmax=vmax,
        mask=metric_data.isna(),
        ax=ax,
    )

    ax.set_xlabel("Plasmid")
    ax.set_ylabel("Sample")
    ax.set_title(metric_name)

    # Rotate x labels
    plt.xticks(rotation=45, ha="right")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved metric heatmap to {output_path}")


def generate_coverage_boxplot_matplotlib(
    boxplot_data: pd.DataFrame,
    output_path: Path | str,
) -> None:
    """Generate boxplot for CoverageOutsideINSERT distribution using matplotlib.

    Args:
        boxplot_data: Pivoted DataFrame (Sample x Plasmid) with coverage values
        output_path: Output PNG file path
    """
    if boxplot_data.empty:
        logger.warning("Empty DataFrame provided to generate_coverage_boxplot_matplotlib")
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="center",
            va="center",
            fontsize=12,
        )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        _close_and_collect(fig)
        return

    # Melt the data for boxplot
    melted = boxplot_data.melt(var_name="Plasmid", value_name="Value")
    melted["Value"] = pd.to_numeric(melted["Value"], errors="coerce")

    fig, ax = plt.subplots(figsize=(10, 8))

    # Create boxplot grouped by Plasmid
    boxplot_data_grouped = melted.groupby("Plasmid")["Value"].apply(list)
    plasmid_labels = list(boxplot_data_grouped.index)
    boxplot_values = [list(vals) for vals in boxplot_data_grouped.values]

    ax.boxplot(
        boxplot_values,
        patch_artist=True,
    )

    ax.set_xlabel("Plasmid")
    ax.set_ylabel("Coverage Outside Insert")
    ax.set_title("CoverageOutsideINSERT Distribution by Plasmid")

    # Set x-axis labels
    ax.set_xticklabels(plasmid_labels)

    # Rotate x labels if many plasmids
    if len(plasmid_labels) > 5:
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved coverage boxplot to {output_path}")


def generate_verdict_bar_matplotlib(
    aggregated_df: pd.DataFrame,
    output_path: Path | str,
    title: str,
    xlabel: str,
) -> None:
    """Generate stacked bar chart for verdict categories using matplotlib.

    Args:
        aggregated_df: DataFrame with verdict categories as columns, grouped items as index
        output_path: Output PNG file path
        title: Plot title
        xlabel: X-axis label (e.g., "Plasmid" or "Sample")
    """
    from .colors import VERDICT_COLORS

    if aggregated_df.empty:
        logger.warning("Empty DataFrame provided to generate_verdict_bar_matplotlib")
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="center",
            va="center",
            fontsize=12,
        )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        _close_and_collect(fig)
        return

    # Define verdict categories in display order
    verdict_categories = ["Clean", "Unclear", "Contaminated", "No Result"]

    # Filter to only categories present in the data
    available_categories = [cat for cat in verdict_categories if cat in aggregated_df.columns]

    # Determine figure width based on number of bars
    figwidth = max(10, len(aggregated_df) * 0.5)
    fig, ax = plt.subplots(figsize=(figwidth, 6))

    # Get x positions for bars
    x_positions = np.arange(len(aggregated_df))
    bar_width = 0.8

    # Iteratively stack bars
    bottom = np.zeros(len(aggregated_df))

    for category in available_categories:
        counts = aggregated_df[category].to_numpy()
        color = VERDICT_COLORS[category]

        ax.bar(
            x_positions,
            counts,
            bar_width,
            bottom=bottom,
            label=category,
            color=color,
            alpha=0.9,
        )

        # Update bottom for next stack layer
        bottom = bottom + counts

    # Set labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_title(title)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(aggregated_df.index)

    # Rotate x-axis labels if more than 5 bars
    if len(aggregated_df) > 5:
        plt.xticks(rotation=45, ha="right")

    # Add legend
    ax.legend()

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved verdict bar chart to {output_path}")
