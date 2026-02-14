"""Matplotlib backend for static PNG plot generation."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd

from .colors import ASSIGNMENT_COLORS

logger = logging.getLogger(__name__)


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
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

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
        plt.close(fig)
        return

    sns.set_theme(style="whitegrid")

    fig, ax = plt.subplots(figsize=(width / 100, height / 100))

    # Create boxplot with custom colors
    categories = ["Plasmid", "Human", "Tied"]
    colors = [ASSIGNMENT_COLORS.get(cat, "#999999") for cat in categories]

    # Prepare data for boxplot
    data_list = [
        reads_df[reads_df["AssignedTo"] == cat]["PlasmidScore"].tolist()
        for cat in categories
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
    plt.close(fig)
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
    import matplotlib.pyplot as plt
    import seaborn as sns

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
        plt.close(fig)
        return

    sns.set_theme(style="whitegrid")

    fig, ax = plt.subplots(figsize=(width / 100, height / 100))

    # Plot each category with its color
    categories = ["Plasmid", "Human", "Tied"]
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
    plt.close(fig)
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
    import matplotlib.pyplot as plt
    import seaborn as sns

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
        plt.close(fig)
        return

    sns.set_theme(style="whitegrid")

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
    plt.close(fig)
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
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    if boxplot_data.empty:
        logger.warning(
            "Empty DataFrame provided to generate_summary_boxplot_matplotlib"
        )
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
        plt.close(fig)
        return

    sns.set_theme(style="whitegrid")

    # Melt the data for boxplot
    melted = boxplot_data.melt(var_name="Plasmid", value_name="Ratio")

    # Apply log10 transformation
    melted["LogRatio"] = np.log10(melted["Ratio"] + log_offset)

    # Height=800px equivalent (~8 inches at 100dpi)
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create boxplot by Plasmid
    plasmids = melted["Plasmid"].unique().tolist()
    data_list = [
        melted[melted["Plasmid"] == p]["LogRatio"].tolist() for p in plasmids
    ]

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
    plt.close(fig)
    logger.info(f"Saved summary boxplot to {output_path}")
