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


def _style_axes(ax: matplotlib.axes.Axes) -> None:
    """Apply clean publication-grade spine and grid styling to an Axes."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#cbd5e1")
    ax.spines["bottom"].set_color("#cbd5e1")
    ax.tick_params(colors="#334155", labelsize=9.5)
    ax.xaxis.label.set_color("#1e293b")
    ax.yaxis.label.set_color("#1e293b")
    ax.title.set_color("#0f172a")


def generate_boxplot_matplotlib(
    reads_df: pd.DataFrame,
    output_path: Path | str,
    title: str | None = None,
    width: int = 1100,
    height: int = 440,
    total_reads: int | None = None,
) -> None:
    """Generate single-sample score distribution boxplot using matplotlib."""
    n_total = total_reads if total_reads is not None else len(reads_df)
    if title is None:
        title = f"Plasmid Alignment Score by Read Assignment (Total Reads: {n_total:,})"
    elif total_reads is not None and "(Total Reads:" not in title:
        title = f"{title} (Total Reads: {n_total:,})"

    if reads_df.empty:
        logger.warning("Empty DataFrame provided to generate_boxplot_matplotlib")
        fig, ax = plt.subplots(figsize=(width / 100, height / 100))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=220, bbox_inches="tight")
        _close_and_collect(fig)
        return

    fig, ax = plt.subplots(figsize=(width / 100, height / 100))
    _style_axes(ax)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, color="#f1f5f9", linewidth=0.9)

    present_cats = set(reads_df["AssignedTo"].dropna().unique())
    categories = ["Plasmid", "Human", "Tied"] + [
        c for c in ("Backbone_Only", "Ambiguous") if c in present_cats
    ]
    colors = [ASSIGNMENT_COLORS.get(cat, "#999999") for cat in categories]

    data_list = [
        reads_df[reads_df["AssignedTo"] == cat]["PlasmidScore"].dropna().tolist()
        for cat in categories
    ]

    ax.axhline(0, color="#94a3b8", linestyle=":", linewidth=1.0, zorder=1)

    bp = ax.boxplot(
        data_list,
        patch_artist=True,
        widths=0.45,
        showfliers=False,
        zorder=2,
        medianprops={"color": "#0f172a", "linewidth": 1.8},
        whiskerprops={"color": "#64748b", "linewidth": 1.2},
        capprops={"color": "#64748b", "linewidth": 1.2},
    )

    ax.set_xticks(range(1, len(categories) + 1))
    ax.set_xticklabels(categories, fontsize=10)

    for patch, color in zip(bp["boxes"], colors, strict=True):
        patch.set_facecolor(color)
        patch.set_edgecolor(color)
        patch.set_linewidth(1.4)
        patch.set_alpha(0.32)

    rng = np.random.default_rng(42)
    for i, _cat in enumerate(categories):
        cat_data = data_list[i]
        if len(cat_data) > 0:
            x = rng.normal(i + 1, 0.055, size=len(cat_data))
            ax.scatter(
                x,
                cat_data,
                alpha=0.45,
                s=9,
                color=colors[i],
                edgecolors="none",
                zorder=3,
            )

    ax.set_xlabel("Read Assignment Category", fontsize=10.5, labelpad=6)
    ax.set_ylabel("Plasmid Alignment Score", fontsize=10.5, labelpad=6)
    ax.set_title(title, fontsize=10, fontweight="medium", color="#475569", pad=8)

    plt.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved boxplot to {output_path}")


def generate_scatter_matplotlib(
    reads_df: pd.DataFrame,
    output_path: Path | str,
    title: str | None = None,
    width: int = 1100,
    height: int = 440,
    total_reads: int | None = None,
) -> None:
    """Generate single-sample score scatter plot using matplotlib."""
    n_total = total_reads if total_reads is not None else len(reads_df)
    if title is None:
        title = f"Plasmid vs. Human Alignment Score per Read (Total Reads: {n_total:,})"
    elif total_reads is not None and "(Total Reads:" not in title:
        title = f"{title} (Total Reads: {n_total:,})"

    if reads_df.empty:
        logger.warning("Empty DataFrame provided to generate_scatter_matplotlib")
        fig, ax = plt.subplots(figsize=(width / 100, height / 100))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=220, bbox_inches="tight")
        _close_and_collect(fig)
        return

    fig, ax = plt.subplots(figsize=(width / 100, height / 100))
    _style_axes(ax)
    ax.grid(True, color="#f1f5f9", linewidth=0.9)

    categories = ["Plasmid", "Human", "Tied", "Backbone_Only", "Ambiguous"]
    for cat in categories:
        cat_data = reads_df[reads_df["AssignedTo"] == cat].dropna(
            subset=["PlasmidScore", "HumanScore"]
        )
        if not cat_data.empty:
            ax.scatter(
                cat_data["PlasmidScore"],
                cat_data["HumanScore"],
                alpha=0.5,
                s=11,
                color=ASSIGNMENT_COLORS.get(cat, "#999999"),
                edgecolors="none",
                label=cat,
                zorder=3,
            )

    p_vals = reads_df["PlasmidScore"].dropna()
    h_vals = reads_df["HumanScore"].dropna()
    s_min = min(
        0.0,
        float(p_vals.min()) if not p_vals.empty else 0.0,
        float(h_vals.min()) if not h_vals.empty else 0.0,
    )
    s_max = max(
        1.0,
        float(p_vals.max()) if not p_vals.empty else 1.0,
        float(h_vals.max()) if not h_vals.empty else 1.0,
    )
    ax.plot(
        [s_min, s_max],
        [s_min, s_max],
        linestyle="--",
        color="#64748b",
        linewidth=1.4,
        label="Equal Score (y = x)",
        zorder=2,
    )

    ax.set_xlabel("Plasmid Alignment Score", fontsize=10.5, labelpad=6)
    ax.set_ylabel("Human Alignment Score", fontsize=10.5, labelpad=6)
    ax.set_title(title, fontsize=10, fontweight="medium", color="#475569", pad=8)
    ax.legend(frameon=True, facecolor="white", edgecolor="#e2e8f0", fontsize=9, loc="upper left")

    plt.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved scatter plot to {output_path}")


def generate_heatmap_matplotlib(
    ratio_data: pd.DataFrame,
    output_path: Path | str,
    threshold: float,
    plot_config: dict[str, Any],
) -> None:
    """Generate summary contamination ratio heatmap using matplotlib."""
    if ratio_data.empty:
        logger.warning("Empty DataFrame provided to generate_heatmap_matplotlib")
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=220, bbox_inches="tight")
        _close_and_collect(fig)
        return

    ratio_data = ratio_data.apply(pd.to_numeric, errors="coerce")

    n_samples, n_plasmids = ratio_data.shape
    fig_width = 8.5
    fig_height = max(4.8, min(8.0, 2.2 + n_samples * 0.45))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    cmap = LinearSegmentedColormap.from_list(
        "contam_ratio", ["#e0f2fe", "#fde68a", "#fca5a5", "#dc2626"], N=256
    )
    show_annot = n_samples <= 15 and n_plasmids <= 12

    sns.heatmap(
        ratio_data,
        cmap=cmap,
        annot=show_annot,
        fmt=".2f",
        annot_kws={"size": 8.5, "weight": "semibold"},
        cbar=True,
        cbar_kws={"shrink": 0.82, "label": "Ratio"},
        linewidths=1.5,
        linecolor="white",
        ax=ax,
    )

    ax.set_xlabel("Plasmid", fontsize=10, labelpad=6)
    ax.set_ylabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8.5)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8.5)

    plt.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved heatmap to {output_path}")


def generate_summary_boxplot_matplotlib(
    boxplot_data: pd.DataFrame,
    output_path: Path | str,
    log_offset: float = 0.001,
    marker_style: str = "outliers",
) -> None:
    """Generate summary contamination ratio boxplot (log scale) using matplotlib."""
    if boxplot_data.empty:
        logger.warning("Empty DataFrame provided to generate_summary_boxplot_matplotlib")
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=220, bbox_inches="tight")
        _close_and_collect(fig)
        return

    melted = boxplot_data.melt(var_name="Plasmid", value_name="Ratio")
    melted["Ratio"] = pd.to_numeric(melted["Ratio"], errors="coerce")
    melted["LogRatio"] = np.log10(melted["Ratio"] + log_offset)

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    _style_axes(ax)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, color="#f1f5f9", linewidth=0.9)

    plasmids = melted["Plasmid"].unique().tolist()
    data_list = [melted[melted["Plasmid"] == p]["LogRatio"].dropna().tolist() for p in plasmids]

    low_log = float(np.log10(0.6 + log_offset))
    high_log = float(np.log10(1.0 + log_offset))
    ax.axhspan(
        low_log,
        high_log,
        color="#f59e0b",
        alpha=0.12,
        label="Unclear (0.6-1.0)",
        zorder=1,
    )
    ax.axhline(
        high_log,
        color="#dc2626",
        linestyle="--",
        linewidth=1.2,
        label="Contaminated (>1.0)",
        zorder=1,
    )

    bp = ax.boxplot(
        data_list,
        patch_artist=True,
        widths=0.45,
        showfliers=False,
        zorder=2,
        medianprops={"color": "#1e3a8a", "linewidth": 1.8},
        whiskerprops={"color": "#64748b", "linewidth": 1.2},
        capprops={"color": "#64748b", "linewidth": 1.2},
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("#2563eb")
        patch.set_edgecolor("#1d4ed8")
        patch.set_alpha(0.22)

    rng = np.random.default_rng(42)
    for i, vals in enumerate(data_list):
        if vals:
            x = rng.normal(i + 1, 0.05, size=len(vals))
            ax.scatter(
                x,
                vals,
                color="#2563eb",
                alpha=0.8,
                s=26,
                zorder=3,
                edgecolors="white",
                linewidth=0.5,
            )

    ax.margins(y=0.14)
    ax.set_xticks(range(1, len(plasmids) + 1))
    ax.set_xticklabels(plasmids, rotation=30, ha="right", fontsize=8.5)
    ax.set_xlabel("Plasmid", fontsize=10, labelpad=6)
    ax.set_ylabel("log10(Contamination Ratio)", fontsize=10, labelpad=6)
    ax.legend(
        frameon=False,
        fontsize=8.5,
        loc="lower right",
        bbox_to_anchor=(1.0, 1.01),
    )

    plt.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved summary boxplot to {output_path}")


def generate_metric_heatmap_matplotlib(
    metric_data: pd.DataFrame,
    output_path: Path | str,
    metric_name: str,
    vmax: float | None = None,
) -> None:
    """Generate heatmap for any metric using matplotlib."""
    if metric_data.empty:
        logger.warning(
            f"Empty DataFrame provided to generate_metric_heatmap_matplotlib for {metric_name}"
        )
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=220, bbox_inches="tight")
        _close_and_collect(fig)
        return

    metric_data = metric_data.apply(pd.to_numeric, errors="coerce")
    plot_data = metric_data.clip(upper=vmax) if vmax is not None else metric_data

    cmap = LinearSegmentedColormap.from_list(
        "coverage", ["#f0f9ff", "#7dd3fc", "#0284c7", "#be123c"], N=256
    )

    n_samples, n_plasmids = metric_data.shape
    figwidth = 8.5
    figheight = max(4.8, min(8.0, 2.2 + n_samples * 0.45))

    fig, ax = plt.subplots(figsize=(figwidth, figheight))
    show_annot = n_samples <= 15 and n_plasmids <= 12
    fmt = ".0f" if "Mismatches" in metric_name else ".2f"

    sns.heatmap(
        plot_data,
        annot=show_annot,
        fmt=fmt,
        annot_kws={"size": 8.5, "weight": "semibold"},
        cbar=True,
        cbar_kws={"shrink": 0.82},
        cmap=cmap,
        vmin=0,
        vmax=vmax,
        linewidths=1.5,
        linecolor="white",
        mask=metric_data.isna(),
        ax=ax,
    )

    ax.set_xlabel("Plasmid", fontsize=10, labelpad=6)
    ax.set_ylabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8.5)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8.5)

    plt.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved metric heatmap to {output_path}")


def generate_coverage_boxplot_matplotlib(
    boxplot_data: pd.DataFrame,
    output_path: Path | str,
) -> None:
    """Generate boxplot for CoverageOutsideINSERT distribution using matplotlib."""
    if boxplot_data.empty:
        logger.warning("Empty DataFrame provided to generate_coverage_boxplot_matplotlib")
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=220, bbox_inches="tight")
        _close_and_collect(fig)
        return

    melted = boxplot_data.melt(var_name="Plasmid", value_name="Value")
    melted["Value"] = pd.to_numeric(melted["Value"], errors="coerce")

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    _style_axes(ax)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, color="#f1f5f9", linewidth=0.9)

    plasmid_labels: list[str] = []
    boxplot_values: list[list[float]] = []
    for plasmid, group in melted.groupby("Plasmid", sort=False):
        plasmid_labels.append(str(plasmid))
        boxplot_values.append(group["Value"].dropna().tolist())

    bp = ax.boxplot(
        boxplot_values,
        patch_artist=True,
        widths=0.45,
        showfliers=False,
        zorder=2,
        medianprops={"color": "#0c4a6e", "linewidth": 1.8},
        whiskerprops={"color": "#64748b", "linewidth": 1.2},
        capprops={"color": "#64748b", "linewidth": 1.2},
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("#0284c7")
        patch.set_edgecolor("#0369a1")
        patch.set_alpha(0.22)

    rng = np.random.default_rng(42)
    for i, vals in enumerate(boxplot_values):
        if vals:
            x = rng.normal(i + 1, 0.05, size=len(vals))
            ax.scatter(
                x,
                vals,
                color="#0284c7",
                alpha=0.8,
                s=26,
                zorder=3,
                edgecolors="white",
                linewidth=0.5,
            )

    ax.margins(y=0.14)
    ax.set_xlabel("Plasmid", fontsize=10, labelpad=6)
    ax.set_ylabel("Coverage Outside Insert", fontsize=10, labelpad=6)
    ax.set_xticks(range(1, len(plasmid_labels) + 1))
    ax.set_xticklabels(plasmid_labels, rotation=30, ha="right", fontsize=8.5)

    plt.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved coverage boxplot to {output_path}")


def generate_verdict_bar_matplotlib(
    aggregated_df: pd.DataFrame,
    output_path: Path | str,
    title: str,
    xlabel: str,
) -> None:
    """Generate stacked bar chart for verdict categories using matplotlib."""
    from .colors import VERDICT_COLORS

    if aggregated_df.empty:
        logger.warning("Empty DataFrame provided to generate_verdict_bar_matplotlib")
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=220, bbox_inches="tight")
        _close_and_collect(fig)
        return

    verdict_categories = ["Clean", "Unclear", "Contaminated", "No Result"]
    available_categories = [cat for cat in verdict_categories if cat in aggregated_df.columns]

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    _style_axes(ax)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, color="#f1f5f9", linewidth=0.9)

    x_positions = np.arange(len(aggregated_df))
    bar_width = 0.5
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
            alpha=0.92,
            edgecolor="white",
            linewidth=0.8,
        )
        bottom = bottom + counts

    ax.margins(y=0.14)
    ax.set_xlabel(xlabel, fontsize=10, labelpad=6)
    ax.set_ylabel("Count", fontsize=10, labelpad=6)
    from matplotlib.ticker import MaxNLocator

    ax.set_xticks(x_positions)
    ax.set_xticklabels(aggregated_df.index, rotation=30, ha="right", fontsize=8.5)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(
        frameon=False,
        fontsize=8.5,
        ncol=len(available_categories),
        loc="lower right",
        bbox_to_anchor=(1.0, 1.01),
    )

    plt.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    _close_and_collect(fig)
    logger.info(f"Saved verdict bar chart to {output_path}")
