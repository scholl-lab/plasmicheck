"""Coverage metrics computation for plasmid regions.

This module provides functions to compute per-base depth and per-region
coverage statistics (mean, median, breadth, uniformity) from BAM alignments.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import numpy.typing as npt
import pysam

from plasmicheck.config import get_config

_cfg = get_config()


def get_depth_array(bam_path: str, contig: str, start: int, end: int) -> npt.NDArray[Any]:
    """Extract per-base depth array for a genomic region.

    Args:
        bam_path: Path to BAM alignment file
        contig: Reference contig name
        start: 0-based inclusive start position
        end: 0-based exclusive end position (half-open interval)

    Returns:
        Numpy array of per-base depths. Empty array if start >= end.

    Note:
        Uses pysam half-open interval convention [start, end).
        Insert region boundaries from cDNA_positions.txt are inclusive,
        so caller must pass (insert_start, insert_end + 1).
    """
    if start >= end:
        return np.array([])

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # count_coverage returns tuple of 4 arrays (A, C, G, T counts)
        cov = bam.count_coverage(
            contig=contig,
            start=start,
            stop=end,
            quality_threshold=0,  # Include all bases (match alignment behavior)
            read_callback="all",
        )
        # Sum A, C, G, T counts to get total depth per position
        depth: npt.NDArray[Any] = (
            np.array(cov[0]) + np.array(cov[1]) + np.array(cov[2]) + np.array(cov[3])
        )

    return depth


def compute_coverage_stats(
    depth_array: npt.NDArray[Any], breadth_thresholds: list[int]
) -> dict[str, float]:
    """Compute coverage statistics from per-base depth array.

    Args:
        depth_array: Numpy array of per-base depths
        breadth_thresholds: List of depth thresholds for breadth calculation (e.g., [1, 5, 10])

    Returns:
        Dictionary with keys:
        - mean_depth: Mean depth across all bases
        - median_depth: Median depth across all bases
        - breadth_1x: Fraction of bases with depth >= 1 (always computed)
        - breadth_Nx: Fraction of bases with depth >= N (for each threshold in breadth_thresholds)
        - cv: Coefficient of variation (std/mean, 0.0 if mean is 0 or array is empty)

    All values are Python floats (not numpy scalars).
    """
    if len(depth_array) == 0:
        stats = {
            "mean_depth": 0.0,
            "median_depth": 0.0,
            "breadth_1x": 0.0,
            "cv": 0.0,
        }
        for threshold in breadth_thresholds:
            if threshold != 1:
                stats[f"breadth_{threshold}x"] = 0.0
        return stats

    mean_val = float(np.mean(depth_array))
    median_val = float(np.median(depth_array))

    stats = {
        "mean_depth": mean_val,
        "median_depth": median_val,
    }

    # Always compute breadth at 1x
    breadth_1x = float(np.sum(depth_array >= 1)) / len(depth_array)
    stats["breadth_1x"] = breadth_1x

    # Breadth at each configured threshold (skip 1 if already computed)
    for threshold in breadth_thresholds:
        if threshold != 1:
            breadth = float(np.sum(depth_array >= threshold)) / len(depth_array)
            stats[f"breadth_{threshold}x"] = breadth

    # Coefficient of variation (check for zero mean)
    if mean_val == 0:
        stats["cv"] = 0.0
    else:
        stats["cv"] = float(np.std(depth_array, ddof=1)) / mean_val

    return stats


def compute_region_coverage_metrics(
    plasmid_bam: str,
    insert_region: tuple[int, int] | None,
    breadth_thresholds: list[int],
) -> tuple[dict[str, dict[str, float]], bool]:
    """Compute coverage metrics for insert and backbone regions.

    Args:
        plasmid_bam: Path to plasmid alignment BAM file
        insert_region: (start, end) inclusive boundaries of insert region, or None for fallback
        breadth_thresholds: List of depth thresholds for breadth calculation

    Returns:
        Tuple of:
        - Dictionary with "insert" and "backbone" keys, each containing coverage stats
        - Boolean indicating fallback mode (True if insert_region was None)

    Fallback behavior:
        When insert_region is None, computes metrics for whole plasmid as "insert",
        and returns empty backbone metrics (all 0.0).
    """
    with pysam.AlignmentFile(plasmid_bam, "rb") as bam:
        contig = bam.get_reference_name(0)
        plasmid_length = bam.lengths[0]

    if insert_region is None:
        # Fallback: whole plasmid as insert, empty backbone
        logging.debug("Computing coverage for whole plasmid (insert region not defined)")
        whole_depth = get_depth_array(plasmid_bam, contig, 0, plasmid_length)
        backbone_stats = {
            "mean_depth": 0.0,
            "median_depth": 0.0,
            "cv": 0.0,
        }
        for threshold in breadth_thresholds:
            backbone_stats[f"breadth_{threshold}x"] = 0.0
        return (
            {
                "insert": compute_coverage_stats(whole_depth, breadth_thresholds),
                "backbone": backbone_stats,
            },
            True,
        )

    # Insert region (inclusive boundaries -> half-open for pysam)
    insert_depth = get_depth_array(plasmid_bam, contig, insert_region[0], insert_region[1] + 1)
    logging.debug(
        f"Computing coverage for insert [{insert_region[0]}, {insert_region[1]}] "
        f"and backbone regions"
    )

    # Backbone regions (before insert + after insert)
    backbone_parts = []
    if insert_region[0] > 0:
        backbone_parts.append(get_depth_array(plasmid_bam, contig, 0, insert_region[0]))
    if insert_region[1] + 1 < plasmid_length:
        backbone_parts.append(
            get_depth_array(plasmid_bam, contig, insert_region[1] + 1, plasmid_length)
        )

    backbone_depth = np.concatenate(backbone_parts) if backbone_parts else np.array([])

    return (
        {
            "insert": compute_coverage_stats(insert_depth, breadth_thresholds),
            "backbone": compute_coverage_stats(backbone_depth, breadth_thresholds),
        },
        False,
    )
