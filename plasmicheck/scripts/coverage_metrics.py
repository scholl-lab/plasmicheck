"""Coverage metrics computation for plasmid regions.

This module provides functions to compute per-base depth and per-region
coverage statistics (mean, median, breadth, uniformity) from BAM alignments.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt
import pysam

from plasmicheck.config import get_config

if TYPE_CHECKING:
    from .compare_alignments import InsertRegionData

_cfg = get_config()


def get_depth_array(
    bam_path: str, contig: str, start: int, end: int, read_callback: Any = "all"
) -> npt.NDArray[Any]:
    """Extract per-base depth array for a genomic region.

    Args:
        bam_path: Path to BAM alignment file
        contig: Reference contig name
        start: 0-based inclusive start position
        end: 0-based exclusive end position (half-open interval)
        read_callback: Filter callback or "all"/"nofilter" passed to pysam.count_coverage

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
            read_callback=read_callback,
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

    # Coefficient of variation (check for zero mean and single element)
    if mean_val == 0 or len(depth_array) == 1:
        stats["cv"] = 0.0
    else:
        stats["cv"] = float(np.std(depth_array, ddof=1)) / mean_val

    return stats


def compute_region_coverage_metrics(
    plasmid_bam: str,
    insert_region: tuple[int, int] | InsertRegionData | None,
    breadth_thresholds: list[int],
    plasmid_depths: npt.NDArray[Any] | None = None,
) -> tuple[dict[str, dict[str, float]], bool]:
    """Compute coverage metrics for insert and backbone regions.

    Args:
        plasmid_bam: Path to plasmid alignment BAM file
        insert_region: (start, end) inclusive boundaries of insert region, or None for fallback
        breadth_thresholds: List of depth thresholds for breadth calculation
        plasmid_depths: Optional precomputed whole-plasmid depth array

    Returns:
        Tuple of:
        - Dictionary with "insert" and "backbone" keys, each containing coverage stats
        - Boolean indicating fallback mode (True if insert_region was None)

    Fallback behavior:
        When insert_region is None, computes metrics for whole plasmid as "insert",
        and returns empty backbone metrics (all 0.0).
    """
    if plasmid_depths is not None:
        contig = ""
        plasmid_length = len(plasmid_depths)
    else:
        with pysam.AlignmentFile(plasmid_bam, "rb") as bam:
            contig = bam.get_reference_name(0)
            plasmid_length = bam.lengths[0]

    if insert_region is None:
        # Fallback: whole plasmid as insert, empty backbone
        logging.debug("Computing coverage for whole plasmid (insert region not defined)")
        whole_depth = (
            plasmid_depths
            if plasmid_depths is not None
            else get_depth_array(plasmid_bam, contig, 0, plasmid_length)
        )
        return (
            {
                "insert": compute_coverage_stats(whole_depth, breadth_thresholds),
                "backbone": compute_coverage_stats(
                    np.array([], dtype=np.int64), breadth_thresholds
                ),
            },
            True,
        )

    insert_start, insert_end = insert_region[0], insert_region[1]
    is_wrapped = insert_start > insert_end
    logging.debug(
        f"Computing coverage for insert [{insert_start}, {insert_end}] and backbone regions"
    )

    if plasmid_depths is not None:
        if is_wrapped:
            insert_depth = np.concatenate(
                [plasmid_depths[insert_start:], plasmid_depths[: insert_end + 1]]
            )
            backbone_depth = plasmid_depths[insert_end + 1 : insert_start]
        else:
            insert_depth = plasmid_depths[insert_start : insert_end + 1]
            backbone_parts = []
            if insert_start > 0:
                backbone_parts.append(plasmid_depths[:insert_start])
            if insert_end + 1 < len(plasmid_depths):
                backbone_parts.append(plasmid_depths[insert_end + 1 :])
            backbone_depth = (
                np.concatenate(backbone_parts) if backbone_parts else np.array([], dtype=np.int64)
            )
    else:
        if is_wrapped:
            insert_depth = np.concatenate(
                [
                    get_depth_array(plasmid_bam, contig, insert_start, plasmid_length),
                    get_depth_array(plasmid_bam, contig, 0, insert_end + 1),
                ]
            )
            backbone_depth = get_depth_array(plasmid_bam, contig, insert_end + 1, insert_start)
        else:
            insert_depth = get_depth_array(plasmid_bam, contig, insert_start, insert_end + 1)
            backbone_parts = []
            if insert_start > 0:
                backbone_parts.append(get_depth_array(plasmid_bam, contig, 0, insert_start))
            if insert_end + 1 < plasmid_length:
                backbone_parts.append(
                    get_depth_array(plasmid_bam, contig, insert_end + 1, plasmid_length)
                )
            backbone_depth = (
                np.concatenate(backbone_parts) if backbone_parts else np.array([], dtype=np.int64)
            )

    return (
        {
            "insert": compute_coverage_stats(insert_depth, breadth_thresholds),
            "backbone": compute_coverage_stats(backbone_depth, breadth_thresholds),
        },
        False,
    )
