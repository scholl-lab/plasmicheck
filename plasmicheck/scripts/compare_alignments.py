from __future__ import annotations

import itertools
import logging
import os
import subprocess
from collections.abc import Iterator
from typing import IO, Any

import pysam

from plasmicheck.config import get_config

from .utils import add_logging_args, configure_logging_from_args

_cfg = get_config()
MATE_BONUS: int = _cfg["scoring"]["mate_bonus"]
CLIPPING_PENALTY: int = _cfg["scoring"]["clipping_penalty"]
MISMATCH_PENALTY: int = _cfg["scoring"]["mismatch_penalty"]
DEFAULT_THRESHOLD: float = _cfg["default_threshold"]
UNCLEAR_RANGE: dict[str, float] = _cfg["unclear_range"]
SAMTOOLS_THREADS: int = _cfg["alignment"]["samtools_threads"]
FILTER_BACKBONE_ONLY: bool = _cfg.get("filter_backbone_only", True)
SCORE_MARGIN: int = _cfg.get("score_margin", 0)


def calculate_alignment_score(read: Any) -> int:
    """Calculate a custom alignment score based on multiple criteria."""
    if read.is_unmapped:
        return 0

    score: int = read.mapping_quality

    # Penalize for clipping
    for operation, length in read.cigartuples:
        if operation in {4, 5}:  # 4: soft clipping, 5: hard clipping
            score -= length * CLIPPING_PENALTY

    # Penalize for mismatches
    mismatches = read.get_tag("NM") if read.has_tag("NM") else 0
    score -= mismatches * MISMATCH_PENALTY

    # Bonus if mate is mapped
    if not read.mate_is_unmapped:
        score += MATE_BONUS

    return score


def parse_insert_region(cdna_positions_file: str) -> tuple[int, int]:
    """Parse the cDNA_positions.txt file to extract the INSERT_REGION."""
    if not os.path.exists(cdna_positions_file):
        raise FileNotFoundError(
            f"Expected cDNA_positions.txt at {cdna_positions_file} but the file was not found."
        )

    logging.debug(f"Looking for cDNA_positions.txt at: {cdna_positions_file}")

    with open(cdna_positions_file) as f:
        lines = f.readlines()

    cdna_start = int(lines[0].split(": ")[1])
    cdna_end = int(lines[1].split(": ")[1])
    return (cdna_start, cdna_end)


def read_overlaps_insert(read: Any, insert_region: tuple[int, int]) -> bool:
    """Check if read alignment overlaps the insert region.

    Uses pysam half-open interval convention:
    - read.reference_start: 0-based, inclusive
    - read.reference_end: 0-based, exclusive

    insert_region uses inclusive boundaries on both ends (from cDNA_positions.txt).
    """
    if read.is_unmapped:
        return False
    read_start = read.reference_start
    read_end = read.reference_end
    if read_start is None or read_end is None:
        return False
    # No overlap: read ends before insert starts, or read starts after insert ends
    return not (read_end <= insert_region[0] or read_start > insert_region[1])


def calculate_coverage_outside_insert(plasmid_bam: str, insert_region: tuple[int, int]) -> float:
    """Calculate the coverage of reads mapping outside the INSERT_REGION."""
    total_aligned_bases_outside_insert = 0

    with pysam.AlignmentFile(plasmid_bam, "rb") as bamfile:
        plasmid_length = bamfile.lengths[bamfile.get_tid(bamfile.get_reference_name(0))]

        total_reference_bases_outside_insert = plasmid_length - (
            insert_region[1] - insert_region[0] + 1
        )
        if total_reference_bases_outside_insert <= 0:
            raise ValueError("The INSERT_REGION covers the entire reference sequence.")

        for read in bamfile.fetch():
            if read.is_unmapped:
                continue

            read_start = read.reference_start
            read_end = read.reference_end
            if read_start is None or read_end is None:
                continue

            if read_end < insert_region[0] or read_start > insert_region[1]:
                total_aligned_bases_outside_insert += read.query_length
            else:
                if read_start < insert_region[0]:
                    outside_start = read_start
                    outside_end = min(read_end, insert_region[0] - 1)
                    total_aligned_bases_outside_insert += outside_end - outside_start + 1
                if read_end > insert_region[1]:
                    outside_start = max(read_start, insert_region[1] + 1)
                    outside_end = read_end
                    total_aligned_bases_outside_insert += outside_end - outside_start + 1

    coverage_outside_insert: float = (
        total_aligned_bases_outside_insert / total_reference_bases_outside_insert
    )
    return coverage_outside_insert


def count_mismatches_near_insert_end(
    plasmid_bam: str, insert_region: tuple[int, int], window_size: int = 100
) -> dict[str, int]:
    """Count reads with/without mismatches or clipping near the INSERT_REGION boundaries."""
    mismatches_near_insert = {"with_mismatches_or_clipping": 0, "without_mismatches_or_clipping": 0}

    with pysam.AlignmentFile(plasmid_bam, "rb") as bamfile:
        for read in bamfile.fetch():
            if read.is_unmapped:
                continue

            read_start = read.reference_start
            read_end = read.reference_end
            if read_start is None or read_end is None:
                continue

            near_insert_start = (
                read_end >= insert_region[0] - window_size
                and read_end <= insert_region[0] + window_size
            )
            near_insert_end = (
                read_start >= insert_region[1] - window_size
                and read_start <= insert_region[1] + window_size
            )

            if near_insert_start or near_insert_end:
                has_mismatches_or_clipping = False
                cigartuples = read.cigartuples
                if cigartuples is None:
                    continue

                for operation, _length in cigartuples:
                    if operation in {1, 2, 3, 4, 5}:
                        has_mismatches_or_clipping = True
                        break

                if has_mismatches_or_clipping:
                    mismatches_near_insert["with_mismatches_or_clipping"] += 1
                else:
                    mismatches_near_insert["without_mismatches_or_clipping"] += 1

    return mismatches_near_insert


# ---------------------------------------------------------------------------
# Streaming name-sorted merge (O(1) memory)
# ---------------------------------------------------------------------------


def _namesort_bam(input_bam: str, output_bam: str) -> None:
    """Sort a BAM by read name using samtools sort -n."""
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(SAMTOOLS_THREADS), "-o", output_bam, input_bam],
        check=True,
    )


def _iter_reads_by_name(bam_path: str) -> Iterator[tuple[str, list[Any]]]:
    """Yield ``(query_name, [reads])`` groups from a name-sorted BAM."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for name, group in itertools.groupby(bam.fetch(until_eof=True), key=lambda r: r.query_name):
            yield name, list(group)


def _best_read(reads: list[Any]) -> Any:
    """Pick the primary alignment (or first) from a group of reads with the same name."""
    for r in reads:
        if not r.is_secondary and not r.is_supplementary:
            return r
    return reads[0]


def _write_assignment(
    outfile: IO[str],
    query_name: str,
    assigned_to: str,
    plasmid_score: int,
    human_score: int,
    plasmid_read: Any | None,
    human_read: Any | None,
) -> None:
    """Write a single line to the reads_assignment TSV."""
    p_cigar = plasmid_read.cigarstring if plasmid_read and plasmid_read.cigarstring else "NA"
    h_cigar = human_read.cigarstring if human_read and human_read.cigarstring else "NA"
    p_mapq = plasmid_read.mapping_quality if plasmid_read else "NA"
    h_mapq = human_read.mapping_quality if human_read else "NA"
    outfile.write(
        f"{query_name}\t{assigned_to}\t{plasmid_score}\t{human_score}\t{p_cigar}\t{h_cigar}\t{p_mapq}\t{h_mapq}\n"
    )


def _assign(
    plasmid_score: int,
    human_score: int,
    *,
    plasmid_read: Any | None = None,
    insert_region: tuple[int, int] | None = None,
    score_margin: int = 0,
) -> str:
    """Assign read to category with optional backbone and margin filtering.

    Order of checks:
    1. Score margin (if enabled and scores differ) -> Ambiguous
    2. Exact tie -> Tied
    3. Plasmid wins + insert_region available -> check overlap -> Plasmid or Backbone_Only
    4. Plasmid wins + no insert_region -> Plasmid (fallback)
    5. Human wins -> Human
    """
    # Score margin check first (if enabled)
    if score_margin > 0:
        score_diff = abs(plasmid_score - human_score)
        if 0 < score_diff < score_margin:
            return "Ambiguous"

    # Exact tie
    if plasmid_score == human_score:
        return "Tied"

    # Plasmid wins
    if plasmid_score > human_score:
        if insert_region is not None and plasmid_read is not None:
            if read_overlaps_insert(plasmid_read, insert_region):
                return "Plasmid"
            else:
                return "Backbone_Only"
        return "Plasmid"

    # Human wins
    return "Human"


def _streaming_compare(
    plasmid_ns_bam: str,
    human_ns_bam: str,
    outfile: IO[str],
    *,
    insert_region: tuple[int, int] | None = None,
    score_margin: int = 0,
) -> dict[str, int]:
    """Two-pointer merge over name-sorted BAMs.  Returns assigned counts."""
    assigned_counts: dict[str, int] = {
        "Plasmid": 0, "Human": 0, "Tied": 0,
        "Backbone_Only": 0, "Ambiguous": 0,
    }

    plasmid_iter = _iter_reads_by_name(plasmid_ns_bam)
    human_iter = _iter_reads_by_name(human_ns_bam)

    p_item: tuple[str, list[Any]] | None = next(plasmid_iter, None)
    h_item: tuple[str, list[Any]] | None = next(human_iter, None)

    while p_item is not None or h_item is not None:
        p_name = p_item[0] if p_item else None
        h_name = h_item[0] if h_item else None

        if p_name is not None and p_name == h_name:
            # Both BAMs have this read — compare scores
            assert p_item is not None
            assert h_item is not None
            p_read = _best_read(p_item[1])
            h_read = _best_read(h_item[1])
            ps = calculate_alignment_score(p_read)
            hs = calculate_alignment_score(h_read)
            assigned = _assign(ps, hs, plasmid_read=p_read, insert_region=insert_region, score_margin=score_margin)
            assigned_counts[assigned] += 1
            _write_assignment(outfile, p_name, assigned, ps, hs, p_read, h_read)
            p_item = next(plasmid_iter, None)
            h_item = next(human_iter, None)

        elif h_name is None or (p_name is not None and p_name < h_name):
            # Plasmid-only read
            assert p_item is not None
            assert p_name is not None
            p_read = _best_read(p_item[1])
            ps = calculate_alignment_score(p_read)
            assigned = _assign(ps, 0, plasmid_read=p_read, insert_region=insert_region, score_margin=score_margin)
            assigned_counts[assigned] += 1
            _write_assignment(outfile, p_name, assigned, ps, 0, p_read, None)
            p_item = next(plasmid_iter, None)

        else:
            # Human-only read
            assert h_item is not None
            h_read = _best_read(h_item[1])
            hs = calculate_alignment_score(h_read)
            assigned = _assign(0, hs, score_margin=score_margin)
            assigned_counts[assigned] += 1
            _write_assignment(outfile, h_name, assigned, 0, hs, None, h_read)
            h_item = next(human_iter, None)

    return assigned_counts


def compare_alignments(
    plasmid_bam: str, human_bam: str, output_basename: str, threshold: float = DEFAULT_THRESHOLD
) -> None:
    logging.info("Starting comparison of alignments")
    logging.debug(f"Processing BAM files: {plasmid_bam} (plasmid), {human_bam} (human)")

    # Parse insert region with graceful fallback
    cdna_positions_file = os.path.join(os.path.dirname(output_basename), "cDNA_positions.txt")
    insert_region: tuple[int, int] | None = None
    try:
        insert_region = parse_insert_region(cdna_positions_file)
        logging.debug(f"INSERT_REGION extracted from {cdna_positions_file}: {insert_region}")
    except (FileNotFoundError, ValueError, IndexError) as e:
        logging.warning(
            f"Backbone filtering unavailable: {e}. "
            "Falling back to pre-v0.33.0 behavior (no insert-region filtering)."
        )

    coverage_outside_insert = 0.0
    mismatches_near_insert: dict[str, int] = {
        "with_mismatches_or_clipping": 0,
        "without_mismatches_or_clipping": 0,
    }
    if insert_region is not None:
        coverage_outside_insert = calculate_coverage_outside_insert(plasmid_bam, insert_region)
        logging.debug(f"Coverage outside INSERT_REGION: {coverage_outside_insert}")
        mismatches_near_insert = count_mismatches_near_insert_end(plasmid_bam, insert_region)
        logging.debug(f"Mismatches near INSERT_REGION: {mismatches_near_insert}")
    else:
        logging.warning("Skipping coverage/mismatch metrics (insert region unavailable).")

    # Group BAMs by read name into temporary files for streaming merge
    plasmid_ns = plasmid_bam.replace(".bam", ".namesorted.bam")
    human_ns = human_bam.replace(".bam", ".namesorted.bam")
    logging.info("Grouping BAMs by read name for streaming comparison...")
    _namesort_bam(plasmid_bam, plasmid_ns)
    _namesort_bam(human_bam, human_ns)

    logging.debug(f"Writing alignment comparison results to {output_basename}.reads_assignment.tsv")

    with open(f"{output_basename}.reads_assignment.tsv", "w") as outfile:
        outfile.write(
            "ReadID\tAssignedTo\tPlasmidScore\tHumanScore\tPlasmidCIGAR\tHumanCIGAR\tPlasmidMapQ\tHumanMapQ\n"
        )
        assigned_counts = _streaming_compare(
            plasmid_ns, human_ns, outfile,
            insert_region=insert_region,
            score_margin=SCORE_MARGIN,
        )

    # Clean up temporary name-sorted BAMs
    for ns_bam in (plasmid_ns, human_ns):
        try:
            os.remove(ns_bam)
        except OSError:
            logging.warning(f"Could not remove temporary file: {ns_bam}")

    plasmid_count = assigned_counts["Plasmid"]
    human_count = assigned_counts["Human"]
    if FILTER_BACKBONE_ONLY:
        # Exclude Backbone_Only and Ambiguous from ratio (new v0.33.0 behavior)
        ratio = plasmid_count / human_count if human_count != 0 else float("inf")
    else:
        # Include all plasmid-favoring categories (pre-v0.33.0 behavior)
        effective_plasmid = plasmid_count + assigned_counts["Backbone_Only"] + assigned_counts["Ambiguous"]
        ratio = effective_plasmid / human_count if human_count != 0 else float("inf")

    if ratio > threshold:
        verdict = "Sample is contaminated with plasmid DNA"
    elif UNCLEAR_RANGE["lower_bound"] <= ratio <= UNCLEAR_RANGE["upper_bound"]:
        verdict = "Sample contamination status is unclear"
    else:
        verdict = "Sample is not contaminated with plasmid DNA"

    logging.info(f"Comparison completed. Verdict: {verdict}, Ratio: {ratio}")

    with open(f"{output_basename}.summary.tsv", "w") as summary_file:
        summary_file.write("Category\tCount\n")
        for category, count in assigned_counts.items():
            summary_file.write(f"{category}\t{count}\n")
        summary_file.write(f"Verdict\t{verdict}\n")
        summary_file.write(f"Ratio\t{ratio}\n")
        summary_file.write(f"CoverageOutsideINSERT\t{coverage_outside_insert:.4f}\n")
        summary_file.write(f"MismatchesNearINSERT\t{mismatches_near_insert}\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Compare alignments and assign reads to plasmid or human reference"
    )
    parser.add_argument("-p", "--plasmid_bam", help="BAM file for plasmid alignment", required=True)
    parser.add_argument("-m", "--human_bam", help="BAM file for human alignment", required=True)
    parser.add_argument("-o", "--output_basename", help="Basename for output files", required=True)
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})",
    )
    add_logging_args(parser)
    args = parser.parse_args()

    configure_logging_from_args(args)

    compare_alignments(args.plasmid_bam, args.human_bam, args.output_basename, args.threshold)
