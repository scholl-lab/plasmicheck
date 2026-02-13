from __future__ import annotations

import json
import logging
import os
from typing import Any

import pysam

from .utils import setup_logging  # Import the setup_logging function

# Resolve the path to config.json in the parent directory of the current script
config_path: str = os.path.join(os.path.dirname(os.path.dirname(__file__)), "config.json")

# Load configuration from JSON file
with open(config_path) as config_file:
    config: dict[str, Any] = json.load(config_file)

# Configuration variables for scoring
MATE_BONUS: int = config["scoring"]["mate_bonus"]
CLIPPING_PENALTY: int = config["scoring"]["clipping_penalty"]
MISMATCH_PENALTY: int = config["scoring"]["mismatch_penalty"]
DEFAULT_THRESHOLD: float = config["default_threshold"]
UNCLEAR_RANGE: dict[str, float] = config["unclear_range"]


def calculate_alignment_score(read: Any) -> int:
    """
    Calculate a custom alignment score based on multiple criteria.
    """
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
    """
    Parse the cDNA_positions.txt file to extract the INSERT_REGION.
    """
    if not os.path.exists(cdna_positions_file):
        raise FileNotFoundError(
            f"Expected cDNA_positions.txt at {cdna_positions_file} but the file was not found."
        )

    logging.debug(f"Looking for cDNA_positions.txt at: {cdna_positions_file}")

    with open(cdna_positions_file) as f:
        lines = f.readlines()

    cdna_start = int(lines[0].split(": ")[1])
    cdna_end = int(lines[1].split(": ")[1])
    insert_region = (cdna_start, cdna_end)

    return insert_region


def calculate_coverage_outside_insert(plasmid_bam: str, insert_region: tuple[int, int]) -> float:
    """
    Calculate the coverage of reads mapping outside the INSERT_REGION.
    Coverage is defined as the total number of aligned bases outside the INSERT_REGION
    divided by the length of the reference region outside the INSERT_REGION.
    """
    total_aligned_bases_outside_insert = 0

    with pysam.AlignmentFile(plasmid_bam, "rb") as bamfile:
        # Calculate the length of the plasmid reference
        plasmid_length = bamfile.lengths[
            bamfile.get_tid(bamfile.get_reference_name(0))
        ]  # Assumes single reference in BAM

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

            # Check if the read overlaps with the region outside the INSERT_REGION
            if read_end < insert_region[0] or read_start > insert_region[1]:
                # Entire read is outside the INSERT_REGION
                total_aligned_bases_outside_insert += read.query_length
            else:
                # Read overlaps with the INSERT_REGION, split into parts
                if read_start < insert_region[0]:  # Part of the read is before the INSERT_REGION
                    outside_start = read_start
                    outside_end = min(read_end, insert_region[0] - 1)
                    total_aligned_bases_outside_insert += outside_end - outside_start + 1

                if read_end > insert_region[1]:  # Part of the read is after the INSERT_REGION
                    outside_start = max(read_start, insert_region[1] + 1)
                    outside_end = read_end
                    total_aligned_bases_outside_insert += outside_end - outside_start + 1

    # Calculate coverage
    coverage_outside_insert = (
        total_aligned_bases_outside_insert / total_reference_bases_outside_insert
        if total_reference_bases_outside_insert > 0
        else 0
    )

    return coverage_outside_insert


def count_mismatches_near_insert_end(
    plasmid_bam: str, insert_region: tuple[int, int], window_size: int = 100
) -> dict[str, int]:
    """
    Count the number of reads with and without mismatches or clipping near the INSERT_REGION end.

    Parameters:
    - plasmid_bam: Path to the BAM file for plasmid alignment.
    - insert_region: Tuple indicating the (start, end) of the INSERT_REGION.
    - window_size: The number of bases around the insert end to consider for counting mismatches and clipping.

    Returns:
    - A dictionary with counts of reads with and without mismatches/clipping near the INSERT_REGION end.
    """
    mismatches_near_insert = {"with_mismatches_or_clipping": 0, "without_mismatches_or_clipping": 0}

    with pysam.AlignmentFile(plasmid_bam, "rb") as bamfile:
        for read in bamfile.fetch():
            if read.is_unmapped:
                continue

            read_start = read.reference_start
            read_end = read.reference_end

            # Check if the read is near the start of the INSERT_REGION
            near_insert_start = (
                read_end >= insert_region[0] - window_size
                and read_end <= insert_region[0] + window_size
            )

            # Check if the read is near the end of the INSERT_REGION
            near_insert_end = (
                read_start >= insert_region[1] - window_size
                and read_start <= insert_region[1] + window_size
            )

            if near_insert_start or near_insert_end:
                has_mismatches_or_clipping = False

                # Check CIGAR string for mismatches, clipping, or other issues
                for operation, _length in read.cigartuples:
                    if operation in {
                        1,
                        2,
                        3,
                        4,
                        5,
                    }:  # Insertion, Deletion, Skipped region, Soft clipping, Hard clipping
                        has_mismatches_or_clipping = True
                        break

                if has_mismatches_or_clipping:
                    mismatches_near_insert["with_mismatches_or_clipping"] += 1
                else:
                    mismatches_near_insert["without_mismatches_or_clipping"] += 1

    return mismatches_near_insert


def compare_alignments(
    plasmid_bam: str, human_bam: str, output_basename: str, threshold: float = DEFAULT_THRESHOLD
) -> None:
    logging.info("Starting comparison of alignments")

    plasmid_samfile = pysam.AlignmentFile(plasmid_bam, "rb")
    human_samfile = pysam.AlignmentFile(human_bam, "rb")

    logging.debug(f"Processing BAM files: {plasmid_bam} (plasmid), {human_bam} (human)")

    # Create dictionaries to store reads by query name
    plasmid_reads = {read.query_name: read for read in plasmid_samfile.fetch(until_eof=True)}
    human_reads = {read.query_name: read for read in human_samfile.fetch(until_eof=True)}

    assigned_counts = {"Plasmid": 0, "Human": 0, "Tied": 0}

    # Extract the INSERT_REGION from cDNA_positions.txt
    cdna_positions_file = os.path.join(os.path.dirname(output_basename), "cDNA_positions.txt")
    if not os.path.exists(cdna_positions_file):
        raise FileNotFoundError(
            f"Expected cDNA_positions.txt at {cdna_positions_file} but the file was not found."
        )

    insert_region = parse_insert_region(cdna_positions_file)
    logging.debug(f"INSERT_REGION extracted from {cdna_positions_file}: {insert_region}")

    # Calculate the coverage outside the INSERT_REGION
    coverage_outside_insert = calculate_coverage_outside_insert(plasmid_bam, insert_region)
    logging.debug(f"Coverage outside INSERT_REGION: {coverage_outside_insert}")

    # Count mismatches near the INSERT_REGION
    mismatches_near_insert = count_mismatches_near_insert_end(plasmid_bam, insert_region)
    logging.debug(f"Mismatches near INSERT_REGION: {mismatches_near_insert}")

    logging.debug(f"Writing alignment comparison results to {output_basename}.reads_assignment.tsv")

    with open(f"{output_basename}.reads_assignment.tsv", "w") as outfile:
        # Include headers for the new fields (CIGAR strings and mapping qualities)
        outfile.write(
            "ReadID\tAssignedTo\tPlasmidScore\tHumanScore\tPlasmidCIGAR\tHumanCIGAR\tPlasmidMapQ\tHumanMapQ\n"
        )
        for query_name, plasmid_read in plasmid_reads.items():
            plasmid_score = calculate_alignment_score(plasmid_read)
            plasmid_cigar = plasmid_read.cigarstring if plasmid_read.cigarstring else "NA"
            plasmid_mapq = plasmid_read.mapping_quality

            human_read = human_reads.get(query_name)
            human_score = calculate_alignment_score(human_read) if human_read else 0
            human_cigar = human_read.cigarstring if human_read and human_read.cigarstring else "NA"
            human_mapq = human_read.mapping_quality if human_read else "NA"

            if plasmid_score > human_score:
                assigned_to = "Plasmid"
                assigned_counts["Plasmid"] += 1
            elif human_score > plasmid_score:
                assigned_to = "Human"
                assigned_counts["Human"] += 1
            else:
                assigned_to = "Tied"
                assigned_counts["Tied"] += 1

            # Write the extended data to the output file
            outfile.write(
                f"{query_name}\t{assigned_to}\t{plasmid_score}\t{human_score}\t{plasmid_cigar}\t{human_cigar}\t{plasmid_mapq}\t{human_mapq}\n"
            )

    plasmid_count = assigned_counts["Plasmid"]
    human_count = assigned_counts["Human"]
    ratio = plasmid_count / human_count if human_count != 0 else float("inf")

    # Determine classification based on ratio
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
        summary_file.write(
            f"CoverageOutsideINSERT\t{coverage_outside_insert:.4f}\n"
        )  # Added coverage information
        summary_file.write(
            f"MismatchesNearINSERT\t{mismatches_near_insert}\n"
        )  # Added mismatches information


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
    parser.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser.add_argument("--log-file", help="Set the log output file", default=None)
    args = parser.parse_args()

    # Setup logging with the specified log level and file
    setup_logging(log_level=args.log_level.upper(), log_file=args.log_file)

    compare_alignments(args.plasmid_bam, args.human_bam, args.output_basename, args.threshold)
