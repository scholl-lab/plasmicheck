from __future__ import annotations

import json
import logging
import os
from typing import Any

from .utils import run_command, setup_logging  # Import setup_logging and run_command functions

# Load configuration from JSON file
with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "config.json")) as config_file:
    config: dict[str, Any] = json.load(config_file)

MINIMAP2_THREADS: int = config["alignment"]["minimap2_threads"]
SAMTOOLS_THREADS: int = config["alignment"]["samtools_threads"]


def align_reads(
    reference_index: str,
    input_file: str,
    output_bam: str,
    alignment_type: str,
    fastq2: str | None = None,
) -> None:
    logging.info(
        f"Starting alignment of {input_file} against {reference_index} as {alignment_type}"
    )

    if alignment_type not in ["human", "plasmid"]:
        logging.error(f"Invalid alignment type: {alignment_type}")
        raise ValueError("alignment_type must be 'human' or 'plasmid'")

    if input_file.endswith(".bam"):
        command = (
            f"samtools fasta {input_file} | minimap2 -t {MINIMAP2_THREADS} -ax sr {reference_index} - "
            f"| samtools view -@ {SAMTOOLS_THREADS} -h -F 4 - | samtools sort -@ {SAMTOOLS_THREADS} -o {output_bam}"
        )
    elif input_file.endswith(".fastq"):
        if fastq2:
            command = (
                f"minimap2 -t {MINIMAP2_THREADS} -ax sr {reference_index} {input_file} {fastq2} "
                f"| samtools view -@ {SAMTOOLS_THREADS} -h -F 4 - | samtools sort -@ {SAMTOOLS_THREADS} -o {output_bam}"
            )
        else:
            command = (
                f"minimap2 -t {MINIMAP2_THREADS} -ax sr {reference_index} {input_file} "
                f"| samtools view -@ {SAMTOOLS_THREADS} -h -F 4 - | samtools sort -@ {SAMTOOLS_THREADS} -o {output_bam}"
            )
    else:
        logging.error(f"Unsupported input file type: {input_file}")
        raise ValueError(
            "Unsupported input file type. Must be .bam, .fastq, or paired FASTQ files."
        )

    # Run the alignment command with retry logic
    run_command(command)

    # Generate the BAI index with retry logic
    run_command(f"samtools index {output_bam}")

    logging.info(
        f"Alignment completed successfully for {input_file}, output written to {output_bam}"
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Align reads to a reference genome and generate a BAI index"
    )
    parser.add_argument(
        "-r", "--reference_index", help="Minimap2 index for the reference genome", required=True
    )
    parser.add_argument(
        "-i",
        "--input_file",
        help="Input file (BAM, interleaved FASTQ, or first FASTQ file for paired FASTQ)",
        required=True,
    )
    parser.add_argument("-o", "--output_bam", help="Output BAM file for alignment", required=True)
    parser.add_argument(
        "-a", "--alignment_type", help="Type of alignment: 'human' or 'plasmid'", required=True
    )
    parser.add_argument("--fastq2", help="Second FASTQ file for paired FASTQ input", default=None)
    parser.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser.add_argument("--log-file", help="Set the log output file", default=None)
    args = parser.parse_args()

    # Setup logging with the specified log level and file
    setup_logging(log_level=args.log_level.upper(), log_file=args.log_file)

    align_reads(
        args.reference_index, args.input_file, args.output_bam, args.alignment_type, args.fastq2
    )
