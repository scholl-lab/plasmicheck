from __future__ import annotations

import logging
import os
from shlex import quote

from plasmicheck.config import get_config

from .utils import run_command


def align_reads(
    reference_index: str,
    input_file: str,
    output_bam: str,
    alignment_type: str,
    fastq2: str | None = None,
    minimap2_threads: int | None = None,
    samtools_threads: int | None = None,
    samtools_sort_memory: str | None = None,
) -> None:
    # Resolve defaults from config when not provided
    cfg = get_config()
    if minimap2_threads is None:
        minimap2_threads = cfg["alignment"]["minimap2_threads"]
    if samtools_threads is None:
        samtools_threads = cfg["alignment"]["samtools_threads"]
    if samtools_sort_memory is None:
        samtools_sort_memory = cfg["alignment"].get("samtools_sort_memory", "2G")

    logging.info(
        f"Starting alignment of {input_file} against {reference_index} as {alignment_type} "
        f"(minimap2: {minimap2_threads}t, samtools: {samtools_threads}t, sort mem: {samtools_sort_memory})"
    )

    if alignment_type not in ["human", "plasmid"]:
        logging.error(f"Invalid alignment type: {alignment_type}")
        raise ValueError("alignment_type must be 'human' or 'plasmid'")

    fastq_extensions = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    sort_flag = "-N " if alignment_type == "human" else ""

    if input_file.endswith(".bam"):
        command = (
            f"samtools collate -@ {quote(str(samtools_threads))} -u -O {quote(input_file)} "
            f"| samtools fastq -@ {quote(str(samtools_threads))} -F 0x900 - "
            f"| minimap2 -t {quote(str(minimap2_threads))} -ax sr {quote(reference_index)} - "
            f"| samtools view -@ {quote(str(samtools_threads))} -u -F 4 - "
            f"| samtools sort {sort_flag}-@ {quote(str(samtools_threads))} "
            f"-m {quote(str(samtools_sort_memory))} -o {quote(output_bam)}"
        )
    elif input_file.endswith(fastq_extensions):
        if fastq2:
            command = (
                f"minimap2 -t {quote(str(minimap2_threads))} -ax sr "
                f"{quote(reference_index)} {quote(input_file)} {quote(fastq2)} "
                f"| samtools view -@ {quote(str(samtools_threads))} -u -F 4 - "
                f"| samtools sort {sort_flag}-@ {quote(str(samtools_threads))} "
                f"-m {quote(str(samtools_sort_memory))} -o {quote(output_bam)}"
            )
        else:
            command = (
                f"minimap2 -t {quote(str(minimap2_threads))} -ax sr "
                f"{quote(reference_index)} {quote(input_file)} "
                f"| samtools view -@ {quote(str(samtools_threads))} -u -F 4 - "
                f"| samtools sort {sort_flag}-@ {quote(str(samtools_threads))} "
                f"-m {quote(str(samtools_sort_memory))} -o {quote(output_bam)}"
            )
    else:
        logging.error(f"Unsupported input file type: {input_file}")
        raise ValueError(
            "Unsupported input file type. Must be .bam, .fastq, or paired FASTQ files."
        )

    # Run the alignment command with retry logic
    run_command(command)

    # Generate the BAI index for coordinate-sorted plasmid BAMs
    if alignment_type == "plasmid":
        run_command(f"samtools index {quote(output_bam)}")
    elif os.path.exists(f"{output_bam}.bai"):
        os.remove(f"{output_bam}.bai")

    logging.info(
        f"Alignment completed successfully for {input_file}, output written to {output_bam}"
    )
