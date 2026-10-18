from __future__ import annotations

import logging
import os
from shlex import quote

from plasmicheck.config import get_config

from .utils import run_command


def _bam_to_fastq_cmd(
    input_file: str,
    samtools_threads: int,
    regions_bed: str | None = None,
    include_unmapped: bool = False,
) -> str:
    """Build shell command that extracts collated FASTQ records from a BAM file."""
    t_arg = quote(str(samtools_threads))
    bam_arg = quote(input_file)
    if not regions_bed:
        return f"samtools collate -@ {t_arg} -u -O {bam_arg} | samtools fastq -@ {t_arg} -F 0x900 -"
    if not os.path.isfile(regions_bed):
        raise FileNotFoundError(f"Regions BED file not found: {regions_bed}")
    bed_arg = quote(regions_bed)
    if not include_unmapped:
        return (
            f"samtools view -@ {t_arg} -u -F 4 -L {bed_arg} {bam_arg} "
            f"| samtools collate -@ {t_arg} -u -O - "
            f"| samtools fastq -@ {t_arg} -F 0x904 -"
        )
    return (
        f"{{ samtools view -@ {t_arg} -u -F 4 -L {bed_arg} {bam_arg} "
        f"| samtools collate -@ {t_arg} -u -O - "
        f"| samtools fastq -@ {t_arg} -F 0x904 -; "
        f"samtools view -@ {t_arg} -u -f 4 {bam_arg} "
        f"| samtools collate -@ {t_arg} -u -O - "
        f"| samtools fastq -@ {t_arg} -f 4 -F 0x900 -; }}"
    )


def align_reads(
    reference_index: str,
    input_file: str,
    output_bam: str,
    alignment_type: str,
    fastq2: str | None = None,
    minimap2_threads: int | None = None,
    samtools_threads: int | None = None,
    samtools_sort_memory: str | None = None,
    regions_bed: str | None = None,
    include_unmapped: bool = False,
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
        extract_cmd = _bam_to_fastq_cmd(
            input_file, samtools_threads, regions_bed=regions_bed, include_unmapped=include_unmapped
        )
        command = (
            f"{extract_cmd} "
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
