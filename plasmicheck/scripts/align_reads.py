from __future__ import annotations

import logging

from plasmicheck.config import get_config

from .utils import add_logging_args, configure_logging_from_args, run_command


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

    if input_file.endswith(".bam"):
        command = (
            f"samtools fasta {input_file} | minimap2 -t {minimap2_threads} -ax sr {reference_index} - "
            f"| samtools view -@ {samtools_threads} -h -F 4 - | samtools sort -@ {samtools_threads} -m {samtools_sort_memory} -o {output_bam}"
        )
    elif input_file.endswith(fastq_extensions):
        if fastq2:
            command = (
                f"minimap2 -t {minimap2_threads} -ax sr {reference_index} {input_file} {fastq2} "
                f"| samtools view -@ {samtools_threads} -h -F 4 - | samtools sort -@ {samtools_threads} -m {samtools_sort_memory} -o {output_bam}"
            )
        else:
            command = (
                f"minimap2 -t {minimap2_threads} -ax sr {reference_index} {input_file} "
                f"| samtools view -@ {samtools_threads} -h -F 4 - | samtools sort -@ {samtools_threads} -m {samtools_sort_memory} -o {output_bam}"
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
    add_logging_args(parser)
    args = parser.parse_args()

    configure_logging_from_args(args)

    align_reads(
        args.reference_index, args.input_file, args.output_bam, args.alignment_type, args.fastq2
    )
