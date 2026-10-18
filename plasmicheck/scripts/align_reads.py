from __future__ import annotations

import logging
import os
from shlex import quote

from plasmicheck.config import get_config

from .utils import run_command


def parse_bed_regions(bed_path: str) -> list[str]:
    if not os.path.isfile(bed_path):
        s_arg = bed_path.strip()
        if (
            s_arg
            and not s_arg.endswith(".bed")
            and ":" in s_arg
            and "-" in s_arg
            and "/" not in s_arg
        ):
            return [r.strip() for r in s_arg.split(",") if r.strip()]
        raise FileNotFoundError(f"Regions BED file not found: {bed_path}")
    regions: list[str] = []
    with open(bed_path, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(("#", "track", "browser")):
                continue
            parts = s.split("\t")
            if len(parts) < 3:
                parts = s.split()
            chrom, start, end = parts[:3]
            regions.append(f"{chrom}:{max(1, int(start) + 1)}-{int(end)}")
    if not regions:
        raise ValueError(f"No valid BED intervals found in {bed_path}")
    return regions


def has_bam_index(bam_path: str) -> bool:
    stem = bam_path[:-4] if bam_path.endswith(".bam") else bam_path
    return any(
        os.path.exists(p)
        for p in (f"{bam_path}.bai", f"{stem}.bai", f"{bam_path}.csi", f"{stem}.csi")
    )


def bam_to_fastq_cmd(
    input_bam: str,
    samtools_threads: int,
    regions: str | list[str] | None = None,
    *,
    regions_bed: str | None = None,
    include_unmapped: bool = True,
) -> str:
    eff_regions = regions if regions is not None else regions_bed
    reg_list = (
        parse_bed_regions(eff_regions) if isinstance(eff_regions, str) else list(eff_regions or [])
    )
    t = quote(str(samtools_threads))
    if reg_list:
        if not has_bam_index(input_bam):
            raise FileNotFoundError(
                f"BAM index (.bai/.csi) required for region extraction: {input_bam}"
            )
        targets = [*reg_list, "*"] if include_unmapped else reg_list
        reg_args = " ".join(quote(r) for r in dict.fromkeys(targets))
        return (
            f"samtools view -@ {t} -u -M {quote(input_bam)} {reg_args} "
            f"| samtools collate -@ {t} -u -O - "
            f"| samtools fastq -@ {t} -n -F 0x900 -"
        )
    return f"samtools collate -@ {t} -u -O {quote(input_bam)} | samtools fastq -@ {t} -n -F 0x900 -"


_bam_to_fastq_cmd = bam_to_fastq_cmd


def align_reads(
    reference_index: str,
    input_file: str,
    output_bam: str,
    alignment_type: str,
    fastq2: str | None = None,
    minimap2_threads: int | None = None,
    samtools_threads: int | None = None,
    samtools_sort_memory: str | None = None,
    regions: str | list[str] | None = None,
    include_unmapped: bool = True,
    regions_bed: str | None = None,
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
    eff_regions = regions if regions is not None else regions_bed

    if input_file.endswith(".bam"):
        command = (
            f"{bam_to_fastq_cmd(input_file, samtools_threads, regions=eff_regions, include_unmapped=include_unmapped)} "
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
