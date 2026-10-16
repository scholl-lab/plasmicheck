from __future__ import annotations

import itertools
import logging
import os
from shlex import quote
from typing import Any

import pysam

from plasmicheck.config import get_config

from .utils import run_command

_cfg = get_config()
PADDING_DEFAULT: int = _cfg["padding"]
FASTA_EXTENSIONS: list[str] = _cfg["alignment"]["fasta_extensions"]
MINIMAP2_THREADS: int = _cfg["alignment"]["minimap2_threads"]
SAMTOOLS_THREADS: int = _cfg["alignment"]["samtools_threads"]
THREADS_MINIMAP2: int = MINIMAP2_THREADS
THREADS_SAMTOOLS: int = SAMTOOLS_THREADS


def find_fasta_file(base_name: str) -> str:
    logging.info(f"Finding FASTA file with base name: {base_name}")
    for ext in FASTA_EXTENSIONS:
        fasta_file = base_name + ext
        if os.path.isfile(fasta_file):
            logging.debug(f"Found FASTA file: {fasta_file}")
            return fasta_file
    logging.error(
        f"No FASTA file found for base name {base_name} with extensions {FASTA_EXTENSIONS}"
    )
    raise FileNotFoundError(
        f"No FASTA file found for base name {base_name} with extensions {FASTA_EXTENSIONS}"
    )


def spliced_alignment(
    human_index: str,
    plasmid_fasta: str,
    output_bam: str,
    padding: int = PADDING_DEFAULT,
    threads_minimap2: int | None = None,
    threads_samtools: int | None = None,
) -> list[tuple[str, int, int]]:
    mm_threads = threads_minimap2 if threads_minimap2 is not None else THREADS_MINIMAP2
    st_threads = threads_samtools if threads_samtools is not None else THREADS_SAMTOOLS
    logging.info(
        f"Performing spliced alignment with {plasmid_fasta} to {human_index} and writing to {output_bam}"
    )
    try:
        run_command(
            f"minimap2 -t {quote(str(mm_threads))} -ax splice "
            f"{quote(human_index)} {quote(plasmid_fasta)} "
            f"| samtools view -@ {quote(str(st_threads))} -h -F 4 - "
            f"| samtools sort -@ {quote(str(st_threads))} -o {quote(output_bam)}"
        )
    except Exception as e:
        logging.error(f"Failed to perform spliced alignment: {e}")
        raise

    if not os.path.isfile(output_bam):
        logging.error(f"Spliced alignment output BAM file not found: {output_bam}")
        raise FileNotFoundError(f"Expected BAM file at {output_bam} but it was not found.")

    logging.info(f"Indexing BAM file: {output_bam}")
    try:
        run_command(f"samtools index {quote(output_bam)}")
    except Exception as e:
        logging.error(f"Failed to index BAM file: {e}")
        raise

    logging.info(f"Extracting spanned human reference regions with padding: {padding}")
    bamfile: Any = pysam.AlignmentFile(output_bam, "rb")
    valid_reads = [
        read for read in bamfile.fetch() if not read.is_unmapped and not read.is_secondary
    ]
    bamfile.close()

    if not valid_reads:
        raise ValueError(f"No human alignment regions found in {output_bam}")

    # Target chromosome filtering: if any non-secondary mapped record has splice junctions,
    # determine target chromosome as the one with largest total spliced alignment reference span
    spliced_reads = [r for r in valid_reads if any(op == 3 for op, _ in (r.cigartuples or []))]
    if spliced_reads:
        chrom_spans: dict[str, int] = {}
        for r in spliced_reads:
            span = (r.reference_end or 0) - (r.reference_start or 0)
            chrom_spans[r.reference_name] = chrom_spans.get(r.reference_name, 0) + span
        target_chrom = max(chrom_spans, key=chrom_spans.get)  # type: ignore[arg-type]
        candidate_reads = [r for r in valid_reads if r.reference_name == target_chrom]
    else:
        candidate_reads = valid_reads

    spanned_regions: list[tuple[str, int, int]] = []
    for read in candidate_reads:
        if read.reference_start is None or read.reference_end is None:
            continue
        start = max(0, read.reference_start - padding)
        end = read.reference_end + padding
        spanned_regions.append((read.reference_name, start, end))

    logging.debug(f"Spanned regions: {spanned_regions}")

    # Merge overlapping/adjacent intervals per chromosome
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(set(spanned_regions)):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def extract_human_reference(
    human_fasta: str, spanned_regions: list[tuple[str, int, int]], output_fasta: str
) -> None:
    logging.info(f"Extracting human reference regions from {human_fasta} to {output_fasta}")

    # Use pysam to open the FASTA file
    fasta: Any = pysam.FastaFile(human_fasta)

    with open(output_fasta, "w") as output_handle:
        for region in spanned_regions:
            chrom: str = region[0]
            start: int = region[1]
            end: int = region[2]
            logging.debug(f"Fetching region {chrom}:{start}-{end}")
            sequence: str = fasta.fetch(chrom, start, end)

            # Write the fetched sequence to the output FASTA file
            output_handle.write(f">{chrom}:{start}-{end}\n{sequence}\n")

    fasta.close()
    logging.info(f"Human reference regions have been extracted to {output_fasta}")


def extract_plasmid_cdna_positions(plasmid_fasta: str, bam_file: str, output_file: str) -> None:
    logging.info(
        f"Extracting cDNA start and end positions from {bam_file} and saving to {output_file}"
    )

    bamfile: Any = pysam.AlignmentFile(bam_file, "rb")
    records_info: list[dict[str, Any]] = []

    for read in bamfile.fetch():
        if read.is_unmapped or read.is_secondary:
            continue

        ops = read.cigartuples or []
        lead = sum(n for op, n in itertools.takewhile(lambda t: t[0] in (4, 5), ops))
        trail = sum(n for op, n in itertools.takewhile(lambda t: t[0] in (4, 5), reversed(ops)))
        if read.is_reverse:
            lead, trail = trail, lead
        rlen = sum(n for op, n in ops if op in (0, 1, 4, 5, 7, 8))
        q_start = lead
        q_end = rlen - trail - 1 if rlen > trail else lead
        junctions = sum(1 for op, _ in ops if op == 3)

        records_info.append(
            {
                "read": read,
                "ops": ops,
                "lead": lead,
                "trail": trail,
                "rlen": rlen,
                "q_start": q_start,
                "q_end": q_end,
                "junctions": junctions,
            }
        )

    bamfile.close()

    def _exons_and_introns(ops_list: list[tuple[int, int]]) -> tuple[list[int], list[int]]:
        exons: list[int] = []
        introns: list[int] = []
        curr = 0
        for op, length in ops_list:
            if op in (0, 2, 7, 8):
                curr += length
            elif op == 3:
                exons.append(curr)
                introns.append(length)
                curr = 0
        if curr > 0 or not exons:
            exons.append(curr)
        return exons, introns

    cdna_start: int | None = None
    cdna_end: int | None = None
    junction_count: int | None = None
    exon_count: int | None = None
    exon_lengths: list[int] | None = None
    intron_lengths: list[int] | None = None
    strand: str | None = None
    chromosome: str | None = None
    genomic_start: int | None = None
    genomic_end: int | None = None

    if records_info:
        spliced_records = [r for r in records_info if r["junctions"] > 0]
        if spliced_records:
            chosen = max(
                spliced_records,
                key=lambda r: (
                    r["q_end"] - r["q_start"],
                    r["junctions"],
                    1 if not r["read"].is_supplementary else 0,
                ),
            )
        else:
            non_supp = [r for r in records_info if not r["read"].is_supplementary]
            chosen = non_supp[0] if non_supp else records_info[0]

        exon_lens, intron_lens = _exons_and_introns(chosen["ops"])

        junction_count = len(intron_lens)
        exon_count = len(exon_lens)
        exon_lengths = exon_lens
        intron_lengths = intron_lens

        chosen_read = chosen["read"]
        strand = "-" if chosen_read.is_reverse else "+"
        chromosome = chosen_read.reference_name
        genomic_start = chosen_read.reference_start
        genomic_end = chosen_read.reference_end
        cdna_start = chosen["q_start"]
        cdna_end = chosen["q_end"]

        # Circular origin wrap-around check
        for other in records_info:
            if other is chosen:
                continue
            other_read = other["read"]
            if (
                other_read.reference_name == chosen_read.reference_name
                and other_read.is_reverse == chosen_read.is_reverse
            ):
                start_p = None
                end_p = None
                if chosen["q_start"] <= 5 and other["q_end"] >= other["rlen"] - 6:
                    start_p, end_p = chosen, other
                elif other["q_start"] <= 5 and chosen["q_end"] >= chosen["rlen"] - 6:
                    start_p, end_p = other, chosen

                if start_p is not None and end_p is not None:
                    g1_s, g1_e = chosen_read.reference_start, chosen_read.reference_end
                    g2_s, g2_e = other_read.reference_start, other_read.reference_end
                    if (
                        g1_s is not None
                        and g1_e is not None
                        and g2_s is not None
                        and g2_e is not None
                    ):
                        genomic_gap = max(g1_s, g2_s) - min(g1_e, g2_e)
                        if -50 <= genomic_gap <= 50_000:
                            cdna_start = end_p["q_start"]
                            cdna_end = start_p["q_end"]
                            genomic_start = min(g1_s, g2_s)
                            genomic_end = max(g1_e, g2_e)
                            first_p, second_p = (chosen, other) if g1_s <= g2_s else (other, chosen)
                            ex1, in1 = _exons_and_introns(first_p["ops"])
                            ex2, in2 = _exons_and_introns(second_p["ops"])
                            if genomic_gap <= 5 and ex1 and ex2:
                                exon_lengths = [*ex1[:-1], ex1[-1] + ex2[0], *ex2[1:]]
                                intron_lengths = [*in1, *in2]
                            else:
                                exon_lengths = [*ex1, *ex2]
                                intron_lengths = [*in1, max(0, genomic_gap), *in2]
                            exon_count = len(exon_lengths)
                            junction_count = len(intron_lengths)
                            break

    # Write the output file with all fields (use NA for missing values)
    with open(output_file, "w") as f:
        if cdna_start is not None and cdna_end is not None:
            insert_region: tuple[int, int] = (cdna_start, cdna_end)
            f.write(f"cDNA start position in plasmid: {cdna_start}\n")
            f.write(f"cDNA end position in plasmid: {cdna_end}\n")
            f.write(f"INSERT_REGION: {insert_region}\n")
        else:
            logging.warning("No alignments found. Writing NA values for cDNA positions.")
            f.write("cDNA start position in plasmid: NA\n")
            f.write("cDNA end position in plasmid: NA\n")
            f.write("INSERT_REGION: NA\n")

        # Write junction metadata fields
        f.write(f"SpliceJunctionCount: {junction_count if junction_count is not None else 'NA'}\n")
        f.write(f"ExonCount: {exon_count if exon_count is not None else 'NA'}\n")

        if exon_lengths is not None and len(exon_lengths) > 0:
            f.write(f"ExonLengths: {','.join(map(str, exon_lengths))}\n")
        else:
            f.write("ExonLengths: NA\n")

        if intron_lengths is not None and len(intron_lengths) > 0:
            f.write(f"IntronLengths: {','.join(map(str, intron_lengths))}\n")
        else:
            f.write("IntronLengths: NA\n")

        f.write(f"Strand: {strand if strand is not None else 'NA'}\n")
        f.write(f"Chromosome: {chromosome if chromosome is not None else 'NA'}\n")
        f.write(f"GenomicStart: {genomic_start if genomic_start is not None else 'NA'}\n")
        f.write(f"GenomicEnd: {genomic_end if genomic_end is not None else 'NA'}\n")

    logging.info(f"cDNA positions and INSERT_REGION saved to {output_file}")
