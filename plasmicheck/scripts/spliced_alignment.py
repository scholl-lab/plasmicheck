from __future__ import annotations

import functools
import gzip
import itertools
import logging
import os
import re
from shlex import quote
from typing import Any

import pysam

from plasmicheck.config import get_config

from .utils import run_command

_cfg = get_config()
PADDING_DEFAULT: int = _cfg["padding"]
DEFAULT_PADDING: int = PADDING_DEFAULT
FASTA_EXTENSIONS: list[str] = _cfg["alignment"]["fasta_extensions"]
MINIMAP2_THREADS: int = _cfg["alignment"]["minimap2_threads"]
SAMTOOLS_THREADS: int = _cfg["alignment"]["samtools_threads"]

_TX_RE = re.compile(r'transcript_id[=\s]+"?([^";\s]+)"?')


def _normalize_chrom(chrom: str) -> str:
    c = chrom.strip()
    if c.upper().startswith("NC_"):
        base = c.split(".")[0].upper()
        if base[3:].isdigit():
            num = int(base[3:])
            if 1 <= num <= 22:
                return str(num)
            if num == 23:
                return "X"
            if num == 24:
                return "Y"
            if num == 12920:
                return "MT"
    if c.lower().startswith("chr"):
        c = c[3:]
    return c.upper() if c.upper() in ("X", "Y", "M", "MT") else c


@functools.lru_cache(maxsize=4)
def _load_gtf_transcripts(
    gtf_path: str,
) -> dict[str, list[tuple[str, int, int, int, list[tuple[int, int]]]]]:
    open_fn: Any = gzip.open if gtf_path.endswith(".gz") else open
    tx_exons: dict[tuple[str, str], list[tuple[int, int]]] = {}
    tx_rank: dict[tuple[str, str], int] = {}
    with open_fn(gtf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("exon", "transcript"):
                continue
            m = _TX_RE.search(fields[8])
            if not m:
                continue
            tx_id = m.group(1)
            norm_chr = _normalize_chrom(fields[0])
            key = (norm_chr, tx_id)
            attrs = fields[8]
            rank = 2 if "MANE_Select" in attrs else (1 if "MANE_Plus_Clinical" in attrs else 0)
            if rank > tx_rank.get(key, 0):
                tx_rank[key] = rank
            if fields[2] == "exon":
                tx_exons.setdefault(key, []).append((int(fields[3]), int(fields[4])))

    by_chrom: dict[str, list[tuple[str, int, int, int, list[tuple[int, int]]]]] = {}
    for (norm_chr, tx_id), raw_exons in tx_exons.items():
        if not raw_exons:
            continue
        exons = sorted(set(raw_exons))
        by_chrom.setdefault(norm_chr, []).append(
            (tx_id, tx_rank.get((norm_chr, tx_id), 0), exons[0][0], exons[-1][1], exons)
        )
    return by_chrom


def _cigar_geometry(
    alignments: list[Any],
) -> tuple[list[tuple[int, int]], set[tuple[int, int]], list[tuple[int, int]]]:
    ref_blocks: list[tuple[int, int]] = []
    introns: set[tuple[int, int]] = set()
    indels: list[tuple[int, int]] = []
    for read in alignments:
        if read.reference_start is None:
            continue
        rpos = int(read.reference_start)
        ops = read.cigartuples or []
        for idx, (op, length) in enumerate(ops):
            if op in (0, 2, 7, 8):
                ref_blocks.append((rpos, rpos + length))
                if op == 2:
                    indels.append((rpos + 1, -length))
                rpos += length
            elif op == 3:
                introns.add((rpos + 1, rpos + length))
                rpos += length
            elif op == 1:
                prev_exon = idx > 0 and ops[idx - 1][0] in (0, 2, 7, 8)
                next_exon = idx + 1 < len(ops) and ops[idx + 1][0] in (0, 2, 7, 8)
                if prev_exon and next_exon:
                    indels.append((rpos + 1, length))
    ref_blocks.sort()
    return ref_blocks, introns, indels


def _resolve_gtf_exons(
    gtf_path: str, chrom: str, alignments: list[Any]
) -> tuple[list[int], list[int], list[tuple[int, int]]] | None:
    ref_blocks, bam_introns, indel_events = _cigar_geometry(alignments)
    if not ref_blocks:
        return None
    aln_s1 = ref_blocks[0][0] + 1
    aln_e1 = max(b[1] for b in ref_blocks)
    candidates = _load_gtf_transcripts(gtf_path).get(_normalize_chrom(chrom), [])
    overlapping = [tx for tx in candidates if tx[3] >= aln_s1 and tx[2] <= aln_e1]
    if not overlapping:
        return None

    def _rank_tx(
        tx: tuple[str, int, int, int, list[tuple[int, int]]],
    ) -> tuple[bool, int, int, int, int]:
        _, mane_rank, tx_min, tx_max, exons = tx
        tx_introns = {(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)}
        exact_matches = len(bam_introns & tx_introns)
        sym_diff = len(bam_introns ^ tx_introns)
        overlap_bp = max(0, min(aln_e1, tx_max) - max(aln_s1, tx_min) + 1)
        return (mane_rank > 0, exact_matches, -sym_diff, mane_rank, overlap_bp)

    _, _, tx_min, tx_max, best_exons = max(overlapping, key=_rank_tx)
    valid_blocks = [b for b in ref_blocks if b[1] > tx_min - 1 and b[0] < tx_max]
    if not valid_blocks:
        return None
    eff_s1 = valid_blocks[0][0] + 1
    eff_e1 = valid_blocks[-1][1]

    clipped_exons: list[tuple[int, int]] = []
    gtf_exon_lengths: list[int] = []
    for s, e in best_exons:
        if e < eff_s1 or s > eff_e1:
            continue
        cs = max(s, eff_s1)
        ce = min(e, eff_e1)
        adj = sum(delta for pos, delta in indel_events if cs <= pos <= ce)
        clipped_exons.append((cs, ce))
        gtf_exon_lengths.append((ce - cs + 1) + adj)

    if not clipped_exons:
        return None
    gtf_introns = [
        (clipped_exons[i][1] + 1, clipped_exons[i + 1][0] - 1)
        for i in range(len(clipped_exons) - 1)
    ]
    gtf_intron_lengths = [e - s + 1 for s, e in gtf_introns]
    return gtf_exon_lengths, gtf_intron_lengths, gtf_introns


def _extract_introns_and_variants(
    alignments: list[Any],
) -> tuple[list[tuple[int, int]], list[str]]:
    introns: set[tuple[int, int]] = set()
    raw_variants: list[tuple[int, str, str]] = []

    for read in alignments:
        if read.reference_start is None:
            continue
        rpos = int(read.reference_start)
        for op, length in read.cigartuples or []:
            if op == 3:
                introns.add((rpos + 1, rpos + length))
                rpos += length
            elif op in (0, 2, 7, 8):
                rpos += length

        qseq = getattr(read, "query_sequence", None)
        if isinstance(qseq, str) and qseq and hasattr(read, "has_tag") and read.has_tag("MD"):
            mismatches: list[tuple[int, str, str]] = []
            for qpos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
                if qpos is None or ref_pos is None or ref_base is None:
                    continue
                alt_base = qseq[qpos].upper()
                ref_upper = ref_base.upper()
                if (
                    (ref_base.islower() or ref_upper != alt_base)
                    and ref_upper != "N"
                    and alt_base != "N"
                ):
                    mismatches.append((ref_pos + 1, ref_upper, alt_base))
            mismatches.sort(key=lambda x: x[0])
            merged: list[tuple[int, str, str]] = []
            for pos, r_b, a_b in mismatches:
                if merged and pos == merged[-1][0] + len(merged[-1][1]):
                    prev_pos, prev_r, prev_a = merged[-1]
                    merged[-1] = (prev_pos, prev_r + r_b, prev_a + a_b)
                else:
                    merged.append((pos, r_b, a_b))
            raw_variants.extend(merged)

    raw_variants.sort(key=lambda x: (x[0], x[1], x[2]))
    dedup_variants = list(dict.fromkeys(f"{pos}:{r_b}>{a_b}" for pos, r_b, a_b in raw_variants))
    return sorted(introns), dedup_variants


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
    mm_threads = threads_minimap2 if threads_minimap2 is not None else MINIMAP2_THREADS
    st_threads = threads_samtools if threads_samtools is not None else SAMTOOLS_THREADS
    logging.info(
        f"Performing spliced alignment with {plasmid_fasta} to {human_index} and writing to {output_bam}"
    )
    try:
        run_command(
            f"minimap2 -t {quote(str(mm_threads))} -ax splice --MD "
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


def extract_plasmid_cdna_positions(
    plasmid_fasta: str,
    bam_file: str,
    output_file: str,
    padding: int | str | None = PADDING_DEFAULT,
    gtf_file: str | None = None,
) -> None:
    if isinstance(padding, str) or padding is None:
        if gtf_file is None and isinstance(padding, str):
            gtf_file = padding
        padding = PADDING_DEFAULT

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
            if op in (0, 1, 7, 8):
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
    genomic_introns: list[tuple[int, int]] = []
    construct_variants: list[str] = []
    map_exon_lengths: list[int] = []
    disagreed_gtf_exons: list[int] | None = None

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
        active_reads = [chosen_read]

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
                            active_reads = [first_p["read"], second_p["read"]]
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

        genomic_introns, construct_variants = _extract_introns_and_variants(active_reads)
        map_exon_lengths = list(exon_lengths) if exon_lengths is not None else []

        if gtf_file is not None and chromosome is not None:
            gtf_res = _resolve_gtf_exons(gtf_file, chromosome, active_reads)
            if gtf_res is not None:
                gtf_exons, gtf_intron_lens, gtf_intron_coords = gtf_res
                genomic_introns = gtf_intron_coords
                if gtf_exons != map_exon_lengths:
                    logging.warning(
                        f"GTF exon lengths {gtf_exons} differ from map-derived exon lengths "
                        f"{map_exon_lengths} for {plasmid_fasta}; overriding with GTF."
                    )
                    disagreed_gtf_exons = gtf_exons
                    exon_lengths = gtf_exons
                    intron_lengths = gtf_intron_lens
                    exon_count = len(gtf_exons)
                    junction_count = len(gtf_intron_lens)
                    if gtf_intron_coords:
                        genomic_start = gtf_intron_coords[0][0] - 1 - gtf_exons[0]
                        genomic_end = gtf_intron_coords[-1][1] + gtf_exons[-1]

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

        if disagreed_gtf_exons is not None:
            f.write(
                f"MapExonLengths: {','.join(map(str, map_exon_lengths)) if map_exon_lengths else 'NA'}\n"
            )
            f.write(f"GTFExonLengths: {','.join(map(str, disagreed_gtf_exons))}\n")

        if gtf_file is not None or genomic_introns or construct_variants:
            f.write(
                f"GenomicIntrons: {','.join(f'{s}-{e}' for s, e in genomic_introns) if genomic_introns else 'NA'}\n"
            )
            f.write(
                f"ConstructVariants: {','.join(construct_variants) if construct_variants else 'NA'}\n"
            )

    logging.info(f"cDNA positions and INSERT_REGION saved to {output_file}")
