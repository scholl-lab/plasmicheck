from __future__ import annotations

import contextlib
import hashlib
import itertools
import logging
import os
import secrets
import subprocess
from collections.abc import Iterator
from dataclasses import dataclass
from typing import IO, Any

import numpy as np
import pysam

from plasmicheck.config import get_config

from .coverage_metrics import (
    compute_coverage_stats,
    compute_region_coverage_metrics,
    get_depth_array,
)
from .resistance_gene_detection import DetectedGene, detect_resistance_genes, locate_genbank_file
from .utils import categorize_value

_cfg = get_config()
MATE_BONUS: int = _cfg["scoring"]["mate_bonus"]
CLIPPING_PENALTY: int = _cfg["scoring"]["clipping_penalty"]
MISMATCH_PENALTY: int = _cfg["scoring"]["mismatch_penalty"]
DEFAULT_THRESHOLD: float = _cfg["default_threshold"]
UNCLEAR_RANGE: dict[str, float] = _cfg["unclear_range"]
SAMTOOLS_THREADS: int = _cfg["alignment"]["samtools_threads"]
FILTER_BACKBONE_ONLY: bool = _cfg.get("filter_backbone_only", True)
SCORE_MARGIN: int = _cfg.get("score_margin", 0)
BREADTH_THRESHOLDS: list[int] = _cfg.get("breadth_thresholds", [5])
CUSTOM_RESISTANCE_PATTERNS: list[str] = _cfg.get("resistance_gene_patterns", [])


def _has_origin_sa_support(read: Any, plasmid_length: int, is_start_clip: bool) -> bool:
    """Check if read has supplementary alignment (SA tag) supporting circular origin crossing."""
    if not (read.has_tag("SA") and plasmid_length > 0):
        return False
    sa_tag = read.get_tag("SA")
    if not isinstance(sa_tag, str):
        return False
    read_rname = getattr(read, "reference_name", "") or ""
    read_strand = "-" if getattr(read, "is_reverse", False) else "+"
    for entry in sa_tag.strip(";").split(";"):
        if not entry:
            continue
        parts = entry.split(",")
        if len(parts) < 6:
            continue
        sa_rname, sa_pos_str, sa_strand, _sa_cigar, _sa_mapq, sa_nm_str = parts[:6]
        try:
            sa_pos = int(sa_pos_str)
            sa_nm = int(sa_nm_str)
        except ValueError:
            continue
        if sa_nm > 2:
            continue
        if read_rname and sa_rname != read_rname:
            continue
        if sa_strand != read_strand:
            continue
        # In 1-based SAM coordinates:
        # If primary read has clip at start (r_start == 0), supplementary alignment should be
        # at the opposite end of the plasmid (near plasmid_length).
        # If primary read has clip at end (r_end == plasmid_length), supplementary alignment
        # should be at coordinate 1 (0-based 0).
        if is_start_clip:
            if sa_pos >= max(1, plasmid_length - 1000):
                return True
        else:
            if sa_pos <= 1000:
                return True
    return False


def calculate_alignment_score(read: Any, plasmid_length: int = 0) -> int:
    """Calculate a custom alignment score based on multiple criteria."""
    if getattr(read, "is_unmapped", False):
        return 0

    ops = read.cigartuples or []
    r_start = read.reference_start
    r_end = read.reference_end
    score: int = read.mapping_quality

    # Penalize for clipping (excusing circular origin boundary clips on plasmids only when supported by SA)
    for idx, (operation, length) in enumerate(ops):
        if operation in {4, 5}:  # 4: soft clipping, 5: hard clipping
            if plasmid_length > 0:
                is_start_clip = idx == 0 and r_start == 0
                is_end_clip = idx == len(ops) - 1 and r_end == plasmid_length
                if (is_start_clip or is_end_clip) and _has_origin_sa_support(
                    read, plasmid_length, is_start_clip
                ):
                    continue
            score -= length * CLIPPING_PENALTY

    # Penalize for mismatches
    mismatches = read.get_tag("NM") if read.has_tag("NM") else 0
    score -= mismatches * MISMATCH_PENALTY

    return score


@dataclass
class InsertRegionData:
    """Parsed cDNA position and junction metadata."""

    start: int
    end: int
    junction_count: int | None = None
    exon_count: int | None = None
    exon_lengths: list[int] | None = None
    intron_lengths: list[int] | None = None
    strand: str | None = None
    chromosome: str | None = None
    genomic_start: int | None = None
    genomic_end: int | None = None

    def __iter__(self) -> Iterator[int]:
        """Support tuple unpacking: start, end = parse_insert_region(file)"""
        return iter((self.start, self.end))

    def __getitem__(self, index: int) -> int:
        """Support tuple-like subscripting: insert_region[0], insert_region[1]"""
        if index == 0:
            return self.start
        elif index == 1:
            return self.end
        else:
            raise IndexError(f"InsertRegionData index out of range: {index}")


def parse_insert_region(cdna_positions_file: str) -> InsertRegionData:
    """Parse the cDNA_positions.txt file to extract the INSERT_REGION and junction metadata."""
    if not os.path.exists(cdna_positions_file):
        raise FileNotFoundError(
            f"Expected cDNA_positions.txt at {cdna_positions_file} but the file was not found."
        )

    logging.debug(f"Looking for cDNA_positions.txt at: {cdna_positions_file}")

    with open(cdna_positions_file) as f:
        lines = [line.strip() for line in f.readlines()]

    def _get_field(prefix: str) -> str | None:
        """Helper to extract field value by prefix."""
        for line in lines:
            if line.startswith(prefix):
                value = line.split(": ", 1)[1].strip()
                return None if value == "NA" else value
        return None

    # Parse cDNA start/end (required fields from lines 0-1)
    start_str = _get_field("cDNA start position in plasmid:")
    end_str = _get_field("cDNA end position in plasmid:")

    # NA values mean no spliced alignment was found — insert region is unavailable.
    # Raise ValueError so callers fall back to no-filtering behavior rather than
    # using a bogus (0, 0) region that would misclassify reads as Backbone_Only.
    if start_str is None or end_str is None:
        raise ValueError(
            "Insert region unavailable (cDNA start/end are NA). "
            "Backbone filtering will be disabled."
        )
    cdna_start = int(start_str)
    cdna_end = int(end_str)

    # Parse junction metadata (optional fields, may not exist in old files)
    junction_count_str = _get_field("SpliceJunctionCount:")
    junction_count = int(junction_count_str) if junction_count_str is not None else None

    exon_count_str = _get_field("ExonCount:")
    exon_count = int(exon_count_str) if exon_count_str is not None else None

    exon_lengths_str = _get_field("ExonLengths:")
    exon_lengths = (
        [int(x) for x in exon_lengths_str.split(",")] if exon_lengths_str is not None else None
    )

    intron_lengths_str = _get_field("IntronLengths:")
    intron_lengths = (
        [int(x) for x in intron_lengths_str.split(",")] if intron_lengths_str is not None else None
    )

    strand = _get_field("Strand:")
    chromosome = _get_field("Chromosome:")

    genomic_start_str = _get_field("GenomicStart:")
    genomic_start = int(genomic_start_str) if genomic_start_str is not None else None

    genomic_end_str = _get_field("GenomicEnd:")
    genomic_end = int(genomic_end_str) if genomic_end_str is not None else None

    return InsertRegionData(
        start=cdna_start,
        end=cdna_end,
        junction_count=junction_count,
        exon_count=exon_count,
        exon_lengths=exon_lengths,
        intron_lengths=intron_lengths,
        strand=strand,
        chromosome=chromosome,
        genomic_start=genomic_start,
        genomic_end=genomic_end,
    )


def read_overlaps_insert(read: Any, insert_region: tuple[int, int] | InsertRegionData) -> bool:
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
    if insert_region[0] > insert_region[1]:
        return bool(read_end - 1 >= insert_region[0] or read_start <= insert_region[1])
    # No overlap: read ends before insert starts, or read starts after insert ends
    return not (read_end <= insert_region[0] or read_start > insert_region[1])


def _compute_junction_positions(
    insert_region: tuple[int, int] | InsertRegionData, plasmid_length: int
) -> list[int]:
    """Compute 0-based plasmid coordinates of exon-exon splice junctions within the insert."""
    if not isinstance(insert_region, InsertRegionData) or not insert_region.exon_lengths:
        return []
    exons = list(insert_region.exon_lengths)
    if len(exons) <= 1:
        return []
    if insert_region.strand == "-":
        exons = exons[::-1]
    positions: list[int] = []
    pos = insert_region.start
    for exon_len in exons[:-1]:
        pos = (pos + exon_len) % plasmid_length if plasmid_length > 0 else pos + exon_len
        positions.append(pos)
    return positions


def _scan_plasmid_bam(
    plasmid_bam: str,
    insert_region: tuple[int, int] | InsertRegionData,
    window_size: int = 100,
    min_junction_anchor: int = 15,
) -> tuple[float, dict[str, int], int, int, int, int, int, float]:
    """Single-pass scan of plasmid_bam for outside-insert coverage, boundary mismatches, and splice-junction reads."""
    total_aligned_bases_outside_insert = 0
    mismatches_near_insert = {"with_mismatches_or_clipping": 0, "without_mismatches_or_clipping": 0}
    splice_junction_reads = 0
    splice_junction_events = 0
    covered_junctions: set[int] = set()
    exon_intron_boundary_reads = 0
    vector_insert_junction_reads = 0
    total_query_length = 0
    primary_read_count = 0

    with pysam.AlignmentFile(plasmid_bam, "rb") as bamfile:
        plasmid_length = 0
        if (
            hasattr(bamfile, "lengths")
            and isinstance(bamfile.lengths, (list, tuple))
            and len(bamfile.lengths) > 0
        ):
            plasmid_length = bamfile.lengths[0]
        else:
            plasmid_length = bamfile.lengths[bamfile.get_tid(bamfile.get_reference_name(0))]

        insert_start, insert_end = insert_region[0], insert_region[1]
        is_wrapped = insert_start > insert_end
        if is_wrapped:
            total_reference_bases_outside_insert = insert_start - insert_end - 1
        else:
            total_reference_bases_outside_insert = plasmid_length - (insert_end - insert_start + 1)
        if total_reference_bases_outside_insert <= 0:
            raise ValueError("The INSERT_REGION covers the entire reference sequence.")

        junc_positions = _compute_junction_positions(insert_region, plasmid_length)
        vi_anchor = 10
        boundaries: list[int] = [
            insert_start,
            (insert_end + 1) % plasmid_length if plasmid_length > 0 else insert_end + 1,
        ]

        for read in bamfile.fetch():
            if read.is_unmapped:
                continue

            read_start = read.reference_start
            read_end = read.reference_end
            if read_start is None or read_end is None:
                continue

            cigartuples = read.cigartuples

            # 1. Coverage outside insert (read_end is 0-based exclusive; insert_end is 0-based inclusive)
            if is_wrapped:
                # Backbone is half-open [insert_end + 1, insert_start)
                outside = max(0, min(read_end, insert_start) - max(read_start, insert_end + 1))
                inside = (read_end - read_start) - outside
                total_aligned_bases_outside_insert += outside
            else:
                # Insert is half-open [insert_start, insert_end + 1)
                inside = max(0, min(read_end, insert_end + 1) - max(read_start, insert_start))
                outside = max(0, (read_end - read_start) - inside)
                total_aligned_bases_outside_insert += outside

            is_primary_usable = not (
                read.is_secondary
                or read.is_supplementary
                or getattr(read, "is_duplicate", False) is True
                or getattr(read, "is_qcfail", False) is True
            )
            if is_primary_usable:
                qlen = getattr(read, "query_length", None)
                if qlen is None or qlen <= 0:
                    qlen = read_end - read_start
                if qlen > 0:
                    total_query_length += qlen
                    primary_read_count += 1

            # 2. Vector-Insert Junction reads (NM <= 2, anchor >= 10 on both sides)
            if (
                is_primary_usable
                and cigartuples is not None
                and (int(read.get_tag("NM")) if read.has_tag("NM") else 0) <= 2
                and plasmid_length > 0
            ):
                for b in boundaries:
                    if vi_anchor <= b <= plasmid_length - vi_anchor:
                        if hasattr(read, "get_overlap") and callable(read.get_overlap):
                            if (
                                read.get_overlap(b - vi_anchor, b + vi_anchor) == 2 * vi_anchor
                                and read_start <= b - vi_anchor
                                and read_end >= b + vi_anchor
                            ):
                                vector_insert_junction_reads += 1
                                break
                        elif read_start <= b - vi_anchor and read_end >= b + vi_anchor:
                            vector_insert_junction_reads += 1
                            break

            # 3. Mismatches near insert boundaries (nearest boundary check)
            near_boundary = False
            for b in boundaries:
                if read_start <= b <= read_end:
                    dist = 0
                else:
                    d_s = abs(read_start - b)
                    d_e = abs(read_end - b)
                    if plasmid_length > 0:
                        d_s = min(d_s % plasmid_length, plasmid_length - (d_s % plasmid_length))
                        d_e = min(d_e % plasmid_length, plasmid_length - (d_e % plasmid_length))
                    dist = min(d_s, d_e)
                if dist <= window_size:
                    near_boundary = True
                    break

            if near_boundary and cigartuples is not None:
                has_mismatch_or_clip = any(op in {1, 2, 3, 4, 5} for op, _ in cigartuples) or (
                    read.has_tag("NM") and int(read.get_tag("NM")) > 0
                )
                if has_mismatch_or_clip:
                    mismatches_near_insert["with_mismatches_or_clipping"] += 1
                else:
                    mismatches_near_insert["without_mismatches_or_clipping"] += 1

            # 4. Direct exon-exon splice-junction spanning reads and events
            if (
                junc_positions
                and is_primary_usable
                and cigartuples is not None
                and not any(op in {4, 5} for op, ln in cigartuples if ln > 2)
                and (int(read.get_tag("NM")) if read.has_tag("NM") else 0) <= 2
            ):
                spanned = False
                for jp in junc_positions:
                    if read_start + min_junction_anchor <= jp <= read_end - min_junction_anchor:
                        spanned = True
                        covered_junctions.add(jp)
                        if hasattr(read, "get_overlap") and callable(read.get_overlap):
                            if (
                                read.get_overlap(jp - min_junction_anchor, jp + min_junction_anchor)
                                == 2 * min_junction_anchor
                            ):
                                splice_junction_events += 1
                        else:
                            splice_junction_events += 1
                if spanned:
                    splice_junction_reads += 1

            # 5. Exon-intron boundary clipped reads
            if junc_positions and is_primary_usable and cigartuples is not None:
                leading_clip = 0
                if cigartuples[0][0] == 4:
                    leading_clip = cigartuples[0][1]
                elif len(cigartuples) > 1 and cigartuples[0][0] == 5 and cigartuples[1][0] == 4:
                    leading_clip = cigartuples[1][1]

                trailing_clip = 0
                if cigartuples[-1][0] == 4:
                    trailing_clip = cigartuples[-1][1]
                elif len(cigartuples) > 1 and cigartuples[-1][0] == 5 and cigartuples[-2][0] == 4:
                    trailing_clip = cigartuples[-2][1]

                has_overlap_fn = hasattr(read, "get_overlap") and callable(read.get_overlap)
                if leading_clip >= min_junction_anchor:
                    for jp in junc_positions:
                        anchor_start = max(read_start, jp)
                        if abs(read_start - jp) <= 2 and (
                            read_end - anchor_start >= min_junction_anchor
                            and (
                                not has_overlap_fn
                                or read.get_overlap(
                                    anchor_start, anchor_start + min_junction_anchor
                                )
                                == min_junction_anchor
                            )
                        ):
                            exon_intron_boundary_reads += 1
                            break

                if trailing_clip >= min_junction_anchor:
                    for jp in junc_positions:
                        anchor_end = min(read_end, jp)
                        if abs(read_end - jp) <= 2 and (
                            anchor_end - read_start >= min_junction_anchor
                            and (
                                not has_overlap_fn
                                or read.get_overlap(anchor_end - min_junction_anchor, anchor_end)
                                == min_junction_anchor
                            )
                        ):
                            exon_intron_boundary_reads += 1
                            break

    mean_read_length = total_query_length / primary_read_count if primary_read_count > 0 else 0.0
    coverage_outside_insert: float = (
        total_aligned_bases_outside_insert / total_reference_bases_outside_insert
    )
    return (
        coverage_outside_insert,
        mismatches_near_insert,
        splice_junction_reads,
        len(covered_junctions),
        splice_junction_events,
        exon_intron_boundary_reads,
        vector_insert_junction_reads,
        mean_read_length,
    )


def calculate_coverage_outside_insert(
    plasmid_bam: str, insert_region: tuple[int, int] | InsertRegionData
) -> float:
    """Calculate the coverage of reads mapping outside the INSERT_REGION."""
    return _scan_plasmid_bam(plasmid_bam, insert_region)[0]


def count_mismatches_near_insert_end(
    plasmid_bam: str, insert_region: tuple[int, int] | InsertRegionData, window_size: int = 100
) -> dict[str, int]:
    """Count reads with/without mismatches or clipping near the INSERT_REGION boundaries."""
    return _scan_plasmid_bam(plasmid_bam, insert_region, window_size=window_size)[1]


# ---------------------------------------------------------------------------
# Streaming name-sorted merge (O(1) memory)
# ---------------------------------------------------------------------------


def _namesort_bam(input_bam: str, output_bam: str, threads: int = SAMTOOLS_THREADS) -> None:
    """Sort names lexicographically to match the streaming merge comparison."""
    subprocess.run(
        ["samtools", "sort", "-N", "-@", str(threads), "-o", output_bam, input_bam],
        check=True,
    )


def _iter_reads_by_name(bam_path: str) -> Iterator[tuple[str, list[Any]]]:
    """Yield ``(query_name, [reads])`` groups from a name-sorted BAM."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for name, group in itertools.groupby(bam.fetch(until_eof=True), key=lambda r: r.query_name):
            yield name, list(group)


def _mate_key(read: Any) -> int | None:
    """Return 0 for unpaired, 1 for read1, 2 for read2, or None for invalid/secondary/supplementary."""
    if (
        getattr(read, "is_unmapped", False)
        or getattr(read, "is_secondary", False)
        or getattr(read, "is_supplementary", False)
    ):
        return None
    if not getattr(read, "is_paired", False):
        return 0
    is_r1 = bool(getattr(read, "is_read1", False))
    is_r2 = bool(getattr(read, "is_read2", False))
    if is_r1 == is_r2:
        return None
    return 1 if is_r1 else 2


def _is_pair(*args: Any) -> bool:
    """Return True if arguments represent a valid mapped paired-end pair (read1 and read2)."""
    if len(args) == 1 and isinstance(args[0], (list, tuple)):
        reads = list(args[0])
    elif len(args) == 2:
        reads = list(args)
    else:
        return False
    return len(reads) == 2 and {_mate_key(r) for r in reads} == {1, 2}


def _primary_reads(reads: list[Any]) -> list[Any]:
    """Return primary records sorted by read1 then read2, or first record as fallback."""
    primary = [
        r
        for r in reads
        if not getattr(r, "is_secondary", False) and not getattr(r, "is_supplementary", False)
    ] or reads[:1]
    if _is_pair(primary):
        primary.sort(
            key=lambda r: (not getattr(r, "is_read1", False), getattr(r, "is_read2", False))
        )
    return primary


def _best_read(reads: list[Any]) -> Any:
    """Pick the primary alignment (or first) from a group of reads with the same name."""
    return _primary_reads(reads)[0]


def _scores(
    reads: list[Any],
    length: int = 0,
    span_origin: bool = False,
    mate_bonus: int = MATE_BONUS,
    plasmid_length: int | None = None,
) -> list[int]:
    """Calculate alignment scores for a read group, awarding mate bonus for concordant pairs."""
    plen = plasmid_length if plasmid_length is not None else length
    scores = [calculate_alignment_score(r, plen) for r in reads]
    if not _is_pair(reads):
        return scores

    a, b = reads[0], reads[1]
    if (
        getattr(a, "mate_is_unmapped", False)
        or getattr(b, "mate_is_unmapped", False)
        or a.reference_end is None
        or b.reference_end is None
        or a.reference_start is None
        or b.reference_start is None
        or a.is_reverse == b.is_reverse
    ):
        return scores

    if (
        hasattr(a, "reference_id")
        and hasattr(b, "reference_id")
        and hasattr(a, "next_reference_id")
        and hasattr(b, "next_reference_id")
        and (
            a.reference_id != b.reference_id
            or a.next_reference_id != b.reference_id
            or b.next_reference_id != a.reference_id
        )
    ):
        return scores

    if (
        hasattr(a, "next_reference_start")
        and hasattr(b, "next_reference_start")
        and a.next_reference_start is not None
        and b.next_reference_start is not None
        and (
            a.next_reference_start != b.reference_start
            or b.next_reference_start != a.reference_start
        )
    ):
        return scores

    forward, reverse = (b, a) if a.is_reverse else (a, b)
    shift = plen if (plen > 0 and forward.reference_start > reverse.reference_start) else 0
    if forward.reference_start > reverse.reference_start and not shift:
        return scores

    span = max(forward.reference_end, reverse.reference_end + shift) - forward.reference_start
    if 0 < span <= 1000:
        return [score + mate_bonus for score in scores]
    return scores


def _is_strictly_inside_insert(
    read: Any, insert_region: tuple[int, int] | InsertRegionData, plasmid_length: int = 0
) -> bool:
    """Return True if read is completely within the insert region."""
    if (
        getattr(read, "is_unmapped", False)
        or read.reference_start is None
        or read.reference_end is None
    ):
        return False
    ins_start, ins_end = insert_region[0], insert_region[1]
    r_start, r_end = read.reference_start, read.reference_end
    if ins_start > ins_end:
        overlaps_backbone = not (r_end <= ins_end + 1 or r_start >= ins_start)
        return bool(not overlaps_backbone and read_overlaps_insert(read, insert_region))
    else:
        return bool(r_start >= ins_start and r_end <= ins_end + 1)


def _is_concordant_bridge_pair(
    p_reads: list[Any],
    insert_region: tuple[int, int] | InsertRegionData,
    plasmid_length: int,
) -> bool:
    """Return True if pair has one mate strictly in insert and other strictly in backbone, with span <= 1000."""
    if not _is_pair(p_reads) or plasmid_length <= 0:
        return False
    a, b = p_reads[0], p_reads[1]
    if getattr(a, "is_unmapped", False) or getattr(b, "is_unmapped", False):
        return False
    if (
        a.reference_end is None
        or b.reference_end is None
        or a.reference_start is None
        or b.reference_start is None
    ):
        return False
    if a.is_reverse == b.is_reverse:
        return False

    forward, reverse = (b, a) if a.is_reverse else (a, b)
    fragment_length = (reverse.reference_end - forward.reference_start) % plasmid_length
    f_len = forward.reference_end - forward.reference_start
    r_len = reverse.reference_end - reverse.reference_start
    if not (max(f_len, r_len) <= fragment_length <= 1000):
        return False

    inside_a = _is_strictly_inside_insert(a, insert_region, plasmid_length)
    outside_a = not read_overlaps_insert(a, insert_region)
    inside_b = _is_strictly_inside_insert(b, insert_region, plasmid_length)
    outside_b = not read_overlaps_insert(b, insert_region)
    return bool((inside_a and outside_b) or (inside_b and outside_a))


class StreamingCompareResult(dict[str, int]):
    """Dict subclass holding assigned read category counts and bridge pair tallies."""

    bridge_pairs: int = 0


def _read_id(query_name: str, key: bytes | None) -> str:
    """Keyed hash for the TSV only; the merge still compares the real names."""
    if key is None:
        return query_name
    return hashlib.blake2b(query_name.encode(), digest_size=8, key=key).hexdigest()


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
    insert_region: tuple[int, int] | InsertRegionData | None = None,
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
    insert_region: tuple[int, int] | InsertRegionData | None = None,
    score_margin: int = 0,
    plasmid_length: int = 0,
    anonymize_key: bytes | None = None,
) -> StreamingCompareResult:
    """Two-pointer merge over name-sorted BAMs.  Returns assigned counts."""
    assigned_counts = StreamingCompareResult(
        {
            "Plasmid": 0,
            "Human": 0,
            "Tied": 0,
            "Backbone_Only": 0,
            "Ambiguous": 0,
        }
    )

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
            p_reads = _primary_reads(p_item[1])
            h_reads = _primary_reads(h_item[1])
            p_read = p_reads[0]
            h_read = h_reads[0]
            p_scores = _scores(p_reads, length=plasmid_length, span_origin=True)
            h_scores = _scores(h_reads, length=0)
            ps = sum(p_scores)
            hs = sum(h_scores)
            if (
                ps != hs
                and len(p_reads) == len(h_reads)
                and len(p_reads) in {1, 2}
                and (len(p_reads) == 1 or (_is_pair(p_reads) and _is_pair(h_reads)))
                and all(
                    _mate_key(pr) is not None
                    and _mate_key(pr) == _mate_key(hr)
                    and pr.cigartuples
                    and hr.cigartuples
                    and pr.cigartuples
                    == (hr.cigartuples[::-1] if pr.is_reverse != hr.is_reverse else hr.cigartuples)
                    and (pr.get_tag("NM") if pr.has_tag("NM") else 0)
                    == (hr.get_tag("NM") if hr.has_tag("NM") else 0)
                    and 40 <= pr.mapping_quality <= 60
                    and 40 <= hr.mapping_quality <= 60
                    and p_score - pr.mapping_quality == h_score - hr.mapping_quality
                    for pr, hr, p_score, h_score in zip(
                        p_reads, h_reads, p_scores, h_scores, strict=True
                    )
                )
            ):
                ps = hs = max(ps, hs)
            rep_p_read = next(
                (
                    r
                    for r in p_reads
                    if insert_region is not None and read_overlaps_insert(r, insert_region)
                ),
                p_read,
            )
            assigned = _assign(
                ps,
                hs,
                score_margin=score_margin,
                plasmid_read=rep_p_read,
                insert_region=insert_region,
            )
            assigned_counts[assigned] += 1
            if (
                assigned == "Plasmid"
                and insert_region is not None
                and plasmid_length > 0
                and _is_concordant_bridge_pair(p_reads, insert_region, plasmid_length)
            ):
                assigned_counts.bridge_pairs += 1

            _write_assignment(
                outfile, _read_id(p_name, anonymize_key), assigned, ps, hs, p_read, h_read
            )
            p_item = next(plasmid_iter, None)
            h_item = next(human_iter, None)

        elif h_name is None or (p_name is not None and p_name < h_name):
            # Plasmid-only read
            assert p_item is not None
            assert p_name is not None
            p_reads = _primary_reads(p_item[1])
            p_read = p_reads[0]
            ps = sum(_scores(p_reads, length=plasmid_length, span_origin=True))
            rep_p_read = next(
                (
                    r
                    for r in p_reads
                    if insert_region is not None and read_overlaps_insert(r, insert_region)
                ),
                p_read,
            )
            assigned = _assign(
                max(ps, 0),
                0,
                score_margin=score_margin,
                plasmid_read=rep_p_read,
                insert_region=insert_region,
            )
            assigned_counts[assigned] += 1
            if (
                assigned == "Plasmid"
                and insert_region is not None
                and plasmid_length > 0
                and _is_concordant_bridge_pair(p_reads, insert_region, plasmid_length)
            ):
                assigned_counts.bridge_pairs += 1

            _write_assignment(
                outfile, _read_id(p_name, anonymize_key), assigned, ps, 0, p_read, None
            )
            p_item = next(plasmid_iter, None)

        else:
            # Human-only read
            assert h_item is not None
            h_reads = _primary_reads(h_item[1])
            h_read = h_reads[0]
            hs = sum(_scores(h_reads, length=0))
            assigned = _assign(0, max(hs, 0), score_margin=score_margin)
            assigned_counts[assigned] += 1
            _write_assignment(
                outfile, _read_id(h_name, anonymize_key), assigned, 0, hs, None, h_read
            )
            h_item = next(human_iter, None)

    return assigned_counts


def compare_alignments(
    plasmid_bam: str,
    human_bam: str,
    output_basename: str,
    threshold: float = DEFAULT_THRESHOLD,
    genbank_file: str | None = None,
    cdna_positions_file: str | None = None,
    filter_backbone_only: bool = FILTER_BACKBONE_ONLY,
    threads: int = SAMTOOLS_THREADS,
    anonymize: bool = False,
    anonymize_key: bytes | None = None,
) -> None:
    logging.info("Starting comparison of alignments")
    logging.debug(f"Processing BAM files: {plasmid_bam} (plasmid), {human_bam} (human)")

    # Parse insert region with graceful fallback
    if cdna_positions_file is None:
        # Backward compatibility: derive from output_basename directory
        cdna_positions_file = os.path.join(os.path.dirname(output_basename), "cDNA_positions.txt")
    insert_region: InsertRegionData | None = None
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
    splice_junction_reads = 0
    splice_junctions_covered = 0
    splice_junction_events = 0
    exon_intron_boundary_reads = 0
    vector_insert_junction_reads = 0
    mean_read_length = 0.0
    if insert_region is not None:
        (
            coverage_outside_insert,
            mismatches_near_insert,
            splice_junction_reads,
            splice_junctions_covered,
            splice_junction_events,
            exon_intron_boundary_reads,
            vector_insert_junction_reads,
            mean_read_length,
        ) = _scan_plasmid_bam(plasmid_bam, insert_region)
        logging.debug(f"Coverage outside INSERT_REGION: {coverage_outside_insert}")
        logging.debug(f"Mismatches near INSERT_REGION: {mismatches_near_insert}")
    else:
        logging.warning("Skipping coverage/mismatch metrics (insert region unavailable).")

    # Compute whole-plasmid depths once if possible for sharing across metrics and genes
    plasmid_depths: Any = None
    plasmid_length = 0
    try:
        with pysam.AlignmentFile(plasmid_bam, "rb") as bam:
            if (
                hasattr(bam, "lengths")
                and isinstance(bam.lengths, (list, tuple))
                and len(bam.lengths) > 0
                and bam.lengths[0] > 0
            ):
                contig = bam.get_reference_name(0)
                plasmid_length = int(bam.lengths[0])
                depth_arr = get_depth_array(plasmid_bam, contig, 0, plasmid_length)
                if isinstance(depth_arr, np.ndarray) and len(depth_arr) > 0:
                    plasmid_depths = depth_arr
    except Exception as e:
        logging.warning(f"Could not compute whole-plasmid depth array: {e}")

    # Compute per-region coverage metrics
    coverage_metrics, coverage_fallback = compute_region_coverage_metrics(
        plasmid_bam, insert_region, BREADTH_THRESHOLDS, plasmid_depths=plasmid_depths
    )
    if coverage_fallback:
        logging.warning(
            "Coverage metrics computed for whole plasmid (insert region not defined). "
            "Backbone metrics are 0.0."
        )

    # Resistance gene detection and per-gene coverage
    detected_genes: list[DetectedGene] = []
    gene_coverage: dict[str, dict[str, float]] = {}
    genbank_located = False

    try:
        # Try provided genbank_file first, fallback to locate_genbank_file
        gb_file = genbank_file
        if gb_file is None:
            gb_file = locate_genbank_file(os.path.dirname(output_basename))
            if gb_file is not None:
                genbank_located = True
                logging.debug(f"Located GenBank file: {gb_file}")

        if gb_file is not None:
            logging.info(f"Detecting resistance genes from {gb_file}")
            detected_genes = detect_resistance_genes(
                gb_file, insert_region, CUSTOM_RESISTANCE_PATTERNS
            )
            logging.info(f"Detected {len(detected_genes)} resistance/reporter genes")

            if detected_genes:
                for gene in detected_genes:
                    if plasmid_depths is not None and len(plasmid_depths) > 0:
                        if gene.start > gene.end:
                            depth_array = np.concatenate(
                                [plasmid_depths[gene.start :], plasmid_depths[: gene.end + 1]]
                            )
                        else:
                            depth_array = plasmid_depths[gene.start : gene.end + 1]
                    else:
                        with pysam.AlignmentFile(plasmid_bam, "rb") as bam:
                            c = bam.get_reference_name(0)
                            plen = bam.lengths[0]
                        if gene.start > gene.end:
                            depth_array = np.concatenate(
                                [
                                    get_depth_array(plasmid_bam, c, gene.start, plen),
                                    get_depth_array(plasmid_bam, c, 0, gene.end + 1),
                                ]
                            )
                        else:
                            depth_array = get_depth_array(plasmid_bam, c, gene.start, gene.end + 1)
                    stats = compute_coverage_stats(depth_array, [1])
                    gene_coverage[gene.sanitized_name] = stats
                    logging.debug(
                        f"Gene {gene.name}: mean_depth={stats['mean_depth']:.2f}, "
                        f"breadth={stats['breadth_1x']:.2f}"
                    )
    except Exception as e:
        logging.warning(f"Resistance gene detection failed: {e}. Continuing without gene data.")
        detected_genes = []
        gene_coverage = {}

    # Group BAMs by read name into temporary files for streaming merge
    output_dir = os.path.dirname(output_basename) or "."
    output_name = os.path.basename(output_basename)
    pid = os.getpid()
    plasmid_ns = os.path.join(output_dir, f"{output_name}.plasmid.namesorted.{pid}.bam")
    human_ns = os.path.join(output_dir, f"{output_name}.human.namesorted.{pid}.bam")
    human_already_namesorted = False
    try:
        with pysam.AlignmentFile(human_bam, "rb") as hbam:
            hd = hbam.header.get("HD", {}) if hasattr(hbam.header, "get") else {}
            if isinstance(hd, dict) and hd.get("SO") == "queryname":
                human_already_namesorted = True
    except Exception:
        pass

    created_temp_bams: list[str] = []
    try:
        logging.info("Grouping BAMs by read name for streaming comparison...")
        created_temp_bams.append(plasmid_ns)
        _namesort_bam(plasmid_bam, plasmid_ns, threads=threads)
        if human_already_namesorted:
            human_ns_to_use = human_bam
        else:
            created_temp_bams.append(human_ns)
            _namesort_bam(human_bam, human_ns, threads=threads)
            human_ns_to_use = human_ns

        logging.debug(
            f"Writing alignment comparison results to {output_basename}.reads_assignment.tsv"
        )

        if anonymize_key is None and anonymize:
            anonymize_key = secrets.token_bytes(16)

        with open(f"{output_basename}.reads_assignment.tsv", "w") as outfile:
            outfile.write(
                "ReadID\tAssignedTo\tPlasmidScore\tHumanScore\tPlasmidCIGAR\tHumanCIGAR\tPlasmidMapQ\tHumanMapQ\n"
            )
            assigned_counts = _streaming_compare(
                plasmid_ns,
                human_ns_to_use,
                outfile,
                insert_region=insert_region,
                score_margin=SCORE_MARGIN,
                plasmid_length=plasmid_length,
                anonymize_key=anonymize_key,
            )
    finally:
        for ns_bam in created_temp_bams:
            with contextlib.suppress(OSError):
                os.remove(ns_bam)

    vector_insert_bridge_pairs = getattr(assigned_counts, "bridge_pairs", 0)

    # Calculate EstimatedContaminationFraction
    junc_positions = (
        _compute_junction_positions(insert_region, plasmid_length)
        if insert_region is not None
        else []
    )
    j_count = len(junc_positions)
    insert_metrics = coverage_metrics["insert"]
    d_insert = insert_metrics.get("mean_depth", 0.0)
    if (
        j_count > 0
        and not coverage_fallback
        and insert_region is not None
        and plasmid_depths is not None
        and plasmid_length > 0
    ):
        try:
            with pysam.AlignmentFile(plasmid_bam, "rb") as bam:
                c_name = bam.get_reference_name(0)

            def primary_cb(r: pysam.AlignedSegment) -> bool:
                return not (
                    r.is_unmapped
                    or r.is_secondary
                    or r.is_supplementary
                    or getattr(r, "is_duplicate", False)
                    or getattr(r, "is_qcfail", False)
                )

            ins_s, ins_e = insert_region[0], insert_region[1]
            if ins_s > ins_e:
                prim_depth = np.concatenate(
                    [
                        get_depth_array(
                            plasmid_bam, c_name, ins_s, plasmid_length, read_callback=primary_cb
                        ),
                        get_depth_array(
                            plasmid_bam, c_name, 0, ins_e + 1, read_callback=primary_cb
                        ),
                    ]
                )
            else:
                prim_depth = get_depth_array(
                    plasmid_bam, c_name, ins_s, ins_e + 1, read_callback=primary_cb
                )
            if len(prim_depth) > 0:
                d_insert = float(np.mean(prim_depth))
        except Exception:
            pass
    a = 15
    w = mean_read_length - 2 * a + 1
    if j_count == 0 or coverage_fallback or d_insert <= 0 or w <= 0:
        estimated_contamination_fraction = "NA"
    elif splice_junction_events == 0:
        estimated_contamination_fraction = f"{0.0:.4f}"
    else:
        est = (splice_junction_events * mean_read_length) / (j_count * w * d_insert)
        estimated_contamination_fraction = f"{est:.4f}"

    plasmid_num = (
        assigned_counts["Plasmid"]
        if filter_backbone_only
        else assigned_counts["Plasmid"] + assigned_counts["Backbone_Only"]
    )
    human_count = assigned_counts["Human"]
    ratio = (
        plasmid_num / human_count
        if human_count != 0
        else (float("inf") if plasmid_num > 0 else 0.0)
    )

    category = categorize_value(ratio, threshold, UNCLEAR_RANGE)
    if category == "Contaminated":
        verdict = "Sample is contaminated with plasmid DNA"
    elif category == "Unclear":
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
        summary_file.write(
            f"VectorInsertJunctionReads\t{vector_insert_junction_reads if insert_region is not None else 'NA'}\n"
        )
        summary_file.write(
            f"VectorInsertBridgePairs\t{vector_insert_bridge_pairs if insert_region is not None else 'NA'}\n"
        )

        # Coverage metrics (new in Phase 9)
        insert_metrics = coverage_metrics["insert"]
        backbone_metrics = coverage_metrics["backbone"]

        # Insert metrics
        summary_file.write(f"MeanDepthInsert\t{insert_metrics['mean_depth']:.2f}\n")
        summary_file.write(f"MedianDepthInsert\t{insert_metrics['median_depth']:.2f}\n")
        summary_file.write(f"BreadthInsert\t{insert_metrics['breadth_1x']:.2f}\n")
        for depth_threshold in BREADTH_THRESHOLDS:
            if depth_threshold != 1:
                summary_file.write(
                    f"BreadthInsert_{depth_threshold}x\t{insert_metrics[f'breadth_{depth_threshold}x']:.2f}\n"
                )
        summary_file.write(f"CoverageCV_Insert\t{insert_metrics['cv']:.2f}\n")

        # Backbone metrics
        summary_file.write(f"MeanDepthBackbone\t{backbone_metrics['mean_depth']:.2f}\n")
        summary_file.write(f"MedianDepthBackbone\t{backbone_metrics['median_depth']:.2f}\n")
        summary_file.write(f"BreadthBackbone\t{backbone_metrics['breadth_1x']:.2f}\n")
        for depth_threshold in BREADTH_THRESHOLDS:
            if depth_threshold != 1:
                summary_file.write(
                    f"BreadthBackbone_{depth_threshold}x\t{backbone_metrics[f'breadth_{depth_threshold}x']:.2f}\n"
                )
        summary_file.write(f"CoverageCV_Backbone\t{backbone_metrics['cv']:.2f}\n")

        # Fallback indicator
        summary_file.write(f"CoverageFallback\t{coverage_fallback}\n")

        # Resistance gene detection results
        if genbank_file is not None or genbank_located:
            summary_file.write("ResistanceGenesAvailable\ttrue\n")
            if gene_coverage:
                gene_names = [g.name for g in detected_genes]
                summary_file.write(f"ResistanceGenes\t{','.join(gene_names)}\n")
                for gene in detected_genes:
                    s = gene.sanitized_name
                    metrics = gene_coverage[s]
                    summary_file.write(f"ResGene_{s}_MeanDepth\t{metrics['mean_depth']:.2f}\n")
                    summary_file.write(f"ResGene_{s}_Breadth\t{metrics['breadth_1x']:.2f}\n")
                    summary_file.write(f"ResGene_{s}_Start\t{gene.start}\n")
                    summary_file.write(f"ResGene_{s}_End\t{gene.end}\n")
                    summary_file.write(f"ResGene_{s}_Region\t{gene.region}\n")
            else:
                summary_file.write("ResistanceGenes\tNone\n")
        else:
            summary_file.write("ResistanceGenesAvailable\tfalse\n")

        # Junction metadata (Phase 15)
        if insert_region is not None:
            jc = insert_region.junction_count
            ec = insert_region.exon_count
            el = insert_region.exon_lengths
            il = insert_region.intron_lengths

            summary_file.write(f"SpliceJunctionCount\t{jc if jc is not None else 'NA'}\n")
            summary_file.write(f"ExonCount\t{ec if ec is not None else 'NA'}\n")

            if el is not None and len(el) > 0:
                summary_file.write(f"ExonLengths\t{','.join(map(str, el))}\n")
            else:
                summary_file.write("ExonLengths\tNA\n")

            if il is not None and len(il) > 0:
                summary_file.write(f"IntronLengths\t{','.join(map(str, il))}\n")
            else:
                summary_file.write("IntronLengths\tNA\n")
            summary_file.write(f"SpliceJunctionReads\t{splice_junction_reads}\n")
            summary_file.write(f"SpliceJunctionEvents\t{splice_junction_events}\n")
            summary_file.write(f"SpliceJunctionsCovered\t{splice_junctions_covered}\n")
            summary_file.write(f"ExonIntronBoundaryReads\t{exon_intron_boundary_reads}\n")
            summary_file.write(
                f"EstimatedContaminationFraction\t{estimated_contamination_fraction}\n"
            )
        else:
            # Fallback for missing insert_region
            summary_file.write("SpliceJunctionCount\tNA\n")
            summary_file.write("ExonCount\tNA\n")
            summary_file.write("ExonLengths\tNA\n")
            summary_file.write("IntronLengths\tNA\n")
            summary_file.write("SpliceJunctionReads\tNA\n")
            summary_file.write("SpliceJunctionEvents\tNA\n")
            summary_file.write("SpliceJunctionsCovered\tNA\n")
            summary_file.write("ExonIntronBoundaryReads\tNA\n")
            summary_file.write("EstimatedContaminationFraction\tNA\n")
