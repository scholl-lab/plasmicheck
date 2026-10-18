#!/usr/bin/env python3
"""In-silico contamination simulator for PlasmiCheck.

Simulates hybrid-capture enrichment pull-down and sequencing of circular
plasmids and host genomic sequences, supporting both real-BAM spike-in
and de novo locus simulation.
"""

from __future__ import annotations

import argparse
import gzip
import random
import shutil
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, TextIO, cast

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

Interval = tuple[int, int]
Pair = tuple[int, int, bool, str, str]
# Pair: fragment start, length, reverse orientation, R1 sequence, R2 sequence


@contextmanager
def open_text(path: Path, mode: str = "r") -> Iterator[TextIO]:
    """Open a plain text or gzip text file via context manager."""
    if str(path).endswith(".gz"):
        gz_mode = mode if ("t" in mode or "b" in mode) else f"{mode}t"
        with gzip.open(path, gz_mode) as f_gz:
            yield cast(TextIO, f_gz)
    else:
        with open(path, mode) as f_txt:
            yield cast(TextIO, f_txt)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def mutate_seq(seq: str, error_rate: float, rng: random.Random) -> str:
    """Introduce random base substitution errors at the given error rate."""
    if error_rate <= 0.0:
        return seq
    bases = list(seq.upper())
    alts = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG", "N": "ACGT"}
    for i, b in enumerate(bases):
        if rng.random() < error_rate:
            choices = alts.get(b, "ACGT")
            bases[i] = rng.choice(choices)
    return "".join(bases)


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA file into a mapping of record name to sequence."""
    records: dict[str, str] = {}
    current_name: str | None = None
    chunks: list[str] = []
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(chunks)
                    chunks = []
                current_name = stripped[1:].split()[0]
            else:
                chunks.append(stripped)
        if current_name is not None:
            records[current_name] = "".join(chunks)
    return records


def parse_interval_string(s: str) -> list[Interval]:
    """Parse comma-separated intervals, e.g. '200-400,500-700'."""
    intervals: list[Interval] = []
    for part in s.split(","):
        p = part.strip()
        if not p:
            continue
        sep = "-" if "-" in p else ":"
        start_s, end_s = p.split(sep, 1)
        intervals.append((int(start_s), int(end_s)))
    return intervals


def load_cdna_targets(path: Path, plasmid_len: int) -> list[Interval]:
    """Load insert targets from cDNA_positions.txt file using parse_insert_region."""
    from plasmicheck.scripts.compare_alignments import parse_insert_region

    data = parse_insert_region(str(path))
    start = data.start
    end = data.end + 1  # 0-based inclusive -> half-open
    if start > end:
        return [(start, plasmid_len), (0, end)]
    return [(start, end)]


def overlap(
    start: int,
    length: int,
    targets: list[Interval],
    source_length: int,
    *,
    circular: bool,
) -> int:
    """Calculate the base-pair overlap between a fragment and target intervals."""
    if not targets or length <= 0 or source_length <= 0:
        return 0

    norm_targets: list[Interval] = []
    for s, e in targets:
        s_mod = s % source_length
        e_mod = e % source_length if e != source_length else source_length
        if s > e:
            norm_targets.append((s_mod, source_length))
            if e_mod > 0:
                norm_targets.append((0, e_mod))
        elif s < e:
            norm_targets.append((s_mod, min(source_length, e)))

    frag_intervals: list[Interval] = []
    if circular:
        f_start = start % source_length
        if length >= source_length:
            frag_intervals.append((0, source_length))
        elif f_start + length <= source_length:
            frag_intervals.append((f_start, f_start + length))
        else:
            frag_intervals.append((f_start, source_length))
            frag_intervals.append((0, (f_start + length) % source_length))
    else:
        frag_intervals.append((start, min(source_length, start + length)))

    total_bp = 0
    for fs, fe in frag_intervals:
        for ts, te in norm_targets:
            ov = min(fe, te) - max(fs, ts)
            if ov > 0:
                total_bp += ov
    return total_bp


def classify_pair(
    start: int,
    length: int,
    is_reverse: bool,
    read_length: int,
    targets: list[Interval],
    source_length: int,
    *,
    circular: bool,
) -> list[str]:
    """Determine biophysical and boundary flags for a simulated pair."""
    flags: list[str] = []
    if circular and (start + length > source_length):
        flags.append("origin_wrap")

    ov = overlap(start, length, targets, source_length, circular=circular)
    if not targets or ov == 0:
        flags.append("backbone_only")
    elif ov == length:
        flags.append("insert_internal")
    else:
        flags.append("vector_boundary")

    if not targets:
        return flags

    coords = [(start + i) % source_length if circular else (start + i) for i in range(length)]
    m1_coords = set(coords[-read_length:] if is_reverse else coords[:read_length])
    m2_coords = set(coords[:read_length] if is_reverse else coords[-read_length:])
    bridge_coords = set(coords[read_length:-read_length]) if length > 2 * read_length else set()

    norm: list[Interval] = []
    for s, e in targets:
        if s > e and circular:
            norm.append((s, source_length))
            norm.append((0, e))
        else:
            norm.append((s, e))

    vector_bounds = {norm[0][0], norm[-1][1]}
    splice_bounds = {norm[i][1] for i in range(len(norm) - 1)}

    def spans(c_set: set[int], bound: int) -> bool:
        prev = (bound - 1) % source_length if circular else bound - 1
        return prev in c_set and bound in c_set

    for sb in sorted(splice_bounds):
        if spans(m1_coords, sb):
            flags.append("r1_splice_junction")
        if spans(m2_coords, sb):
            flags.append("r2_splice_junction")
        if spans(bridge_coords, sb) and not spans(m1_coords, sb) and not spans(m2_coords, sb):
            flags.append("fragment_bridge")
            flags.append("splice_junction_bridge")

    for vb in sorted(vector_bounds):
        if spans(m1_coords, vb):
            flags.append("r1_vector_boundary")
        if spans(m2_coords, vb):
            flags.append("r2_vector_boundary")
        if spans(bridge_coords, vb) and not spans(m1_coords, vb) and not spans(m2_coords, vb):
            flags.append("fragment_bridge")
            flags.append("vector_boundary_bridge")

    return list(dict.fromkeys(flags))


def simulate_pairs(
    sequence: str,
    targets: list[Interval],
    count: int,
    *,
    circular: bool,
    read_length: int,
    fragment_mean: float,
    fragment_sd: float,
    min_overlap: int,
    epsilon: float,
    error_rate: float,
    rng: random.Random,
) -> Iterator[Pair]:
    """Simulate paired-end reads from a template sequence using biophysical retention."""
    if count < 0:
        raise ValueError(f"count must be non-negative, got {count}")
    if count == 0:
        return

    source_length = len(sequence)
    if source_length == 0:
        raise ValueError("sequence must not be empty")
    if read_length <= 0:
        raise ValueError(f"read_length must be positive, got {read_length}")
    if read_length > source_length:
        raise ValueError(f"read_length ({read_length}) exceeds sequence length ({source_length})")
    if fragment_mean < read_length:
        raise ValueError(f"fragment_mean ({fragment_mean}) must be >= read_length ({read_length})")
    if fragment_sd < 0:
        raise ValueError(f"fragment_sd must be non-negative, got {fragment_sd}")
    if not (0.0 <= epsilon <= 1.0):
        raise ValueError(f"epsilon must be in [0.0, 1.0], got {epsilon}")
    if not (0.0 <= error_rate <= 1.0):
        raise ValueError(f"error_rate must be in [0.0, 1.0], got {error_rate}")
    if min_overlap <= 0:
        raise ValueError(f"min_overlap must be positive, got {min_overlap}")

    if epsilon <= 0.0:
        if not targets:
            raise ValueError("Impossible capture setting: epsilon is 0.0 and targets is empty")
        max_ov = overlap(0, source_length, targets, source_length, circular=circular)
        if max_ov < min_overlap:
            raise ValueError(
                f"Impossible capture setting: max overlap ({max_ov}) < min_overlap ({min_overlap})"
            )

    accepted = 0
    attempts = 0
    max_attempts = max(100_000, count * 10_000)

    while accepted < count:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError(
                f"Exceeded maximum attempts ({max_attempts}) without yielding {count} pairs."
            )

        flen = round(rng.gauss(fragment_mean, fragment_sd))
        if flen < read_length:
            continue
        if flen > source_length:
            flen = source_length

        if circular:
            start = rng.randint(0, source_length - 1)
        else:
            max_start = source_length - flen
            if max_start < 0:
                continue
            start = rng.randint(0, max_start)

        ov = overlap(start, flen, targets, source_length, circular=circular)
        p_keep = 1.0 if ov >= min_overlap else epsilon
        if p_keep <= 0.0 or (p_keep < 1.0 and rng.random() >= p_keep):
            continue

        if circular:
            if start + flen <= source_length:
                frag_seq = sequence[start : start + flen]
            else:
                frag_seq = sequence[start:] + sequence[: flen - (source_length - start)]
        else:
            frag_seq = sequence[start : start + flen]

        is_rev = rng.random() < 0.5
        oriented = reverse_complement(frag_seq) if is_rev else frag_seq
        r1 = oriented[:read_length]
        r2 = reverse_complement(oriented[-read_length:])

        if error_rate > 0.0:
            r1 = mutate_seq(r1, error_rate, rng)
            r2 = mutate_seq(r2, error_rate, rng)

        yield (start, flen, is_rev, r1, r2)
        accepted += 1


def extract_background(
    bam: Path,
    r1: Path,
    r2: Path,
    *,
    reserved_prefix: str = "pcsim:",
) -> int:
    """Extract host primary paired reads from BAM using samtools collate | fastq."""
    bam = Path(bam)
    r1 = Path(r1)
    r2 = Path(r2)

    if not bam.exists():
        raise FileNotFoundError(f"Input BAM file not found: {bam}")
    if bam.resolve() in (r1.resolve(), r2.resolve()):
        raise ValueError("Output FASTQ path cannot match input BAM path")

    r1.parent.mkdir(parents=True, exist_ok=True)
    r2.parent.mkdir(parents=True, exist_ok=True)

    def cleanup() -> None:
        for p in (r1, r2):
            if p.exists():
                p.unlink()

    if shutil.which("samtools") is not None:
        collate_cmd = ["samtools", "collate", "-O", "-u", str(bam)]
        fastq_cmd = [
            "samtools",
            "fastq",
            "-n",
            "-F",
            "0x900",
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-0",
            "/dev/null",
            "-s",
            "/dev/null",
            "-",
        ]

        p1 = subprocess.Popen(collate_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(fastq_cmd, stdin=p1.stdout, stderr=subprocess.PIPE)
        if p1.stdout is not None:
            p1.stdout.close()
        _, err2 = p2.communicate()
        _, err1 = p1.communicate()

        if p1.returncode != 0:
            cleanup()
            raise RuntimeError(f"samtools collate failed: {err1.decode(errors='replace')}")
        if p2.returncode != 0:
            cleanup()
            raise RuntimeError(f"samtools fastq failed: {err2.decode(errors='replace')}")
    else:
        import pysam

        mates: dict[str, dict[int, tuple[str, str]]] = {}
        with pysam.AlignmentFile(str(bam), "rb") as bam_in:
            for rec in bam_in.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or not rec.is_paired:
                    continue
                if rec.is_read1 == rec.is_read2 or not rec.query_name or not rec.query_sequence:
                    continue
                seq = (
                    reverse_complement(rec.query_sequence) if rec.is_reverse else rec.query_sequence
                )
                raw_qual: str | None = getattr(rec, "qual", None)
                quals = raw_qual or ("I" * len(seq))
                if rec.is_reverse and raw_qual:
                    quals = raw_qual[::-1]
                mate_num = 1 if rec.is_read1 else 2
                mates.setdefault(rec.query_name, {})[mate_num] = (seq, quals)
        with open_text(r1, "w") as out1, open_text(r2, "w") as out2:
            for qname, pair_dict in mates.items():
                if 1 in pair_dict and 2 in pair_dict:
                    s1, q1 = pair_dict[1]
                    s2, q2 = pair_dict[2]
                    out1.write(f"@{qname}\n{s1}\n+\n{q1}\n")
                    out2.write(f"@{qname}\n{s2}\n+\n{q2}\n")

    retained_pairs = 0
    try:
        with open_text(r1) as f1, open_text(r2) as f2:
            while True:
                h1 = f1.readline()
                h2 = f2.readline()
                if not h1 and not h2:
                    break
                if not h1 or not h2:
                    raise ValueError("Mismatched FASTQ length between R1 and R2")

                for _ in range(3):
                    f1.readline()
                    f2.readline()

                name1 = h1.strip().split()[0].lstrip("@")
                name2 = h2.strip().split()[0].lstrip("@")
                base1 = name1[:-2] if name1.endswith(("/1", "/2")) else name1
                base2 = name2[:-2] if name2.endswith(("/1", "/2")) else name2

                if base1 != base2:
                    raise ValueError(f"Mismatched read pair names: {name1} vs {name2}")
                if name1.startswith(reserved_prefix) or name2.startswith(reserved_prefix):
                    raise ValueError(
                        f"QNAME collision: background read starts with '{reserved_prefix}': {name1}"
                    )
                retained_pairs += 1
    except Exception:
        cleanup()
        raise

    return retained_pairs


def write_pairs(
    pairs: Iterable[Pair],
    r1: TextIO,
    r2: TextIO,
    truth: TextIO,
    *,
    prefix: str = "pcsim:",
    source: str = "plasmid",
    targets: list[Interval] | None = None,
    source_length: int | None = None,
    circular: bool = False,
    read_length: int | None = None,
) -> int:
    """Write simulated pairs to FASTQ streams and record ground truth."""
    try:
        if truth.tell() == 0:
            truth.write("read_id\tsource\tstart\tlength\tis_reverse\toverlap\tflags\n")
    except (OSError, AttributeError):
        pass

    count = 0
    clean_prefix = prefix if prefix.endswith((":", "_")) else f"{prefix}_"

    for idx, (start, length, is_rev, r1_seq, r2_seq) in enumerate(pairs, 1):
        qname = f"{clean_prefix}{source}_{idx:06d}"
        rlen = read_length if read_length is not None else len(r1_seq)

        flags_str = "none"
        ov = 0
        if targets is not None and source_length is not None:
            ov = overlap(start, length, targets, source_length, circular=circular)
            flag_list = classify_pair(
                start, length, is_rev, rlen, targets, source_length, circular=circular
            )
            flags_str = ";".join(flag_list) if flag_list else "none"

        r1.write(f"@{qname}/1\n{r1_seq}\n+\n{'I' * len(r1_seq)}\n")
        r2.write(f"@{qname}/2\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")
        truth.write(f"{qname}\t{source}\t{start}\t{length}\t{int(is_rev)}\t{ov}\t{flags_str}\n")
        count += 1

    return count


def resolve_targets(args: argparse.Namespace, plasmid_len: int) -> list[Interval]:
    """Resolve target intervals from CLI arguments."""
    if args.insert_exons:
        return parse_interval_string(args.insert_exons)
    if args.cdna_positions:
        return load_cdna_targets(args.cdna_positions, plasmid_len)
    if args.insert_start is not None and args.insert_end is not None:
        s, e = args.insert_start, args.insert_end
        if s > e:
            return [(s, plasmid_len), (0, e)]
        return [(s, e)]
    return [(0, plasmid_len)]


def build_parser() -> argparse.ArgumentParser:
    """Build CLI parser for simulate_contamination."""
    parser = argparse.ArgumentParser(
        description="In-silico contamination simulator for PlasmiCheck."
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    def add_common_args(p: argparse.ArgumentParser) -> None:
        p.add_argument("--plasmid-fasta", type=Path, required=True, help="Plasmid FASTA")
        p.add_argument("--cdna-positions", type=Path, default=None, help="cDNA_positions.txt file")
        p.add_argument("--insert-start", type=int, default=None, help="Insert start (0-based)")
        p.add_argument("--insert-end", type=int, default=None, help="Insert end (0-based)")
        p.add_argument("--insert-exons", type=str, default=None, help="Exons: '200-400,500-700'")
        p.add_argument("--out-prefix", type=str, default="simulated", help="Output prefix")
        p.add_argument("--out-r1", type=Path, default=None, help="Output R1 FASTQ(.gz)")
        p.add_argument("--out-r2", type=Path, default=None, help="Output R2 FASTQ(.gz)")
        p.add_argument("--out-truth", type=Path, default=None, help="Output truth TSV")
        p.add_argument("--read-length", type=int, default=150, help="Read length in bp")
        p.add_argument("--fragment-mean", type=float, default=300.0, help="Fragment mean length")
        p.add_argument("--fragment-sd", type=float, default=50.0, help="Fragment length SD")
        p.add_argument("--min-overlap", type=int, default=45, help="Min overlap bp")
        p.add_argument("--epsilon", type=float, default=0.01, help="Backbone retention probability")
        p.add_argument("--error-rate", type=float, default=0.005, help="Sequencing error rate")
        p.add_argument("--seed", type=int, default=42, help="Random seed")
        p.add_argument("--prefix", type=str, default="pcsim:", help="Read name prefix")

    sp = subparsers.add_parser("spike-in", help="Extract background from BAM and spike in plasmid")
    sp.add_argument("--bam", type=Path, required=True, help="Input background BAM")
    sp.add_argument("--plasmid-count", type=int, default=50, help="Number of plasmid pairs")
    sp.add_argument("--ratio", type=float, default=None, help="Spike-in ratio to background")
    add_common_args(sp)

    dp = subparsers.add_parser("de-novo", help="Simulate host locus and plasmid de novo")
    dp.add_argument("--host-fasta", type=Path, required=True, help="Host sequence FASTA")
    dp.add_argument("--host-count", type=int, default=200, help="Number of host pairs")
    dp.add_argument("--host-exons", type=str, default=None, help="Host exons: '100-300,500-700'")
    dp.add_argument("--plasmid-count", type=int, default=50, help="Number of plasmid pairs")
    add_common_args(dp)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Execute simulator CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)

    out_prefix = args.out_prefix
    r1_path: Path = args.out_r1 or Path(f"{out_prefix}_R1.fastq.gz")
    r2_path: Path = args.out_r2 or Path(f"{out_prefix}_R2.fastq.gz")
    truth_path: Path = args.out_truth or Path(f"{out_prefix}_truth.tsv")

    if args.mode == "spike-in":
        if args.bam.resolve() in (r1_path.resolve(), r2_path.resolve()):
            raise ValueError("Output FASTQ path cannot match input BAM path")

        retained_bg = extract_background(args.bam, r1_path, r2_path, reserved_prefix=args.prefix)
        n_plasmid = (
            round(retained_bg * args.ratio) if args.ratio is not None else args.plasmid_count
        )

        if n_plasmid > 0:
            plasmid_recs = read_fasta(args.plasmid_fasta)
            plasmid_seq = next(iter(plasmid_recs.values()))
            targets = resolve_targets(args, len(plasmid_seq))
            rng = random.Random(args.seed)
            pairs = simulate_pairs(
                plasmid_seq,
                targets,
                n_plasmid,
                circular=True,
                read_length=args.read_length,
                fragment_mean=args.fragment_mean,
                fragment_sd=args.fragment_sd,
                min_overlap=args.min_overlap,
                epsilon=args.epsilon,
                error_rate=args.error_rate,
                rng=rng,
            )
            truth_path.parent.mkdir(parents=True, exist_ok=True)
            with (
                open_text(r1_path, "a") as f1,
                open_text(r2_path, "a") as f2,
                open(truth_path, "w") as ft,
            ):
                write_pairs(
                    pairs,
                    f1,
                    f2,
                    ft,
                    prefix=args.prefix,
                    source="plasmid",
                    targets=targets,
                    source_length=len(plasmid_seq),
                    circular=True,
                    read_length=args.read_length,
                )

    elif args.mode == "de-novo":
        host_recs = read_fasta(args.host_fasta)
        host_seq = next(iter(host_recs.values()))
        plasmid_recs = read_fasta(args.plasmid_fasta)
        plasmid_seq = next(iter(plasmid_recs.values()))

        host_targets = (
            parse_interval_string(args.host_exons) if args.host_exons else [(0, len(host_seq))]
        )
        plasmid_targets = resolve_targets(args, len(plasmid_seq))

        rng = random.Random(args.seed)
        host_pairs = simulate_pairs(
            host_seq,
            host_targets,
            args.host_count,
            circular=False,
            read_length=args.read_length,
            fragment_mean=args.fragment_mean,
            fragment_sd=args.fragment_sd,
            min_overlap=args.min_overlap,
            epsilon=args.epsilon,
            error_rate=args.error_rate,
            rng=rng,
        )
        plasmid_pairs = simulate_pairs(
            plasmid_seq,
            plasmid_targets,
            args.plasmid_count,
            circular=True,
            read_length=args.read_length,
            fragment_mean=args.fragment_mean,
            fragment_sd=args.fragment_sd,
            min_overlap=args.min_overlap,
            epsilon=args.epsilon,
            error_rate=args.error_rate,
            rng=rng,
        )

        r1_path.parent.mkdir(parents=True, exist_ok=True)
        r2_path.parent.mkdir(parents=True, exist_ok=True)
        truth_path.parent.mkdir(parents=True, exist_ok=True)

        with (
            open_text(r1_path, "w") as f1,
            open_text(r2_path, "w") as f2,
            open(truth_path, "w") as ft,
        ):
            write_pairs(
                host_pairs,
                f1,
                f2,
                ft,
                prefix=f"{args.prefix}host_",
                source="host",
                targets=host_targets,
                source_length=len(host_seq),
                circular=False,
                read_length=args.read_length,
            )
            write_pairs(
                plasmid_pairs,
                f1,
                f2,
                ft,
                prefix=f"{args.prefix}plasmid_",
                source="plasmid",
                targets=plasmid_targets,
                source_length=len(plasmid_seq),
                circular=True,
                read_length=args.read_length,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
