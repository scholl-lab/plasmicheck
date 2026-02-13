#!/usr/bin/env python3
"""Generate deterministic synthetic test data for PlasmiCheck integration tests.

This script creates a minimal, self-contained test dataset (~50 KB total) that
exercises the full pipeline in <10 seconds on CI.  All sequences are synthetic
(no real genomic data) and deterministically seeded for reproducibility.

Usage
-----
    python generate_test_data.py          # writes files to the same directory
    python generate_test_data.py --outdir /tmp/synth

Design
------
- **human_ref.fasta** — 2 fake chromosomes (~2 000 bp each).
- **plasmid.gb** — GenBank file with a backbone + cDNA insert copied from
  ``chr_fake1[200:700]`` so that spliced alignment can discover it.
- **contaminated_R1/R2.fastq** — 200 read pairs, 60 % from plasmid, 40 % from
  human reference, with 1 % per-base mutation rate for realism.
- **not_contaminated_R1/R2.fastq** — 200 read pairs, 100 % from human ref.
- Read length: 150 bp, insert size ~300 bp (paired-end).
"""

from __future__ import annotations

import argparse
import random
from datetime import datetime, timezone
from pathlib import Path

SEED = 42
READ_LEN = 150
INSERT_SIZE = 300
N_READS = 200
CONTAM_FRAC = 0.60  # fraction of plasmid-origin reads in the contaminated sample
MUTATION_RATE = 0.01

# Region of chr_fake1 that is the "cDNA insert" in the plasmid
INSERT_START = 200
INSERT_END = 700


def _random_seq(length: int, rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", k=length))


def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    bases = list(seq)
    for i in range(len(bases)):
        if rng.random() < rate:
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def _quality_string(length: int) -> str:
    """Return a constant high-quality string (Phred ~30)."""
    return "I" * length


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def _write_genbank(path: Path, name: str, seq: str, insert_start: int, insert_end: int) -> None:
    """Write a minimal GenBank file with a CDS feature marking the insert."""
    today = datetime.now(tz=timezone.utc).strftime("%d-%b-%Y").upper()
    locus_len = len(seq)

    # GenBank fixed-width format
    lines: list[str] = [
        f"LOCUS       {name:<16} {locus_len} bp    DNA     circular SYN {today}",
        f"DEFINITION  Synthetic test plasmid {name}.",
        f"ACCESSION   {name}",
        f"VERSION     {name}.1",
        "KEYWORDS    .",
        "SOURCE      synthetic construct",
        "  ORGANISM  synthetic construct",
        "            other sequences; artificial sequences.",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{locus_len}",
        '                     /organism="synthetic construct"',
        '                     /mol_type="other DNA"',
        f"     CDS             {insert_start + 1}..{insert_end}",
        '                     /gene="test_insert"',
        '                     /note="cDNA insert from chr_fake1"',
        "ORIGIN",
    ]

    # Sequence block in GenBank format (60 bases per line, numbered)
    for i in range(0, len(seq), 60):
        chunk = seq[i : i + 60]
        # Split into groups of 10
        groups = [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
        lines.append(f"{i + 1:>9} {' '.join(groups)}")
    lines.append("//")

    path.write_text("\n".join(lines) + "\n")


def _generate_read_pair(source_seq: str, rng: random.Random) -> tuple[str, str]:
    """Generate a paired-end read pair from *source_seq*."""
    max_start = len(source_seq) - INSERT_SIZE
    if max_start < 0:
        max_start = 0
    start = rng.randint(0, max_start)
    fragment = source_seq[start : start + INSERT_SIZE]
    if len(fragment) < INSERT_SIZE:
        fragment = fragment + _random_seq(INSERT_SIZE - len(fragment), rng)

    r1 = fragment[:READ_LEN]
    r2 = _reverse_complement(fragment[-READ_LEN:])

    # Apply mutations for realism
    r1 = _mutate(r1, MUTATION_RATE, rng)
    r2 = _mutate(r2, MUTATION_RATE, rng)
    return r1, r2


def _write_fastq(
    path_r1: Path,
    path_r2: Path,
    read_pairs: list[tuple[str, str, str]],
) -> None:
    """Write paired FASTQ files.  *read_pairs* is [(name, r1_seq, r2_seq), ...]."""
    with open(path_r1, "w") as fh1, open(path_r2, "w") as fh2:
        for name, r1, r2 in read_pairs:
            qual = _quality_string(READ_LEN)
            fh1.write(f"@{name}/1\n{r1}\n+\n{qual}\n")
            fh2.write(f"@{name}/2\n{r2}\n+\n{qual}\n")


def main(outdir: Path | None = None) -> None:
    rng = random.Random(SEED)

    if outdir is None:
        outdir = Path(__file__).resolve().parent
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Generate human reference ---
    chr1_seq = _random_seq(2000, rng)
    chr2_seq = _random_seq(2000, rng)
    _write_fasta(
        outdir / "human_ref.fasta",
        [("chr_fake1", chr1_seq), ("chr_fake2", chr2_seq)],
    )

    # --- Generate plasmid (GenBank) ---
    # Backbone = random, insert = copy of chr_fake1[200:700]
    backbone_5 = _random_seq(400, rng)  # upstream of insert
    insert_seq = chr1_seq[INSERT_START:INSERT_END]  # 500 bp
    backbone_3 = _random_seq(400, rng)  # downstream of insert
    plasmid_seq = backbone_5 + insert_seq + backbone_3
    plasmid_insert_start = len(backbone_5)
    plasmid_insert_end = plasmid_insert_start + len(insert_seq)

    _write_genbank(
        outdir / "plasmid.gb",
        "pTEST_SYNTH",
        plasmid_seq,
        plasmid_insert_start,
        plasmid_insert_end,
    )

    # --- Generate contaminated reads ---
    n_plasmid = int(N_READS * CONTAM_FRAC)
    n_human = N_READS - n_plasmid
    pairs: list[tuple[str, str, str]] = []
    for i in range(n_plasmid):
        r1, r2 = _generate_read_pair(plasmid_seq, rng)
        pairs.append((f"contam_plasmid_{i:04d}", r1, r2))
    for i in range(n_human):
        source = chr1_seq if rng.random() < 0.5 else chr2_seq
        r1, r2 = _generate_read_pair(source, rng)
        pairs.append((f"contam_human_{i:04d}", r1, r2))
    rng.shuffle(pairs)
    _write_fastq(
        outdir / "contaminated_R1.fastq",
        outdir / "contaminated_R2.fastq",
        pairs,
    )

    # --- Generate clean reads ---
    pairs_clean: list[tuple[str, str, str]] = []
    for i in range(N_READS):
        source = chr1_seq if rng.random() < 0.5 else chr2_seq
        r1, r2 = _generate_read_pair(source, rng)
        pairs_clean.append((f"clean_human_{i:04d}", r1, r2))
    _write_fastq(
        outdir / "not_contaminated_R1.fastq",
        outdir / "not_contaminated_R2.fastq",
        pairs_clean,
    )

    # --- Expected verdict files ---
    expected = outdir / "expected"
    expected.mkdir(exist_ok=True)
    (expected / "contaminated.verdict.txt").write_text("contaminated\n")
    (expected / "not_contaminated.verdict.txt").write_text("not contaminated\n")

    print(f"Synthetic test data generated in {outdir}")
    for f in sorted(outdir.rglob("*")):
        if f.is_file():
            size = f.stat().st_size
            print(f"  {f.relative_to(outdir)}  ({size:,} bytes)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic test data for PlasmiCheck")
    parser.add_argument("--outdir", type=Path, default=None, help="Output directory")
    args = parser.parse_args()
    main(args.outdir)
