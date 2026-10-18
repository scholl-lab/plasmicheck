#!/usr/bin/env python3
"""Tests for in-silico contamination simulator."""

from __future__ import annotations

import gzip
import random
import subprocess
import sys
from pathlib import Path

import pytest

# Add scripts directory to sys.path so simulate_contamination can be imported as a top-level module
SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from simulate_contamination import (  # noqa: E402
    classify_pair,
    extract_background,
    main,
    mutate_seq,
    overlap,
    reverse_complement,
    simulate_pairs,
    write_pairs,
)


def _make_bam_from_sam(sam_text: str, out_bam: Path) -> None:
    """Helper to compile a SAM text into a BAM file using samtools view."""
    sam_path = out_bam.with_suffix(".sam")
    sam_path.write_text(sam_text)
    subprocess.run(["samtools", "view", "-b", "-o", str(out_bam), str(sam_path)], check=True)
    sam_path.unlink()


def test_deterministic_seeding() -> None:
    """Simulating with the same seed produces identical pairs; different seeds differ."""
    seq = "ACGT" * 250  # 1000 bp
    targets = [(200, 400)]

    pairs1 = list(
        simulate_pairs(
            seq,
            targets,
            20,
            circular=True,
            read_length=50,
            fragment_mean=150,
            fragment_sd=20,
            min_overlap=30,
            epsilon=0.1,
            error_rate=0.01,
            rng=random.Random(12345),
        )
    )
    pairs2 = list(
        simulate_pairs(
            seq,
            targets,
            20,
            circular=True,
            read_length=50,
            fragment_mean=150,
            fragment_sd=20,
            min_overlap=30,
            epsilon=0.1,
            error_rate=0.01,
            rng=random.Random(12345),
        )
    )
    pairs_other = list(
        simulate_pairs(
            seq,
            targets,
            20,
            circular=True,
            read_length=50,
            fragment_mean=150,
            fragment_sd=20,
            min_overlap=30,
            epsilon=0.1,
            error_rate=0.01,
            rng=random.Random(99999),
        )
    )

    assert len(pairs1) == 20
    assert pairs1 == pairs2
    assert pairs1 != pairs_other


def test_step_retention_epsilon_zero() -> None:
    """When epsilon=0.0, every yielded fragment must satisfy overlap >= min_overlap."""
    seq = "ACGT" * 250  # 1000 bp
    targets = [(200, 500)]  # 300 bp target
    min_overlap = 45

    pairs = list(
        simulate_pairs(
            seq,
            targets,
            30,
            circular=True,
            read_length=50,
            fragment_mean=150,
            fragment_sd=20,
            min_overlap=min_overlap,
            epsilon=0.0,
            error_rate=0.0,
            rng=random.Random(42),
        )
    )

    assert len(pairs) == 30
    for start, length, _, _, _ in pairs:
        ov = overlap(start, length, targets, len(seq), circular=True)
        assert ov >= min_overlap, f"Fragment at {start} len {length} has ov {ov} < {min_overlap}"


def test_step_retention_epsilon_one() -> None:
    """When epsilon=1.0, backbone fragments (overlap=0) are retained alongside on-target."""
    seq = "ACGT" * 500  # 2000 bp
    targets = [(200, 300)]  # small 100 bp target in 2000 bp plasmid

    pairs = list(
        simulate_pairs(
            seq,
            targets,
            50,
            circular=True,
            read_length=50,
            fragment_mean=120,
            fragment_sd=10,
            min_overlap=45,
            epsilon=1.0,
            error_rate=0.0,
            rng=random.Random(42),
        )
    )

    overlaps = [overlap(s, flen, targets, len(seq), circular=True) for s, flen, _, _, _ in pairs]
    # At epsilon=1.0 on a 100bp target in 2000bp, most fragments will be backbone-only (overlap == 0)
    assert any(ov == 0 for ov in overlaps)


def test_impossible_parameter_rejection() -> None:
    """Simulator raises ValueError on impossible or invalid parameters."""
    seq = "ACGT" * 50  # 200 bp

    # Negative count
    with pytest.raises(ValueError, match="count must be non-negative"):
        list(
            simulate_pairs(
                seq,
                [(10, 50)],
                -1,
                circular=True,
                read_length=50,
                fragment_mean=100,
                fragment_sd=10,
                min_overlap=20,
                epsilon=0.1,
                error_rate=0.0,
                rng=random.Random(1),
            )
        )

    # Empty sequence
    with pytest.raises(ValueError, match="sequence must not be empty"):
        list(
            simulate_pairs(
                "",
                [(10, 50)],
                5,
                circular=True,
                read_length=50,
                fragment_mean=100,
                fragment_sd=10,
                min_overlap=20,
                epsilon=0.1,
                error_rate=0.0,
                rng=random.Random(1),
            )
        )

    # Read length > sequence length
    with pytest.raises(ValueError, match=r"read_length .* exceeds sequence length"):
        list(
            simulate_pairs(
                seq,
                [(10, 50)],
                5,
                circular=True,
                read_length=250,
                fragment_mean=300,
                fragment_sd=10,
                min_overlap=20,
                epsilon=0.1,
                error_rate=0.0,
                rng=random.Random(1),
            )
        )

    # fragment_mean < read_length
    with pytest.raises(ValueError, match=r"fragment_mean .* must be >= read_length"):
        list(
            simulate_pairs(
                seq,
                [(10, 50)],
                5,
                circular=True,
                read_length=100,
                fragment_mean=80,
                fragment_sd=10,
                min_overlap=20,
                epsilon=0.1,
                error_rate=0.0,
                rng=random.Random(1),
            )
        )

    # Invalid epsilon
    with pytest.raises(ValueError, match="epsilon must be in"):
        list(
            simulate_pairs(
                seq,
                [(10, 50)],
                5,
                circular=True,
                read_length=50,
                fragment_mean=100,
                fragment_sd=10,
                min_overlap=20,
                epsilon=1.5,
                error_rate=0.0,
                rng=random.Random(1),
            )
        )

    # Invalid error_rate
    with pytest.raises(ValueError, match="error_rate must be in"):
        list(
            simulate_pairs(
                seq,
                [(10, 50)],
                5,
                circular=True,
                read_length=50,
                fragment_mean=100,
                fragment_sd=10,
                min_overlap=20,
                epsilon=0.1,
                error_rate=-0.1,
                rng=random.Random(1),
            )
        )

    # Invalid min_overlap <= 0
    with pytest.raises(ValueError, match="min_overlap must be positive"):
        list(
            simulate_pairs(
                seq,
                [(10, 50)],
                5,
                circular=True,
                read_length=50,
                fragment_mean=100,
                fragment_sd=10,
                min_overlap=0,
                epsilon=0.1,
                error_rate=0.0,
                rng=random.Random(1),
            )
        )

    # Impossible capture setting: epsilon=0.0 with empty targets
    with pytest.raises(
        ValueError, match=r"Impossible capture setting: epsilon is 0\.0 and targets is empty"
    ):
        list(
            simulate_pairs(
                seq,
                [],
                5,
                circular=True,
                read_length=50,
                fragment_mean=100,
                fragment_sd=10,
                min_overlap=20,
                epsilon=0.0,
                error_rate=0.0,
                rng=random.Random(1),
            )
        )

    # Impossible capture setting: epsilon=0.0 but max overlap < min_overlap
    with pytest.raises(
        ValueError, match=r"Impossible capture setting: max overlap .* < min_overlap"
    ):
        list(
            simulate_pairs(
                seq,
                [(10, 20)],
                5,
                circular=True,
                read_length=50,
                fragment_mean=100,
                fragment_sd=10,
                min_overlap=30,
                epsilon=0.0,
                error_rate=0.0,
                rng=random.Random(1),
            )
        )


def test_circular_origin_wrapping() -> None:
    """Fragments wrapping around circular origin (position 0) reconstruct correct sequences."""
    # Plasmid of length 500: head is 100 'A', body is 300 'C', tail is 100 'G'
    seq = ("A" * 100) + ("C" * 300) + ("G" * 100)
    source_len = len(seq)
    targets = [(0, 500)]

    pairs = list(
        simulate_pairs(
            seq,
            targets,
            100,
            circular=True,
            read_length=40,
            fragment_mean=120,
            fragment_sd=20,
            min_overlap=20,
            epsilon=1.0,
            error_rate=0.0,
            rng=random.Random(42),
        )
    )

    wrapping_pairs = [p for p in pairs if p[0] + p[1] > source_len]
    assert len(wrapping_pairs) > 0, "Expected at least one origin-wrapping pair"

    for start, flen, is_rev, r1, r2 in wrapping_pairs:
        # Expected fragment sequence
        expected_frag = seq[start:] + seq[: flen - (source_len - start)]
        assert len(expected_frag) == flen
        oriented = reverse_complement(expected_frag) if is_rev else expected_frag
        assert r1 == oriented[:40]
        assert r2 == reverse_complement(oriented[-40:])

        # classify_pair identifies origin_wrap
        flags = classify_pair(start, flen, is_rev, 40, targets, source_len, circular=True)
        assert "origin_wrap" in flags


def test_reverse_complement_pairing() -> None:
    """Both forward and reverse fragment orientations produce inward-facing paired reads."""
    seq = "A" * 100 + "C" * 100 + "G" * 100 + "T" * 100  # 400 bp
    targets = [(0, 400)]

    pairs = list(
        simulate_pairs(
            seq,
            targets,
            50,
            circular=False,
            read_length=30,
            fragment_mean=100,
            fragment_sd=0,
            min_overlap=20,
            epsilon=1.0,
            error_rate=0.0,
            rng=random.Random(101),
        )
    )

    has_forward = any(not is_rev for _, _, is_rev, _, _ in pairs)
    has_reverse = any(is_rev for _, _, is_rev, _, _ in pairs)
    assert has_forward and has_reverse, "Expected both forward and reverse oriented fragments"

    for start, flen, is_rev, r1, r2 in pairs:
        frag = seq[start : start + flen]
        if not is_rev:
            assert r1 == frag[:30]
            assert r2 == reverse_complement(frag[-30:])
        else:
            rc_frag = reverse_complement(frag)
            assert r1 == rc_frag[:30]
            assert r2 == reverse_complement(rc_frag[-30:])


def test_error_injection() -> None:
    """Mutate_seq introduces substitutions at requested rate with valid bases."""
    seq = "A" * 1000
    rng = random.Random(42)

    # 0 error rate
    assert mutate_seq(seq, 0.0, rng) == seq

    # 15% error rate
    mutated = mutate_seq(seq, 0.15, rng)
    diffs = sum(1 for a, b in zip(seq, mutated, strict=True) if a != b)
    # Check within reasonable binomial interval for N=1000, p=0.15 (150 +/- 30)
    assert 110 <= diffs <= 190
    assert all(c in "ACGT" for c in mutated)


def test_multi_exon_boundary_truth() -> None:
    """Classify_pair flags splice junctions, vector boundaries, and fragment bridges accurately."""
    # Synthetic plasmid of 1000 bp with 3-exon insert:
    # Exon 1: 200..350 (150 bp)
    # Exon 2: 350..500 (150 bp)
    # Exon 3: 500..650 (150 bp)
    # Splice junctions: 350, 500
    # Vector-insert boundaries: 200, 650
    targets = [(200, 350), (350, 500), (500, 650)]
    source_len = 1000
    read_len = 100

    # 1. Fragment where Mate 1 directly spans splice junction 350:
    # starts at 300, length 250 -> Mate 1 spans 300..400 (contains 349 and 350)
    flags1 = classify_pair(300, 250, False, read_len, targets, source_len, circular=True)
    assert "r1_splice_junction" in flags1
    assert "insert_internal" in flags1

    # 2. Fragment where unsequenced bridge spans junction 350:
    # starts at 220, length 260 -> Mate 1 is 220..320, Mate 2 is 380..480.
    # Bridge is 320..380, spanning 350!
    flags2 = classify_pair(220, 260, False, read_len, targets, source_len, circular=True)
    assert "fragment_bridge" in flags2
    assert "splice_junction_bridge" in flags2
    assert "r1_splice_junction" not in flags2
    assert "r2_splice_junction" not in flags2

    # 3. Fragment spanning vector-insert boundary 200 in Mate 1:
    # starts at 150, length 200 -> Mate 1 is 150..250 (contains 199 and 200)
    flags3 = classify_pair(150, 200, False, read_len, targets, source_len, circular=True)
    assert "vector_boundary" in flags3
    assert "r1_vector_boundary" in flags3

    # 4. Backbone-only fragment: starts at 750, length 150 -> 750..900 (overlap 0)
    flags4 = classify_pair(750, 150, False, read_len, targets, source_len, circular=True)
    assert "backbone_only" in flags4
    assert "insert_internal" not in flags4
    assert "vector_boundary" not in flags4


def test_write_pairs_qname_prefixing_and_truth(tmp_path: Path) -> None:
    """Write_pairs prefixes read names with pcsim: and writes truth TSV."""
    seq = "ACGT" * 250  # 1000 bp
    targets = [(200, 500)]
    pairs = list(
        simulate_pairs(
            seq,
            targets,
            5,
            circular=True,
            read_length=50,
            fragment_mean=150,
            fragment_sd=10,
            min_overlap=30,
            epsilon=0.0,
            error_rate=0.0,
            rng=random.Random(42),
        )
    )

    r1_file = tmp_path / "test_R1.fastq"
    r2_file = tmp_path / "test_R2.fastq"
    truth_file = tmp_path / "truth.tsv"

    with open(r1_file, "w") as f1, open(r2_file, "w") as f2, open(truth_file, "w") as ft:
        written = write_pairs(
            pairs,
            f1,
            f2,
            ft,
            prefix="pcsim:",
            source="plasmid",
            targets=targets,
            source_length=len(seq),
            circular=True,
            read_length=50,
        )

    assert written == 5

    r1_lines = r1_file.read_text().splitlines()
    r2_lines = r2_file.read_text().splitlines()
    assert len(r1_lines) == 20  # 5 records * 4 lines
    assert len(r2_lines) == 20
    assert r1_lines[0].startswith("@pcsim:plasmid_000001/1")
    assert r2_lines[0].startswith("@pcsim:plasmid_000001/2")

    truth_lines = truth_file.read_text().splitlines()
    assert truth_lines[0] == "read_id\tsource\tstart\tlength\tis_reverse\toverlap\tflags"
    assert len(truth_lines) == 6  # header + 5 pairs
    first_record = truth_lines[1].split("\t")
    assert first_record[0] == "pcsim:plasmid_000001"
    assert first_record[1] == "plasmid"


def test_extract_background_orphan_disposal(tmp_path: Path) -> None:
    """Extract_background retains complete paired reads and discards orphans."""
    bam = tmp_path / "input.bam"
    r1 = tmp_path / "out_R1.fastq.gz"
    r2 = tmp_path / "out_R2.fastq.gz"

    seq20 = "A" * 20
    qual20 = "I" * 20
    sam_text = (
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:2000\n"
        # Complete pair 1
        f"read1\t99\tchr1\t100\t60\t20M\t=\t200\t120\t{seq20}\t{qual20}\n"
        f"read1\t147\tchr1\t200\t60\t20M\t=\t100\t-120\t{seq20}\t{qual20}\n"
        # Complete pair 2
        f"read2\t99\tchr1\t300\t60\t20M\t=\t400\t120\t{seq20}\t{qual20}\n"
        f"read2\t147\tchr1\t400\t60\t20M\t=\t300\t-120\t{seq20}\t{qual20}\n"
        # Orphan read (marked paired with mapped mate, but mate is absent from stream)
        f"orphan\t99\tchr1\t500\t60\t20M\t=\t600\t120\t{seq20}\t{qual20}\n"
    )
    _make_bam_from_sam(sam_text, bam)

    retained = extract_background(bam, r1, r2, reserved_prefix="pcsim:")
    assert retained == 2

    with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
        r1_names = [line.strip() for i, line in enumerate(f1) if i % 4 == 0]
        r2_names = [line.strip() for i, line in enumerate(f2) if i % 4 == 0]

    assert r1_names == r2_names
    assert sorted(r1_names) == ["@read1", "@read2"]


def test_extract_background_qname_collision_rejection(tmp_path: Path) -> None:
    """Extract_background raises ValueError on background reads starting with reserved prefix."""
    bam = tmp_path / "collision.bam"
    r1 = tmp_path / "out_R1.fastq.gz"
    r2 = tmp_path / "out_R2.fastq.gz"

    seq20 = "A" * 20
    qual20 = "I" * 20
    sam_text = (
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:2000\n"
        f"pcsim:host_read\t99\tchr1\t100\t60\t20M\t=\t200\t120\t{seq20}\t{qual20}\n"
        f"pcsim:host_read\t147\tchr1\t200\t60\t20M\t=\t100\t-120\t{seq20}\t{qual20}\n"
    )
    _make_bam_from_sam(sam_text, bam)

    with pytest.raises(ValueError, match="QNAME collision: background read starts with 'pcsim:'"):
        extract_background(bam, r1, r2, reserved_prefix="pcsim:")


def test_extract_background_validation(tmp_path: Path) -> None:
    """Extract_background validates existence and prevents input overwrite."""
    missing = tmp_path / "nonexistent.bam"
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"

    with pytest.raises(FileNotFoundError, match="Input BAM file not found"):
        extract_background(missing, r1, r2)

    fake_bam = tmp_path / "fake.bam"
    fake_bam.write_text("fake")
    with pytest.raises(ValueError, match="Output FASTQ path cannot match input BAM path"):
        extract_background(fake_bam, fake_bam, r2)


def test_valid_multi_member_fastq_gz_output(tmp_path: Path) -> None:
    """Appending plasmid reads to gzipped background FASTQ produces valid multi-member gzip."""
    bam = tmp_path / "test.bam"
    r1 = tmp_path / "out_R1.fastq.gz"
    r2 = tmp_path / "out_R2.fastq.gz"
    truth = tmp_path / "truth.tsv"

    seq20 = "A" * 20
    qual20 = "I" * 20
    sam_text = (
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:2000\n"
        f"bg_read1\t99\tchr1\t100\t60\t20M\t=\t200\t120\t{seq20}\t{qual20}\n"
        f"bg_read1\t147\tchr1\t200\t60\t20M\t=\t100\t-120\t{seq20}\t{qual20}\n"
    )
    _make_bam_from_sam(sam_text, bam)

    # 1. Extract background
    retained = extract_background(bam, r1, r2, reserved_prefix="pcsim:")
    assert retained == 1

    # 2. Spike in plasmid reads via gzip append mode "at"
    plasmid_seq = "C" * 500
    plasmid_pairs = list(
        simulate_pairs(
            plasmid_seq,
            [(100, 300)],
            3,
            circular=True,
            read_length=20,
            fragment_mean=80,
            fragment_sd=5,
            min_overlap=15,
            epsilon=0.0,
            error_rate=0.0,
            rng=random.Random(42),
        )
    )

    with gzip.open(r1, "at") as f1, gzip.open(r2, "at") as f2, open(truth, "w") as ft:
        write_pairs(plasmid_pairs, f1, f2, ft, prefix="pcsim:", source="plasmid")

    # 3. Verify seamless decompression of multi-member gzip
    with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
        r1_names = [line.strip() for i, line in enumerate(f1) if i % 4 == 0]
        r2_names = [line.strip() for i, line in enumerate(f2) if i % 4 == 0]

    assert len(r1_names) == 4  # 1 background + 3 plasmid
    assert len(r2_names) == 4
    assert r1_names[0] == "@bg_read1"
    assert r1_names[1] == "@pcsim:plasmid_000001/1"
    assert r1_names[2] == "@pcsim:plasmid_000002/1"
    assert r1_names[3] == "@pcsim:plasmid_000003/1"


def test_cli_spike_in_and_de_novo(tmp_path: Path) -> None:
    """CLI subcommands spike-in and de-novo execute successfully end-to-end."""
    bam = tmp_path / "host.bam"
    plasmid_fa = tmp_path / "plasmid.fasta"
    host_fa = tmp_path / "host.fasta"

    plasmid_fa.write_text(">plasmid1\n" + "ACGT" * 250 + "\n")
    host_fa.write_text(">chr1\n" + "TGCA" * 250 + "\n")

    seq20 = "A" * 20
    qual20 = "I" * 20
    sam_text = (
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:2000\n"
        f"h1\t99\tchr1\t100\t60\t20M\t=\t200\t120\t{seq20}\t{qual20}\n"
        f"h1\t147\tchr1\t200\t60\t20M\t=\t100\t-120\t{seq20}\t{qual20}\n"
    )
    _make_bam_from_sam(sam_text, bam)

    # Test spike-in CLI
    spike_prefix = tmp_path / "spike"
    rc_spike = main(
        [
            "spike-in",
            "--bam",
            str(bam),
            "--plasmid-fasta",
            str(plasmid_fa),
            "--insert-start",
            "200",
            "--insert-end",
            "500",
            "--plasmid-count",
            "5",
            "--read-length",
            "20",
            "--fragment-mean",
            "60",
            "--fragment-sd",
            "5",
            "--min-overlap",
            "10",
            "--out-prefix",
            str(spike_prefix),
        ]
    )
    assert rc_spike == 0
    assert (tmp_path / "spike_R1.fastq.gz").exists()
    assert (tmp_path / "spike_truth.tsv").exists()

    # Test de-novo CLI
    denovo_prefix = tmp_path / "denovo"
    rc_denovo = main(
        [
            "de-novo",
            "--host-fasta",
            str(host_fa),
            "--plasmid-fasta",
            str(plasmid_fa),
            "--host-exons",
            "100-300,500-700",
            "--insert-exons",
            "200-400,600-800",
            "--host-count",
            "10",
            "--plasmid-count",
            "5",
            "--read-length",
            "20",
            "--fragment-mean",
            "60",
            "--fragment-sd",
            "5",
            "--min-overlap",
            "10",
            "--out-prefix",
            str(denovo_prefix),
        ]
    )
    assert rc_denovo == 0
    assert (tmp_path / "denovo_R1.fastq.gz").exists()
    assert (tmp_path / "denovo_truth.tsv").exists()
