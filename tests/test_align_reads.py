"""Exercise shell argument handling without installed alignment tools."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pysam
import pytest

from plasmicheck.scripts.align_reads import (
    align_reads,
    bam_to_fastq_cmd,
    has_bam_index,
    parse_bed_regions,
)
from plasmicheck.scripts.spliced_alignment import spliced_alignment
from plasmicheck.scripts.utils import run_command


def _run_with_fake_tools(command: str) -> subprocess.CompletedProcess[str]:
    return run_command(
        'minimap2() { printf \'%s\\n\' "$@" > "$COMMAND_LOG/minimap2"; }; '
        'samtools() { printf \'%s\\n\' "$@" > "$COMMAND_LOG/samtools-$1"; }; ' + command
    )


@pytest.fixture
def command_log(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("COMMAND_LOG", str(tmp_path))
    monkeypatch.setattr("plasmicheck.scripts.align_reads.run_command", _run_with_fake_tools)
    monkeypatch.setattr("plasmicheck.scripts.spliced_alignment.run_command", _run_with_fake_tools)
    monkeypatch.setattr(
        "plasmicheck.scripts.utils._cfg", {"retry_settings": {"retries": 1, "delay": 0}}
    )
    return tmp_path


@pytest.mark.unit
@pytest.mark.parametrize("input_kind", ["bam", "fastq", "paired"])
@pytest.mark.parametrize("suffix", [" $(printf injected > injected)", " quote's ; &"])
def test_alignment_preserves_shell_metacharacters(
    command_log: Path, input_kind: str, suffix: str
) -> None:
    reference = f"reference{suffix}.mmi"
    input_file = f"reads{suffix}.{'bam' if input_kind == 'bam' else 'fastq'}"
    fastq2 = f"reverse{suffix}.fastq" if input_kind == "paired" else None
    output_bam = f"output{suffix}.bam"
    memory = f"2G{suffix}"

    align_reads(reference, input_file, output_bam, "plasmid", fastq2, 2, 1, memory)

    minimap_args = (command_log / "minimap2").read_text().splitlines()
    expected_reads = ["-"] if input_kind == "bam" else [input_file]
    if fastq2:
        expected_reads.append(fastq2)
    assert minimap_args == ["-t", "2", "-ax", "sr", reference, *expected_reads]
    if input_kind == "bam":
        assert (command_log / "samtools-collate").read_text().splitlines() == [
            "collate",
            "-@",
            "1",
            "-u",
            "-O",
            input_file,
        ]
        assert (command_log / "samtools-fastq").read_text().splitlines() == [
            "fastq",
            "-@",
            "1",
            "-n",
            "-F",
            "0x900",
            "-",
        ]
    assert (command_log / "samtools-sort").read_text().splitlines() == [
        "sort",
        "-@",
        "1",
        "-m",
        memory,
        "-o",
        output_bam,
    ]
    assert (command_log / "samtools-index").read_text().splitlines() == ["index", output_bam]
    assert not (command_log / "injected").exists()


@pytest.mark.unit
def test_human_alignment_sorts_by_name_without_index(command_log: Path) -> None:
    align_reads("human.mmi", "reads.fastq", "human_out.bam", "human", None, 2, 1, "2G")
    assert (command_log / "samtools-sort").read_text().splitlines() == [
        "sort",
        "-N",
        "-@",
        "1",
        "-m",
        "2G",
        "-o",
        "human_out.bam",
    ]
    assert not (command_log / "samtools-index").exists()


@pytest.mark.unit
def test_spliced_alignment_preserves_shell_metacharacters(command_log: Path) -> None:
    suffix = " quote's ; $(printf injected > injected)"
    reference = f"human{suffix}.mmi"
    plasmid = f"plasmid{suffix}.fasta"
    output = f"spliced{suffix}.bam"
    with pysam.AlignmentFile(output, "wb", header={"SQ": [{"SN": "ref", "LN": 2000}]}) as bam:
        seg = pysam.AlignedSegment()
        seg.query_name = "plasmid"
        seg.reference_id = 0
        seg.reference_start = 100
        seg.cigartuples = [(0, 50)]
        bam.write(seg)
    pysam.index(output)

    assert spliced_alignment(reference, plasmid, output) == [("ref", 0, 1150)]

    assert (command_log / "minimap2").read_text().splitlines()[-2:] == [reference, plasmid]
    assert (command_log / "samtools-sort").read_text().splitlines()[-2:] == ["-o", output]
    assert (command_log / "samtools-index").read_text().splitlines() == ["index", output]
    assert not (command_log / "injected").exists()


@pytest.mark.unit
def test_parse_bed_regions(tmp_path: Path) -> None:
    bed = tmp_path / "regions.bed"
    bed.write_text(
        "# comment line\n"
        "track name=test\n"
        "browser position chr1:1-100\n"
        "\n"
        "chr1\t0\t100\tgeneA\n"
        "chr11\t128890000\t128920000\n"
        "chr2 -5 50\n"
    )
    assert parse_bed_regions(str(bed)) == [
        "chr1:1-100",
        "chr11:128890001-128920000",
        "chr2:1-50",
    ]

    empty_bed = tmp_path / "empty.bed"
    empty_bed.write_text("# only comments\n\ntrack name=empty\n")
    with pytest.raises(ValueError, match="No valid BED intervals found"):
        parse_bed_regions(str(empty_bed))


@pytest.mark.unit
@pytest.mark.parametrize(
    "index_name", ["sample.bam.bai", "sample.bai", "sample.bam.csi", "sample.csi"]
)
def test_has_bam_index(tmp_path: Path, index_name: str) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    assert not has_bam_index(str(bam))
    (tmp_path / index_name).write_bytes(b"")
    assert has_bam_index(str(bam))


@pytest.mark.unit
def test_bam_to_fastq_cmd(tmp_path: Path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")

    assert bam_to_fastq_cmd(str(bam), 4) == (
        f"samtools collate -@ 4 -u -O {bam} | samtools fastq -@ 4 -n -F 0x900 -"
    )

    with pytest.raises(FileNotFoundError, match="BAM index"):
        bam_to_fastq_cmd(str(bam), 4, regions=["chr1:101-200"])

    (tmp_path / "sample.bam.bai").write_bytes(b"")
    assert bam_to_fastq_cmd(
        str(bam), 2, regions=["chr1:101-200", "chr1:101-200", "chr2:1-50", "*"]
    ) == (
        f"samtools view -@ 2 -u -M {bam} chr1:101-200 chr2:1-50 '*' "
        f"| samtools collate -@ 2 -u -O - "
        f"| samtools fastq -@ 2 -n -F 0x900 -"
    )

    bed = tmp_path / "targets.bed"
    bed.write_text("chr1\t100\t200\n")
    assert bam_to_fastq_cmd(str(bam), 2, regions=str(bed)) == (
        f"samtools view -@ 2 -u -M {bam} chr1:101-200 '*' "
        f"| samtools collate -@ 2 -u -O - "
        f"| samtools fastq -@ 2 -n -F 0x900 -"
    )


@pytest.mark.unit
def test_align_reads_with_regions_on_indexed_bam(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    (tmp_path / "sample.bam.bai").write_bytes(b"")
    commands: list[str] = []
    monkeypatch.setattr("plasmicheck.scripts.align_reads.run_command", commands.append)

    align_reads(
        "ref.mmi",
        str(bam),
        "out.bam",
        "plasmid",
        minimap2_threads=2,
        samtools_threads=1,
        samtools_sort_memory="2G",
        regions=["chr1:101-200"],
        include_unmapped=True,
    )

    assert len(commands) == 2
    assert commands[0].startswith(
        f"samtools view -@ 1 -u -M {bam} chr1:101-200 '*' "
        f"| samtools collate -@ 1 -u -O - "
        f"| samtools fastq -@ 1 -n -F 0x900 - "
        f"| minimap2 -t 2 -ax sr ref.mmi -"
    )
    assert commands[1] == "samtools index out.bam"
