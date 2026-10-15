"""Exercise shell argument handling without installed alignment tools."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pysam
import pytest

from plasmicheck.scripts.align_reads import align_reads
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
        assert (command_log / "samtools-fasta").read_text().splitlines() == ["fasta", input_file]
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
def test_spliced_alignment_preserves_shell_metacharacters(command_log: Path) -> None:
    suffix = " quote's ; $(printf injected > injected)"
    reference = f"human{suffix}.mmi"
    plasmid = f"plasmid{suffix}.fasta"
    output = f"spliced{suffix}.bam"
    with pysam.AlignmentFile(output, "wb", header={"SQ": [{"SN": "ref", "LN": 1000}]}):
        pass
    pysam.index(output)

    assert spliced_alignment(reference, plasmid, output) == []

    assert (command_log / "minimap2").read_text().splitlines()[-2:] == [reference, plasmid]
    assert (command_log / "samtools-sort").read_text().splitlines()[-2:] == ["-o", output]
    assert (command_log / "samtools-index").read_text().splitlines() == ["index", output]
    assert not (command_log / "injected").exists()
