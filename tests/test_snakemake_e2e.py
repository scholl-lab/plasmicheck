"""End-to-end Snakemake cohort workflow run on the synthetic fixture (#129).

Runs the real Snakemake CLI (prepare -> scan per sample -> summary) on two BAMs and
checks the verdicts. Guards the {sample} wildcard against matching rule markers.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

ROOT = Path(__file__).resolve().parents[1]
SYN = ROOT / "tests" / "data" / "synthetic"


def _tool(name: str) -> str:
    path = shutil.which(name)
    assert path, f"{name} is required for the Snakemake end-to-end test"
    return path


def test_cohort_workflow_end_to_end(tmp_path: Path) -> None:
    minimap2, samtools, snakemake = _tool("minimap2"), _tool("samtools"), _tool("snakemake")
    bams = {}
    for sample, prefix in (("SAMPLE_A", "contaminated"), ("SAMPLE_B", "not_contaminated")):
        bam = tmp_path / f"{prefix}.bam"
        aln = subprocess.run(
            [
                minimap2,
                "-ax",
                "sr",
                str(SYN / "human_ref.fasta"),
                str(SYN / f"{prefix}_R1.fastq"),
                str(SYN / f"{prefix}_R2.fastq"),
            ],
            check=True,
            capture_output=True,
        )
        subprocess.run([samtools, "sort", "-o", str(bam), "-"], input=aln.stdout, check=True)
        subprocess.run([samtools, "index", str(bam)], check=True)
        bams[sample] = bam
    shutil.copy(SYN / "human_ref.fasta", tmp_path / "human_ref.fasta")
    (tmp_path / "samples.tsv").write_text(
        "sample\tbam\n" + "".join(f"{s}\t{b}\n" for s, b in bams.items())
    )
    out = tmp_path / "output"
    (tmp_path / "config.yaml").write_text(
        f'human_fasta: "{tmp_path / "human_ref.fasta"}"\n'
        f'samples: "{tmp_path / "samples.tsv"}"\n'
        f'constructs:\n  - "{SYN / "plasmid.gb"}"\n'
        f'output_dir: "{out}"\nsummary_dir: "{out / "summary"}"\n'
        f'log_dir: "{tmp_path / "logs"}"\ncache_dir: "{tmp_path / "cache"}"\n'
        "include_unmapped: true\n"
    )
    res = subprocess.run(
        [
            snakemake,
            "-s",
            str(ROOT / "snakemake" / "plasmicheck.smk"),
            "--configfile",
            str(tmp_path / "config.yaml"),
            "-c",
            "2",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr[-4000:]
    verdicts = (out / "summary" / "data" / "verdict_table.tsv").read_text()
    row_a = next(ln for ln in verdicts.splitlines() if "SAMPLE_A" in ln)
    row_b = next(ln for ln in verdicts.splitlines() if "SAMPLE_B" in ln)
    assert "is contaminated" in row_a and "not contaminated" not in row_a
    assert "not contaminated" in row_b
