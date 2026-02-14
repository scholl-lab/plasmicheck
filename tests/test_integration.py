"""Integration tests using synthetic test data.

These tests exercise the real pipeline end-to-end against tiny synthetic
FASTA/FASTQ/GenBank files (~50 KB total).  They require ``minimap2`` and
``samtools`` on PATH and are marked ``@pytest.mark.integration``.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

# Skip entire module when external tools are unavailable.
pytestmark = pytest.mark.integration

_HAS_MINIMAP2 = shutil.which("minimap2") is not None
_HAS_SAMTOOLS = shutil.which("samtools") is not None
_SKIP_REASON = "minimap2 and/or samtools not found on PATH"


@pytest.fixture
def output_dir(tmp_path: Path) -> Path:
    d = tmp_path / "output"
    d.mkdir()
    return d


@pytest.mark.skipif(not (_HAS_MINIMAP2 and _HAS_SAMTOOLS), reason=_SKIP_REASON)
class TestConvertGenBank:
    def test_convert_plasmid_to_fasta(self, synthetic_data_dir: Path, tmp_path: Path) -> None:
        from plasmicheck.scripts.convert_plasmidfile_to_fasta import convert_plasmidfile_to_fasta

        gb = str(synthetic_data_dir / "plasmid.gb")
        out = str(tmp_path / "plasmid.fasta")
        convert_plasmidfile_to_fasta(gb, out, "genbank", shift_bases=500, generate_shifted=False)
        result = Path(out)
        assert result.exists()
        content = result.read_text()
        assert content.startswith(">")
        assert len(content) > 100


@pytest.mark.skipif(not (_HAS_MINIMAP2 and _HAS_SAMTOOLS), reason=_SKIP_REASON)
class TestIndexCreation:
    def test_create_indexes_for_human_ref(self, synthetic_data_dir: Path, tmp_path: Path) -> None:
        import shutil as _shutil

        from plasmicheck.scripts.create_indexes import create_indexes

        # Copy to tmp so we don't pollute test data dir
        fasta = tmp_path / "human_ref.fasta"
        _shutil.copy(synthetic_data_dir / "human_ref.fasta", fasta)

        create_indexes(str(fasta), overwrite=True)

        mmi = tmp_path / "human_ref.mmi"
        fai = tmp_path / "human_ref.fasta.fai"
        assert mmi.exists(), ".mmi index not created"
        assert fai.exists(), ".fai index not created"


@pytest.mark.skipif(not (_HAS_MINIMAP2 and _HAS_SAMTOOLS), reason=_SKIP_REASON)
class TestPipelineContaminated:
    def test_pipeline_detects_contamination(
        self, synthetic_data_dir: Path, output_dir: Path
    ) -> None:
        from plasmicheck.scripts.run_pipeline import run_pipeline

        human_fasta = str(synthetic_data_dir / "human_ref.fasta")
        plasmid_gb = str(synthetic_data_dir / "plasmid.gb")
        r1 = str(synthetic_data_dir / "contaminated_R1.fastq")
        r2 = str(synthetic_data_dir / "contaminated_R2.fastq")

        run_pipeline(
            human_fasta=human_fasta,
            plasmid_files=plasmid_gb,
            output_folder=str(output_dir),
            keep_intermediate=True,
            overwrite=True,
            threshold=0.8,
            sequencing_files_r1=r1,
            sequencing_files_r2=r2,
        )

        # Find the summary file
        summaries = list(output_dir.rglob("*.summary.tsv"))
        assert len(summaries) >= 1, "No summary.tsv found"
        content = summaries[0].read_text()
        assert "Verdict" in content
        # With 60% plasmid reads and threshold 0.8, expect contamination
        assert "contaminated" in content.lower()


@pytest.mark.skipif(not (_HAS_MINIMAP2 and _HAS_SAMTOOLS), reason=_SKIP_REASON)
class TestPipelineClean:
    def test_pipeline_detects_no_contamination(
        self, synthetic_data_dir: Path, output_dir: Path
    ) -> None:
        from plasmicheck.scripts.run_pipeline import run_pipeline

        human_fasta = str(synthetic_data_dir / "human_ref.fasta")
        plasmid_gb = str(synthetic_data_dir / "plasmid.gb")
        r1 = str(synthetic_data_dir / "not_contaminated_R1.fastq")
        r2 = str(synthetic_data_dir / "not_contaminated_R2.fastq")

        run_pipeline(
            human_fasta=human_fasta,
            plasmid_files=plasmid_gb,
            output_folder=str(output_dir),
            keep_intermediate=True,
            overwrite=True,
            threshold=0.8,
            sequencing_files_r1=r1,
            sequencing_files_r2=r2,
        )

        summaries = list(output_dir.rglob("*.summary.tsv"))
        assert len(summaries) >= 1, "No summary.tsv found"
        content = summaries[0].read_text()
        assert "Verdict" in content
        assert "not contaminated" in content.lower()
