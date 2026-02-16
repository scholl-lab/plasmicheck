"""Unit tests for plasmicheck.scripts.spliced_alignment."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

from plasmicheck.scripts.spliced_alignment import extract_plasmid_cdna_positions, find_fasta_file


class TestFindFastaFile:
    @pytest.mark.unit
    def test_finds_fasta_extension(self, tmp_path: Path) -> None:
        base = tmp_path / "reference"
        (tmp_path / "reference.fasta").write_text(">seq\nATCG")
        result = find_fasta_file(str(base))
        assert result.endswith(".fasta")

    @pytest.mark.unit
    def test_finds_fa_extension(self, tmp_path: Path) -> None:
        base = tmp_path / "reference"
        (tmp_path / "reference.fa").write_text(">seq\nATCG")
        result = find_fasta_file(str(base))
        assert result.endswith(".fa")

    @pytest.mark.unit
    def test_no_matching_file(self, tmp_path: Path) -> None:
        base = tmp_path / "reference"
        with pytest.raises(FileNotFoundError, match="No FASTA file found"):
            find_fasta_file(str(base))

    @pytest.mark.unit
    def test_prefers_first_extension(self, tmp_path: Path) -> None:
        """When multiple extensions exist, returns the first match from FASTA_EXTENSIONS."""
        base = tmp_path / "reference"
        # .fasta comes before .fa in the extension list
        (tmp_path / "reference.fasta").write_text(">seq1\nATCG")
        (tmp_path / "reference.fa").write_text(">seq2\nGCTA")
        result = find_fasta_file(str(base))
        assert result.endswith(".fasta")


class TestExtractJunctionMetadata:
    """Test junction metadata extraction from BAM CIGAR strings."""

    def _create_mock_read(
        self,
        *,
        is_unmapped: bool = False,
        is_secondary: bool = False,
        is_supplementary: bool = False,
        cigartuples: list[tuple[int, int]] | None = None,
        is_reverse: bool = False,
        reference_name: str = "chr1",
        reference_start: int = 1000,
        reference_end: int = 2000,
    ) -> Any:
        """Create a mock pysam read with specified attributes."""
        read = MagicMock()
        read.is_unmapped = is_unmapped
        read.is_secondary = is_secondary
        read.is_supplementary = is_supplementary
        read.cigartuples = cigartuples
        read.is_reverse = is_reverse
        read.reference_name = reference_name
        read.reference_start = reference_start
        read.reference_end = reference_end
        return read

    @pytest.mark.unit
    def test_two_exons_one_junction(self, tmp_path: Path) -> None:
        """Test extraction of junction metadata from BAM with 2 exons and 1 junction."""
        output_file = tmp_path / "cDNA_positions.txt"

        # Mock BAM with primary alignment: 50M (exon1) 100N (intron) 50M (exon2)
        mock_read = self._create_mock_read(
            cigartuples=[(0, 50), (3, 100), (0, 50)],
            reference_name="chr1",
            reference_start=1000,
            reference_end=1200,
            is_reverse=False,
        )

        mock_bamfile = MagicMock()
        mock_bamfile.fetch.return_value = [mock_read]
        mock_bamfile.__enter__.return_value = mock_bamfile
        mock_bamfile.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bamfile):
            extract_plasmid_cdna_positions("plasmid.fa", "test.bam", str(output_file))

        content = output_file.read_text()
        assert "SpliceJunctionCount: 1" in content
        assert "ExonCount: 2" in content
        assert "ExonLengths: 50,50" in content
        assert "IntronLengths: 100" in content
        assert "Strand: +" in content
        assert "Chromosome: chr1" in content
        assert "GenomicStart: 1000" in content
        assert "GenomicEnd: 1200" in content

    @pytest.mark.unit
    def test_single_exon_no_junctions(self, tmp_path: Path) -> None:
        """Test extraction with single exon (no junctions)."""
        output_file = tmp_path / "cDNA_positions.txt"

        # Mock BAM with primary alignment: 100M (single exon, no junctions)
        mock_read = self._create_mock_read(
            cigartuples=[(0, 100)],
            reference_name="chr2",
            reference_start=5000,
            reference_end=5100,
            is_reverse=True,
        )

        mock_bamfile = MagicMock()
        mock_bamfile.fetch.return_value = [mock_read]
        mock_bamfile.__enter__.return_value = mock_bamfile
        mock_bamfile.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bamfile):
            extract_plasmid_cdna_positions("plasmid.fa", "test.bam", str(output_file))

        content = output_file.read_text()
        assert "SpliceJunctionCount: 0" in content
        assert "ExonCount: 1" in content
        assert "ExonLengths: 100" in content
        assert "IntronLengths: NA" in content
        assert "Strand: -" in content
        assert "Chromosome: chr2" in content
        assert "GenomicStart: 5000" in content
        assert "GenomicEnd: 5100" in content

    @pytest.mark.unit
    def test_secondary_supplementary_alignments_skipped(self, tmp_path: Path) -> None:
        """Test that secondary/supplementary alignments are skipped for junction counting."""
        output_file = tmp_path / "cDNA_positions.txt"

        # Primary alignment with 1 junction
        primary = self._create_mock_read(
            cigartuples=[(0, 50), (3, 100), (0, 50)],
            reference_name="chr1",
            reference_start=1000,
            reference_end=1200,
        )

        # Secondary alignment with different junction structure (should be ignored)
        secondary = self._create_mock_read(
            is_secondary=True,
            cigartuples=[(0, 30), (3, 200), (0, 30), (3, 150), (0, 30)],
            reference_name="chr1",
            reference_start=2000,
            reference_end=2440,
        )

        # Supplementary alignment (should be ignored)
        supplementary = self._create_mock_read(
            is_supplementary=True,
            cigartuples=[(0, 100)],
            reference_name="chr2",
            reference_start=3000,
            reference_end=3100,
        )

        mock_bamfile = MagicMock()
        mock_bamfile.fetch.return_value = [secondary, supplementary, primary]
        mock_bamfile.__enter__.return_value = mock_bamfile
        mock_bamfile.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bamfile):
            extract_plasmid_cdna_positions("plasmid.fa", "test.bam", str(output_file))

        content = output_file.read_text()
        # Should use primary alignment metadata (1 junction, 2 exons)
        assert "SpliceJunctionCount: 1" in content
        assert "ExonCount: 2" in content
        assert "ExonLengths: 50,50" in content
        assert "IntronLengths: 100" in content

    @pytest.mark.unit
    def test_strand_and_coordinates_extraction(self, tmp_path: Path) -> None:
        """Test correct extraction of strand and genomic coordinates."""
        output_file = tmp_path / "cDNA_positions.txt"

        mock_read = self._create_mock_read(
            cigartuples=[(0, 75), (3, 250), (0, 75)],
            reference_name="chr3",
            reference_start=10000,
            reference_end=10400,
            is_reverse=True,
        )

        mock_bamfile = MagicMock()
        mock_bamfile.fetch.return_value = [mock_read]
        mock_bamfile.__enter__.return_value = mock_bamfile
        mock_bamfile.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bamfile):
            extract_plasmid_cdna_positions("plasmid.fa", "test.bam", str(output_file))

        content = output_file.read_text()
        assert "Strand: -" in content
        assert "Chromosome: chr3" in content
        assert "GenomicStart: 10000" in content
        assert "GenomicEnd: 10400" in content

    @pytest.mark.unit
    def test_file_format_all_fields_present(self, tmp_path: Path) -> None:
        """Test that output file contains all 11 expected fields."""
        output_file = tmp_path / "cDNA_positions.txt"

        mock_read = self._create_mock_read(
            cigartuples=[(0, 50), (3, 100), (0, 50)],
        )

        mock_bamfile = MagicMock()
        mock_bamfile.fetch.return_value = [mock_read]
        mock_bamfile.__enter__.return_value = mock_bamfile
        mock_bamfile.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bamfile):
            extract_plasmid_cdna_positions("plasmid.fa", "test.bam", str(output_file))

        lines = output_file.read_text().strip().split("\n")
        assert len(lines) == 11

        # Verify all expected keys are present
        expected_keys = [
            "cDNA start position in plasmid:",
            "cDNA end position in plasmid:",
            "INSERT_REGION:",
            "SpliceJunctionCount:",
            "ExonCount:",
            "ExonLengths:",
            "IntronLengths:",
            "Strand:",
            "Chromosome:",
            "GenomicStart:",
            "GenomicEnd:",
        ]
        for key in expected_keys:
            assert any(key in line for line in lines), f"Missing key: {key}"

    @pytest.mark.unit
    def test_unmapped_plasmid_produces_na_values(self, tmp_path: Path) -> None:
        """Test that unmapped plasmids produce NA for all junction fields."""
        output_file = tmp_path / "cDNA_positions.txt"

        # Mock BAM with only unmapped reads
        mock_unmapped = self._create_mock_read(is_unmapped=True)

        mock_bamfile = MagicMock()
        mock_bamfile.fetch.return_value = [mock_unmapped]
        mock_bamfile.__enter__.return_value = mock_bamfile
        mock_bamfile.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bamfile):
            extract_plasmid_cdna_positions("plasmid.fa", "test.bam", str(output_file))

        content = output_file.read_text()
        # All fields should be NA
        assert "cDNA start position in plasmid: NA" in content
        assert "cDNA end position in plasmid: NA" in content
        assert "INSERT_REGION: NA" in content
        assert "SpliceJunctionCount: NA" in content
        assert "ExonCount: NA" in content
        assert "ExonLengths: NA" in content
        assert "IntronLengths: NA" in content
        assert "Strand: NA" in content
        assert "Chromosome: NA" in content
        assert "GenomicStart: NA" in content
        assert "GenomicEnd: NA" in content
