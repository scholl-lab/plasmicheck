"""Unit tests for plasmicheck.scripts.spliced_alignment."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pysam
import pytest

from plasmicheck.scripts.spliced_alignment import (
    extract_plasmid_cdna_positions,
    find_fasta_file,
    spliced_alignment,
)


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

    @pytest.mark.unit
    def test_minus_strand_cdna_positions(self, tmp_path: Path) -> None:
        """Test minus-strand cDNA positions have flipped lead and trail clips."""
        output_file = tmp_path / "cDNA_positions.txt"
        # Query len: 30 + 100 + 10 = 140. Reverse strand: lead becomes 10, trail becomes 30.
        mock_read = self._create_mock_read(
            cigartuples=[(4, 30), (0, 100), (4, 10)],
            reference_name="chr1",
            reference_start=1000,
            reference_end=1100,
            is_reverse=True,
        )
        mock_bamfile = MagicMock()
        mock_bamfile.fetch.return_value = [mock_read]
        with patch("pysam.AlignmentFile", return_value=mock_bamfile):
            extract_plasmid_cdna_positions("plasmid.fa", "test.bam", str(output_file))

        content = output_file.read_text()
        assert "cDNA start position in plasmid: 10" in content
        assert "cDNA end position in plasmid: 109" in content
        assert "INSERT_REGION: (10, 109)" in content
        assert "Strand: -" in content

    @pytest.mark.unit
    def test_exon_lengths_indel_accumulation(self, tmp_path: Path) -> None:
        """Test that D and =/X indels within an exon are accumulated into that exon."""
        output_file = tmp_path / "cDNA_positions.txt"
        # 50M 2D 48M 100N 60M -> exon 1 is 50+2+48 = 100, exon 2 is 60
        mock_read = self._create_mock_read(
            cigartuples=[(0, 50), (2, 2), (0, 48), (3, 100), (0, 60)],
            reference_name="chr1",
            reference_start=1000,
            reference_end=1260,
        )
        mock_bamfile = MagicMock()
        mock_bamfile.fetch.return_value = [mock_read]
        with patch("pysam.AlignmentFile", return_value=mock_bamfile):
            extract_plasmid_cdna_positions("plasmid.fa", "test.bam", str(output_file))

        content = output_file.read_text()
        assert "SpliceJunctionCount: 1" in content
        assert "ExonCount: 2" in content
        assert "ExonLengths: 100,60" in content
        assert "IntronLengths: 100" in content

    @pytest.mark.unit
    def test_origin_wrapping_cdna_positions(self, tmp_path: Path) -> None:
        """Test circular origin wrap-around support produces cdna_start > cdna_end."""
        output_file = tmp_path / "cDNA_positions.txt"
        # Start piece: 0..1142 on plasmid (1143 bp), ref 1120..2262
        start_piece = self._create_mock_read(
            cigartuples=[(0, 1143), (4, 11189)],
            reference_name="chr1",
            reference_start=1120,
            reference_end=2262,
            is_reverse=False,
        )
        # End piece: 12210..12331 on plasmid (122 bp), ref 1000..1121
        end_piece = self._create_mock_read(
            cigartuples=[(4, 12210), (0, 122)],
            reference_name="chr1",
            reference_start=1000,
            reference_end=1121,
            is_reverse=False,
        )
        mock_bamfile = MagicMock()
        mock_bamfile.fetch.return_value = [start_piece, end_piece]
        with patch("pysam.AlignmentFile", return_value=mock_bamfile):
            extract_plasmid_cdna_positions("plasmid.fa", "test.bam", str(output_file))

        content = output_file.read_text()
        assert "cDNA start position in plasmid: 12210" in content
        assert "cDNA end position in plasmid: 1142" in content
        assert "INSERT_REGION: (12210, 1142)" in content
        assert "GenomicStart: 1000" in content
        assert "GenomicEnd: 2262" in content


class TestSplicedAlignmentFunction:
    @pytest.mark.unit
    def test_empty_region_value_error(self, tmp_path: Path) -> None:
        """BAM with no mapped alignments should raise ValueError."""
        bam_path = tmp_path / "empty.bam"
        with pysam.AlignmentFile(str(bam_path), "wb", header={"SQ": [{"SN": "chr1", "LN": 1000}]}):
            pass
        pysam.index(str(bam_path))

        with (
            patch("plasmicheck.scripts.spliced_alignment.run_command"),
            pytest.raises(ValueError, match="No human alignment regions found"),
        ):
            spliced_alignment("human.mmi", "plasmid.fa", str(bam_path))

    @pytest.mark.unit
    def test_target_chromosome_filtering_and_interval_merging(self, tmp_path: Path) -> None:
        """Target chromosome filtering selects chrom with spliced span and merges overlapping intervals."""
        bam_path = tmp_path / "aligned.bam"
        header = {"SQ": [{"SN": "chr1", "LN": 100000}, {"SN": "chr17", "LN": 100000}]}
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
            # Spliced read 1 on chr1: 1000..1500 (span 500)
            a1 = pysam.AlignedSegment()
            a1.query_name = "p"
            a1.reference_id = 0
            a1.reference_start = 1000
            a1.cigartuples = [(0, 100), (3, 300), (0, 100)]
            out.write(a1)
            # Spliced read 2 on chr1: 1200..1800 (overlaps with a1)
            a2 = pysam.AlignedSegment()
            a2.query_name = "p"
            a2.reference_id = 0
            a2.reference_start = 1200
            a2.cigartuples = [(0, 100), (3, 400), (0, 100)]
            out.write(a2)
            # Unspliced read on chr17 (e.g. vector backbone hit)
            a3 = pysam.AlignedSegment()
            a3.query_name = "p"
            a3.reference_id = 1
            a3.reference_start = 5000
            a3.cigartuples = [(0, 50)]
            out.write(a3)

        pysam.index(str(bam_path))

        with patch("plasmicheck.scripts.spliced_alignment.run_command"):
            merged = spliced_alignment("human.mmi", "plasmid.fa", str(bam_path), padding=100)

        # chr17 should be excluded; chr1 intervals should be merged
        assert len(merged) == 1
        assert merged[0][0] == "chr1"
        assert merged[0][1] == 900  # 1000 - 100
        assert merged[0][2] == 1900  # 1800 + 100

    @pytest.mark.unit
    def test_threads_parameters_passed_to_commands(self, tmp_path: Path) -> None:
        """Custom thread parameters are passed to minimap2 and samtools commands."""
        bam_path = tmp_path / "aligned.bam"
        header = {"SQ": [{"SN": "chr1", "LN": 100000}]}
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
            a = pysam.AlignedSegment()
            a.query_name = "p"
            a.reference_id = 0
            a.reference_start = 1000
            a.cigartuples = [(0, 100)]
            out.write(a)
        pysam.index(str(bam_path))

        with patch("plasmicheck.scripts.spliced_alignment.run_command") as mock_cmd:
            spliced_alignment(
                "human.mmi", "plasmid.fa", str(bam_path), threads_minimap2=12, threads_samtools=6
            )
            cmd = mock_cmd.call_args_list[0][0][0]
            assert "minimap2 -t 12" in cmd
            assert "samtools view -@ 6" in cmd
            assert "samtools sort -@ 6" in cmd
