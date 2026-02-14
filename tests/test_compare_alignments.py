"""Unit tests for plasmicheck.scripts.compare_alignments."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from plasmicheck.scripts.compare_alignments import (
    _assign,
    _best_read,
    _collate_bam,
    _name_group_bam,
    _resort_supplementary,
    _streaming_compare,
    _write_assignment,
    calculate_alignment_score,
    parse_insert_region,
)


class TestCalculateAlignmentScore:
    @pytest.mark.unit
    def test_unmapped_read_returns_zero(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(is_unmapped=True)
        assert calculate_alignment_score(read) == 0

    @pytest.mark.unit
    def test_simple_mapped_read(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, no clipping, 0 mismatches, mate mapped -> 60 + 10 = 70
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(0, 100)],  # 100M
            nm_tag=0,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 70

    @pytest.mark.unit
    def test_with_soft_clipping(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, 10bp soft clip, 0 mismatches, mate mapped -> 60 - 10 + 10 = 60
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(4, 10), (0, 90)],  # 10S90M
            nm_tag=0,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 60

    @pytest.mark.unit
    def test_with_mismatches(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, no clipping, NM=5, mate mapped -> 60 - 5 + 10 = 65
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=5,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 65

    @pytest.mark.unit
    def test_combined_penalties(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, 10bp soft clip, NM=5, mate mapped -> 60 - 10 - 5 + 10 = 55
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(4, 10), (0, 90)],
            nm_tag=5,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 55

    @pytest.mark.unit
    def test_mate_unmapped_no_bonus(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, no clipping, 0 mismatches, mate UNmapped -> 60 + 0 = 60
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=0,
            mate_is_unmapped=True,
        )
        assert calculate_alignment_score(read) == 60

    @pytest.mark.unit
    def test_hard_clipping(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, 20bp hard clip, 0 mismatches, mate mapped -> 60 - 20 + 10 = 50
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(5, 20), (0, 80)],  # 20H80M
            nm_tag=0,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 50

    @pytest.mark.unit
    def test_no_nm_tag(self, mock_pysam_read: Any) -> None:
        # No NM tag -> 0 mismatches assumed
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=None,  # no NM tag
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 70


class TestParseInsertRegion:
    @pytest.mark.unit
    def test_valid_file(self, cdna_positions_file: Path) -> None:
        result = parse_insert_region(str(cdna_positions_file))
        assert result == (100, 500)

    @pytest.mark.unit
    def test_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match=r"cDNA_positions\.txt"):
            parse_insert_region(str(tmp_path / "nonexistent.txt"))

    @pytest.mark.unit
    @pytest.mark.parametrize("start,end", [(0, 1000), (50, 200), (1, 99999)])
    def test_various_positions(self, tmp_path: Path, start: int, end: int) -> None:
        f = tmp_path / "cDNA_positions.txt"
        f.write_text(
            f"cDNA start position in plasmid: {start}\n"
            f"cDNA end position in plasmid: {end}\n"
            f"INSERT_REGION: ({start}, {end})\n"
        )
        assert parse_insert_region(str(f)) == (start, end)


class TestAssign:
    @pytest.mark.unit
    def test_plasmid_wins(self) -> None:
        assert _assign(70, 30) == "Plasmid"

    @pytest.mark.unit
    def test_human_wins(self) -> None:
        assert _assign(30, 70) == "Human"

    @pytest.mark.unit
    def test_tied(self) -> None:
        assert _assign(50, 50) == "Tied"


class TestBestRead:
    @pytest.mark.unit
    def test_picks_primary(self, mock_pysam_read: Any) -> None:
        primary = mock_pysam_read(query_name="r1")
        primary.is_secondary = False
        primary.is_supplementary = False
        secondary = mock_pysam_read(query_name="r1")
        secondary.is_secondary = True
        secondary.is_supplementary = False
        result = _best_read([secondary, primary])
        assert result is primary

    @pytest.mark.unit
    def test_falls_back_to_first(self, mock_pysam_read: Any) -> None:
        sec1 = mock_pysam_read(query_name="r1")
        sec1.is_secondary = True
        sec1.is_supplementary = False
        sec2 = mock_pysam_read(query_name="r1")
        sec2.is_secondary = True
        sec2.is_supplementary = False
        result = _best_read([sec1, sec2])
        assert result is sec1


class TestWriteAssignment:
    @pytest.mark.unit
    def test_full_data(self, mock_pysam_read: Any) -> None:
        buf = io.StringIO()
        p_read = mock_pysam_read(cigarstring="100M", mapping_quality=60)
        h_read = mock_pysam_read(cigarstring="50M50S", mapping_quality=30)
        _write_assignment(buf, "read1", "Plasmid", 70, 30, p_read, h_read)
        line = buf.getvalue()
        assert line.startswith("read1\tPlasmid\t70\t30\t100M\t50M50S\t60\t30\n")

    @pytest.mark.unit
    def test_missing_human_read(self, mock_pysam_read: Any) -> None:
        buf = io.StringIO()
        p_read = mock_pysam_read(cigarstring="100M", mapping_quality=60)
        _write_assignment(buf, "read1", "Plasmid", 70, 0, p_read, None)
        line = buf.getvalue()
        assert "NA" in line  # human CIGAR/MAPQ should be NA


class TestStreamingCompare:
    @pytest.mark.unit
    def test_matched_reads(self, mock_pysam_read: Any) -> None:
        """Both BAMs have the same reads — should compare and assign."""
        p_read = mock_pysam_read(query_name="r1", mapping_quality=60)
        p_read.is_secondary = False
        p_read.is_supplementary = False
        h_read = mock_pysam_read(query_name="r1", mapping_quality=30)
        h_read.is_secondary = False
        h_read.is_supplementary = False

        with (
            patch(
                "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
                side_effect=[iter([("r1", [p_read])]), iter([("r1", [h_read])])],
            ),
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
            assert counts["Plasmid"] == 1
            assert counts["Human"] == 0

    @pytest.mark.unit
    def test_plasmid_only_read(self, mock_pysam_read: Any) -> None:
        """Read in plasmid BAM only."""
        p_read = mock_pysam_read(query_name="r1", mapping_quality=60)
        p_read.is_secondary = False
        p_read.is_supplementary = False

        with (
            patch(
                "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
                side_effect=[iter([("r1", [p_read])]), iter([])],
            ),
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
            assert counts["Plasmid"] == 1

    @pytest.mark.unit
    def test_human_only_read(self, mock_pysam_read: Any) -> None:
        """Read in human BAM only."""
        h_read = mock_pysam_read(query_name="r1", mapping_quality=60)
        h_read.is_secondary = False
        h_read.is_supplementary = False

        with (
            patch(
                "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
                side_effect=[iter([]), iter([("r1", [h_read])])],
            ),
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
            assert counts["Human"] == 1


class TestCollateBam:
    @pytest.mark.unit
    def test_collate_bam_calls_samtools_collate(self, tmp_path: Path) -> None:
        """Verify _collate_bam calls samtools collate with correct arguments."""
        input_bam = str(tmp_path / "input.bam")
        output_bam = str(tmp_path / "output.bam")

        with (
            patch("subprocess.run") as mock_run,
            patch("plasmicheck.scripts.compare_alignments._resort_supplementary") as mock_resort,
        ):
            mock_run.return_value = None
            _collate_bam(input_bam, output_bam, threads=4)

            # Verify collate was called
            assert mock_run.call_count == 1
            call_args = mock_run.call_args[0][0]
            assert call_args[0] == "samtools"
            assert call_args[1] == "collate"
            assert "-@" in call_args
            assert "4" in call_args
            assert "-o" in call_args
            assert input_bam in call_args
            # Verify no -f (fast mode) flag
            assert "-f" not in call_args

            # Verify _resort_supplementary was called
            assert mock_resort.call_count == 1

    @pytest.mark.unit
    def test_collate_bam_fallback_on_error(self, tmp_path: Path) -> None:
        """Verify fallback to sort -n when collate fails with CalledProcessError."""
        import subprocess

        input_bam = str(tmp_path / "input.bam")
        output_bam = str(tmp_path / "output.bam")

        with (
            patch("subprocess.run") as mock_run,
            patch("plasmicheck.scripts.compare_alignments._namesort_bam_fallback") as mock_fallback,
        ):
            # First call (collate) fails, second call (fallback) succeeds
            mock_run.side_effect = [
                subprocess.CalledProcessError(1, "samtools collate"),
                None,  # fallback call succeeds
            ]

            _collate_bam(input_bam, output_bam)

            # Verify fallback was called
            assert mock_fallback.call_count == 1
            mock_fallback.assert_called_once_with(input_bam, output_bam)

    @pytest.mark.unit
    def test_collate_bam_fallback_on_missing_samtools(self, tmp_path: Path) -> None:
        """Verify fallback to sort -n when samtools is not found."""
        input_bam = str(tmp_path / "input.bam")
        output_bam = str(tmp_path / "output.bam")

        with (
            patch("subprocess.run") as mock_run,
            patch("plasmicheck.scripts.compare_alignments._namesort_bam_fallback") as mock_fallback,
        ):
            # First call (collate) fails with FileNotFoundError
            mock_run.side_effect = [FileNotFoundError("samtools not found"), None]

            _collate_bam(input_bam, output_bam)

            # Verify fallback was called
            assert mock_fallback.call_count == 1
            mock_fallback.assert_called_once_with(input_bam, output_bam)

    @pytest.mark.unit
    def test_collate_temp_file_cleanup(self, tmp_path: Path) -> None:
        """Verify temp file is cleaned up after collate."""
        import tempfile
        from unittest.mock import MagicMock

        input_bam = str(tmp_path / "input.bam")
        output_bam = str(tmp_path / "output.bam")

        # Create a real temp file to track
        with tempfile.NamedTemporaryFile(
            prefix="collate_", suffix=".bam", delete=False
        ) as real_temp:
            temp_path = real_temp.name

        # Mock NamedTemporaryFile to return our tracked temp file
        mock_temp = MagicMock()
        mock_temp.name = temp_path
        mock_temp.__enter__ = MagicMock(return_value=mock_temp)
        mock_temp.__exit__ = MagicMock(return_value=None)

        with (
            patch(
                "plasmicheck.scripts.compare_alignments.tempfile.NamedTemporaryFile",
                return_value=mock_temp,
            ),
            patch("subprocess.run") as mock_run,
            patch("plasmicheck.scripts.compare_alignments._resort_supplementary"),
        ):
            mock_run.return_value = None

            _collate_bam(input_bam, output_bam)

            # Verify temp file was cleaned up
            assert not Path(temp_path).exists()


class TestResortSupplementary:
    @pytest.mark.unit
    def test_resort_supplementary_ordering(self, tmp_path: Path) -> None:
        """Verify supplementary reads are re-ordered: primary, supplementary, secondary."""
        try:
            import pysam
        except ImportError:
            pytest.skip("pysam not available")

        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"

        # Create a small BAM with reads in wrong order
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # Create reads with same query name but wrong order
            # Secondary first
            read1 = pysam.AlignedSegment()
            read1.query_name = "read1"
            read1.flag = 256  # secondary
            read1.reference_id = 0
            read1.reference_start = 100
            read1.mapping_quality = 30
            read1.cigar = [(0, 50)]  # 50M
            read1.query_sequence = "A" * 50
            read1.query_qualities = [30] * 50

            # Primary second
            read2 = pysam.AlignedSegment()
            read2.query_name = "read1"
            read2.flag = 0  # primary
            read2.reference_id = 0
            read2.reference_start = 200
            read2.mapping_quality = 60
            read2.cigar = [(0, 50)]
            read2.query_sequence = "A" * 50
            read2.query_qualities = [30] * 50

            # Supplementary third
            read3 = pysam.AlignedSegment()
            read3.query_name = "read1"
            read3.flag = 2048  # supplementary
            read3.reference_id = 0
            read3.reference_start = 300
            read3.mapping_quality = 40
            read3.cigar = [(0, 50)]
            read3.query_sequence = "A" * 50
            read3.query_qualities = [30] * 50

            # Write in wrong order: secondary, primary, supplementary
            outf.write(read1)
            outf.write(read2)
            outf.write(read3)

        # Re-sort supplementary
        _resort_supplementary(str(input_bam), str(output_bam))

        # Read output and verify order
        with pysam.AlignmentFile(str(output_bam), "rb") as inf:
            reads = list(inf.fetch(until_eof=True))
            assert len(reads) == 3
            # Order should be: primary (flag 0), supplementary (flag 2048), secondary (flag 256)
            assert reads[0].flag == 0  # primary first
            assert reads[1].flag == 2048  # supplementary second
            assert reads[2].flag == 256  # secondary third


class TestNameGroupBam:
    @pytest.mark.unit
    def test_name_group_bam_uses_collate_by_default(self, tmp_path: Path) -> None:
        """Verify _name_group_bam calls _collate_bam by default."""
        input_bam = str(tmp_path / "input.bam")
        output_bam = str(tmp_path / "output.bam")

        with (
            patch("plasmicheck.scripts.compare_alignments._collate_bam") as mock_collate,
            patch("plasmicheck.scripts.compare_alignments._namesort_bam_fallback") as mock_fallback,
        ):
            mock_collate.return_value = None

            _name_group_bam(input_bam, output_bam)

            # Verify collate was called
            assert mock_collate.call_count == 1
            mock_collate.assert_called_once_with(input_bam, output_bam)

            # Verify fallback was NOT called
            assert mock_fallback.call_count == 0
