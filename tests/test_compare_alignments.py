"""Unit tests for plasmicheck.scripts.compare_alignments."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from plasmicheck.scripts.compare_alignments import (
    InsertRegionData,
    _assign,
    _best_read,
    _streaming_compare,
    _write_assignment,
    calculate_alignment_score,
    parse_insert_region,
    read_overlaps_insert,
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
        # Should return InsertRegionData
        assert isinstance(result, InsertRegionData)
        assert result.start == 100
        assert result.end == 500
        # Tuple unpacking should still work
        start, end = result
        assert (start, end) == (100, 500)

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
        result = parse_insert_region(str(f))
        assert isinstance(result, InsertRegionData)
        # Tuple unpacking should still work
        s, e = result
        assert (s, e) == (start, end)

    @pytest.mark.unit
    def test_new_format_fields(self, tmp_path: Path) -> None:
        """Test parsing of new 11-line format with junction metadata."""
        f = tmp_path / "cDNA_positions.txt"
        f.write_text(
            "cDNA start position in plasmid: 100\n"
            "cDNA end position in plasmid: 500\n"
            "INSERT_REGION: (100, 500)\n"
            "SpliceJunctionCount: 2\n"
            "ExonCount: 3\n"
            "ExonLengths: 50,60,70\n"
            "IntronLengths: 100,150\n"
            "Strand: +\n"
            "Chromosome: chr1\n"
            "GenomicStart: 1000\n"
            "GenomicEnd: 1530\n"
        )
        result = parse_insert_region(str(f))
        assert isinstance(result, InsertRegionData)
        assert result.start == 100
        assert result.end == 500
        assert result.junction_count == 2
        assert result.exon_count == 3
        assert result.exon_lengths == [50, 60, 70]
        assert result.intron_lengths == [100, 150]
        assert result.strand == "+"
        assert result.chromosome == "chr1"
        assert result.genomic_start == 1000
        assert result.genomic_end == 1530

    @pytest.mark.unit
    def test_old_format_backward_compat(self, tmp_path: Path) -> None:
        """Test backward compatibility with old 3-line format."""
        f = tmp_path / "cDNA_positions.txt"
        f.write_text(
            "cDNA start position in plasmid: 200\n"
            "cDNA end position in plasmid: 800\n"
            "INSERT_REGION: (200, 800)\n"
        )
        result = parse_insert_region(str(f))
        assert isinstance(result, InsertRegionData)
        assert result.start == 200
        assert result.end == 800
        # New fields should be None (not present in old format)
        assert result.junction_count is None
        assert result.exon_count is None
        assert result.exon_lengths is None
        assert result.intron_lengths is None
        assert result.strand is None
        assert result.chromosome is None
        assert result.genomic_start is None
        assert result.genomic_end is None

    @pytest.mark.unit
    def test_na_values_raises_valueerror(self, tmp_path: Path) -> None:
        """Test that NA cDNA start/end raises ValueError (insert region unavailable)."""
        f = tmp_path / "cDNA_positions.txt"
        f.write_text(
            "cDNA start position in plasmid: NA\n"
            "cDNA end position in plasmid: NA\n"
            "INSERT_REGION: NA\n"
            "SpliceJunctionCount: NA\n"
            "ExonCount: NA\n"
            "ExonLengths: NA\n"
            "IntronLengths: NA\n"
            "Strand: NA\n"
            "Chromosome: NA\n"
            "GenomicStart: NA\n"
            "GenomicEnd: NA\n"
        )
        with pytest.raises(ValueError, match="Insert region unavailable"):
            parse_insert_region(str(f))


class TestInsertRegionDataIterator:
    """Test InsertRegionData.__iter__ for backward compatibility."""

    @pytest.mark.unit
    def test_iter_produces_start_end_tuple(self) -> None:
        """Test that __iter__ returns (start, end) for tuple unpacking."""
        data = InsertRegionData(start=100, end=500, junction_count=2)
        result = list(data)
        assert result == [100, 500]

    @pytest.mark.unit
    def test_tuple_unpacking(self) -> None:
        """Test that tuple unpacking works as expected."""
        data = InsertRegionData(start=100, end=500, junction_count=2)
        start, end = data
        assert start == 100
        assert end == 500


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


class TestReadOverlapsInsert:
    """Test the read_overlaps_insert() function with boundary cases."""

    @pytest.mark.unit
    def test_read_entirely_within_insert(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=200, reference_end=400)
        assert read_overlaps_insert(read, (100, 500)) is True

    @pytest.mark.unit
    def test_read_entirely_outside_before(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=0, reference_end=50)
        assert read_overlaps_insert(read, (100, 500)) is False

    @pytest.mark.unit
    def test_read_entirely_outside_after(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=600, reference_end=700)
        assert read_overlaps_insert(read, (100, 500)) is False

    @pytest.mark.unit
    def test_read_overlaps_insert_start(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=50, reference_end=150)
        assert read_overlaps_insert(read, (100, 500)) is True

    @pytest.mark.unit
    def test_read_overlaps_insert_end(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=450, reference_end=550)
        assert read_overlaps_insert(read, (100, 500)) is True

    @pytest.mark.unit
    def test_read_exactly_at_insert_start_boundary(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=100, reference_end=200)
        assert read_overlaps_insert(read, (100, 500)) is True

    @pytest.mark.unit
    def test_read_exactly_at_insert_end_boundary(self, mock_pysam_read: Any) -> None:
        # reference_end=501 is exclusive, so last base is 500 which IS in insert
        read = mock_pysam_read(reference_start=400, reference_end=501)
        assert read_overlaps_insert(read, (100, 500)) is True

    @pytest.mark.unit
    def test_read_ends_exactly_at_insert_start(self, mock_pysam_read: Any) -> None:
        # reference_end=100 is exclusive, so last base is 99, BEFORE insert
        read = mock_pysam_read(reference_start=0, reference_end=100)
        assert read_overlaps_insert(read, (100, 500)) is False

    @pytest.mark.unit
    def test_read_starts_exactly_after_insert_end(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=501, reference_end=600)
        assert read_overlaps_insert(read, (100, 500)) is False

    @pytest.mark.unit
    def test_unmapped_read(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(is_unmapped=True)
        assert read_overlaps_insert(read, (100, 500)) is False

    @pytest.mark.unit
    def test_none_reference_positions(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=None, reference_end=None)
        assert read_overlaps_insert(read, (100, 500)) is False


class TestAssignExtended:
    """Test the extended _assign() function with 5 categories."""

    @pytest.mark.unit
    def test_plasmid_wins_no_filtering(self) -> None:
        # Backward compat: no insert_region
        assert _assign(70, 30) == "Plasmid"

    @pytest.mark.unit
    def test_human_wins_no_filtering(self) -> None:
        assert _assign(30, 70) == "Human"

    @pytest.mark.unit
    def test_tied_no_filtering(self) -> None:
        assert _assign(50, 50) == "Tied"

    @pytest.mark.unit
    def test_plasmid_wins_overlaps_insert(self, mock_pysam_read: Any) -> None:
        read_in_insert = mock_pysam_read(reference_start=200, reference_end=400)
        assert _assign(70, 30, plasmid_read=read_in_insert, insert_region=(100, 500)) == "Plasmid"

    @pytest.mark.unit
    def test_plasmid_wins_outside_insert(self, mock_pysam_read: Any) -> None:
        read_outside = mock_pysam_read(reference_start=0, reference_end=50)
        assert (
            _assign(70, 0, plasmid_read=read_outside, insert_region=(100, 500)) == "Backbone_Only"
        )

    @pytest.mark.unit
    def test_plasmid_only_read_outside_insert(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=0, reference_end=50)
        assert _assign(55, 0, plasmid_read=read, insert_region=(100, 500)) == "Backbone_Only"

    @pytest.mark.unit
    def test_plasmid_only_read_overlaps_insert(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=200, reference_end=400)
        assert _assign(55, 0, plasmid_read=read, insert_region=(100, 500)) == "Plasmid"

    @pytest.mark.unit
    def test_score_margin_makes_ambiguous(self) -> None:
        # diff=2 < margin=5
        assert _assign(52, 50, score_margin=5) == "Ambiguous"

    @pytest.mark.unit
    def test_score_margin_not_triggered_when_diff_exceeds(self) -> None:
        # diff=40 >= margin=5
        assert _assign(70, 30, score_margin=5) == "Plasmid"

    @pytest.mark.unit
    def test_score_margin_zero_disabled(self) -> None:
        # margin=0 disabled, diff=1 enough
        assert _assign(51, 50, score_margin=0) == "Plasmid"

    @pytest.mark.unit
    def test_tied_not_ambiguous_with_margin(self) -> None:
        # Exact tie is always Tied, never Ambiguous
        assert _assign(50, 50, score_margin=10) == "Tied"

    @pytest.mark.unit
    def test_human_not_affected_by_insert_region(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(reference_start=200, reference_end=400)
        # Human wins regardless of insert region
        assert _assign(30, 70, plasmid_read=read, insert_region=(100, 500)) == "Human"


class TestStreamingCompareWithFiltering:
    """Test _streaming_compare() with insert_region parameter."""

    @pytest.mark.unit
    def test_plasmid_only_read_classified_backbone_only(self, mock_pysam_read: Any) -> None:
        # Plasmid-only read outside insert region
        p_read = mock_pysam_read(
            query_name="r1",
            mapping_quality=60,
            reference_start=0,
            reference_end=50,
        )
        p_read.is_secondary = False
        p_read.is_supplementary = False

        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("r1", [p_read])]), iter([])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf, insert_region=(100, 500))
            assert counts["Backbone_Only"] == 1
            assert counts["Plasmid"] == 0

    @pytest.mark.unit
    def test_plasmid_only_read_classified_plasmid_when_overlaps_insert(
        self, mock_pysam_read: Any
    ) -> None:
        # Plasmid-only read inside insert region
        p_read = mock_pysam_read(
            query_name="r1",
            mapping_quality=60,
            reference_start=200,
            reference_end=400,
        )
        p_read.is_secondary = False
        p_read.is_supplementary = False

        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("r1", [p_read])]), iter([])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf, insert_region=(100, 500))
            assert counts["Plasmid"] == 1
            assert counts["Backbone_Only"] == 0

    @pytest.mark.unit
    def test_no_insert_region_falls_back_to_old_behavior(self, mock_pysam_read: Any) -> None:
        # insert_region=None -> all plasmid-only reads are Plasmid
        p_read = mock_pysam_read(
            query_name="r1",
            mapping_quality=60,
            reference_start=0,
            reference_end=50,
        )
        p_read.is_secondary = False
        p_read.is_supplementary = False

        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("r1", [p_read])]), iter([])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf, insert_region=None)
            assert counts["Plasmid"] == 1
            assert counts["Backbone_Only"] == 0

    @pytest.mark.unit
    def test_score_margin_creates_ambiguous_in_stream(self, mock_pysam_read: Any) -> None:
        # Both reads with small score diff
        p_read = mock_pysam_read(
            query_name="r1",
            mapping_quality=52,
            cigartuples=[(0, 100)],
            nm_tag=0,
            mate_is_unmapped=False,
        )
        p_read.is_secondary = False
        p_read.is_supplementary = False

        h_read = mock_pysam_read(
            query_name="r1",
            mapping_quality=50,
            cigartuples=[(0, 100)],
            nm_tag=0,
            mate_is_unmapped=False,
        )
        h_read.is_secondary = False
        h_read.is_supplementary = False

        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("r1", [p_read])]), iter([("r1", [h_read])])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf, score_margin=5)
            # ps=62, hs=60, diff=2 < margin=5 -> Ambiguous
            assert counts["Ambiguous"] == 1


class TestBackwardCompatibility:
    """Test backward compatibility with existing code."""

    @pytest.mark.unit
    def test_existing_assign_api_unchanged(self) -> None:
        # Positional args still work
        assert _assign(70, 30) == "Plasmid"
        assert _assign(30, 70) == "Human"
        assert _assign(50, 50) == "Tied"

    @pytest.mark.unit
    def test_existing_streaming_compare_api_unchanged(self, mock_pysam_read: Any) -> None:
        # _streaming_compare without insert_region works identically
        p_read = mock_pysam_read(query_name="r1", mapping_quality=60)
        p_read.is_secondary = False
        p_read.is_supplementary = False

        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("r1", [p_read])]), iter([])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
            # Old behavior: plasmid-only -> Plasmid
            assert counts["Plasmid"] == 1
            # New categories initialized but zero
            assert counts["Backbone_Only"] == 0
            assert counts["Ambiguous"] == 0


class TestResistanceGeneIntegration:
    """Test resistance gene detection integration in compare_alignments."""

    @pytest.mark.unit
    def test_summary_tsv_includes_resistance_genes(self, tmp_path: Path) -> None:
        """Verify summary.tsv includes ResGene_ rows when GenBank file is provided."""
        from unittest.mock import MagicMock

        from plasmicheck.scripts.compare_alignments import compare_alignments
        from plasmicheck.scripts.resistance_gene_detection import DetectedGene

        # Mock BAM files
        plasmid_bam = tmp_path / "plasmid.bam"
        human_bam = tmp_path / "human.bam"
        output_base = tmp_path / "output"
        cdna_file = tmp_path / "cDNA_positions.txt"

        # Create cDNA positions file
        cdna_file.write_text(
            "cDNA start position in plasmid: 100\ncDNA end position in plasmid: 500\n"
        )

        # Mock GenBank file
        genbank_file = str(tmp_path / "plasmid.gb")

        # Create mock detected genes
        mock_genes = [
            DetectedGene(
                name="AmpR",
                sanitized_name="AmpR",
                start=1000,
                end=1800,
                region="backbone",
                qualifier_source="product",
            )
        ]

        # Mock coverage stats
        mock_stats = {
            "mean_depth": 42.5,
            "median_depth": 40.0,
            "breadth_1x": 0.95,
            "cv": 0.15,
        }

        with (
            patch(
                "plasmicheck.scripts.compare_alignments.detect_resistance_genes",
                return_value=mock_genes,
            ),
            patch("plasmicheck.scripts.compare_alignments.get_depth_array", return_value=[]),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_coverage_stats",
                return_value=mock_stats,
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics"
            ) as mock_cov,
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch("plasmicheck.scripts.compare_alignments._streaming_compare") as mock_stream,
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("pysam.AlignmentFile") as mock_bam,
        ):
            # Setup coverage metrics mock
            mock_cov.return_value = (
                {
                    "insert": {
                        "mean_depth": 10.0,
                        "median_depth": 9.0,
                        "breadth_1x": 0.8,
                        "breadth_5x": 0.7,
                        "cv": 0.2,
                    },
                    "backbone": {
                        "mean_depth": 5.0,
                        "median_depth": 4.0,
                        "breadth_1x": 0.6,
                        "breadth_5x": 0.5,
                        "cv": 0.3,
                    },
                },
                False,
            )

            # Setup streaming compare mock
            mock_stream.return_value = {
                "Plasmid": 100,
                "Human": 10,
                "Tied": 5,
                "Backbone_Only": 2,
                "Ambiguous": 1,
            }

            # Setup BAM mock
            mock_bam_instance = MagicMock()
            mock_bam_instance.get_reference_name.return_value = "plasmid"
            mock_bam_instance.lengths = [10000]  # Mock plasmid length
            mock_bam_instance.get_tid.return_value = 0
            mock_bam.return_value.__enter__.return_value = mock_bam_instance

            # Run compare_alignments
            compare_alignments(
                str(plasmid_bam),
                str(human_bam),
                str(output_base),
                genbank_file=genbank_file,
            )

            # Read summary.tsv
            summary_file = output_base.parent / (output_base.name + ".summary.tsv")
            assert summary_file.exists()
            content = summary_file.read_text()

            # Verify resistance gene rows
            assert "ResistanceGenesAvailable\ttrue" in content
            assert "ResistanceGenes\tAmpR" in content
            assert "ResGene_AmpR_MeanDepth\t42.50" in content
            assert "ResGene_AmpR_Breadth\t0.95" in content
            assert "ResGene_AmpR_Start\t1000" in content
            assert "ResGene_AmpR_End\t1800" in content
            assert "ResGene_AmpR_Region\tbackbone" in content

    @pytest.mark.unit
    def test_summary_tsv_no_genbank(self, tmp_path: Path) -> None:
        """Verify summary.tsv shows ResistanceGenesAvailable: false when no GenBank."""
        from unittest.mock import MagicMock

        from plasmicheck.scripts.compare_alignments import compare_alignments

        plasmid_bam = tmp_path / "plasmid.bam"
        human_bam = tmp_path / "human.bam"
        output_base = tmp_path / "output"
        cdna_file = tmp_path / "cDNA_positions.txt"

        cdna_file.write_text(
            "cDNA start position in plasmid: 100\ncDNA end position in plasmid: 500\n"
        )

        with (
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics"
            ) as mock_cov,
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch("plasmicheck.scripts.compare_alignments._streaming_compare") as mock_stream,
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("pysam.AlignmentFile") as mock_bam,
        ):
            mock_cov.return_value = (
                {
                    "insert": {
                        "mean_depth": 10.0,
                        "median_depth": 9.0,
                        "breadth_1x": 0.8,
                        "breadth_5x": 0.7,
                        "cv": 0.2,
                    },
                    "backbone": {
                        "mean_depth": 5.0,
                        "median_depth": 4.0,
                        "breadth_1x": 0.6,
                        "breadth_5x": 0.5,
                        "cv": 0.3,
                    },
                },
                False,
            )
            mock_stream.return_value = {
                "Plasmid": 100,
                "Human": 10,
                "Tied": 5,
                "Backbone_Only": 2,
                "Ambiguous": 1,
            }

            mock_bam_instance = MagicMock()
            mock_bam_instance.lengths = [10000]
            mock_bam_instance.get_tid.return_value = 0
            mock_bam.return_value.__enter__.return_value = mock_bam_instance

            compare_alignments(
                str(plasmid_bam), str(human_bam), str(output_base), genbank_file=None
            )

            summary_file = output_base.parent / (output_base.name + ".summary.tsv")
            assert summary_file.exists()
            content = summary_file.read_text()

            assert "ResistanceGenesAvailable\tfalse" in content
            assert "ResGene_" not in content

    @pytest.mark.unit
    def test_summary_tsv_genbank_no_genes(self, tmp_path: Path) -> None:
        """Verify ResistanceGenes: None when GenBank provided but no genes detected."""
        from unittest.mock import MagicMock

        from plasmicheck.scripts.compare_alignments import compare_alignments

        plasmid_bam = tmp_path / "plasmid.bam"
        human_bam = tmp_path / "human.bam"
        output_base = tmp_path / "output"
        cdna_file = tmp_path / "cDNA_positions.txt"
        genbank_file = str(tmp_path / "plasmid.gb")

        cdna_file.write_text(
            "cDNA start position in plasmid: 100\ncDNA end position in plasmid: 500\n"
        )

        with (
            patch(
                "plasmicheck.scripts.compare_alignments.detect_resistance_genes", return_value=[]
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics"
            ) as mock_cov,
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch("plasmicheck.scripts.compare_alignments._streaming_compare") as mock_stream,
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("pysam.AlignmentFile") as mock_bam,
        ):
            mock_cov.return_value = (
                {
                    "insert": {
                        "mean_depth": 10.0,
                        "median_depth": 9.0,
                        "breadth_1x": 0.8,
                        "breadth_5x": 0.7,
                        "cv": 0.2,
                    },
                    "backbone": {
                        "mean_depth": 5.0,
                        "median_depth": 4.0,
                        "breadth_1x": 0.6,
                        "breadth_5x": 0.5,
                        "cv": 0.3,
                    },
                },
                False,
            )
            mock_stream.return_value = {
                "Plasmid": 100,
                "Human": 10,
                "Tied": 5,
                "Backbone_Only": 2,
                "Ambiguous": 1,
            }

            mock_bam_instance = MagicMock()
            mock_bam_instance.lengths = [10000]
            mock_bam_instance.get_tid.return_value = 0
            mock_bam.return_value.__enter__.return_value = mock_bam_instance

            compare_alignments(
                str(plasmid_bam),
                str(human_bam),
                str(output_base),
                genbank_file=genbank_file,
            )

            summary_file = output_base.parent / (output_base.name + ".summary.tsv")
            assert summary_file.exists()
            content = summary_file.read_text()

            assert "ResistanceGenesAvailable\ttrue" in content
            assert "ResistanceGenes\tNone" in content
            assert "ResGene_" not in content


class TestJunctionMetadataInSummary:
    """Test junction metadata fields in summary.tsv (Phase 15)."""

    @pytest.mark.unit
    def test_summary_includes_junction_fields(self, tmp_path: Path) -> None:
        """Verify summary.tsv includes junction metadata from InsertRegionData."""
        from unittest.mock import MagicMock

        from plasmicheck.scripts.compare_alignments import compare_alignments

        plasmid_bam = tmp_path / "plasmid.bam"
        human_bam = tmp_path / "human.bam"
        output_base = tmp_path / "output"
        cdna_file = tmp_path / "cDNA_positions.txt"

        # Create new-format cDNA positions file with junction metadata
        cdna_file.write_text(
            "cDNA start position in plasmid: 100\n"
            "cDNA end position in plasmid: 500\n"
            "INSERT_REGION: (100, 500)\n"
            "SpliceJunctionCount: 2\n"
            "ExonCount: 3\n"
            "ExonLengths: 50,60,70\n"
            "IntronLengths: 100,150\n"
            "Strand: +\n"
            "Chromosome: chr1\n"
            "GenomicStart: 1000\n"
            "GenomicEnd: 1530\n"
        )

        with (
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics"
            ) as mock_cov,
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch("plasmicheck.scripts.compare_alignments._streaming_compare") as mock_stream,
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("pysam.AlignmentFile") as mock_bam,
        ):
            mock_cov.return_value = (
                {
                    "insert": {
                        "mean_depth": 10.0,
                        "median_depth": 9.0,
                        "breadth_1x": 0.8,
                        "breadth_5x": 0.7,
                        "cv": 0.2,
                    },
                    "backbone": {
                        "mean_depth": 5.0,
                        "median_depth": 4.0,
                        "breadth_1x": 0.6,
                        "breadth_5x": 0.5,
                        "cv": 0.3,
                    },
                },
                False,
            )
            mock_stream.return_value = {
                "Plasmid": 100,
                "Human": 10,
                "Tied": 5,
                "Backbone_Only": 2,
                "Ambiguous": 1,
            }

            mock_bam_instance = MagicMock()
            mock_bam_instance.lengths = [10000]
            mock_bam_instance.get_tid.return_value = 0
            mock_bam.return_value.__enter__.return_value = mock_bam_instance

            compare_alignments(
                str(plasmid_bam),
                str(human_bam),
                str(output_base),
                genbank_file=None,
                cdna_positions_file=str(cdna_file),
            )

            summary_file = output_base.parent / (output_base.name + ".summary.tsv")
            assert summary_file.exists()
            content = summary_file.read_text()

            # Verify junction metadata rows
            assert "SpliceJunctionCount\t2" in content
            assert "ExonCount\t3" in content
            assert "ExonLengths\t50,60,70" in content
            assert "IntronLengths\t100,150" in content

    @pytest.mark.unit
    def test_summary_junction_fields_na_when_old_format(self, tmp_path: Path) -> None:
        """Verify junction fields are NA when using old 3-line cDNA_positions.txt."""
        from unittest.mock import MagicMock

        from plasmicheck.scripts.compare_alignments import compare_alignments

        plasmid_bam = tmp_path / "plasmid.bam"
        human_bam = tmp_path / "human.bam"
        output_base = tmp_path / "output"
        cdna_file = tmp_path / "cDNA_positions.txt"

        # Create old-format cDNA positions file (3 lines only)
        cdna_file.write_text(
            "cDNA start position in plasmid: 100\n"
            "cDNA end position in plasmid: 500\n"
            "INSERT_REGION: (100, 500)\n"
        )

        with (
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics"
            ) as mock_cov,
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch("plasmicheck.scripts.compare_alignments._streaming_compare") as mock_stream,
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("pysam.AlignmentFile") as mock_bam,
        ):
            mock_cov.return_value = (
                {
                    "insert": {
                        "mean_depth": 10.0,
                        "median_depth": 9.0,
                        "breadth_1x": 0.8,
                        "breadth_5x": 0.7,
                        "cv": 0.2,
                    },
                    "backbone": {
                        "mean_depth": 5.0,
                        "median_depth": 4.0,
                        "breadth_1x": 0.6,
                        "breadth_5x": 0.5,
                        "cv": 0.3,
                    },
                },
                False,
            )
            mock_stream.return_value = {
                "Plasmid": 100,
                "Human": 10,
                "Tied": 5,
                "Backbone_Only": 2,
                "Ambiguous": 1,
            }

            mock_bam_instance = MagicMock()
            mock_bam_instance.lengths = [10000]
            mock_bam_instance.get_tid.return_value = 0
            mock_bam.return_value.__enter__.return_value = mock_bam_instance

            compare_alignments(
                str(plasmid_bam),
                str(human_bam),
                str(output_base),
                genbank_file=None,
                cdna_positions_file=str(cdna_file),
            )

            summary_file = output_base.parent / (output_base.name + ".summary.tsv")
            assert summary_file.exists()
            content = summary_file.read_text()

            # Verify junction metadata rows are NA (old format)
            assert "SpliceJunctionCount\tNA" in content
            assert "ExonCount\tNA" in content
            assert "ExonLengths\tNA" in content
            assert "IntronLengths\tNA" in content
