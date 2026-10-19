"""Unit tests for plasmicheck.scripts.compare_alignments."""

from __future__ import annotations

import io
import shutil
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pysam
import pytest

from plasmicheck.scripts.compare_alignments import (
    DEFAULT_THRESHOLD,
    InsertRegionData,
    _assign,
    _best_read,
    _has_origin_sa_support,
    _is_concordant_bridge_pair,
    _is_pair,
    _is_strictly_inside_insert,
    _mate_key,
    _namesort_bam,
    _read_id,
    _scan_plasmid_bam,
    _scores,
    _streaming_compare,
    _write_assignment,
    calculate_alignment_score,
    compare_alignments,
    parse_insert_region,
    read_overlaps_insert,
)


@pytest.fixture
def mock_pysam_read() -> Any:
    """Factory fixture for creating mock pysam aligned reads with is_paired support."""

    def _make_read(
        *,
        is_unmapped: bool = False,
        mapping_quality: int = 60,
        cigartuples: list[tuple[int, int]] | None = None,
        mate_is_unmapped: bool = False,
        is_paired: bool = True,
        is_read1: bool = True,
        is_read2: bool = False,
        is_secondary: bool = False,
        is_supplementary: bool = False,
        is_reverse: bool = False,
        reference_id: int = 0,
        next_reference_id: int = 0,
        reference_start: int | None = 0,
        reference_end: int | None = 100,
        next_reference_start: int | None = None,
        nm_tag: int | None = 0,
        query_name: str = "test_read",
        cigarstring: str | None = "100M",
        query_length: int = 100,
        reference_name: str = "plasmid",
        tags: dict[str, Any] | None = None,
    ) -> MagicMock:
        read = MagicMock()
        read.is_unmapped = is_unmapped
        read.mapping_quality = mapping_quality
        read.cigartuples = cigartuples if cigartuples is not None else [(0, 100)]
        read.mate_is_unmapped = mate_is_unmapped
        read.is_paired = is_paired
        read.is_read1 = is_read1 if is_paired else False
        read.is_read2 = is_read2 if is_paired else False
        read.is_secondary = is_secondary
        read.is_supplementary = is_supplementary
        read.is_reverse = is_reverse
        read.reference_id = reference_id
        read.next_reference_id = next_reference_id
        read.reference_start = reference_start
        read.reference_end = reference_end
        read.next_reference_start = (
            next_reference_start if next_reference_start is not None else reference_start
        )
        read.query_name = query_name
        read.cigarstring = cigarstring
        read.query_length = query_length
        read.reference_name = reference_name

        tag_dict: dict[str, Any] = {}
        if nm_tag is not None:
            tag_dict["NM"] = nm_tag
        if tags is not None:
            tag_dict.update(tags)

        read.has_tag.side_effect = lambda t: t in tag_dict
        read.get_tag.side_effect = lambda t: tag_dict.get(t, 0)

        def _get_overlap(start: int, end: int) -> int:
            if reference_start is None or reference_end is None:
                return 0
            s = max(reference_start, start)
            e = min(reference_end, end)
            return max(0, e - s)

        read.get_overlap.side_effect = _get_overlap

        return read

    return _make_read


class TestNameSortedMerge:
    @pytest.mark.unit
    def test_sort_order_matches_string_comparison(self) -> None:
        with patch("plasmicheck.scripts.compare_alignments.subprocess.run") as run:
            _namesort_bam("input.bam", "sorted.bam")

        command = run.call_args.args[0]
        assert command[:3] == ["samtools", "sort", "-N"]

    @pytest.mark.integration
    @pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not found on PATH")
    @pytest.mark.parametrize(
        ("plasmid_names", "human_names", "expected"),
        [
            (
                ["read2", "read10"],
                ["read10"],
                [("read10", "Tied"), ("read2", "Plasmid")],
            ),
            (
                ["read10"],
                ["read2", "read10"],
                [("read10", "Tied"), ("read2", "Human")],
            ),
            (
                ["read1", "read01", "read001", "read2", "read10"],
                ["read001", "read1", "read10"],
                [
                    ("read001", "Tied"),
                    ("read01", "Plasmid"),
                    ("read1", "Tied"),
                    ("read10", "Tied"),
                    ("read2", "Plasmid"),
                ],
            ),
        ],
    )
    def test_each_name_is_compared_once(
        self,
        tmp_path: Path,
        plasmid_names: list[str],
        human_names: list[str],
        expected: list[tuple[str, str]],
    ) -> None:
        for label, names in (("plasmid", plasmid_names), ("human", human_names)):
            source = tmp_path / f"{label}.bam"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "ref", "LN": 1000}]}
            with pysam.AlignmentFile(str(source), "wb", header=header) as bam:
                for name in names:
                    read = pysam.AlignedSegment()
                    read.query_name = name
                    read.query_sequence = "A" * 100
                    read.reference_id = 0
                    read.reference_start = 0
                    read.mapping_quality = 60
                    read.cigarstring = "100M"
                    bam.write(read)
            _namesort_bam(str(source), str(tmp_path / f"{label}.sorted.bam"))

        output = io.StringIO()
        counts = _streaming_compare(
            str(tmp_path / "plasmid.sorted.bam"), str(tmp_path / "human.sorted.bam"), output
        )

        assignments = [tuple(line.split("\t")[:2]) for line in output.getvalue().splitlines()]
        assert assignments == expected
        assert sum(counts.values()) == len(expected)


class TestCalculateAlignmentScore:
    @pytest.mark.unit
    def test_unmapped_read_returns_zero(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(is_unmapped=True)
        assert calculate_alignment_score(read) == 0

    @pytest.mark.unit
    def test_simple_mapped_read(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, no clipping, 0 mismatches -> 60
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(0, 100)],  # 100M
            nm_tag=0,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 60

    @pytest.mark.unit
    def test_with_soft_clipping(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, 10bp soft clip, 0 mismatches -> 60 - 10 = 50
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(4, 10), (0, 90)],  # 10S90M
            nm_tag=0,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 50

    @pytest.mark.unit
    def test_with_mismatches(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, no clipping, NM=5 -> 60 - 5 = 55
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=5,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 55

    @pytest.mark.unit
    def test_combined_penalties(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, 10bp soft clip, NM=5 -> 60 - 10 - 5 = 45
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(4, 10), (0, 90)],
            nm_tag=5,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 45

    @pytest.mark.unit
    def test_mate_unmapped_no_bonus(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, no clipping, 0 mismatches -> 60
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=0,
            mate_is_unmapped=True,
        )
        assert calculate_alignment_score(read) == 60

    @pytest.mark.unit
    def test_hard_clipping(self, mock_pysam_read: Any) -> None:
        # MAPQ=60, 20bp hard clip, 0 mismatches -> 60 - 20 = 40
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(5, 20), (0, 80)],  # 20H80M
            nm_tag=0,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 40

    @pytest.mark.unit
    def test_no_nm_tag(self, mock_pysam_read: Any) -> None:
        # No NM tag -> 0 mismatches assumed
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=None,  # no NM tag
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 60

    @pytest.mark.unit
    def test_unpaired_read_no_mate_bonus(self, mock_pysam_read: Any) -> None:
        read = mock_pysam_read(
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=0,
            is_paired=False,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(read) == 60


class TestScores:
    """Test _scores() fragment-level mate bonus and pair concordance."""

    @pytest.mark.unit
    def test_concordant_pair_gets_mate_bonus(self, mock_pysam_read: Any) -> None:
        r1 = mock_pysam_read(
            mapping_quality=60,
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            reference_start=100,
            reference_end=200,
            next_reference_start=300,
        )
        r2 = mock_pysam_read(
            mapping_quality=60,
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            reference_start=300,
            reference_end=400,
            next_reference_start=100,
        )
        # 60 + 10 = 70 for each mate
        assert _scores([r1, r2]) == [70, 70]

    @pytest.mark.unit
    def test_distant_mates_no_mate_bonus(self, mock_pysam_read: Any) -> None:
        # Span > 1000 bp across intron on human
        r1 = mock_pysam_read(
            mapping_quality=60,
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            reference_start=100,
            reference_end=200,
            next_reference_start=5000,
        )
        r2 = mock_pysam_read(
            mapping_quality=60,
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            reference_start=5000,
            reference_end=5100,
            next_reference_start=100,
        )
        assert _scores([r1, r2], length=0) == [60, 60]

    @pytest.mark.unit
    def test_circular_origin_crossing_pair_gets_mate_bonus(self, mock_pysam_read: Any) -> None:
        # Plasmid length 12332, forward mate near end, reverse mate near start
        r1 = mock_pysam_read(
            mapping_quality=60,
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            reference_start=12182,
            reference_end=12283,
            next_reference_start=50,
        )
        r2 = mock_pysam_read(
            mapping_quality=60,
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            reference_start=50,
            reference_end=151,
            next_reference_start=12182,
        )
        # Circular span: (151 + 12332) - 12182 = 301 <= 1000 -> mate bonus awarded!
        assert _scores([r1, r2], length=12332, span_origin=True) == [70, 70]


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
    @pytest.mark.parametrize("extra_side", ["plasmid", "human"])
    def test_numeric_names_in_lexicographic_order(
        self, mock_pysam_read: Any, extra_side: str
    ) -> None:
        shared = mock_pysam_read(query_name="read10", mapping_quality=60)
        extra = mock_pysam_read(query_name="read2", mapping_quality=60)
        plasmid = [("read10", [shared])]
        human = [("read10", [shared])]
        (plasmid if extra_side == "plasmid" else human).append(("read2", [extra]))

        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter(plasmid), iter(human)],
        ):
            output = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", output)

        assignments = [tuple(line.split("\t")[:2]) for line in output.getvalue().splitlines()]
        assert assignments == [("read10", "Tied"), ("read2", extra_side.title())]
        assert sum(counts.values()) == 2

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

    @pytest.mark.unit
    def test_paired_read_group_scores_summed(self, mock_pysam_read: Any) -> None:
        """Paired read group: both primary mates' scores are summed."""
        # Plasmid mates: 60 + 10 mate bonus = 70 each -> sum = 140
        p_r1 = mock_pysam_read(
            query_name="pair1",
            mapping_quality=60,
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            reference_start=0,
            reference_end=100,
            next_reference_start=200,
            mate_is_unmapped=False,
        )
        p_r2 = mock_pysam_read(
            query_name="pair1",
            mapping_quality=60,
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            reference_start=200,
            reference_end=300,
            next_reference_start=0,
            mate_is_unmapped=False,
        )
        # Human mates: 30 + 10 mate bonus = 40 each -> sum = 80
        h_r1 = mock_pysam_read(
            query_name="pair1",
            mapping_quality=30,
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            reference_start=0,
            reference_end=100,
            next_reference_start=200,
            mate_is_unmapped=False,
        )
        h_r2 = mock_pysam_read(
            query_name="pair1",
            mapping_quality=30,
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            reference_start=200,
            reference_end=300,
            next_reference_start=0,
            mate_is_unmapped=False,
        )

        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("pair1", [p_r1, p_r2])]), iter([("pair1", [h_r1, h_r2])])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
            assert counts["Plasmid"] == 1
            fields = buf.getvalue().splitlines()[0].split("\t")
            assert fields[0] == "pair1"
            assert fields[1] == "Plasmid"
            assert fields[2] == "140"
            assert fields[3] == "80"

    @pytest.mark.unit
    def test_plasmid_only_negative_score_tied(self, mock_pysam_read: Any) -> None:
        """One-sided plasmid read with negative score becomes Tied, raw score in TSV."""
        p_read = mock_pysam_read(
            query_name="r_neg",
            mapping_quality=0,
            cigartuples=[(4, 50), (0, 50)],
            is_paired=False,
            mate_is_unmapped=True,
        )

        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("r_neg", [p_read])]), iter([])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
            assert counts["Tied"] == 1
            assert counts["Plasmid"] == 0
            fields = buf.getvalue().splitlines()[0].split("\t")
            assert fields[1] == "Tied"
            assert fields[2] == "-50"
            assert fields[3] == "0"

    @pytest.mark.unit
    def test_human_only_negative_score_tied(self, mock_pysam_read: Any) -> None:
        """One-sided human read with negative score becomes Tied, raw score in TSV."""
        h_read = mock_pysam_read(
            query_name="r_neg_h",
            mapping_quality=0,
            cigartuples=[(4, 63), (0, 37)],
            is_paired=False,
            mate_is_unmapped=True,
        )

        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([]), iter([("r_neg_h", [h_read])])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
            assert counts["Tied"] == 1
            assert counts["Human"] == 0
            fields = buf.getvalue().splitlines()[0].split("\t")
            assert fields[1] == "Tied"
            assert fields[2] == "0"
            assert fields[3] == "-63"


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

    @pytest.mark.unit
    def test_origin_wrapping_insert(self, mock_pysam_read: Any) -> None:
        """Origin-wrapping insert region (start > end), e.g. (800, 100)."""
        # Read in first half [850, 950]
        r1 = mock_pysam_read(reference_start=850, reference_end=950)
        assert read_overlaps_insert(r1, (800, 100)) is True

        # Read in second half [10, 50]
        r2 = mock_pysam_read(reference_start=10, reference_end=50)
        assert read_overlaps_insert(r2, (800, 100)) is True

        # Read overlapping start boundary [750, 850]
        r3 = mock_pysam_read(reference_start=750, reference_end=850)
        assert read_overlaps_insert(r3, (800, 100)) is True

        # Read strictly in backbone [200, 300]
        r_backbone = mock_pysam_read(reference_start=200, reference_end=300)
        assert read_overlaps_insert(r_backbone, (800, 100)) is False


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
        # Both reads with small score diff due to alignment penalties
        p_read = mock_pysam_read(
            query_name="r1",
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=1,
            mate_is_unmapped=False,
        )
        p_read.is_secondary = False
        p_read.is_supplementary = False

        h_read = mock_pysam_read(
            query_name="r1",
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=3,
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
            # ps=59, hs=57, diff=2 < margin=5 -> Ambiguous
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


class TestVerdictAndRatio:
    """Test ratio calculation edge cases (0/0) and verdict determination (#100)."""

    @pytest.mark.unit
    def test_zero_zero_counts_produce_zero_ratio_not_contaminated(self, tmp_path: Path) -> None:
        """0 plasmid reads and 0 human reads (0/0) produces ratio 0.0 and not contaminated."""
        output_base = tmp_path / "zero_counts"
        dummy_metrics = {
            "insert": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
            "backbone": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
        }

        with (
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch(
                "plasmicheck.scripts.compare_alignments._streaming_compare",
                return_value={
                    "Plasmid": 0,
                    "Human": 0,
                    "Tied": 0,
                    "Backbone_Only": 0,
                    "Ambiguous": 0,
                },
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics",
                return_value=(dummy_metrics, False),
            ),
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch("pysam.AlignmentFile"),
        ):
            compare_alignments("p.bam", "h.bam", str(output_base))

        content = (tmp_path / "zero_counts.summary.tsv").read_text()
        assert "Verdict\tSample is not contaminated with plasmid DNA" in content
        assert "Ratio\t0.0" in content

    @pytest.mark.unit
    def test_backbone_only_with_zero_human_produces_zero_ratio_when_filtered(
        self, tmp_path: Path
    ) -> None:
        """When filter_backbone_only=True, Backbone_Only reads with 0 human produce ratio 0.0."""
        output_base = tmp_path / "backbone_counts"
        dummy_metrics = {
            "insert": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
            "backbone": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
        }

        with (
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch(
                "plasmicheck.scripts.compare_alignments._streaming_compare",
                return_value={
                    "Plasmid": 0,
                    "Human": 0,
                    "Tied": 0,
                    "Backbone_Only": 15,
                    "Ambiguous": 0,
                },
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics",
                return_value=(dummy_metrics, False),
            ),
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch("pysam.AlignmentFile"),
        ):
            compare_alignments("p.bam", "h.bam", str(output_base), filter_backbone_only=True)

        content = (tmp_path / "backbone_counts.summary.tsv").read_text()
        assert "Verdict\tSample is not contaminated with plasmid DNA" in content
        assert "Ratio\t0.0" in content

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "plasmid_reads,expected_ratio", [(70, 0.7), (80, 0.8), (90, 0.9), (95, 0.95), (100, 1.0)]
    )
    def test_ratio_point_nine_default_threshold_unclear(
        self, tmp_path: Path, plasmid_reads: int, expected_ratio: float
    ) -> None:
        """Ratios in [0.6, 1.0] with default threshold (0.8) fall in UNCLEAR_RANGE -> unclear."""
        output_base = tmp_path / f"unclear_counts_{plasmid_reads}"
        dummy_metrics = {
            "insert": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
            "backbone": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
        }

        with (
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch(
                "plasmicheck.scripts.compare_alignments._streaming_compare",
                return_value={
                    "Plasmid": plasmid_reads,
                    "Human": 100,
                    "Tied": 0,
                    "Backbone_Only": 0,
                    "Ambiguous": 0,
                },
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics",
                return_value=(dummy_metrics, False),
            ),
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch("pysam.AlignmentFile"),
        ):
            compare_alignments("p.bam", "h.bam", str(output_base), threshold=DEFAULT_THRESHOLD)

        content = (tmp_path / f"unclear_counts_{plasmid_reads}.summary.tsv").read_text()
        assert "Verdict\tSample contamination status is unclear" in content
        assert f"Ratio\t{expected_ratio}" in content

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "threshold,plasmid_reads,expected_verdict",
        [
            (2.0, 150, "Sample contamination status is unclear"),
            (2.0, 210, "Sample is contaminated with plasmid DNA"),
            (0.4, 40, "Sample is not contaminated with plasmid DNA"),
            (0.4, 50, "Sample is contaminated with plasmid DNA"),
        ],
    )
    def test_custom_threshold_verdict_boundaries(
        self, tmp_path: Path, threshold: float, plasmid_reads: int, expected_verdict: str
    ) -> None:
        """Custom threshold extends unclear upper bound or acts as strict boundary below lower_bound."""
        output_base = tmp_path / f"custom_thresh_{threshold}_{plasmid_reads}"
        dummy_metrics = {
            "insert": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
            "backbone": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
        }
        with (
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch(
                "plasmicheck.scripts.compare_alignments._streaming_compare",
                return_value={
                    "Plasmid": plasmid_reads,
                    "Human": 100,
                    "Tied": 0,
                    "Backbone_Only": 0,
                    "Ambiguous": 0,
                },
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics",
                return_value=(dummy_metrics, False),
            ),
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch("pysam.AlignmentFile"),
        ):
            compare_alignments("p.bam", "h.bam", str(output_base), threshold=threshold)

        content = (tmp_path / f"custom_thresh_{threshold}_{plasmid_reads}.summary.tsv").read_text()
        assert f"Verdict\t{expected_verdict}" in content


class TestAccuracyDiscoveriesV038:
    """Tests for v0.38.0 accuracy improvements (boundary off-by-one, circular origin, MAPQ tie, splice junctions)."""

    @pytest.mark.unit
    def test_coverage_outside_insert_no_false_positive_at_end_boundary(
        self, mock_pysam_read: Any
    ) -> None:
        """Read ending at insert_end + 1 (0-based exclusive) must have 0.0 outside-insert coverage."""
        from plasmicheck.scripts.compare_alignments import calculate_coverage_outside_insert

        # insert_region is (100, 500) inclusive -> [100, 501) half-open
        read = mock_pysam_read(reference_start=400, reference_end=501, cigartuples=[(0, 101)])
        mock_bam = MagicMock()
        mock_bam.lengths = [1000]
        mock_bam.fetch.return_value = [read]
        mock_bam.__enter__.return_value = mock_bam
        mock_bam.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bam):
            cov = calculate_coverage_outside_insert("p.bam", (100, 500))
        assert cov == 0.0

    @pytest.mark.unit
    def test_circular_origin_boundary_clip_excused(self, mock_pysam_read: Any) -> None:
        """Soft-clip at circular origin is excused only when supported by valid SA alignment."""
        # Unsupported clip without SA tag is penalized: 24 - 35 = -11
        r_unsupported = mock_pysam_read(
            mapping_quality=24,
            cigartuples=[(4, 35), (0, 66)],
            reference_start=0,
            reference_end=66,
            nm_tag=0,
            is_paired=True,
            mate_is_unmapped=False,
        )
        assert calculate_alignment_score(r_unsupported, 0) == -11
        assert calculate_alignment_score(r_unsupported, 12332) == -11

        # Valid SA tag at opposite end of plasmid excuses the origin clip
        r_supported = mock_pysam_read(
            mapping_quality=24,
            cigartuples=[(4, 35), (0, 66)],
            reference_start=0,
            reference_end=66,
            nm_tag=0,
            is_paired=True,
            mate_is_unmapped=False,
            tags={"SA": "plasmid,12297,+,35M66S,60,0;"},
        )
        assert calculate_alignment_score(r_supported, 12332) == 24

    @pytest.mark.unit
    def test_mapq_asymmetry_identical_cigar_nm_ties(self, mock_pysam_read: Any) -> None:
        """Identical CIGAR and NM with MAPQ=60 on plasmid vs 48..59 on human ties instead of false Plasmid win."""
        p_read = mock_pysam_read(
            query_name="r_sym", mapping_quality=60, cigartuples=[(0, 101)], nm_tag=0
        )
        h_read = mock_pysam_read(
            query_name="r_sym", mapping_quality=53, cigartuples=[(0, 101)], nm_tag=0
        )
        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("r_sym", [p_read])]), iter([("r_sym", [h_read])])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
        assert counts["Tied"] == 1
        assert counts["Plasmid"] == 0

    @pytest.mark.unit
    def test_splice_junction_reads_counted_on_both_strands(self, mock_pysam_read: Any) -> None:
        """_scan_plasmid_bam counts unclipped reads spanning cDNA exon-exon splice junctions."""
        from plasmicheck.scripts.compare_alignments import _scan_plasmid_bam

        # Minus-strand insert starting at 100 with exons [200, 100] -> reversed to [100, 200] -> junction at 200
        ir = InsertRegionData(
            start=100,
            end=399,
            junction_count=1,
            exon_count=2,
            exon_lengths=[200, 100],
            intron_lengths=[500],
            strand="-",
        )
        spanning_read = mock_pysam_read(
            reference_start=160,
            reference_end=260,
            cigartuples=[(0, 100)],
        )
        clipped_read = mock_pysam_read(
            reference_start=160,
            reference_end=260,
            cigartuples=[(4, 20), (0, 80)],
        )
        mock_bam = MagicMock()
        mock_bam.lengths = [5000]
        mock_bam.fetch.return_value = [spanning_read, clipped_read]
        mock_bam.__enter__.return_value = mock_bam
        mock_bam.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bam):
            cov, _mm, sj_reads, sj_cov, *_rest = _scan_plasmid_bam("p.bam", ir)
        assert cov == 0.0
        assert sj_reads == 1
        assert sj_cov == 1

    @pytest.mark.unit
    def test_mapq_tie_breaker_reverse_strand_cigar_inversion(self, mock_pysam_read: Any) -> None:
        """MAPQ tie-breaker reverses CIGAR when read orientations differ (pr.is_reverse != hr.is_reverse)."""
        p_read = mock_pysam_read(
            query_name="r_inv",
            mapping_quality=60,
            cigartuples=[(4, 10), (0, 90)],
            nm_tag=0,
            is_reverse=False,
            is_paired=False,
        )
        # Human read is reverse strand with reversed CIGAR -> should tie
        h_read_rev = mock_pysam_read(
            query_name="r_inv",
            mapping_quality=50,
            cigartuples=[(0, 90), (4, 10)],
            nm_tag=0,
            is_reverse=True,
            is_paired=False,
        )
        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("r_inv", [p_read])]), iter([("r_inv", [h_read_rev])])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
        assert counts["Tied"] == 1
        assert counts["Plasmid"] == 0

        # When reverse-strand CIGAR is NOT inverted, tie-breaker does not fire
        h_read_non_inv = mock_pysam_read(
            query_name="r_inv",
            mapping_quality=50,
            cigartuples=[(4, 10), (0, 90)],
            nm_tag=0,
            is_reverse=True,
            is_paired=False,
        )
        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("r_inv", [p_read])]), iter([("r_inv", [h_read_non_inv])])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
        # Plasmid score 50 > human score 40
        assert counts["Plasmid"] == 1
        assert counts["Tied"] == 0

    @pytest.mark.unit
    def test_mapq_tie_breaker_mate_bonus_asymmetry_preserves_plasmid_win(
        self, mock_pysam_read: Any
    ) -> None:
        """MAPQ tie-breaker does not trigger when score differences arise from mate bonuses."""
        pr1 = mock_pysam_read(
            query_name="pair_mate",
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            reference_start=100,
            reference_end=200,
            next_reference_start=300,
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=0,
        )
        pr2 = mock_pysam_read(
            query_name="pair_mate",
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            reference_start=300,
            reference_end=400,
            next_reference_start=100,
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=0,
        )
        # Human reads have identical CIGAR, NM, and MAPQ, but span > 1000 (discordant -> no mate bonus)
        hr1 = mock_pysam_read(
            query_name="pair_mate",
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            reference_start=100,
            reference_end=200,
            next_reference_start=5000,
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=0,
        )
        hr2 = mock_pysam_read(
            query_name="pair_mate",
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            reference_start=5000,
            reference_end=5100,
            next_reference_start=100,
            mapping_quality=60,
            cigartuples=[(0, 100)],
            nm_tag=0,
        )
        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[iter([("pair_mate", [pr1, pr2])]), iter([("pair_mate", [hr1, hr2])])],
        ):
            buf = io.StringIO()
            counts = _streaming_compare("p.bam", "h.bam", buf)
        # Plasmid pair received mate bonus; tie-breaker guard preserves Plasmid win
        assert counts["Plasmid"] == 1
        assert counts["Tied"] == 0

    @pytest.mark.unit
    def test_has_origin_sa_support_validation(self, mock_pysam_read: Any) -> None:
        """_has_origin_sa_support validates strand, contig, NM <= 2, and opposite-end position."""
        r_start_clip = mock_pysam_read(
            reference_name="pContig",
            reference_start=0,
            reference_end=66,
            is_reverse=False,
            tags={"SA": "pContig,12297,+,35M66S,60,0;"},
        )
        assert _has_origin_sa_support(r_start_clip, 12332, is_start_clip=True) is True

        # Wrong strand
        r_wrong_strand = mock_pysam_read(
            reference_name="pContig",
            reference_start=0,
            reference_end=66,
            is_reverse=False,
            tags={"SA": "pContig,12297,-,35M66S,60,0;"},
        )
        assert _has_origin_sa_support(r_wrong_strand, 12332, is_start_clip=True) is False

        # Wrong contig
        r_wrong_contig = mock_pysam_read(
            reference_name="pContig",
            reference_start=0,
            reference_end=66,
            is_reverse=False,
            tags={"SA": "otherContig,12297,+,35M66S,60,0;"},
        )
        assert _has_origin_sa_support(r_wrong_contig, 12332, is_start_clip=True) is False

        # NM > 2
        r_high_nm = mock_pysam_read(
            reference_name="pContig",
            reference_start=0,
            reference_end=66,
            is_reverse=False,
            tags={"SA": "pContig,12297,+,35M66S,60,3;"},
        )
        assert _has_origin_sa_support(r_high_nm, 12332, is_start_clip=True) is False

        # Position not near opposite end
        r_wrong_pos = mock_pysam_read(
            reference_name="pContig",
            reference_start=0,
            reference_end=66,
            is_reverse=False,
            tags={"SA": "pContig,5000,+,35M66S,60,0;"},
        )
        assert _has_origin_sa_support(r_wrong_pos, 12332, is_start_clip=True) is False

        # End clip with SA near start
        r_end_clip = mock_pysam_read(
            reference_name="pContig",
            reference_start=12266,
            reference_end=12332,
            is_reverse=False,
            tags={"SA": "pContig,50,+,35M66S,60,0;"},
        )
        assert _has_origin_sa_support(r_end_clip, 12332, is_start_clip=False) is True

    @pytest.mark.unit
    def test_vector_insert_junction_reads_in_scan(self, mock_pysam_read: Any) -> None:
        """_scan_plasmid_bam counts reads overlapping vector-insert boundaries with >=10bp anchors."""
        ir = (100, 500)
        # Read 1: crosses boundary 100 with anchor >= 10 and NM=0
        r1 = mock_pysam_read(
            reference_start=85,
            reference_end=125,
            cigartuples=[(0, 40)],
            nm_tag=0,
        )
        # Read 2: crosses boundary 501 with anchor >= 10 and NM=1
        r2 = mock_pysam_read(
            reference_start=490,
            reference_end=530,
            cigartuples=[(0, 40)],
            nm_tag=1,
        )
        # Read 3: crosses boundary 100 but NM=3 (> 2) -> not counted
        r3 = mock_pysam_read(
            reference_start=85,
            reference_end=125,
            cigartuples=[(0, 40)],
            nm_tag=3,
        )
        # Read 4: strictly inside insert -> not crossing boundary
        r4 = mock_pysam_read(
            reference_start=200,
            reference_end=300,
            cigartuples=[(0, 100)],
            nm_tag=0,
        )
        mock_bam = MagicMock()
        mock_bam.lengths = [1000]
        mock_bam.fetch.return_value = [r1, r2, r3, r4]
        mock_bam.__enter__.return_value = mock_bam
        mock_bam.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bam):
            *_, vi_reads, _mean_len = _scan_plasmid_bam("p.bam", ir)
        assert vi_reads == 2

    @pytest.mark.unit
    def test_wrapped_insert_scan_inside_calculation(self, mock_pysam_read: Any) -> None:
        """_scan_plasmid_bam calculates inside/outside correctly for origin-wrapping insert."""
        # Insert wraps around origin: [800, 1000) + [0, 201). Backbone is [201, 800) = 599 bp.
        ir = (800, 200)
        # Read inside insert [850, 950] -> outside = 0, inside = 100
        r1 = mock_pysam_read(reference_start=850, reference_end=950, cigartuples=[(0, 100)])
        # Read straddling boundary [750, 850] -> outside = [750, 800) = 50, inside = 50
        r2 = mock_pysam_read(reference_start=750, reference_end=850, cigartuples=[(0, 100)])
        # Read strictly in backbone [300, 400] -> outside = 100, inside = 0
        r3 = mock_pysam_read(reference_start=300, reference_end=400, cigartuples=[(0, 100)])

        mock_bam = MagicMock()
        mock_bam.lengths = [1000]
        mock_bam.fetch.return_value = [r1, r2, r3]
        mock_bam.__enter__.return_value = mock_bam
        mock_bam.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bam):
            cov, *_rest = _scan_plasmid_bam("p.bam", ir)
        # Total outside bases = 0 + 50 + 100 = 150. Reference backbone bases = 599.
        assert cov == pytest.approx(150 / 599)

    @pytest.mark.unit
    def test_concordant_bridge_pairs_and_streaming(self, mock_pysam_read: Any) -> None:
        """Concordant bridge pairs (one mate in insert, one in backbone, span <= 1000) are tallied."""
        ir = (100, 500)
        # Inward facing pair: mate 1 strictly in insert (200..300), mate 2 strictly in backbone (600..700)
        pr1 = mock_pysam_read(
            query_name="bridge1",
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            reference_start=200,
            reference_end=300,
            mapping_quality=60,
            cigartuples=[(0, 100)],
        )
        pr2 = mock_pysam_read(
            query_name="bridge1",
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            reference_start=600,
            reference_end=700,
            mapping_quality=60,
            cigartuples=[(0, 100)],
        )
        assert _mate_key(pr1) == 1
        assert _mate_key(pr2) == 2
        assert _is_pair([pr1, pr2]) is True
        assert _is_strictly_inside_insert(pr1, ir, 1000) is True
        assert _is_strictly_inside_insert(pr2, ir, 1000) is False
        assert _is_concordant_bridge_pair([pr1, pr2], ir, 1000) is True

        # Non-bridge pair: both mates strictly inside insert
        pr_in1 = mock_pysam_read(
            query_name="in_pair",
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            reference_start=150,
            reference_end=250,
        )
        pr_in2 = mock_pysam_read(
            query_name="in_pair",
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            reference_start=300,
            reference_end=400,
        )
        assert _is_concordant_bridge_pair([pr_in1, pr_in2], ir, 1000) is False

        # Streaming compare matched branch counts bridge pair
        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[
                iter([("bridge1", [pr1, pr2])]),
                iter([]),
            ],
        ):
            buf = io.StringIO()
            counts_p_only = _streaming_compare(
                "p.bam", "h.bam", buf, insert_region=ir, plasmid_length=1000
            )
        assert counts_p_only["Plasmid"] == 1
        assert counts_p_only.bridge_pairs == 1

        # Streaming compare matched branch where plasmid wins over human
        hr1 = mock_pysam_read(
            query_name="bridge1",
            is_paired=True,
            is_read1=True,
            is_read2=False,
            is_reverse=False,
            mapping_quality=20,
        )
        hr2 = mock_pysam_read(
            query_name="bridge1",
            is_paired=True,
            is_read1=False,
            is_read2=True,
            is_reverse=True,
            mapping_quality=20,
        )
        with patch(
            "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
            side_effect=[
                iter([("bridge1", [pr1, pr2])]),
                iter([("bridge1", [hr1, hr2])]),
            ],
        ):
            buf = io.StringIO()
            counts_matched = _streaming_compare(
                "p.bam", "h.bam", buf, insert_region=ir, plasmid_length=1000
            )
        assert counts_matched["Plasmid"] == 1
        assert counts_matched.bridge_pairs == 1

    @pytest.mark.unit
    def test_mismatches_near_insert_nearest_boundary_circular(self, mock_pysam_read: Any) -> None:
        """_scan_plasmid_bam nearest boundary distance calculates circular shortest distance."""
        ir = (100, 950)  # Boundaries: 100 and 951 on a 1000bp plasmid
        # Read at [5, 45] is 54 bp from boundary 951 across plasmid origin (1000 - 951 + 5 = 54 <= 100)
        read_circular_near = mock_pysam_read(
            reference_start=5,
            reference_end=45,
            cigartuples=[(0, 40)],
            nm_tag=1,
        )
        # Read at [500, 600] is ~400 bp away from both boundaries
        read_far = mock_pysam_read(
            reference_start=500,
            reference_end=600,
            cigartuples=[(0, 100)],
            nm_tag=1,
        )
        mock_bam = MagicMock()
        mock_bam.lengths = [1000]
        mock_bam.fetch.return_value = [read_circular_near, read_far]
        mock_bam.__enter__.return_value = mock_bam
        mock_bam.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bam):
            _, mm, *_rest = _scan_plasmid_bam("p.bam", ir, window_size=100)
        assert mm["with_mismatches_or_clipping"] == 1
        assert mm["without_mismatches_or_clipping"] == 0

    @pytest.mark.unit
    def test_splice_junction_events_and_exon_intron_boundaries(self, mock_pysam_read: Any) -> None:
        """_scan_plasmid_bam counts multiple junction events per read and exon-intron boundary clips."""
        # 3 exons -> 2 junctions at 200 and 300
        ir = InsertRegionData(
            start=100,
            end=399,
            junction_count=2,
            exon_count=3,
            exon_lengths=[100, 100, 100],
            strand="+",
        )
        # Read 1: Spans both junctions (150 to 350)
        r_multi = mock_pysam_read(
            reference_start=150,
            reference_end=350,
            cigartuples=[(0, 200)],
            nm_tag=0,
            query_length=200,
        )
        # Read 2: Leading clip of 20 bp starting at junction 200
        r_lead_clip = mock_pysam_read(
            reference_start=201,
            reference_end=300,
            cigartuples=[(4, 20), (0, 99)],
            nm_tag=0,
            query_length=119,
        )
        # Read 3: Trailing clip of 20 bp ending at junction 300
        r_trail_clip = mock_pysam_read(
            reference_start=100,
            reference_end=299,
            cigartuples=[(0, 199), (4, 20)],
            nm_tag=0,
            query_length=219,
        )
        mock_bam = MagicMock()
        mock_bam.lengths = [5000]
        mock_bam.fetch.return_value = [r_multi, r_lead_clip, r_trail_clip]
        mock_bam.__enter__.return_value = mock_bam
        mock_bam.__exit__.return_value = None

        with patch("pysam.AlignmentFile", return_value=mock_bam):
            (
                _cov,
                _mm,
                sj_reads,
                sj_covered,
                sj_events,
                ei_boundary_reads,
                _vi_reads,
                mean_read_len,
            ) = _scan_plasmid_bam("p.bam", ir)

        assert sj_reads == 1
        assert sj_events == 2
        assert sj_covered == 2
        assert ei_boundary_reads == 2
        assert mean_read_len == pytest.approx((200 + 119 + 219) / 3)

    @pytest.mark.unit
    def test_estimated_contamination_fraction_in_summary(self, tmp_path: Path) -> None:
        """compare_alignments calculates EstimatedContaminationFraction correctly and handles fallbacks."""
        output_base = tmp_path / "est_test"
        cdna_file = tmp_path / "cDNA_positions.txt"
        cdna_file.write_text(
            "cDNA start position in plasmid: 100\n"
            "cDNA end position in plasmid: 500\n"
            "INSERT_REGION: (100, 500)\n"
            "SpliceJunctionCount: 1\n"
            "ExonCount: 2\n"
            "ExonLengths: 200,201\n"
            "IntronLengths: 100\n"
            "Strand: +\n"
        )
        dummy_metrics = {
            "insert": {
                "mean_depth": 50.0,
                "median_depth": 50.0,
                "breadth_1x": 1.0,
                "breadth_5x": 1.0,
                "cv": 0.1,
            },
            "backbone": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
        }

        empty_mm: dict[str, int] = {}
        # Case 1: Valid calculation
        # J_count = 1, events = 5, mean_read_length = 150.0, d_insert = 50.0
        # W = 150 - 30 + 1 = 121
        # est = (5 * 150) / (1 * 121 * 50) = 750 / 6050 = 0.1240
        mock_scan_results: tuple[float, dict[str, int], int, int, int, int, int, float] = (
            0.0,
            empty_mm,
            5,
            1,
            5,
            0,
            0,
            150.0,
        )
        with (
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch(
                "plasmicheck.scripts.compare_alignments._streaming_compare",
                return_value={
                    "Plasmid": 10,
                    "Human": 10,
                    "Tied": 0,
                    "Backbone_Only": 0,
                    "Ambiguous": 0,
                },
            ),
            patch(
                "plasmicheck.scripts.compare_alignments._scan_plasmid_bam",
                return_value=mock_scan_results,
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics",
                return_value=(dummy_metrics, False),
            ),
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch("pysam.AlignmentFile"),
        ):
            compare_alignments(
                "p.bam", "h.bam", str(output_base), cdna_positions_file=str(cdna_file)
            )

        content = (tmp_path / "est_test.summary.tsv").read_text()
        assert "EstimatedContaminationFraction\t0.1240" in content
        assert "SpliceJunctionEvents\t5" in content

        # Case 2: Zero events -> 0.0000
        mock_zero_events: tuple[float, dict[str, int], int, int, int, int, int, float] = (
            0.0,
            empty_mm,
            0,
            0,
            0,
            0,
            0,
            150.0,
        )
        with (
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch(
                "plasmicheck.scripts.compare_alignments._streaming_compare",
                return_value={
                    "Plasmid": 10,
                    "Human": 10,
                    "Tied": 0,
                    "Backbone_Only": 0,
                    "Ambiguous": 0,
                },
            ),
            patch(
                "plasmicheck.scripts.compare_alignments._scan_plasmid_bam",
                return_value=mock_zero_events,
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics",
                return_value=(dummy_metrics, False),
            ),
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch("pysam.AlignmentFile"),
        ):
            compare_alignments(
                "p.bam", "h.bam", str(output_base), cdna_positions_file=str(cdna_file)
            )

        content = (tmp_path / "est_test.summary.tsv").read_text()
        assert "EstimatedContaminationFraction\t0.0000" in content

        # Case 3: Fallback when coverage_fallback=True -> NA
        with (
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch(
                "plasmicheck.scripts.compare_alignments._streaming_compare",
                return_value={
                    "Plasmid": 10,
                    "Human": 10,
                    "Tied": 0,
                    "Backbone_Only": 0,
                    "Ambiguous": 0,
                },
            ),
            patch(
                "plasmicheck.scripts.compare_alignments._scan_plasmid_bam",
                return_value=mock_scan_results,
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics",
                return_value=(dummy_metrics, True),
            ),
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch("pysam.AlignmentFile"),
        ):
            compare_alignments(
                "p.bam", "h.bam", str(output_base), cdna_positions_file=str(cdna_file)
            )

        content = (tmp_path / "est_test.summary.tsv").read_text()
        assert "EstimatedContaminationFraction\tNA" in content

    @pytest.mark.unit
    def test_temp_bams_naming_and_cleanup_on_error(self, tmp_path: Path) -> None:
        """Temp BAMs are created in output dir with PID and cleaned up in finally on error."""
        import os

        sub_dir = tmp_path / "sub_dir"
        sub_dir.mkdir()
        output_base = sub_dir / "sample"

        pid = os.getpid()
        plasmid_temp = sub_dir / f"sample.plasmid.namesorted.{pid}.bam"
        human_temp = sub_dir / f"sample.human.namesorted.{pid}.bam"

        def _fake_namesort(in_bam: str, out_bam: str, threads: int = 1) -> None:
            Path(out_bam).write_text("fake bam")

        with (
            pytest.raises(RuntimeError, match="Simulated merge failure"),
            patch(
                "plasmicheck.scripts.compare_alignments._namesort_bam", side_effect=_fake_namesort
            ),
            patch(
                "plasmicheck.scripts.compare_alignments._streaming_compare",
                side_effect=RuntimeError("Simulated merge failure"),
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics",
                return_value=({"insert": {}, "backbone": {}}, False),
            ),
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch("pysam.AlignmentFile"),
        ):
            compare_alignments("p.bam", "h.bam", str(output_base))

        # Both temp BAMs must be removed by finally block
        assert not plasmid_temp.exists()
        assert not human_temp.exists()


@pytest.mark.unit
def test_read_id_hashing() -> None:
    """Test keyed hash helper for anonymization."""
    assert _read_id("read1", None) == "read1"
    key1 = b"k" * 16
    key2 = b"j" * 16
    hashed1 = _read_id("read1", key1)
    hashed2 = _read_id("read1", key2)
    assert hashed1 == _read_id("read1", key1)
    assert hashed1 != hashed2
    assert len(hashed1) == 16  # 8 bytes hex


@pytest.mark.unit
def test_streaming_compare_anonymize(mock_pysam_read: Any) -> None:
    """Test _streaming_compare hashes ReadID when anonymize_key is provided."""
    p_read = mock_pysam_read(query_name="read1", mapping_quality=60)
    p_read.is_secondary = False
    p_read.is_supplementary = False
    h_read = mock_pysam_read(query_name="read1", mapping_quality=30)
    h_read.is_secondary = False
    h_read.is_supplementary = False

    with patch(
        "plasmicheck.scripts.compare_alignments._iter_reads_by_name",
        side_effect=[
            iter([("read1", [p_read])]),
            iter([("read1", [h_read])]),
            iter([("read1", [p_read])]),
            iter([("read1", [h_read])]),
        ],
    ):
        buf_plain = io.StringIO()
        plain = _streaming_compare("p.bam", "h.bam", buf_plain)
        buf_anon = io.StringIO()
        anon = _streaming_compare("p.bam", "h.bam", buf_anon, anonymize_key=b"k" * 16)

    assert plain == anon
    assert "read1" in buf_plain.getvalue()
    assert "read1" not in buf_anon.getvalue()
    hashed = _read_id("read1", b"k" * 16)
    assert hashed in buf_anon.getvalue()


class TestDepthNormalisedPlasmidExcess:
    """Tests for Issue #126 depth-normalised plasmid excess fraction helpers."""

    @pytest.mark.unit
    def test_compute_f_excess_clean_and_excess(self) -> None:
        from plasmicheck.scripts.compare_alignments import (
            compute_f_excess,
            compute_f_junction,
            compute_f_plasmid,
        )

        # Clean sample where dp_locus / dp_sample_median <= r_cohort -> f_excess == 0.0
        assert compute_f_excess(55.0, 100.0, r_cohort=0.60) == 0.0
        assert compute_f_excess(60.0, 100.0, r_cohort=0.60) == 0.0

        # Invalid/zero/non-finite inputs -> 0.0
        assert compute_f_excess(0.0, 100.0, r_cohort=0.60) == 0.0
        assert compute_f_excess(100.0, 0.0, r_cohort=0.60) == 0.0
        assert compute_f_excess(100.0, 100.0, r_cohort=0.0) == 0.0
        assert compute_f_excess(float("nan"), 100.0, r_cohort=0.60) == 0.0
        assert compute_f_excess(float("inf"), 100.0, r_cohort=0.60) == 0.0

        # Single-exon construct (0 splice junctions) with 2x locus depth excess (120 vs expected 60)
        f_ex = compute_f_excess(120.0, 100.0, r_cohort=0.60)
        assert f_ex == pytest.approx(0.5)
        f_junc = compute_f_junction(
            0.0, 120.0, plasmid_ratio=0.1, is_contaminated=False, is_unclear=False
        )
        assert f_junc == 0.0
        f_pl = compute_f_plasmid(f_ex, f_junc)
        assert f_pl > 0.0
        assert f_pl == pytest.approx(f_ex)

    @pytest.mark.unit
    def test_compute_f_junction_and_f_plasmid(self) -> None:
        from plasmicheck.scripts.compare_alignments import compute_f_junction, compute_f_plasmid

        # Junction reads driven
        fj = compute_f_junction(
            35.5, 101.0, plasmid_ratio=0.5, is_contaminated=False, read_len=101, min_overlap=15
        )
        # scale = 101 / 71; 35.5 * (101/71) / 101 = 0.5
        assert fj == pytest.approx(0.5)

        # Ratio fallback when contaminated vs unclear
        fj_contam = compute_f_junction(
            0.0, 100.0, plasmid_ratio=3.0, is_contaminated=True, is_unclear=False
        )
        assert fj_contam == pytest.approx(3.0 / 4.0)

        fj_unclear = compute_f_junction(
            0.0, 100.0, plasmid_ratio=3.0, is_contaminated=False, is_unclear=True
        )
        assert fj_unclear == pytest.approx(0.5 * (3.0 / 4.0))

        assert compute_f_plasmid(0.2, 0.75) == pytest.approx(0.75)
        assert compute_f_plasmid(1.5, 0.3) == pytest.approx(1.0)

    @pytest.mark.unit
    def test_compute_bed_median_depth(self, tmp_path: Path) -> None:
        import numpy as np

        from plasmicheck.scripts.compare_alignments import compute_bed_median_depth

        bed_path = tmp_path / "targets.bed"
        bed_path.write_text("# header comment\nchr11\t100\t200\nchr16\t300\t400\nchr3\t500\t600\n")
        bam_path = tmp_path / "sample.bam"
        bam_path.write_text("dummy")

        mock_af = MagicMock()
        mock_af.references = ("11", "16", "3")
        # Return depths 80, 100, 140 across the 3 100-bp intervals -> median 100.0
        mock_af.count_coverage.side_effect = [
            (np.full(100, 80), np.zeros(100), np.zeros(100), np.zeros(100)),
            (np.full(100, 100), np.zeros(100), np.zeros(100), np.zeros(100)),
            (np.full(100, 140), np.zeros(100), np.zeros(100), np.zeros(100)),
        ]
        mock_af.__enter__.return_value = mock_af
        mock_af.__exit__.return_value = None

        with (
            patch("pysam.AlignmentFile", return_value=mock_af),
            patch(
                "plasmicheck.scripts.compare_alignments.subprocess.run",
                side_effect=OSError("force pysam fallback"),
            ),
        ):
            med = compute_bed_median_depth(str(bam_path), str(bed_path))
        assert med == pytest.approx(100.0)
        for call in mock_af.count_coverage.call_args_list:
            assert call.kwargs.get("quality_threshold") == 0

    @pytest.mark.unit
    def test_compute_bed_median_depth_samtools_bedcov_and_chr_mismatch(
        self, tmp_path: Path
    ) -> None:
        import shutil
        import subprocess

        from plasmicheck.scripts.compare_alignments import compute_bed_median_depth

        bed_path = tmp_path / "targets.bed"
        bed_path.write_text(
            "# header comment\nchr11\t100\t200\nchr16\t300\t500\nchr3\t500\t600\nchrUnk\t0\t100\n"
        )
        bam_path = tmp_path / "sample.bam"
        bam_path.write_text("dummy")

        mock_af = MagicMock()
        mock_af.references = ("11", "16", "3")
        mock_af.__enter__.return_value = mock_af
        mock_af.__exit__.return_value = None

        captured_bed_lines: list[str] = []

        def fake_run(cmd: list[str], **_kwargs: Any) -> subprocess.CompletedProcess[str]:
            assert cmd[:3] == ["samtools", "bedcov", "-j"]
            captured_bed_lines.extend(Path(cmd[3]).read_text().splitlines())
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout="11\t100\t200\t8000\n16\t300\t500\t20000\n3\t500\t600\t14000\n",
                stderr="",
            )

        with (
            patch("pysam.AlignmentFile", return_value=mock_af),
            patch("plasmicheck.scripts.compare_alignments.subprocess.run", side_effect=fake_run),
        ):
            med = compute_bed_median_depth(str(bam_path), str(bed_path))

        assert med == pytest.approx(100.0)
        assert captured_bed_lines == ["11\t100\t200", "16\t300\t500", "3\t500\t600"]
        mock_af.count_coverage.assert_not_called()

        # Verify identical results between samtools bedcov and pysam fallback on real covered loci
        if shutil.which("samtools"):
            real_bed = tmp_path / "real_targets.bed"
            # Mix chr-prefixed, unprefixed, and absent contig names across covered target loci
            real_bed.write_text(
                "chr11\t128911000\t128912000\n16\t1210000\t1211000\nchr3\t184357000\t184358000\nchr99\t0\t1000\n"
            )
            for bam_rel in (
                "tests/data/real/contaminated.renamed.bam",
                "tests/data/real/not_contaminated.renamed.bam",
            ):
                if Path(bam_rel).exists():
                    st_val = compute_bed_median_depth(bam_rel, str(real_bed))
                    with patch(
                        "plasmicheck.scripts.compare_alignments.subprocess.run",
                        side_effect=OSError("force pysam fallback"),
                    ):
                        py_val = compute_bed_median_depth(bam_rel, str(real_bed))
                    assert st_val == pytest.approx(py_val)
                    assert st_val > 200.0


class TestFlagVcfVariantsEngine:
    """Tests for Issue #127 VCF / TSV variant flagging engine."""

    @pytest.fixture
    def apa_summary_rows(self) -> list[dict[str, Any]]:
        return [
            # APA18-T: contaminated with KCNJ5 G151R and L168R constructs
            {
                "sample": "APA18-T",
                "plasmid": "pcDNA5_FRT_TO_KCNJ5_G151R",
                "chrom": "chr11",
                "genomic_start": 128910000,
                "genomic_end": 128913000,
                "introns": [],
                "construct_variants": ["128911724:G>A"],
                "verdict": "Sample is contaminated with plasmid DNA",
                "ratio": 12.0,
                "spliced_reads": 0,
                "f_plasmid": 0.60,
                "f_excess": 0.60,
            },
            {
                "sample": "APA18-T",
                "plasmid": "pcDNA5_FRT_TO_KCNJ5_L168R",
                "chrom": "chr11",
                "genomic_start": 128910000,
                "genomic_end": 128913000,
                "introns": [],
                "construct_variants": ["128911776:T>G"],
                "verdict": "Sample is contaminated with plasmid DNA",
                "ratio": 12.0,
                "spliced_reads": 0,
                "f_plasmid": 0.60,
                "f_excess": 0.60,
            },
            # APA19-T: contaminated with KCNJ5 G151R, L168R, CLCN2 p.Tyr17His MNV, and CACNA1H
            {
                "sample": "APA19-T",
                "plasmid": "pcDNA5_FRT_TO_KCNJ5_G151R",
                "chrom": "chr11",
                "genomic_start": 128910000,
                "genomic_end": 128913000,
                "introns": [],
                "construct_variants": ["128911724:G>A"],
                "verdict": "Sample is contaminated with plasmid DNA",
                "ratio": 8.0,
                "spliced_reads": 0,
                "f_plasmid": 0.50,
                "f_excess": 0.50,
            },
            {
                "sample": "APA19-T",
                "plasmid": "pcDNA5_FRT_TO_KCNJ5_L168R",
                "chrom": "chr11",
                "genomic_start": 128910000,
                "genomic_end": 128913000,
                "introns": [],
                "construct_variants": ["128911776:T>G"],
                "verdict": "Sample is contaminated with plasmid DNA",
                "ratio": 8.0,
                "spliced_reads": 0,
                "f_plasmid": 0.50,
                "f_excess": 0.50,
            },
            {
                "sample": "APA19-T",
                "plasmid": "pcDNA5_FRT_TO_CLCN2_Y17H",
                "chrom": "chr3",
                "genomic_start": 184350000,
                "genomic_end": 184365000,
                "introns": [(184353147, 184353240), (184361300, 184361400)],
                "construct_variants": ["184361431:AC>GT"],
                "verdict": "Sample is contaminated with plasmid DNA",
                "ratio": 15.0,
                "spliced_reads": 42,
                "f_plasmid": 0.70,
                "f_excess": 0.65,
            },
            # APA41-T: contaminated with KCNJ5 L168R and CACNA1H spliced construct
            {
                "sample": "APA41-T",
                "plasmid": "pcDNA5_FRT_TO_KCNJ5_L168R",
                "chrom": "chr11",
                "genomic_start": 128910000,
                "genomic_end": 128913000,
                "introns": [],
                "construct_variants": ["128911776:T>G"],
                "verdict": "Sample is contaminated with plasmid DNA",
                "ratio": 5.0,
                "spliced_reads": 0,
                "f_plasmid": 0.40,
                "f_excess": 0.40,
            },
            {
                "sample": "APA41-T",
                "plasmid": "pcDNA5_FRT_TO_CACNA1H",
                "chrom": "chr16",
                "genomic_start": 1200000,
                "genomic_end": 1270000,
                "introns": [(1207119, 1207220), (1210654, 1210750), (1265000, 1265500)],
                "construct_variants": ["1208076:C>G"],
                "verdict": "Sample is contaminated with plasmid DNA",
                "ratio": 25.0,
                "spliced_reads": 80,
                "f_plasmid": 0.80,
                "f_excess": 0.75,
            },
            # APA23-N (contaminated normal) vs APA23-TS (clean tumour)
            {
                "sample": "APA23-N",
                "plasmid": "pcDNA5_FRT_TO_KCNJ5_L168R",
                "chrom": "chr11",
                "genomic_start": 128910000,
                "genomic_end": 128913000,
                "introns": [],
                "construct_variants": ["128911776:T>G"],
                "verdict": "Sample is contaminated with plasmid DNA",
                "ratio": 53.5,
                "spliced_reads": 0,
                "f_plasmid": 0.85,
                "f_excess": 0.80,
            },
            {
                "sample": "APA23-TS",
                "plasmid": "pcDNA5_FRT_TO_KCNJ5_L168R",
                "chrom": "chr11",
                "genomic_start": 128910000,
                "genomic_end": 128913000,
                "introns": [],
                "construct_variants": ["128911776:T>G"],
                "verdict": "Sample is not contaminated with plasmid DNA",
                "ratio": 0.0,
                "spliced_reads": 0,
                "f_plasmid": 0.0,
                "f_excess": 0.0,
            },
        ]

    @pytest.mark.unit
    def test_flag_vcf_variants_tsv_fixture(
        self, tmp_path: Path, apa_summary_rows: list[dict[str, Any]]
    ) -> None:
        from plasmicheck.scripts.compare_alignments import flag_vcf_variants

        tsv_in = tmp_path / "variants.tsv"
        tsv_out = tmp_path / "flagged_variants.tsv"
        rows = [
            "sample\tchr\tpos\tref\talt\tdp\tad\tvaf\tref_ctx_in\talt_ctx_in\tgene",
            # APA18_To: KCNJ5 G151R (AD=49) + L168R (AD=15) -> construct_allele, cooc=yes, ratio=15/49=0.306122
            "APA18_To\tchr11\t128911724\tG\tA\t200\t49\t0.245\t-\t-\tKCNJ5",
            "APA18_To\tchr11\t128911776\tT\tG\t200\t15\t0.075\t-\t-\tKCNJ5",
            # APA18_TvsN: KCNJ5 G151R + L168R
            "APA18_TvsN\tchr11\t128911724\tG\tA\t200\t49\t0.245\t-\t-\tKCNJ5",
            "APA18_TvsN\tchr11\t128911776\tT\tG\t200\t15\t0.075\t-\t-\tKCNJ5",
            # APA19_To: KCNJ5 G151R + L168R and CLCN2 MNV (chr3:184361431 AC>GT)
            "APA19_To\tchr11\t128911724\tG\tA\t180\t46\t0.255\t-\t-\tKCNJ5",
            "APA19_To\tchr11\t128911776\tT\tG\t180\t136\t0.755\t-\t-\tKCNJ5",
            "APA19_To\tchr3\t184361431\tAC\tGT\t150\t60\t0.400\t-\t-\tCLCN2",
            # APA41_To: KCNJ5 L168R construct_allele
            "APA41_To\tchr11\t128911776\tT\tG\t160\t32\t0.200\t-\t-\tKCNJ5",
            # APA23-N (contaminated normal -> construct_allele) vs APA23_To (clean tumour -> none!)
            "APA23-N\tchr11\t128911776\tT\tG\t140\t70\t0.500\t-\t-\tKCNJ5",
            "APA23_To\tchr11\t128911776\tT\tG\t150\t35\t0.233\t-\t-\tKCNJ5",
            # CACNA1H 15-bp truncated TSV ref deletion with 4-bp left-norm shift (intron 1207119-1207220, pos=1207114 -> del_start=1207115)
            "APA41_To\tchr16\t1207114\tACGTACGTACGTACG\tA\t120\t40\t0.333\t-\t-\tCACNA1H",
            # CACNA1H full intron deletion with 4-bp left-norm shift (intron 1210654-1210750, len=97, pos=1210649)
            f"APA41_To\tchr16\t1210649\t{'A' * 98}\tA\t120\t45\t0.375\t-\t-\tCACNA1H",
            # CACNA1H splice-boundary insertion via alt_ctx_in != '-'
            "APA41_To\tchr16\t1265005\tC\tCAGGT\t110\t30\t0.272\t-\tAGGT\tCACNA1H",
            # CLCN2 splice-boundary insertion within +-5 bp of intron boundary (184353147-1 = 184353146)
            "APA19_To\tchr3\t184353143\tG\tGACCT\t130\t35\t0.269\t-\t-\tCLCN2",
        ]
        tsv_in.write_text("\n".join(rows) + "\n")

        flagged = flag_vcf_variants(str(tsv_in), str(tsv_out), apa_summary_rows)
        assert len(flagged) == 14
        assert tsv_out.exists()

        # APA18_To G151R & L168R: construct_allele, cooc=yes, ratio=0.306122
        assert flagged[0]["plasmid_artefact"] == "construct_allele"
        assert flagged[0]["construct_allele_cooccurrence"] == "yes"
        assert flagged[0]["construct_allele_ad_ratio"] == "0.306122"
        assert float(flagged[0]["plasmid_explainable_vaf"]) == pytest.approx(0.30)
        assert flagged[1]["plasmid_artefact"] == "construct_allele"
        assert flagged[1]["construct_allele_cooccurrence"] == "yes"
        assert flagged[1]["construct_allele_ad_ratio"] == "0.306122"

        # APA18_TvsN G151R & L168R
        assert flagged[2]["plasmid_artefact"] == "construct_allele"
        assert flagged[3]["plasmid_artefact"] == "construct_allele"

        # APA19_To G151R, L168R, and CLCN2 MNV (AC>GT)
        assert flagged[4]["plasmid_artefact"] == "construct_allele"
        assert flagged[4]["construct_allele_cooccurrence"] == "yes"
        assert flagged[4]["construct_allele_ad_ratio"] == f"{46 / 136:.6f}"
        assert flagged[5]["plasmid_artefact"] == "construct_allele"
        assert flagged[6]["plasmid_artefact"] == "construct_allele"
        assert flagged[6]["matching_constructs"] == "pcDNA5_FRT_TO_CLCN2_Y17H"

        # APA41_To single KCNJ5 L168R: construct_allele, cooc=no
        assert flagged[7]["plasmid_artefact"] == "construct_allele"
        assert flagged[7]["construct_allele_cooccurrence"] == "no"
        assert flagged[7]["construct_allele_ad_ratio"] == "NA"

        # APA23-N is construct_allele; APA23_To MUST be "none" (not conflated with APA23-N)
        assert flagged[8]["plasmid_artefact"] == "construct_allele"
        assert flagged[9]["plasmid_artefact"] == "none"
        assert flagged[9]["plasmid_flag"] == "none"
        assert float(flagged[9]["f_plasmid"]) == 0.0
        assert float(flagged[9]["plasmid_explainable_vaf"]) == 0.0

        # CACNA1H and CLCN2 cdna_junction_artefact deletions and insertions
        for idx in (10, 11, 12, 13):
            assert flagged[idx]["plasmid_artefact"] == "cdna_junction_artefact"
            assert float(flagged[idx]["plasmid_explainable_vaf"]) > 0.0

    @pytest.mark.unit
    def test_flag_vcf_variants_vcf_fixture(
        self, tmp_path: Path, apa_summary_rows: list[dict[str, Any]]
    ) -> None:
        from plasmicheck.scripts.compare_alignments import flag_vcf_variants

        vcf_in = tmp_path / "variants.vcf"
        tsv_out = tmp_path / "vcf_flagged.tsv"
        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr11,length=135086622>\n"
            "##contig=<ID=chr16,length=90338345>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAPA18_To\tAPA23_To\tAPA41_To\n"
            "chr11\t128911724\t.\tG\tA\t.\tPASS\t.\tGT:DP:AD\t0/1:200:151,49\t0/0:150:150,0\t0/0:160:160,0\n"
            "chr11\t128911776\t.\tT\tG\t.\tPASS\t.\tGT:DP:AD\t0/1:200:185,15\t0/1:150:115,35\t0/1:160:128,32\n"
            f"chr16\t1210649\t.\t{'A' * 98}\tA\t.\tPASS\t.\tGT:DP:AD\t0/0:100:100,0\t0/0:100:100,0\t0/1:120:75,45\n"
        )
        vcf_in.write_text(vcf_content)

        flagged = flag_vcf_variants(str(vcf_in), str(tsv_out), apa_summary_rows)
        by_sample_pos = {(r["sample"], int(r["pos"])): r for r in flagged}

        # 0/0 calls skipped; 5 non-ref calls emitted
        assert len(flagged) == 5

        # APA18_To has both G151R and L168R with exact AD ratio 15/49 = 0.306122
        r_18_g151r = by_sample_pos[("APA18_To", 128911724)]
        assert r_18_g151r["plasmid_artefact"] == "construct_allele"
        assert r_18_g151r["construct_allele_cooccurrence"] == "yes"
        assert r_18_g151r["construct_allele_ad_ratio"] == "0.306122"

        r_18_l168r = by_sample_pos[("APA18_To", 128911776)]
        assert r_18_l168r["plasmid_artefact"] == "construct_allele"
        assert r_18_l168r["construct_allele_ad_ratio"] == "0.306122"

        # APA23_To somatic L168R is "none"
        r_23_l168r = by_sample_pos[("APA23_To", 128911776)]
        assert r_23_l168r["plasmid_artefact"] == "none"

        # APA41_To has L168R (construct_allele) and CACNA1H intron deletion (cdna_junction_artefact)
        assert by_sample_pos[("APA41_To", 128911776)]["plasmid_artefact"] == "construct_allele"
        assert by_sample_pos[("APA41_To", 1210649)]["plasmid_artefact"] == "cdna_junction_artefact"

    @pytest.mark.unit
    def test_compare_alignments_emits_depth_excess_and_flagged_variants(
        self, tmp_path: Path
    ) -> None:
        from plasmicheck.scripts.compare_alignments import compare_alignments

        output_base = tmp_path / "APA18-T"
        cdna_file = tmp_path / "cDNA_positions.txt"
        cdna_file.write_text(
            "cDNA start position in plasmid: 100\n"
            "cDNA end position in plasmid: 500\n"
            "INSERT_REGION: (100, 500)\n"
            "SpliceJunctionCount: 0\n"
            "ExonCount: 1\n"
            "ExonLengths: 401\n"
            "IntronLengths: NA\n"
            "Strand: +\n"
            "Chromosome: chr11\n"
            "GenomicStart: 128910000\n"
            "GenomicEnd: 128913000\n"
            "GenomicIntrons: NA\n"
            "ConstructVariants: 128911724:G>A\n"
        )
        vcf_file = tmp_path / "input_variants.tsv"
        vcf_file.write_text(
            "sample\tchr\tpos\tref\talt\tdp\tad\tvaf\n"
            "APA18_To\tchr11\t128911724\tG\tA\t200\t49\t0.245\n"
        )
        dummy_metrics = {
            "insert": {
                "mean_depth": 120.0,
                "median_depth": 120.0,
                "breadth_1x": 1.0,
                "breadth_5x": 1.0,
                "cv": 0.1,
            },
            "backbone": {
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "breadth_1x": 0.0,
                "breadth_5x": 0.0,
                "cv": 0.0,
            },
        }
        empty_mm: dict[str, int] = {}
        mock_scan_results: tuple[float, dict[str, int], int, int, int, int, int, float] = (
            0.0,
            empty_mm,
            0,
            0,
            0,
            0,
            0,
            101.0,
        )
        with (
            patch("plasmicheck.scripts.compare_alignments._namesort_bam"),
            patch(
                "plasmicheck.scripts.compare_alignments._streaming_compare",
                return_value={
                    "Plasmid": 50,
                    "Human": 10,
                    "Tied": 0,
                    "Backbone_Only": 0,
                    "Ambiguous": 0,
                },
            ),
            patch(
                "plasmicheck.scripts.compare_alignments._scan_plasmid_bam",
                return_value=mock_scan_results,
            ),
            patch(
                "plasmicheck.scripts.compare_alignments.compute_region_coverage_metrics",
                return_value=(dummy_metrics, False),
            ),
            patch("plasmicheck.scripts.compare_alignments.os.remove"),
            patch("plasmicheck.scripts.compare_alignments.locate_genbank_file", return_value=None),
            patch("pysam.AlignmentFile"),
        ):
            compare_alignments(
                "p.bam",
                "h.bam",
                str(output_base),
                cdna_positions_file=str(cdna_file),
                vcf_file=str(vcf_file),
                dp_sample_median=100.0,
                dp_locus=120.0,
                r_cohort=0.60,
            )
            with patch(
                "plasmicheck.scripts.compare_alignments.compute_bed_median_depth"
            ) as mock_bed_dp:
                compare_alignments(
                    "p.bam",
                    "h.bam",
                    str(tmp_path / "APA18-T-zero"),
                    cdna_positions_file=str(cdna_file),
                    target_bed="targets.bed",
                    dp_sample_median=0.0,
                    sample_bam="s.bam",
                )
            mock_bed_dp.assert_not_called()

        summary_content = (tmp_path / "APA18-T.summary.tsv").read_text()
        assert "dp_sample_median\t100.000000" in summary_content
        assert "dp_locus\t120.000000" in summary_content
        assert "r_cohort\t0.600000" in summary_content
        assert "f_excess\t0.500000" in summary_content
        assert "ConstructVariants\t128911724:G>A" in summary_content

        flagged_tsv = tmp_path / "APA18-T.flagged_variants.tsv"
        assert flagged_tsv.exists()
        flagged_text = flagged_tsv.read_text()
        assert "construct_allele" in flagged_text
