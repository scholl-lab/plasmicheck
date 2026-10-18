"""Unit tests for plasmicheck.scripts.generate_report."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from plasmicheck.scripts.generate_report import (
    downsample_data,
    extract_verdict_from_summary,
    generate_plots,
    load_data,
    main,
)


@pytest.mark.unit
def test_sample_plots_disable_cloud_sharing(
    tmp_path: Path, sample_reads_assignment_df: pd.DataFrame
) -> None:
    boxplot, _, scatter, _ = generate_plots(sample_reads_assignment_df, str(tmp_path))
    html_outputs = [boxplot, scatter, *(p.read_text() for p in tmp_path.rglob("*.html"))]
    assert len(html_outputs) == 4
    for html in html_outputs:
        assert '"showSendToCloud": false' in html


class TestDownsampleData:
    @pytest.mark.unit
    def test_under_limit_unchanged(self, sample_reads_assignment_df: pd.DataFrame) -> None:
        df = sample_reads_assignment_df  # 5 rows
        result, downsampled = downsample_data(df, 100)
        assert not downsampled
        assert len(result) == len(df)

    @pytest.mark.unit
    def test_over_limit_sampled(self) -> None:
        df = pd.DataFrame({"A": range(100)})
        result, downsampled = downsample_data(df, 10)
        assert downsampled
        assert len(result) == 10

    @pytest.mark.unit
    def test_exactly_at_limit(self) -> None:
        df = pd.DataFrame({"A": range(50)})
        result, downsampled = downsample_data(df, 50)
        assert not downsampled
        assert len(result) == 50

    @pytest.mark.unit
    def test_reproducible(self) -> None:
        df = pd.DataFrame({"A": range(100)})
        r1, _ = downsample_data(df, 10)
        r2, _ = downsample_data(df, 10)
        pd.testing.assert_frame_equal(r1.reset_index(drop=True), r2.reset_index(drop=True))


class TestExtractVerdictFromSummary:
    @pytest.mark.unit
    def test_has_verdict_row(self, sample_summary_df: pd.DataFrame) -> None:
        verdict = extract_verdict_from_summary(sample_summary_df)
        assert "contaminated" in verdict.lower()

    @pytest.mark.unit
    def test_no_verdict_row(self) -> None:
        df = pd.DataFrame(
            {
                "Category": ["Plasmid", "Human"],
                "Count": [100, 50],
            }
        )
        verdict = extract_verdict_from_summary(df)
        assert verdict == "Verdict not found in summary file"

    @pytest.mark.unit
    def test_multiple_verdict_rows(self) -> None:
        df = pd.DataFrame(
            {
                "Category": ["Verdict", "Verdict"],
                "Count": ["first verdict", "second verdict"],
            }
        )
        verdict = extract_verdict_from_summary(df)
        assert verdict == "first verdict"


class TestLoadData:
    @pytest.mark.unit
    def test_valid_tsv_files(self, tmp_path: Path) -> None:
        reads_file = tmp_path / "reads.tsv"
        summary_file = tmp_path / "summary.tsv"

        reads_file.write_text(
            "ReadID\tAssignedTo\tPlasmidScore\tHumanScore\tPlasmidCIGAR\tHumanCIGAR\tPlasmidMapQ\tHumanMapQ\n"
            "read1\tPlasmid\t70\t30\t100M\t50M50S\t60\t20\n"
        )
        summary_file.write_text(
            "Category\tCount\nPlasmid\t100\nHuman\t50\nVerdict\tSample is contaminated\n"
        )

        reads_df, summary_df = load_data(str(reads_file), str(summary_file))
        assert "AssignedTo" in reads_df.columns
        assert "ReadID" not in reads_df.columns
        assert "Category" in summary_df.columns
        assert len(reads_df) == 1
        assert len(summary_df) == 3


@pytest.mark.unit
def test_report_autoescape_escapes_script(tmp_path: Path) -> None:
    """Verify Jinja2 template autoescaping prevents HTML injection."""
    reads_file = tmp_path / "reads.tsv"
    summary_file = tmp_path / "summary.tsv"
    output_folder = tmp_path / "output"
    output_folder.mkdir()

    reads_file.write_text("AssignedTo\tPlasmidScore\tHumanScore\nPlasmid\t70\t30\n")
    summary_file.write_text(
        "Category\tCount\nPlasmid\t100\nHuman\t50\nVerdict\tSample <script>alert(1)</script> is clean\n"
    )

    with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
        mock_plots.return_value = ("<div>box</div>", None, "<div>scatter</div>", None)
        main(
            str(reads_file),
            str(summary_file),
            str(output_folder),
            command_line="<script>bad_cmd</script>",
        )

    report_path = output_folder / "report_interactive.html"
    assert report_path.exists()
    content = report_path.read_text()
    assert "<script>alert(1)</script>" not in content
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in content
    assert "<script>bad_cmd</script>" not in content
    assert "&lt;script&gt;bad_cmd&lt;/script&gt;" in content


@pytest.mark.unit
def test_main_uses_summary_ratio_and_handles_nan_score(tmp_path: Path) -> None:
    """Verify summary.tsv Ratio is respected and NaN scores are tolerated."""
    reads_file = tmp_path / "reads.tsv"
    summary_file = tmp_path / "summary.tsv"
    output_folder = tmp_path / "output"
    output_folder.mkdir()

    reads_file.write_text("AssignedTo\tPlasmidScore\tHumanScore\nPlasmid\t\t30\nHuman\t50\t\n")
    summary_file.write_text(
        "Category\tCount\nPlasmid\t1\nHuman\t1\nVerdict\tSample is contaminated\nRatio\t4.5\n"
    )

    with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
        mock_plots.return_value = ("<div>box</div>", None, "<div>scatter</div>", None)
        main(
            str(reads_file),
            str(summary_file),
            str(output_folder),
        )

    report_path = output_folder / "report_interactive.html"
    assert report_path.exists()
    content = report_path.read_text()
    assert "4.500" in content


class TestCoverageMetricsInReport:
    """Tests for coverage metrics display in HTML reports (Phase 9 Plan 02)."""

    @pytest.mark.unit
    def test_coverage_metrics_in_report(
        self,
        tmp_path: Path,
        sample_reads_assignment_df: pd.DataFrame,
        sample_summary_df_with_coverage: pd.DataFrame,
    ) -> None:
        """Verify coverage metrics are passed through to the template."""
        reads_file = tmp_path / "reads_assignment.tsv"
        summary_file = tmp_path / "summary.tsv"
        output_folder = tmp_path / "output"
        output_folder.mkdir()

        sample_reads_assignment_df.to_csv(reads_file, sep="\t", index=False)
        sample_summary_df_with_coverage.to_csv(summary_file, sep="\t", index=False)

        # Mock plot generation to avoid plotly/kaleido dependencies
        with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
            mock_plots.return_value = (
                "<div>box</div>",
                None,
                "<div>scatter</div>",
                None,
            )

            main(
                str(reads_file),
                str(summary_file),
                str(output_folder),
                command_line="test command",
            )

        # Verify interactive report was created
        report_path = output_folder / "report_interactive.html"
        assert report_path.exists()
        html_content = report_path.read_text()

        # Verify Coverage Analysis card is present
        assert "Coverage Analysis" in html_content

        # Verify at least one coverage metric value appears in the report
        assert "45.25" in html_content  # MeanDepthInsert from fixture

    @pytest.mark.unit
    def test_coverage_metrics_missing_graceful(
        self,
        tmp_path: Path,
        sample_reads_assignment_df: pd.DataFrame,
        sample_summary_df: pd.DataFrame,
    ) -> None:
        """Verify report renders without error when summary.tsv has no coverage rows."""
        reads_file = tmp_path / "reads_assignment.tsv"
        summary_file = tmp_path / "summary.tsv"
        output_folder = tmp_path / "output"
        output_folder.mkdir()

        sample_reads_assignment_df.to_csv(reads_file, sep="\t", index=False)
        sample_summary_df.to_csv(summary_file, sep="\t", index=False)

        # Mock plot generation
        with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
            mock_plots.return_value = (
                "<div>box</div>",
                None,
                "<div>scatter</div>",
                None,
            )

            main(
                str(reads_file),
                str(summary_file),
                str(output_folder),
                command_line="test command",
            )

        # Verify report was created and contains Coverage Analysis card
        report_path = output_folder / "report_interactive.html"
        assert report_path.exists()
        html_content = report_path.read_text()

        assert "Coverage Analysis" in html_content
        # Default values should be shown (0.00)
        assert "0.00" in html_content

    @pytest.mark.unit
    def test_coverage_fallback_warning(
        self,
        tmp_path: Path,
        sample_reads_assignment_df: pd.DataFrame,
    ) -> None:
        """Verify fallback warning appears when CoverageFallback is True."""
        reads_file = tmp_path / "reads_assignment.tsv"
        summary_file = tmp_path / "summary.tsv"
        output_folder = tmp_path / "output"
        output_folder.mkdir()

        sample_reads_assignment_df.to_csv(reads_file, sep="\t", index=False)

        # Create summary with CoverageFallback=True
        summary_df_fallback = pd.DataFrame(
            {
                "Category": [
                    "Plasmid",
                    "Human",
                    "Tied",
                    "Backbone_Only",
                    "Ambiguous",
                    "Verdict",
                    "Ratio",
                    "CoverageOutsideINSERT",
                    "MismatchesNearINSERT",
                    "MeanDepthInsert",
                    "MedianDepthInsert",
                    "BreadthInsert",
                    "BreadthInsert_5x",
                    "CoverageCV_Insert",
                    "MeanDepthBackbone",
                    "MedianDepthBackbone",
                    "BreadthBackbone",
                    "BreadthBackbone_5x",
                    "CoverageCV_Backbone",
                    "CoverageFallback",
                ],
                "Count": [
                    100,
                    50,
                    10,
                    15,
                    5,
                    "Sample is contaminated with plasmid DNA",
                    2.0,
                    0.1234,
                    "{'with_mismatches_or_clipping': 5, 'without_mismatches_or_clipping': 10}",
                    45.25,
                    42.00,
                    0.95,
                    0.85,
                    0.22,
                    0.00,
                    0.00,
                    0.00,
                    0.00,
                    0.00,
                    "True",  # CoverageFallback is True
                ],
            }
        )
        summary_df_fallback.to_csv(summary_file, sep="\t", index=False)

        # Mock plot generation
        with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
            mock_plots.return_value = (
                "<div>box</div>",
                None,
                "<div>scatter</div>",
                None,
            )

            main(
                str(reads_file),
                str(summary_file),
                str(output_folder),
                command_line="test command",
            )

        # Verify fallback warning appears
        report_path = output_folder / "report_interactive.html"
        assert report_path.exists()
        html_content = report_path.read_text()

        assert "Insert region not defined" in html_content


@pytest.mark.unit
def test_truncate_list_filter_semantics() -> None:
    """Verify _truncate_list_filter semantics for empty, missing, and truncated items."""
    from plasmicheck.scripts.generate_report import _truncate_list_filter

    assert _truncate_list_filter("") == ""
    assert _truncate_list_filter("NA") == "NA"
    assert _truncate_list_filter("nan") == "nan"
    assert _truncate_list_filter("1,2,3,4,5") == "1,2,3,4,5"
    assert _truncate_list_filter("1,2,3,4,5,6") == "1,2,3,4,5,..."


@pytest.mark.unit
def test_compute_gauge_position_bounds() -> None:
    """Verify _compute_gauge_position under standard and zero/inverted unclear intervals."""
    from plasmicheck.scripts.generate_report import _compute_gauge_position

    # Standard bounds: 0.6 to 1.0
    pos_low = _compute_gauge_position(0.3, 0.6, 1.0)
    assert 3.0 <= pos_low < 33.34
    pos_mid = _compute_gauge_position(0.8, 0.6, 1.0)
    assert 33.33 <= pos_mid <= 66.67
    pos_high = _compute_gauge_position(5.0, 0.6, 1.0)
    assert 66.67 <= pos_high <= 95.0

    # Inverted / zero-width unclear interval (e.g. upper <= lower)
    pos_inverted = _compute_gauge_position(0.3, 0.6, 0.4)
    assert 3.0 <= pos_inverted <= 95.0


@pytest.mark.unit
def test_report_effective_threshold_from_summary_tsv(tmp_path: Path) -> None:
    """Verify main() extracts threshold from summary.tsv Threshold row and passes it."""
    from plasmicheck.scripts.generate_report import main

    reads_file = tmp_path / "reads.tsv"
    summary_file = tmp_path / "summary.tsv"
    output_folder = tmp_path / "output"
    output_folder.mkdir()

    reads_file.write_text("AssignedTo\tPlasmidScore\tHumanScore\nPlasmid\t70\t30\n")
    summary_file.write_text(
        "Category\tCount\nPlasmid\t100\nHuman\t50\nRatio\t0.9\nThreshold\t0.95\nVerdict\tSample contamination status is unclear\n"
    )

    with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
        mock_plots.return_value = ("<div>box</div>", None, "<div>scatter</div>", None)
        main(str(reads_file), str(summary_file), str(output_folder))

        # Check threshold passed to generate_plots was 0.95
        assert mock_plots.call_count == 1
        call_kwargs = mock_plots.call_args[1]
        assert call_kwargs["threshold"] == 0.95
        assert call_kwargs["unclear"]["upper_bound"] == 0.95

    report_path = output_folder / "report_interactive.html"
    assert report_path.exists()
    content = report_path.read_text()
    assert "0.95" in content


@pytest.mark.unit
def test_report_custom_threshold_reclassifies_verdict(tmp_path: Path) -> None:
    """Verify custom threshold reclassifies display verdict from stored ratio."""
    from plasmicheck.scripts.generate_report import main

    reads_file = tmp_path / "reads.tsv"
    summary_file = tmp_path / "summary.tsv"
    output_folder = tmp_path / "output"
    output_folder.mkdir()

    reads_file.write_text("AssignedTo\tPlasmidScore\tHumanScore\nPlasmid\t70\t30\n")
    # Stored verdict is Contaminated (based on default 0.8), but Ratio is 0.9
    summary_file.write_text(
        "Category\tCount\nPlasmid\t100\nHuman\t50\nRatio\t0.9\nVerdict\tSample is contaminated with plasmid DNA\n"
    )

    with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
        mock_plots.return_value = ("<div>box</div>", None, "<div>scatter</div>", None)
        # Pass threshold=1.2: with 0.6..1.2, ratio 0.9 becomes Unclear
        main(str(reads_file), str(summary_file), str(output_folder), threshold=1.2)

    report_path = output_folder / "report_interactive.html"
    assert report_path.exists()
    content = report_path.read_text()
    # Display verdict should be reclassified to Unclear
    assert "Unclear" in content


@pytest.mark.unit
def test_report_anonymize_paths(tmp_path: Path) -> None:
    """Verify --anonymize strips directory paths from report metadata."""
    from plasmicheck.scripts.generate_report import main

    reads_file = tmp_path / "reads.tsv"
    summary_file = tmp_path / "summary.tsv"
    output_folder = tmp_path / "output"
    output_folder.mkdir()

    reads_file.write_text("AssignedTo\tPlasmidScore\tHumanScore\nPlasmid\t70\t30\n")
    summary_file.write_text(
        "Category\tCount\nPlasmid\t100\nHuman\t50\nRatio\t2.0\nVerdict\tSample is contaminated with plasmid DNA\n"
    )

    with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
        mock_plots.return_value = ("<div>box</div>", None, "<div>scatter</div>", None)
        main(
            str(reads_file),
            str(summary_file),
            str(output_folder),
            human_fasta="/secret/patient/data/human.fa",
            plasmid_gb="/private/path/plasmid.gb",
            sequencing_file="/confidential/patient123/reads.bam",
            command_line="plasmicheck pipeline -hf /secret/patient/data/human.fa -pf /private/path/plasmid.gb",
            anonymize=True,
        )

    report_path = output_folder / "report_interactive.html"
    assert report_path.exists()
    content = report_path.read_text()
    assert "/secret/patient/data" not in content
    assert "/private/path" not in content
    assert "/confidential/patient123" not in content
    assert "human.fa" in content
    assert "plasmid.gb" in content
    assert "reads.bam" in content
    assert "--anonymize was enabled" in content


@pytest.mark.unit
def test_load_verified_igv_b64_integrity(tmp_path: Path) -> None:
    """Verify vendored igv.min.js.gz SHA-256 checks succeed and reject tampered bytes."""
    from plasmicheck.scripts.generate_report import _load_verified_igv_b64

    b64 = _load_verified_igv_b64()
    assert len(b64) > 400_000

    tampered = tmp_path / "tampered.js.gz"
    tampered.write_bytes(b"corrupted-gzip-payload")
    with (
        patch("plasmicheck.scripts.generate_report.get_resource_path", return_value=tampered),
        pytest.raises(RuntimeError, match="compressed SHA-256 mismatch"),
    ):
        _load_verified_igv_b64()


@pytest.mark.unit
def test_build_igv_context_and_embed_in_report(tmp_path: Path) -> None:
    """Verify _build_igv_context embeds offline igv.js, BED tracks, BAM reads, and loci."""
    import pysam

    from plasmicheck.scripts.generate_report import _build_igv_context, main

    base_dir = tmp_path / "sample_out"
    align_dir = base_dir / "alignments"
    ref_dir = base_dir / "references"
    data_dir = base_dir / "data"
    reports_dir = base_dir / "reports"
    for d in (align_dir, ref_dir, data_dir, reports_dir):
        d.mkdir(parents=True)

    # Write plasmid FASTA (1000 bp)
    fasta_path = ref_dir / "pTEST.fasta"
    fasta_path.write_text(">pTEST\n" + ("ACGT" * 250) + "\n")

    # Write cDNA_positions.txt with 2 exons and human genomic coordinates
    cdna_path = data_dir / "cDNA_positions.txt"
    cdna_path.write_text(
        "cDNA start position in plasmid: 200\n"
        "cDNA end position in plasmid: 500\n"
        "SpliceJunctionCount: 1\n"
        "ExonCount: 2\n"
        "ExonLengths: 150,150\n"
        "IntronLengths: 1000\n"
        "Strand: +\n"
        "Chromosome: chr16\n"
        "GenomicStart: 1153737\n"
        "GenomicEnd: 1220991\n"
    )

    # Write igv_session.xml
    (base_dir / "igv_session.xml").write_text("<Session/>")

    # Write real small BAM with 5 mapped reads (including 5' VI, J1 EE at 350, and J1 EI clipped at 350)
    bam_path = align_dir / "plasmid_alignment.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "pTEST", "LN": 1000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
        for i in range(4):
            seg = pysam.AlignedSegment(out_bam.header)
            seg.query_name = f"read_{i}"
            seg.query_sequence = "ACGT" * 25
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = (
                180 + i * 45
            )  # 180..280 (5' VI), 225..325, 270..370 (J1 EE), 315..415 (J1 EE)
            seg.mapping_quality = 60
            seg.cigartuples = [(0, 100)]
            seg.query_qualities = pysam.qualitystring_to_array("I" * 100)
            out_bam.write(seg)
        # Add 5th read soft-clipped at J1 (350) -> Human_ExonIntron (EI)
        ei_seg = pysam.AlignedSegment(out_bam.header)
        ei_seg.query_name = "read_ei"
        ei_seg.query_sequence = "ACGT" * 25
        ei_seg.flag = 0
        ei_seg.reference_id = 0
        ei_seg.reference_start = 350
        ei_seg.mapping_quality = 60
        ei_seg.cigartuples = [(4, 25), (0, 75)]  # 25S75M starting at 350
        ei_seg.query_qualities = pysam.qualitystring_to_array("I" * 100)
        out_bam.write(ei_seg)

    ctx = _build_igv_context(str(base_dir), reports_dir=str(reports_dir), max_reads=3)
    assert ctx["igv_embedded"] is True
    assert ctx["igv_downsampled"] is True
    assert ctx["igv_read_count"] == 3
    assert ctx["igv_total_primary_reads"] == 5
    assert ctx["igv_splice_junction_count"] == 1
    assert ctx["igv_human_genomic_locus"] == "chr16:1,153,737-1,220,991 (+)"
    assert ctx["igv_session_relpath"] == "../igv_session.xml"
    assert len(ctx["igv_loci"]) == 5  # Full, Insert, 5', 3', Splice Junction 1
    j1_locus = ctx["igv_loci"][4]
    assert j1_locus["ee_count"] == 2
    assert j1_locus["ei_count"] == 1
    assert "Spliced human intron: 1,000 bp" in j1_locus["note"]
    assert "const IGV_GZ_B64 =" in ctx["igv_payload"]
    assert "const IGV_BROWSER_OPTIONS =" in ctx["igv_payload"]
    assert '"format": "bedgraph"' in ctx["igv_payload"]
    assert '"colorBy": "tag:PC"' in ctx["igv_payload"]

    # Also render full interactive + static report with #121/#122 metrics
    reads_file = data_dir / "reads.tsv"
    summary_file = data_dir / "summary.tsv"
    reads_file.write_text("AssignedTo\tPlasmidScore\tHumanScore\nPlasmid\t90\t20\n")
    summary_file.write_text(
        "Category\tCount\n"
        "Plasmid\t80\n"
        "Human\t20\n"
        "Ratio\t4.0\n"
        "Verdict\tSample is contaminated with plasmid DNA\n"
        "SpliceJunctionCount\t1\n"
        "ExonCount\t2\n"
        "SpliceJunctionEvents\t42\n"
        "VectorInsertJunctionReads\t11\n"
        "EstimatedContaminationFraction\t0.8250\n"
    )

    with patch("plasmicheck.scripts.generate_report.generate_plots") as mock_plots:
        png_box = reports_dir / "box.png"
        png_scat = reports_dir / "scat.png"
        png_box.write_bytes(b"\x89PNG\r\n\x1a\n")
        png_scat.write_bytes(b"\x89PNG\r\n\x1a\n")
        mock_plots.return_value = (
            "<div>box</div>",
            str(png_box),
            "<div>scatter</div>",
            str(png_scat),
        )
        main(
            str(reads_file),
            str(summary_file),
            str(base_dir),
            reports_dir=str(reports_dir),
            static_report=True,
        )

    interactive_html = (reports_dir / "report_interactive.html").read_text()
    assert "Est. Contamination:" in interactive_html
    assert "82.5% (0.8250)" in interactive_html
    assert "chr16:1,153,737-1,220,991 (+)" in interactive_html
    assert "Splice-Junction Events:" in interactive_html
    assert "Vector Junctions:" in interactive_html
    assert "igvLociTable" in interactive_html
    assert "DecompressionStream('gzip')" in interactive_html
    assert "../igv_session.xml" in interactive_html

    static_html = (reports_dir / "report_non_interactive.html").read_text()
    assert (
        "Interactive IGV alignment browser is only active in report_interactive.html."
        in static_html
    )
    assert "const IGV_GZ_B64 =" not in static_html
