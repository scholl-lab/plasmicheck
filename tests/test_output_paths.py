"""Unit tests for OutputPaths dataclass."""

from __future__ import annotations

from pathlib import Path

import pytest

from plasmicheck.output_paths import OutputPaths
from plasmicheck.scripts.generate_summary_reports import (
    _extract_sample_plasmid,
    read_compare_outputs,
)


@pytest.mark.unit
def test_output_paths_frozen() -> None:
    """Test that OutputPaths is frozen (immutable after creation)."""
    paths = OutputPaths(
        output_root=Path("/tmp/output"),
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    with pytest.raises(AttributeError, match="cannot assign to field"):
        paths.sample_name = "new_sample"  # type: ignore[misc]


@pytest.mark.unit
def test_output_paths_hashable() -> None:
    """Test that OutputPaths is hashable (frozen implies this)."""
    paths = OutputPaths(
        output_root=Path("/tmp/output"),
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    # Should be hashable and usable in sets/dicts
    paths_set = {paths}
    assert paths in paths_set
    paths_dict = {paths: "value"}
    assert paths_dict[paths] == "value"


@pytest.mark.unit
def test_base_path_construction() -> None:
    """Test base path construction."""
    paths = OutputPaths(
        output_root=Path("/data/output"),
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    assert paths.base == Path("/data/output/sample1/plasmid1")


@pytest.mark.unit
def test_subdirectory_properties() -> None:
    """Test all subdirectory properties."""
    paths = OutputPaths(
        output_root=Path("/data/output"),
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    base = Path("/data/output/sample1/plasmid1")
    assert paths.alignments_dir == base / "alignments"
    assert paths.references_dir == base / "references"
    assert paths.data_dir == base / "data"
    assert paths.reports_dir == base / "reports"


@pytest.mark.unit
def test_alignment_file_paths() -> None:
    """Test alignment file path properties have correct parent directories."""
    paths = OutputPaths(
        output_root=Path("/data/output"),
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    alignments_dir = Path("/data/output/sample1/plasmid1/alignments")

    assert paths.plasmid_bam == alignments_dir / "plasmid_alignment.bam"
    assert paths.plasmid_bam.parent == alignments_dir

    assert paths.spliced_human_bam == alignments_dir / "spliced_human_alignment.bam"
    assert paths.spliced_human_bam.parent == alignments_dir

    assert paths.spliced_bam == alignments_dir / "spliced_alignment.bam"
    assert paths.spliced_bam.parent == alignments_dir


@pytest.mark.unit
def test_reference_file_paths() -> None:
    """Test reference file path properties have correct parent directories."""
    paths = OutputPaths(
        output_root=Path("/data/output"),
        sample_name="sample1",
        plasmid_name="pUC19",
    )

    references_dir = Path("/data/output/sample1/pUC19/references")

    assert paths.plasmid_fasta == references_dir / "pUC19.fasta"
    assert paths.plasmid_fasta.parent == references_dir

    assert paths.plasmid_genbank == references_dir / "pUC19.gb"
    assert paths.plasmid_genbank.parent == references_dir

    assert paths.plasmid_index == references_dir / "pUC19.mmi"
    assert paths.plasmid_index.parent == references_dir

    assert paths.spliced_fasta == references_dir / "spliced_reference.fasta"
    assert paths.spliced_fasta.parent == references_dir

    assert paths.spliced_index == references_dir / "spliced_reference.mmi"
    assert paths.spliced_index.parent == references_dir


@pytest.mark.unit
def test_data_file_paths() -> None:
    """Test data file path properties have correct parent directories."""
    paths = OutputPaths(
        output_root=Path("/data/output"),
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    data_dir = Path("/data/output/sample1/plasmid1/data")

    # comparison_prefix is a string (used as prefix for TSV files)
    assert paths.comparison_prefix == str(data_dir / "comparison_result")
    assert isinstance(paths.comparison_prefix, str)

    assert paths.comparison_reads == data_dir / "comparison_result.reads_assignment.tsv"
    assert paths.comparison_reads.parent == data_dir

    assert paths.comparison_summary == data_dir / "comparison_result.summary.tsv"
    assert paths.comparison_summary.parent == data_dir

    assert paths.cdna_positions == data_dir / "cDNA_positions.txt"
    assert paths.cdna_positions.parent == data_dir


@pytest.mark.unit
def test_report_file_paths() -> None:
    """Test report file path properties have correct parent directories."""
    paths = OutputPaths(
        output_root=Path("/data/output"),
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    reports_dir = Path("/data/output/sample1/plasmid1/reports")

    assert paths.report_interactive == reports_dir / "report_interactive.html"
    assert paths.report_interactive.parent == reports_dir

    assert paths.report_static == reports_dir / "report_non_interactive.html"
    assert paths.report_static.parent == reports_dir


@pytest.mark.unit
def test_base_level_files() -> None:
    """Test files that live at base level (not in subdirectories)."""
    paths = OutputPaths(
        output_root=Path("/data/output"),
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    base = Path("/data/output/sample1/plasmid1")

    assert paths.md5sum_file == base / "md5sum.txt"
    assert paths.md5sum_file.parent == base


@pytest.mark.unit
def test_ensure_directories_creates_all_subdirs(tmp_path: Path) -> None:
    """Test ensure_directories creates all 5 subdirectories."""
    paths = OutputPaths(
        output_root=tmp_path,
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    # Initially, directories don't exist
    assert not paths.alignments_dir.exists()
    assert not paths.references_dir.exists()
    assert not paths.data_dir.exists()
    assert not paths.reports_dir.exists()
    # Create all directories
    paths.ensure_directories()

    # Now all should exist
    assert paths.alignments_dir.exists()
    assert paths.alignments_dir.is_dir()

    assert paths.references_dir.exists()
    assert paths.references_dir.is_dir()

    assert paths.data_dir.exists()
    assert paths.data_dir.is_dir()

    assert paths.reports_dir.exists()
    assert paths.reports_dir.is_dir()


@pytest.mark.unit
def test_ensure_directories_idempotent(tmp_path: Path) -> None:
    """Test ensure_directories can be called multiple times safely."""
    paths = OutputPaths(
        output_root=tmp_path,
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    # Call multiple times - should not raise
    paths.ensure_directories()
    paths.ensure_directories()
    paths.ensure_directories()

    # Directories should still exist
    assert paths.alignments_dir.exists()
    assert paths.references_dir.exists()
    assert paths.data_dir.exists()
    assert paths.reports_dir.exists()


@pytest.mark.unit
def test_plasmid_name_in_file_names() -> None:
    """Test that plasmid_name is correctly used in file names."""
    paths = OutputPaths(
        output_root=Path("/output"),
        sample_name="sample1",
        plasmid_name="pUC19_modified",
    )

    # Plasmid name should appear in reference file names (in filename itself)
    assert paths.plasmid_fasta.name == "pUC19_modified.fasta"
    assert paths.plasmid_genbank.name == "pUC19_modified.gb"
    assert paths.plasmid_index.name == "pUC19_modified.mmi"

    # Other files don't have plasmid name in filename (only in parent path)
    assert paths.spliced_fasta.name == "spliced_reference.fasta"
    assert paths.comparison_reads.name == "comparison_result.reads_assignment.tsv"


@pytest.mark.unit
def test_all_paths_return_path_objects() -> None:
    """Test that all path properties return Path objects (except comparison_prefix)."""
    paths = OutputPaths(
        output_root=Path("/output"),
        sample_name="sample1",
        plasmid_name="plasmid1",
    )

    # Directory properties
    assert isinstance(paths.base, Path)
    assert isinstance(paths.alignments_dir, Path)
    assert isinstance(paths.references_dir, Path)
    assert isinstance(paths.data_dir, Path)
    assert isinstance(paths.reports_dir, Path)

    # Alignment file properties
    assert isinstance(paths.plasmid_bam, Path)
    assert isinstance(paths.spliced_human_bam, Path)
    assert isinstance(paths.spliced_bam, Path)

    # Reference file properties
    assert isinstance(paths.plasmid_fasta, Path)
    assert isinstance(paths.plasmid_genbank, Path)
    assert isinstance(paths.plasmid_index, Path)
    assert isinstance(paths.spliced_fasta, Path)
    assert isinstance(paths.spliced_index, Path)

    # Data file properties (except comparison_prefix which is str)
    assert isinstance(paths.comparison_prefix, str)  # Exception: used as string prefix
    assert isinstance(paths.comparison_reads, Path)
    assert isinstance(paths.comparison_summary, Path)
    assert isinstance(paths.cdna_positions, Path)

    # Report file properties
    assert isinstance(paths.report_interactive, Path)
    assert isinstance(paths.report_static, Path)

    # Base-level files
    assert isinstance(paths.md5sum_file, Path)


@pytest.mark.unit
def test_extract_sample_plasmid_integration_with_output_paths(tmp_path: Path) -> None:
    """Test that OutputPaths directory structure is compatible with _extract_sample_plasmid."""
    # Create OutputPaths
    paths = OutputPaths(
        output_root=tmp_path,
        sample_name="sampleA",
        plasmid_name="plasmidX",
    )

    # Create directories
    paths.ensure_directories()

    # Write dummy comparison files
    paths.comparison_summary.write_text("Category\tMetric\tCount\nVerdict\tContaminated\t1\n")
    paths.comparison_reads.write_text("ReadName\tAssignedTo\nread1\tPlasmid\n")

    # Extract sample/plasmid from the paths
    sample_from_summary, plasmid_from_summary = _extract_sample_plasmid(
        str(paths.comparison_summary)
    )
    sample_from_reads, plasmid_from_reads = _extract_sample_plasmid(str(paths.comparison_reads))

    # Should match the original values
    assert sample_from_summary == "sampleA"
    assert plasmid_from_summary == "plasmidX"
    assert sample_from_reads == "sampleA"
    assert plasmid_from_reads == "plasmidX"


@pytest.mark.unit
def test_output_paths_with_read_compare_outputs(tmp_path: Path) -> None:
    """Integration test: OutputPaths creates structure, read_compare_outputs reads it."""
    # Create two sample/plasmid combinations using OutputPaths
    paths1 = OutputPaths(
        output_root=tmp_path,
        sample_name="sampleA",
        plasmid_name="plasmidX",
    )
    paths2 = OutputPaths(
        output_root=tmp_path,
        sample_name="sampleB",
        plasmid_name="plasmidY",
    )

    # Create directories
    paths1.ensure_directories()
    paths2.ensure_directories()

    # Write realistic TSV files for sample1/plasmid1
    paths1.comparison_reads.write_text(
        "ReadName\tAssignedTo\tPlasmidScore\tHumanScore\n"
        "read1\tPlasmid\t50\t30\n"
        "read2\tPlasmid\t60\t40\n"
    )
    # Note: read_compare_outputs renames "Count" to "Value"
    paths1.comparison_summary.write_text(
        "Category\tMetric\tCount\n"
        "Verdict\tVerdict\tContaminated\n"
        "Ratio\tRatio\t1.5\n"
        "Plasmid\tPlasmid\t2\n"
        "Human\tHuman\t0\n"
    )

    # Write realistic TSV files for sample2/plasmid2
    paths2.comparison_reads.write_text(
        "ReadName\tAssignedTo\tPlasmidScore\tHumanScore\n"
        "read3\tHuman\t20\t50\n"
        "read4\tHuman\t25\t55\n"
    )
    paths2.comparison_summary.write_text(
        "Category\tMetric\tCount\n"
        "Verdict\tVerdict\tNotContaminated\n"
        "Ratio\tRatio\t0.4\n"
        "Plasmid\tPlasmid\t0\n"
        "Human\tHuman\t2\n"
    )

    # Read outputs with read_compare_outputs
    reads_df, summary_df = read_compare_outputs(str(tmp_path))

    # Verify reads DataFrame is empty (#109.1: read_compare_outputs skips .reads_assignment.tsv)
    assert reads_df.empty

    # Verify summary DataFrame
    assert len(summary_df) == 8  # 4 rows from sample1 + 4 from sample2
    assert set(summary_df["Sample"]) == {"sampleA", "sampleB"}
    assert set(summary_df["Plasmid"]) == {"plasmidX", "plasmidY"}

    # Check verdict for each sample
    verdict_df = summary_df[summary_df["Category"] == "Verdict"]
    assert len(verdict_df) == 2
    sample1_verdict = verdict_df[verdict_df["Sample"] == "sampleA"]["Value"].iloc[0]
    sample2_verdict = verdict_df[verdict_df["Sample"] == "sampleB"]["Value"].iloc[0]
    assert sample1_verdict == "Contaminated"
    assert sample2_verdict == "NotContaminated"
