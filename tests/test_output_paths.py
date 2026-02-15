"""Unit tests for OutputPaths dataclass."""

from __future__ import annotations

from pathlib import Path

import pytest

from plasmicheck.output_paths import OutputPaths


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
    assert paths.plots_dir == base / "plots"


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
    assert not paths.plots_dir.exists()

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

    assert paths.plots_dir.exists()
    assert paths.plots_dir.is_dir()


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
    assert paths.plots_dir.exists()


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
    assert isinstance(paths.plots_dir, Path)

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
