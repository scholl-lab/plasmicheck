"""Centralized output path construction for PlasmiCheck pipeline.

This module provides a frozen dataclass that defines the directory structure
and file paths for all pipeline outputs, eliminating scattered path construction
across the codebase.

Output structure:
    output_root/sample_name/plasmid_name/
    ├── alignments/          # BAM alignment files
    ├── references/          # FASTA and index files
    ├── data/                # TSV results and positions
    ├── reports/             # HTML reports, plots, and assets
    │   ├── assets/          # plotly.min.js
    │   └── plots/           # Standalone viewable plot HTMLs
    └── md5sum.txt           # Checksums (base level)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class OutputPaths:
    """Centralized path construction for pipeline outputs.

    All pipeline components should construct paths via this class.
    Provides type safety and single source of truth for output structure.

    Args:
        output_root: Base output directory (e.g., "output")
        sample_name: Sanitized sample identifier
        plasmid_name: Sanitized plasmid identifier
    """

    output_root: Path
    sample_name: str
    plasmid_name: str

    # -------------------------------------------------------------------------
    # Directory properties
    # -------------------------------------------------------------------------

    @property
    def base(self) -> Path:
        """Sample/plasmid base directory."""
        return self.output_root / self.sample_name / self.plasmid_name

    @property
    def alignments_dir(self) -> Path:
        """Alignment files directory."""
        return self.base / "alignments"

    @property
    def references_dir(self) -> Path:
        """Reference sequences and indexes directory."""
        return self.base / "references"

    @property
    def data_dir(self) -> Path:
        """Analysis data and results directory."""
        return self.base / "data"

    @property
    def reports_dir(self) -> Path:
        """HTML reports directory."""
        return self.base / "reports"

    # -------------------------------------------------------------------------
    # Alignment file properties
    # -------------------------------------------------------------------------

    @property
    def plasmid_bam(self) -> Path:
        """Plasmid alignment BAM file."""
        return self.alignments_dir / "plasmid_alignment.bam"

    @property
    def spliced_human_bam(self) -> Path:
        """Spliced human reference alignment BAM file."""
        return self.alignments_dir / "spliced_human_alignment.bam"

    @property
    def spliced_bam(self) -> Path:
        """Spliced alignment BAM file (plasmid cDNA to human genome)."""
        return self.alignments_dir / "spliced_alignment.bam"

    # -------------------------------------------------------------------------
    # Reference file properties
    # -------------------------------------------------------------------------

    @property
    def plasmid_fasta(self) -> Path:
        """Plasmid reference FASTA file."""
        return self.references_dir / f"{self.plasmid_name}.fasta"

    @property
    def plasmid_genbank(self) -> Path:
        """Plasmid GenBank file (from xDNA conversion or copy)."""
        return self.references_dir / f"{self.plasmid_name}.gb"

    @property
    def plasmid_index(self) -> Path:
        """Plasmid minimap2 index file."""
        return self.references_dir / f"{self.plasmid_name}.mmi"

    @property
    def spliced_fasta(self) -> Path:
        """Spliced human reference FASTA file."""
        return self.references_dir / "spliced_reference.fasta"

    @property
    def spliced_index(self) -> Path:
        """Spliced human reference minimap2 index file."""
        return self.references_dir / "spliced_reference.mmi"

    # -------------------------------------------------------------------------
    # Data file properties
    # -------------------------------------------------------------------------

    @property
    def comparison_prefix(self) -> str:
        """Comparison output prefix (for .reads_assignment.tsv and .summary.tsv).

        Returns string (not Path) because compare_alignments uses it as a prefix.
        """
        return str(self.data_dir / "comparison_result")

    @property
    def comparison_reads(self) -> Path:
        """Per-read assignment TSV file."""
        return self.data_dir / "comparison_result.reads_assignment.tsv"

    @property
    def comparison_summary(self) -> Path:
        """Comparison summary TSV file."""
        return self.data_dir / "comparison_result.summary.tsv"

    @property
    def cdna_positions(self) -> Path:
        """cDNA positions file (insert region boundaries)."""
        return self.data_dir / "cDNA_positions.txt"

    # -------------------------------------------------------------------------
    # Report file properties
    # -------------------------------------------------------------------------

    @property
    def report_interactive(self) -> Path:
        """Interactive HTML report."""
        return self.reports_dir / "report_interactive.html"

    @property
    def report_static(self) -> Path:
        """Non-interactive HTML report."""
        return self.reports_dir / "report_non_interactive.html"

    # -------------------------------------------------------------------------
    # Base-level files
    # -------------------------------------------------------------------------

    @property
    def md5sum_file(self) -> Path:
        """MD5 checksum file (at base level, covers all subdirs)."""
        return self.base / "md5sum.txt"

    @property
    def igv_session(self) -> Path:
        """IGV session file; resource paths inside are relative to base."""
        return self.base / "igv_session.xml"

    # -------------------------------------------------------------------------
    # Methods
    # -------------------------------------------------------------------------

    def ensure_directories(self) -> None:
        """Create all output subdirectories.

        Creates the 4 subdirectories (alignments, references, data, reports)
        with parents=True and exist_ok=True for idempotent directory creation.
        """
        for directory in [
            self.alignments_dir,
            self.references_dir,
            self.data_dir,
            self.reports_dir,
        ]:
            directory.mkdir(parents=True, exist_ok=True)
