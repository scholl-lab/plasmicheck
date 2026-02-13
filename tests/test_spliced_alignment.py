"""Unit tests for plasmicheck.scripts.spliced_alignment."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pysam", reason="pysam not available on this platform")

from plasmicheck.scripts.spliced_alignment import find_fasta_file


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
