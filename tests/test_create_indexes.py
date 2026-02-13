"""Unit tests for plasmicheck.scripts.create_indexes."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from plasmicheck.scripts.create_indexes import create_indexes


class TestCreateIndexes:
    @pytest.mark.unit
    @patch("plasmicheck.scripts.create_indexes.subprocess.run")
    def test_uses_file_lock(self, mock_run: MagicMock, tmp_path: Path) -> None:
        """Verify that FileLock is used during index creation."""
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")
        with patch("plasmicheck.scripts.create_indexes.FileLock") as mock_lock:
            mock_lock.return_value.__enter__ = MagicMock(return_value=None)
            mock_lock.return_value.__exit__ = MagicMock(return_value=False)
            create_indexes(str(fasta), overwrite=True)
            mock_lock.assert_called_once_with(str(tmp_path / "ref.lock"), timeout=600)

    @pytest.mark.unit
    @patch("plasmicheck.scripts.create_indexes.subprocess.run")
    def test_skips_when_indexes_exist(self, mock_run: MagicMock, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")
        # Create pre-existing index files
        (tmp_path / "ref.mmi").write_text("index")
        (tmp_path / "ref.fasta.fai").write_text("index")

        create_indexes(str(fasta), overwrite=False)
        mock_run.assert_not_called()

    @pytest.mark.unit
    @patch("plasmicheck.scripts.create_indexes.subprocess.run")
    def test_fast_path_skips_without_lock(self, mock_run: MagicMock, tmp_path: Path) -> None:
        """Fast-path avoids acquiring the lock when both indexes exist."""
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")
        (tmp_path / "ref.mmi").write_text("index")
        (tmp_path / "ref.fasta.fai").write_text("index")

        with patch("plasmicheck.scripts.create_indexes.FileLock") as mock_lock:
            create_indexes(str(fasta), overwrite=False)
            mock_lock.assert_not_called()

    @pytest.mark.unit
    @patch("plasmicheck.scripts.create_indexes.subprocess.run")
    def test_partial_indexes_recreates_both(self, mock_run: MagicMock, tmp_path: Path) -> None:
        """When only one index file exists, both are recreated."""
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")
        (tmp_path / "ref.mmi").write_text("old")
        # ref.fasta.fai intentionally missing

        create_indexes(str(fasta), overwrite=False)
        assert mock_run.call_count == 2  # both minimap2 + samtools

    @pytest.mark.unit
    @patch("plasmicheck.scripts.create_indexes.subprocess.run")
    def test_overwrite_recreates_indexes(self, mock_run: MagicMock, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")
        (tmp_path / "ref.mmi").write_text("old")
        (tmp_path / "ref.fasta.fai").write_text("old")

        create_indexes(str(fasta), overwrite=True)
        assert mock_run.call_count == 2  # minimap2 + samtools

    @pytest.mark.unit
    @patch("plasmicheck.scripts.create_indexes.subprocess.run")
    def test_creates_indexes_when_none_exist(self, mock_run: MagicMock, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")

        create_indexes(str(fasta), overwrite=False)
        assert mock_run.call_count == 2
