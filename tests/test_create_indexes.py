"""Unit tests for plasmicheck.scripts.create_indexes."""

from __future__ import annotations

import os
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

    @pytest.mark.unit
    @patch("plasmicheck.scripts.create_indexes.subprocess.run")
    def test_threads_parameter_passed_to_minimap2(
        self, mock_run: MagicMock, tmp_path: Path
    ) -> None:
        """threads parameter is forwarded to minimap2."""
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")

        create_indexes(str(fasta), overwrite=True, threads=16)
        minimap2_call = mock_run.call_args_list[0][0][0]
        assert minimap2_call[0] == "minimap2"
        assert minimap2_call[1] == "-t"
        assert minimap2_call[2] == "16"

    @pytest.mark.unit
    @patch("plasmicheck.scripts.create_indexes.subprocess.run")
    def test_empty_fai_recreates_indexes(self, mock_run: MagicMock, tmp_path: Path) -> None:
        """When .fai exists but is 0 bytes, both indexes are recreated."""
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")
        (tmp_path / "ref.mmi").write_text("index")
        (tmp_path / "ref.fasta.fai").touch()  # 0 bytes

        create_indexes(str(fasta), overwrite=False)
        assert mock_run.call_count == 2

    @pytest.mark.unit
    def test_atomic_index_creation_and_replace(self, tmp_path: Path) -> None:
        """Minimap2 and samtools write to temporary files which are atomically replaced."""
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")

        # Fake subprocess.run to simulate tools creating their temp files
        def fake_run(cmd: list[str], check: bool = True) -> None:
            if cmd[0] == "minimap2":
                tmp_out = cmd[-2]
                Path(tmp_out).write_text("fake_mmi_content")
            elif cmd[0] == "samtools":
                assert "-o" in cmd
                tmp_out = cmd[cmd.index("-o") + 1]
                Path(tmp_out).write_text("fake_fai_content")

        with patch("plasmicheck.scripts.create_indexes.subprocess.run", side_effect=fake_run):
            create_indexes(str(fasta), overwrite=True)

        mmi_file = tmp_path / "ref.mmi"
        assert mmi_file.exists()
        assert mmi_file.read_text() == "fake_mmi_content"

        fai_file = tmp_path / "ref.fasta.fai"
        assert fai_file.exists()
        assert fai_file.read_text() == "fake_fai_content"

        # Verify no leftover .tmp files
        tmp_mmi = tmp_path / f"ref.mmi.tmp.{os.getpid()}"
        assert not tmp_mmi.exists()
        tmp_fai = tmp_path / f"ref.fasta.fai.{os.getpid()}.tmp"
        assert not tmp_fai.exists()

    @pytest.mark.unit
    def test_atomic_fai_cleanup_on_error(self, tmp_path: Path) -> None:
        """Temporary .fai file is cleaned up if samtools faidx fails."""
        import subprocess

        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">seq\nATCG")

        def fake_run(cmd: list[str], check: bool = True) -> None:
            if cmd[0] == "minimap2":
                tmp_out = cmd[-2]
                Path(tmp_out).write_text("fake_mmi_content")
            elif cmd[0] == "samtools":
                tmp_out = cmd[cmd.index("-o") + 1]
                Path(tmp_out).write_text("partial_fai_content")
                raise subprocess.CalledProcessError(1, cmd)

        with (
            patch("plasmicheck.scripts.create_indexes.subprocess.run", side_effect=fake_run),
            pytest.raises(subprocess.CalledProcessError),
        ):
            create_indexes(str(fasta), overwrite=True)

        tmp_fai = tmp_path / f"ref.fasta.fai.{os.getpid()}.tmp"
        assert not tmp_fai.exists()
