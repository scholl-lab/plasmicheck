"""Unit tests for plasmicheck.scripts.utils."""

from __future__ import annotations

import hashlib
import os
import subprocess
import tarfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from plasmicheck.scripts.utils import (
    archive_output_folder,
    calculate_md5,
    quality_control,
    run_command,
    sanitize_filename,
    validate_file_existence,
    validate_file_format,
    write_md5sum,
)

# -- sanitize_filename ----------------------------------------------------------


class TestSanitizeFilename:
    @pytest.mark.unit
    def test_alphanumeric_passthrough(self) -> None:
        assert sanitize_filename("my-file_name") == "my-file_name"

    @pytest.mark.unit
    def test_spaces_replaced(self) -> None:
        assert sanitize_filename("my file") == "my_file"

    @pytest.mark.unit
    def test_special_chars_replaced(self) -> None:
        assert sanitize_filename("file@#$.txt") == "file____txt"

    @pytest.mark.unit
    def test_dots_replaced(self) -> None:
        assert sanitize_filename("name.ext") == "name_ext"

    @pytest.mark.unit
    def test_empty_string(self) -> None:
        assert sanitize_filename("") == ""

    @pytest.mark.unit
    def test_already_clean(self) -> None:
        assert sanitize_filename("clean_file-name") == "clean_file-name"


# -- validate_file_existence ----------------------------------------------------


class TestValidateFileExistence:
    @pytest.mark.unit
    def test_all_files_exist(self, tmp_path: Path) -> None:
        f1 = tmp_path / "a.txt"
        f2 = tmp_path / "b.txt"
        f1.write_text("x")
        f2.write_text("y")
        validate_file_existence([str(f1), str(f2)])  # no exception

    @pytest.mark.unit
    def test_one_missing(self, tmp_path: Path) -> None:
        f1 = tmp_path / "exists.txt"
        f1.write_text("x")
        missing = str(tmp_path / "nope.txt")
        with pytest.raises(FileNotFoundError, match=r"nope\.txt"):
            validate_file_existence([str(f1), missing])

    @pytest.mark.unit
    def test_multiple_missing(self, tmp_path: Path) -> None:
        m1 = str(tmp_path / "miss1.txt")
        m2 = str(tmp_path / "miss2.txt")
        with pytest.raises(FileNotFoundError, match=r"miss1\.txt") as exc_info:
            validate_file_existence([m1, m2])
        assert "miss2.txt" in str(exc_info.value)


# -- validate_file_format -------------------------------------------------------


class TestValidateFileFormat:
    @pytest.mark.unit
    def test_valid_fasta(self) -> None:
        validate_file_format(["ref.fasta", "other.fa"], "fasta")

    @pytest.mark.unit
    def test_valid_fastq_gz(self) -> None:
        validate_file_format(["reads.fastq.gz", "reads2.fq.gz"], "fastq")

    @pytest.mark.unit
    def test_invalid_extension(self) -> None:
        with pytest.raises(ValueError, match="Incorrect format"):
            validate_file_format(["reads.xyz"], "fasta")

    @pytest.mark.unit
    def test_unsupported_format_name(self) -> None:
        with pytest.raises(ValueError, match="Unsupported format"):
            validate_file_format(["file.txt"], "nonexistent_format")

    @pytest.mark.unit
    def test_string_format_not_list(self) -> None:
        validate_file_format(["data.bam"], "bam")  # string format, not list

    @pytest.mark.unit
    def test_list_of_formats(self) -> None:
        validate_file_format(["ref.fasta", "reads.fastq"], ["fasta", "fastq"])


# -- quality_control ------------------------------------------------------------


class TestQualityControl:
    @pytest.mark.unit
    def test_valid_files_valid_format(self, tmp_path: Path) -> None:
        f = tmp_path / "test.fasta"
        f.write_text(">seq\nATCG")
        quality_control([str(f)], "fasta")

    @pytest.mark.unit
    def test_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            quality_control([str(tmp_path / "nope.fasta")], "fasta")

    @pytest.mark.unit
    def test_wrong_format(self, tmp_path: Path) -> None:
        f = tmp_path / "test.xyz"
        f.write_text("data")
        with pytest.raises(ValueError):
            quality_control([str(f)], "fasta")

    @pytest.mark.unit
    def test_no_format_check(self, tmp_path: Path) -> None:
        f = tmp_path / "anything.xyz"
        f.write_text("data")
        quality_control([str(f)], None)  # only checks existence


# -- calculate_md5 --------------------------------------------------------------


class TestCalculateMd5:
    @pytest.mark.unit
    def test_known_content(self, tmp_path: Path) -> None:
        f = tmp_path / "test.txt"
        content = b"hello world"
        f.write_bytes(content)
        expected = hashlib.md5(content).hexdigest()
        assert calculate_md5(str(f)) == expected


# -- write_md5sum ---------------------------------------------------------------


class TestWriteMd5sum:
    @pytest.mark.unit
    def test_creates_md5sum_file(self, tmp_path: Path) -> None:
        f = tmp_path / "data.txt"
        f.write_text("test content")
        write_md5sum(str(f), "input", str(tmp_path))

        md5_file = tmp_path / "md5sum.txt"
        assert md5_file.exists()
        lines = md5_file.read_text().strip().split("\n")
        assert len(lines) == 1
        parts = lines[0].split("\t")
        assert parts[0] == "input"
        assert parts[1] == str(f)
        assert len(parts[2]) == 32  # md5 hex length

    @pytest.mark.unit
    def test_appends_to_existing(self, tmp_path: Path) -> None:
        f1 = tmp_path / "a.txt"
        f2 = tmp_path / "b.txt"
        f1.write_text("aaa")
        f2.write_text("bbb")
        write_md5sum(str(f1), "type1", str(tmp_path))
        write_md5sum(str(f2), "type2", str(tmp_path))

        md5_file = tmp_path / "md5sum.txt"
        lines = md5_file.read_text().strip().split("\n")
        assert len(lines) == 2


# -- archive_output_folder ------------------------------------------------------


class TestArchiveOutputFolder:
    @pytest.mark.unit
    def test_creates_tar_gz_default_name(self, tmp_path: Path) -> None:
        # Create folder structure: parent/sample/plasmid/
        sample_dir = tmp_path / "sample1" / "plasmidA"
        sample_dir.mkdir(parents=True)
        (sample_dir / "result.txt").write_text("data")

        archive_path = archive_output_folder(str(sample_dir))
        assert os.path.isfile(archive_path)
        assert archive_path.endswith(".tar.gz")
        assert "sample1" in os.path.basename(archive_path)
        assert "plasmidA" in os.path.basename(archive_path)

    @pytest.mark.unit
    def test_custom_archive_name(self, tmp_path: Path) -> None:
        folder = tmp_path / "out"
        folder.mkdir()
        (folder / "file.txt").write_text("content")

        archive_path = archive_output_folder(str(folder), archive_name="custom.tar.gz")
        assert os.path.isfile(archive_path)
        assert archive_path.endswith("custom.tar.gz")

    @pytest.mark.unit
    def test_archive_contains_files(self, tmp_path: Path) -> None:
        folder = tmp_path / "data_folder"
        folder.mkdir()
        (folder / "a.txt").write_text("aaa")
        (folder / "b.txt").write_text("bbb")

        archive_path = archive_output_folder(str(folder), archive_name="test.tar.gz")
        with tarfile.open(archive_path, "r:gz") as tar:
            names = tar.getnames()
            assert any("a.txt" in n for n in names)
            assert any("b.txt" in n for n in names)


# -- run_command ----------------------------------------------------------------


class TestRunCommand:
    @pytest.mark.unit
    @patch("plasmicheck.scripts.utils.subprocess.run")
    def test_successful_command(self, mock_run: MagicMock) -> None:
        mock_run.return_value = subprocess.CompletedProcess(
            args="echo hello", returncode=0, stdout="hello\n", stderr=""
        )
        result = run_command("echo hello")
        assert result.returncode == 0
        mock_run.assert_called_once()

    @pytest.mark.unit
    @patch("plasmicheck.scripts.utils.time.sleep")
    @patch("plasmicheck.scripts.utils.subprocess.run")
    def test_retries_on_failure(self, mock_run: MagicMock, mock_sleep: MagicMock) -> None:
        # First call fails, second succeeds
        mock_run.side_effect = [
            subprocess.CalledProcessError(1, "cmd", stderr="error"),
            subprocess.CompletedProcess(args="cmd", returncode=0, stdout="ok", stderr=""),
        ]
        result = run_command("cmd")
        assert result.returncode == 0
        assert mock_run.call_count == 2
        mock_sleep.assert_called_once()

    @pytest.mark.unit
    @patch("plasmicheck.scripts.utils.time.sleep")
    @patch("plasmicheck.scripts.utils.subprocess.run")
    def test_exhausted_retries_raises(self, mock_run: MagicMock, mock_sleep: MagicMock) -> None:
        mock_run.side_effect = subprocess.CalledProcessError(1, "cmd", stderr="error")
        with pytest.raises(subprocess.CalledProcessError):
            run_command("cmd")
        assert mock_run.call_count == 3  # default retries from config
