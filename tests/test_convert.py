"""Unit tests for plasmicheck.scripts.convert_plasmidfile_to_fasta."""

from __future__ import annotations

import pytest

from plasmicheck.scripts.convert_plasmidfile_to_fasta import sanitize_filename


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
