"""Unit tests for plasmicheck.scripts.convert_plasmidfile_to_fasta."""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from plasmicheck.scripts.convert_plasmidfile_to_fasta import convert_plasmidfile_to_fasta
from plasmicheck.scripts.utils import sanitize_filename


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


class TestRecordIdOverride:
    @pytest.mark.unit
    def test_xdna_record_id_is_sanitized_filename(self, tmp_path: Path) -> None:
        """xDNA conversion should use sanitized filename as record ID."""
        # Find a real xDNA file
        xdna_file = "reference/plasmid/PIRES2_KCNJ5_eGFP.xdna"
        if not os.path.exists(xdna_file):
            pytest.skip(f"Test file {xdna_file} not found")

        output_file = tmp_path / "output.fasta"
        convert_plasmidfile_to_fasta(
            input_file=xdna_file,
            output_file=str(output_file),
            file_type="xdna",
            overwrite=True,
        )

        # Read output FASTA and verify header uses sanitized filename
        with open(output_file) as f:
            first_line = f.readline()
        expected_id = sanitize_filename(Path(xdna_file).stem)
        assert first_line.startswith(f">{expected_id}")
        # Should NOT contain broken xDNA author ID like [author
        assert "[author" not in first_line

    @pytest.mark.unit
    def test_genbank_record_id_is_sanitized_filename(self, tmp_path: Path) -> None:
        """GenBank conversion should use sanitized filename as record ID."""
        # Find a real GenBank file
        gb_file = "reference/plasmid/pCMV6_CACNA1D_genbank.gb"
        if not os.path.exists(gb_file):
            pytest.skip(f"Test file {gb_file} not found")

        output_file = tmp_path / "output.fasta"
        convert_plasmidfile_to_fasta(
            input_file=gb_file,
            output_file=str(output_file),
            file_type="genbank",
            overwrite=True,
        )

        # Read output FASTA and verify header uses sanitized filename
        with open(output_file) as f:
            first_line = f.readline()
        expected_id = sanitize_filename(Path(gb_file).stem)
        assert first_line.startswith(f">{expected_id}")

    @pytest.mark.unit
    def test_fasta_description_contains_provenance(self, tmp_path: Path) -> None:
        """Converted FASTA should contain provenance information."""
        xdna_file = "reference/plasmid/PIRES2_KCNJ5_eGFP.xdna"
        if not os.path.exists(xdna_file):
            pytest.skip(f"Test file {xdna_file} not found")

        output_file = tmp_path / "output.fasta"
        convert_plasmidfile_to_fasta(
            input_file=xdna_file,
            output_file=str(output_file),
            file_type="xdna",
            overwrite=True,
        )

        # Verify description contains provenance
        with open(output_file) as f:
            first_line = f.readline()
        assert "converted_from=xdna" in first_line


class TestActionableErrors:
    @pytest.mark.unit
    def test_unsupported_file_type_error(self, tmp_path: Path) -> None:
        """Unsupported file_type should raise ValueError with supported formats."""
        input_file = tmp_path / "test.txt"
        input_file.write_text("test")
        output_file = tmp_path / "output.fasta"

        with pytest.raises(ValueError) as exc_info:
            convert_plasmidfile_to_fasta(
                input_file=str(input_file),
                output_file=str(output_file),
                file_type="invalid",
                overwrite=True,
            )
        error_msg = str(exc_info.value)
        # Error message should mention supported formats
        assert "Unsupported file type" in error_msg or "Supported formats" in error_msg
        assert "snapgene" in error_msg.lower()
        assert "xdna" in error_msg.lower()
        assert "genbank" in error_msg.lower()

    @pytest.mark.unit
    def test_conversion_error_message_is_actionable(self, tmp_path: Path) -> None:
        """BioPython errors should be wrapped with actionable messages."""
        input_file = tmp_path / "broken.xdna"
        # Create minimal xDNA header (112 bytes) but with version error
        header = bytearray(112)
        header[0] = 0x00
        header[111] = 0xFF
        input_file.write_bytes(bytes(header))

        output_file = tmp_path / "output.fasta"

        # Mock SeqIO to raise ValueError with "version" in message
        with patch("plasmicheck.scripts.convert_plasmidfile_to_fasta.SeqIO") as mock_seqio:
            mock_seqio.parse.side_effect = ValueError("Unsupported xDNA version 5")

            with pytest.raises(ValueError) as exc_info:
                convert_plasmidfile_to_fasta(
                    input_file=str(input_file),
                    output_file=str(output_file),
                    file_type="xdna",
                    overwrite=True,
                )

            error_msg = str(exc_info.value)
            # Should mention exporting as GenBank
            assert "Export as GenBank" in error_msg
            # Should mention SnapGene or Serial Cloner
            assert "SnapGene" in error_msg or "Serial Cloner" in error_msg

    @pytest.mark.unit
    def test_corrupt_file_error_message(self, tmp_path: Path) -> None:
        """Corrupt file errors should have actionable messages."""
        input_file = tmp_path / "corrupt.gb"
        input_file.write_text("LOCUS test\ntruncated...")

        output_file = tmp_path / "output.fasta"

        # Mock SeqIO to raise ValueError with "corrupt" in message
        with patch("plasmicheck.scripts.convert_plasmidfile_to_fasta.SeqIO") as mock_seqio:
            mock_seqio.parse.side_effect = ValueError("File appears corrupt")

            with pytest.raises(ValueError) as exc_info:
                convert_plasmidfile_to_fasta(
                    input_file=str(input_file),
                    output_file=str(output_file),
                    file_type="genbank",
                    overwrite=True,
                )

            error_msg = str(exc_info.value)
            # Should mention re-exporting from cloning software
            assert "corrupt" in error_msg.lower() or "truncated" in error_msg.lower()
            assert "SnapGene" in error_msg or "Benchling" in error_msg
