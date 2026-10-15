"""Unit tests for plasmicheck.scripts.convert_plasmidfile_to_fasta."""

from __future__ import annotations

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


def _create_synthetic_xdna(
    path: Path, seq: str = "ATCG" * 50, record_id: str = "test_plasmid"
) -> Path:
    from Bio.Seq import Seq
    from Bio.SeqIO.XdnaIO import XdnaWriter
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(Seq(seq), id=record_id, description="synthetic plasmid")
    rec.annotations["molecule_type"] = "DNA"
    with open(path, "wb") as f:
        XdnaWriter(f).write_file([rec])
    return path


class TestRecordIdOverride:
    @pytest.mark.unit
    def test_xdna_record_id_is_sanitized_filename(self, tmp_path: Path) -> None:
        """xDNA conversion should use sanitized filename as record ID."""
        xdna_file = _create_synthetic_xdna(tmp_path / "My Test Plasmid.xdna")
        output_file = tmp_path / "output.fasta"
        convert_plasmidfile_to_fasta(
            input_file=str(xdna_file),
            output_file=str(output_file),
            file_type="xdna",
            overwrite=True,
        )

        with open(output_file) as f:
            first_line = f.readline()
        expected_id = sanitize_filename(xdna_file.stem)
        assert first_line.startswith(f">{expected_id}")
        assert "[author" not in first_line

    @pytest.mark.unit
    def test_genbank_record_id_is_sanitized_filename(self, tmp_path: Path) -> None:
        """GenBank conversion should use sanitized filename as record ID."""
        gb_file = Path("tests/data/synthetic/plasmid.gb")
        output_file = tmp_path / "output.fasta"
        convert_plasmidfile_to_fasta(
            input_file=str(gb_file),
            output_file=str(output_file),
            file_type="genbank",
            overwrite=True,
        )

        with open(output_file) as f:
            first_line = f.readline()
        expected_id = sanitize_filename(gb_file.stem)
        assert first_line.startswith(f">{expected_id}")

    @pytest.mark.unit
    def test_fasta_description_contains_provenance(self, tmp_path: Path) -> None:
        """Converted FASTA should contain provenance information."""
        xdna_file = _create_synthetic_xdna(tmp_path / "synthetic.xdna")
        output_file = tmp_path / "output.fasta"
        convert_plasmidfile_to_fasta(
            input_file=str(xdna_file),
            output_file=str(output_file),
            file_type="xdna",
            overwrite=True,
        )

        with open(output_file) as f:
            first_line = f.readline()
        assert "converted_from=xdna" in first_line

    @pytest.mark.unit
    def test_xdna_origin_wrapping_fallback(self, tmp_path: Path) -> None:
        """Origin-wrapping xDNA file falls back to binary sequence extraction."""
        import struct

        xdna_path = tmp_path / "origin_wrap.xdna"
        seq = b"ATCGATCG" * 20
        header = bytearray(112)
        packed = struct.pack(">BBB25xII60xI12x", 0, 1, 1, len(seq), 0, 0)
        header[: len(packed)] = packed
        xdna_path.write_bytes(header + seq)

        output_file = tmp_path / "output.fasta"
        with patch(
            "Bio.SeqIO.parse",
            side_effect=ValueError("End location (4) must be greater than start (10580)"),
        ):
            convert_plasmidfile_to_fasta(
                input_file=str(xdna_path),
                output_file=str(output_file),
                file_type="xdna",
                overwrite=True,
            )

        with open(output_file) as f:
            lines = f.readlines()
        assert lines[0].startswith(f">{xdna_path.stem}")
        assert "".join(line.strip() for line in lines[1:]) == seq.decode("ascii")

    @pytest.mark.unit
    def test_generate_shifted_deprecation_noop(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """generate_shifted logs a deprecation warning and does not write shifted reference."""
        import logging

        gb_file = Path("tests/data/synthetic/plasmid.gb")
        output_file = tmp_path / "output.fasta"

        with caplog.at_level(logging.WARNING):
            convert_plasmidfile_to_fasta(
                input_file=str(gb_file),
                output_file=str(output_file),
                file_type="genbank",
                generate_shifted=True,
                overwrite=True,
            )

        assert "--generate_shifted is deprecated and ignored" in caplog.text
        content = output_file.read_text()
        assert "_shifted" not in content


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
