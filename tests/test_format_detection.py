"""Unit tests for plasmicheck.scripts.format_detection."""

from __future__ import annotations

from pathlib import Path

import pytest

from plasmicheck.scripts.format_detection import (
    detect_format,
    is_snapgene,
    is_xdna,
    read_magic_bytes,
)


class TestIsSnapgene:
    @pytest.mark.unit
    def test_valid_snapgene_signature(self, tmp_path: Path) -> None:
        """Valid SnapGene magic bytes should be detected."""
        test_file = tmp_path / "test.dna"
        # SnapGene signature: b'\x09\x00\x00\x00\x0ESnapGene'
        test_file.write_bytes(b"\x09\x00\x00\x00\x0eSnapGene" + b"\x00" * 100)
        magic_bytes = read_magic_bytes(str(test_file))
        assert is_snapgene(magic_bytes) is True

    @pytest.mark.unit
    def test_not_snapgene(self, tmp_path: Path) -> None:
        """Random bytes should not be detected as SnapGene."""
        test_file = tmp_path / "random.bin"
        test_file.write_bytes(b"\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f")
        magic_bytes = read_magic_bytes(str(test_file))
        assert is_snapgene(magic_bytes) is False

    @pytest.mark.unit
    def test_empty_file(self, tmp_path: Path) -> None:
        """Empty file should return False without crashing."""
        test_file = tmp_path / "empty.bin"
        test_file.write_bytes(b"")
        magic_bytes = read_magic_bytes(str(test_file))
        assert is_snapgene(magic_bytes) is False


class TestIsXdna:
    @pytest.mark.unit
    def test_valid_xdna_signature(self, tmp_path: Path) -> None:
        """Valid xDNA signature should be detected."""
        test_file = tmp_path / "test.xdna"
        # xDNA: first byte 0x00 (version), byte at index 111 is 0xFF
        header = bytearray(112)
        header[0] = 0x00  # version
        header[111] = 0xFF  # marker
        test_file.write_bytes(bytes(header))
        assert is_xdna(str(test_file)) is True

    @pytest.mark.unit
    def test_not_xdna_wrong_version(self, tmp_path: Path) -> None:
        """Wrong version byte should not be detected as xDNA."""
        test_file = tmp_path / "wrong_version.bin"
        header = bytearray(112)
        header[0] = 0x01  # wrong version
        header[111] = 0xFF
        test_file.write_bytes(bytes(header))
        assert is_xdna(str(test_file)) is False

    @pytest.mark.unit
    def test_not_xdna_missing_marker(self, tmp_path: Path) -> None:
        """Missing marker byte should not be detected as xDNA."""
        test_file = tmp_path / "missing_marker.bin"
        header = bytearray(112)
        header[0] = 0x00
        header[111] = 0x00  # wrong marker
        test_file.write_bytes(bytes(header))
        assert is_xdna(str(test_file)) is False

    @pytest.mark.unit
    def test_file_too_short(self, tmp_path: Path) -> None:
        """File shorter than 112 bytes should return False gracefully."""
        test_file = tmp_path / "short.bin"
        test_file.write_bytes(b"\x00" * 50)
        assert is_xdna(str(test_file)) is False


class TestDetectFormat:
    @pytest.mark.unit
    def test_detect_snapgene(self, tmp_path: Path) -> None:
        """SnapGene file should be detected."""
        test_file = tmp_path / "test.dna"
        test_file.write_bytes(b"\x09\x00\x00\x00\x0eSnapGene" + b"\x00" * 100)
        assert detect_format(str(test_file)) == "snapgene"

    @pytest.mark.unit
    def test_detect_xdna(self, tmp_path: Path) -> None:
        """xDNA file should be detected."""
        test_file = tmp_path / "test.xdna"
        header = bytearray(112)
        header[0] = 0x00
        header[111] = 0xFF
        test_file.write_bytes(bytes(header))
        assert detect_format(str(test_file)) == "xdna"

    @pytest.mark.unit
    def test_detect_genbank(self, tmp_path: Path) -> None:
        """GenBank file should be detected."""
        test_file = tmp_path / "test.gb"
        test_file.write_text("LOCUS       test_sequence\n")
        assert detect_format(str(test_file)) == "genbank"

    @pytest.mark.unit
    def test_detect_fasta(self, tmp_path: Path) -> None:
        """FASTA file should be detected."""
        test_file = tmp_path / "test.fasta"
        test_file.write_text(">test_sequence\nATCG\n")
        assert detect_format(str(test_file)) == "fasta"

    @pytest.mark.unit
    def test_unrecognized_raises(self, tmp_path: Path) -> None:
        """Unrecognized binary file should raise ValueError with actionable message."""
        test_file = tmp_path / "unknown.bin"
        test_file.write_bytes(b"\xff\xff\xff\xff" * 50)
        with pytest.raises(ValueError) as exc_info:
            detect_format(str(test_file))
        assert "Supported formats" in str(exc_info.value)
        assert "SnapGene" in str(exc_info.value)
        assert "xDNA" in str(exc_info.value)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "xdna_file",
        [
            "reference/plasmid/PIRES2_KCNJ5_eGFP.xdna",
            "reference/plasmid/pENTR2B_Dual_CACNA1D.xdna",
            "reference/plasmid/pENTR_2B_Dual_CACNA1H_myc_DDK.xdna",
        ],
    )
    def test_detect_real_xdna_files(self, xdna_file: str) -> None:
        """Real xDNA files should be detected correctly."""
        import os

        if not os.path.exists(xdna_file):
            pytest.skip(f"Test file {xdna_file} not found")
        assert detect_format(xdna_file) == "xdna"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "genbank_file",
        [
            "reference/plasmid/pCMV6-AC-IRES-GFP-Puro_CACNA1H.gb",
            "reference/plasmid/pCMV6_CACNA1D_genbank.gb",
            "reference/plasmid/pCMV6_CACNA1H_MycDDK.gb",
        ],
    )
    def test_detect_real_genbank_files(self, genbank_file: str) -> None:
        """Real GenBank files should be detected correctly."""
        import os

        if not os.path.exists(genbank_file):
            pytest.skip(f"Test file {genbank_file} not found")
        assert detect_format(genbank_file) == "genbank"
