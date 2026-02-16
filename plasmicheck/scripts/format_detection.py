from __future__ import annotations

import logging
from pathlib import Path


def read_magic_bytes(file_path: str, n_bytes: int = 16) -> bytes:
    """Read the first N bytes from a file.

    Args:
        file_path: Path to file to read
        n_bytes: Number of bytes to read (default: 16)

    Returns:
        First n_bytes of the file as bytes
    """
    with open(file_path, "rb") as f:
        return f.read(n_bytes)


def is_snapgene(magic_bytes: bytes) -> bool:
    """Check if file is SnapGene format based on magic bytes.

    SnapGene .dna files start with: b'\\x09\\x00\\x00\\x00\\x0ESnapGene'

    Args:
        magic_bytes: First bytes of file (at least 13 bytes)

    Returns:
        True if SnapGene format detected
    """
    snapgene_signature = b"\x09\x00\x00\x00\x0eSnapGene"
    return magic_bytes[:13] == snapgene_signature


def is_xdna(file_path: str) -> bool:
    """Check if file is xDNA format.

    xDNA files have:
    - First byte is version (typically 0x00)
    - Byte at index 111 is 0xFF

    Args:
        file_path: Path to file to check

    Returns:
        True if xDNA format detected
    """
    with open(file_path, "rb") as f:
        header = f.read(112)

    if len(header) < 112:
        return False

    # First byte is version (typically 0x00)
    # Byte at index 111 (112th byte) should be 0xFF
    return header[0] == 0x00 and header[111] == 0xFF


def detect_format(file_path: str) -> str:
    """Detect file format from magic bytes and content.

    Args:
        file_path: Path to file to analyze

    Returns:
        Format name: "snapgene", "xdna", "genbank", or "fasta"

    Raises:
        ValueError: If format cannot be determined
    """
    # Read magic bytes for binary format detection
    magic_bytes = read_magic_bytes(file_path, n_bytes=16)

    # Check SnapGene first (binary format)
    if is_snapgene(magic_bytes):
        logging.info("Detected format: snapgene")
        return "snapgene"

    # Check xDNA (binary format, needs 112 bytes)
    if is_xdna(file_path):
        logging.info("Detected format: xdna")
        return "xdna"

    # Try text-based formats (GenBank, FASTA)
    try:
        with open(file_path, encoding="utf-8") as f:
            first_line = f.readline().strip()

        if first_line.startswith("LOCUS"):
            logging.info("Detected format: genbank")
            return "genbank"
        elif first_line.startswith(">"):
            logging.info("Detected format: fasta")
            return "fasta"
    except UnicodeDecodeError:
        # Binary file that's not SnapGene or xDNA
        pass

    # Format not recognized
    filename = Path(file_path).name
    raise ValueError(
        f"Unrecognized file format: {filename}. "
        f"Supported formats: SnapGene (.dna), xDNA (.xdna), GenBank (.gb/.gbk), FASTA (.fasta/.fa). "
        f"If unsure, export as GenBank from your cloning software (SnapGene, Serial Cloner, Benchling)."
    )
