from __future__ import annotations

import logging
import os
import shutil
import struct
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plasmicheck.config import get_config

from .utils import sanitize_filename

_cfg = get_config()
DEFAULT_SHIFT_BASES: int = _cfg["shift_bases"]
CONVERSION_CONFIG: dict[str, Any] = _cfg["conversion"]


def convert_plasmidfile_to_fasta(
    input_file: str,
    output_file: str,
    file_type: str,
    shift_bases: int = DEFAULT_SHIFT_BASES,
    generate_shifted: bool = False,
    overwrite: bool = False,
    preserve_genbank: bool = False,
) -> None:
    """Convert a plasmid file to a FASTA file, optionally generating a shifted reference.

    Args:
        preserve_genbank: If True, preserve the GenBank file in the output directory.
            For xDNA files, retains the intermediate .gb file.
            For GenBank files, creates a copy in the output directory.
    """
    logging.info(f"Starting conversion of {input_file} to {output_file} as {file_type} file")

    if os.path.exists(output_file) and not overwrite:
        logging.warning(f"File {output_file} already exists. Skipping conversion.")
        return

    # Determine the sanitized record ID from filename
    sanitized_record_id = sanitize_filename(Path(input_file).stem)

    # Determine the input file type and convert accordingly.
    try:
        if file_type == "genbank":
            # Parse records, override ID, then write
            gb_records: Any = list(SeqIO.parse(input_file, "genbank"))  # type: ignore[no-untyped-call]
            for record in gb_records:
                record.id = sanitized_record_id
                record.name = sanitized_record_id
                record.description = "converted_from=genbank"
            SeqIO.write(gb_records, output_file, "fasta")
            # Preserve GenBank file if requested
            if preserve_genbank:
                gb_output_file: str = os.path.splitext(output_file)[0] + ".gb"
                shutil.copyfile(input_file, gb_output_file)
                logging.info(f"Preserved GenBank file: {gb_output_file}")
        elif file_type == "xdna":
            try:
                xdna_records: Any = list(SeqIO.parse(input_file, "xdna"))  # type: ignore[no-untyped-call]
                for record in xdna_records:
                    record.id = sanitized_record_id
                    record.name = sanitized_record_id
                    record.description = "converted_from=xdna"
                SeqIO.write(xdna_records, output_file, "fasta")
                if preserve_genbank:
                    gb_output_file = os.path.splitext(output_file)[0] + ".gb"
                    SeqIO.write(xdna_records, gb_output_file, "genbank")
                    logging.info(f"Preserved GenBank file: {gb_output_file}")
            except ValueError as ve:
                ve_msg = str(ve).lower()
                if "version" in ve_msg or "corrupt" in ve_msg or "truncated" in ve_msg:
                    raise
                # Fallback: origin-wrapping features cause Biopython XdnaIO ValueError.
                # Extract sequence directly from binary xDNA file (no feature annotations available).
                with open(input_file, "rb") as f:
                    header = f.read(112)
                    if len(header) < 112:
                        raise
                    (version, _st, _top, length, _neg, _com) = struct.unpack(
                        ">BBB25xII60xI12x", header
                    )
                    if version != 0:
                        raise ValueError(f"Unsupported xDNA version {version}") from ve
                    f.seek(112)
                    seq_str = f.read(length).decode("ascii", errors="ignore")
                fallback_record = SeqRecord(Seq(seq_str), id=sanitized_record_id, description="")
                SeqIO.write([fallback_record], output_file, "fasta")
                gb_output_file = os.path.splitext(output_file)[0] + ".gb"
                if os.path.exists(gb_output_file):
                    os.remove(gb_output_file)
        elif file_type == "snapgene":
            # Parse SnapGene file
            sg_records: Any = list(SeqIO.parse(input_file, "snapgene"))  # type: ignore[no-untyped-call]
            for record in sg_records:
                record.id = sanitized_record_id
                record.name = sanitized_record_id
                record.description = "converted_from=snapgene"
            SeqIO.write(sg_records, output_file, "fasta")
            # Preserve as GenBank if requested
            if preserve_genbank:
                gb_output_file = os.path.splitext(output_file)[0] + ".gb"
                # Convert snapgene -> genbank
                SeqIO.write(sg_records, gb_output_file, "genbank")
                logging.info(f"Preserved GenBank file: {gb_output_file}")
        elif file_type == "fasta":
            # Already FASTA — re-write with sanitized record IDs
            fa_records: Any = list(SeqIO.parse(input_file, "fasta"))  # type: ignore[no-untyped-call]
            for record in fa_records:
                record.id = sanitized_record_id
                record.name = sanitized_record_id
                record.description = "converted_from=fasta"
            SeqIO.write(fa_records, output_file, "fasta")
        else:
            raise ValueError(
                f"Unsupported file type: {file_type}. Supported: snapgene, xdna, genbank, fasta"
            )
        logging.info(f"Conversion to FASTA completed successfully for {output_file}")
    except Exception as e:
        error_msg = str(e).lower()
        # Provide actionable error messages based on failure type
        if "version" in error_msg:
            raise ValueError(
                "Unsupported xDNA version. Export as GenBank (.gb) from SnapGene or Serial Cloner, then use plasmicheck convert."
            ) from e
        elif "corrupt" in error_msg or "truncated" in error_msg:
            raise ValueError(
                "File appears corrupt or truncated. Re-export from your cloning software (SnapGene, Serial Cloner, Benchling)."
            ) from e
        else:
            raise ValueError(
                f"Failed to parse {input_file}. Supported formats: SnapGene (.dna), xDNA (.xdna), GenBank (.gb/.gbk), FASTA (.fasta/.fa). "
                f"If unsure, export as GenBank from your cloning software."
            ) from e

    if generate_shifted:
        logging.warning(
            "--generate_shifted is deprecated and ignored; minimap2 handles origin-spanning reads natively."
        )
