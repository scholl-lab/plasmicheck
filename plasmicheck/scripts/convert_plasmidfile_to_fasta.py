from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path
from typing import Any

from Bio import SeqIO

from plasmicheck.config import get_config

from .utils import add_logging_args, configure_logging_from_args, sanitize_filename

_cfg = get_config()
DEFAULT_SHIFT_BASES: int = _cfg["shift_bases"]
CONVERSION_CONFIG: dict[str, Any] = _cfg["conversion"]


def create_sanitized_copy(input_file: str) -> str:
    """Create a sanitized copy of the input file."""
    sanitized_filename: str = sanitize_filename(os.path.basename(input_file))
    sanitized_filepath: str = os.path.join(os.path.dirname(input_file), sanitized_filename)
    shutil.copyfile(input_file, sanitized_filepath)
    logging.debug(f"Created sanitized copy: {sanitized_filepath}")
    return sanitized_filepath


def convert_xdna_to_gb(input_file: str, output_file: str) -> None:
    """Convert xDNA to GenBank using SeqIO."""
    logging.info(f"Converting xDNA file {input_file} to GenBank format")
    records: Any = SeqIO.parse(input_file, "xdna")  # type: ignore[no-untyped-call]
    count: int = SeqIO.write(records, output_file, "genbank")
    logging.info(f"Converted {count} records from xDNA to GenBank")


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
            For GenBank files, creates a sanitized copy in the output directory.
    """
    logging.info(f"Starting conversion of {input_file} to {output_file} as {file_type} file")

    # Create a sanitized copy of the input file
    sanitized_input_file: str = create_sanitized_copy(input_file)

    if os.path.exists(output_file) and not overwrite:
        logging.warning(f"File {output_file} already exists. Skipping conversion.")
        return

    # Determine the sanitized record ID from filename
    sanitized_record_id = sanitize_filename(Path(input_file).stem)

    # Determine the input file type and convert accordingly
    try:
        if file_type == "genbank":
            # Parse records, override ID, then write
            gb_records: Any = list(SeqIO.parse(sanitized_input_file, "genbank"))  # type: ignore[no-untyped-call]
            for record in gb_records:
                record.id = sanitized_record_id
                record.name = sanitized_record_id
                record.description = "converted_from=genbank"
            SeqIO.write(gb_records, output_file, "fasta")
            # Preserve GenBank file if requested
            if preserve_genbank:
                gb_output_file: str = os.path.splitext(output_file)[0] + ".gb"
                shutil.copyfile(sanitized_input_file, gb_output_file)
                logging.info(f"Preserved GenBank file: {gb_output_file}")
        elif file_type == "xdna":
            temp_gb_file: str = os.path.splitext(output_file)[0] + ".gb"
            convert_xdna_to_gb(sanitized_input_file, temp_gb_file)
            # Parse records, override ID, then write
            xdna_records: Any = list(SeqIO.parse(temp_gb_file, "genbank"))  # type: ignore[no-untyped-call]
            for record in xdna_records:
                record.id = sanitized_record_id
                record.name = sanitized_record_id
                record.description = "converted_from=xdna"
            SeqIO.write(xdna_records, output_file, "fasta")
            # Only remove temp file if not preserving GenBank
            if not preserve_genbank:
                os.remove(temp_gb_file)
            else:
                logging.info(f"Preserved GenBank file: {temp_gb_file}")
        elif file_type == "snapgene":
            # Parse SnapGene file
            sg_records: Any = list(SeqIO.parse(sanitized_input_file, "snapgene"))  # type: ignore[no-untyped-call]
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
        logging.info(f"Generating shifted reference for {output_file}")
        try:
            # Read the converted FASTA file
            records: list[Any] = list(SeqIO.parse(output_file, "fasta"))  # type: ignore[no-untyped-call]

            with open(output_file, "a") as output_handle:
                for record in records:
                    # Generate the shifted sequence
                    shifted_seq: Any = record.seq[shift_bases:] + record.seq[:shift_bases]
                    shifted_record: Any = record[:]
                    shifted_record.id = f"{record.id}_shifted"
                    shifted_record.description = (
                        f"{record.description} (shifted by {shift_bases} bases)"
                    )
                    shifted_record.seq = shifted_seq
                    # Write the shifted records to the same file
                    SeqIO.write(shifted_record, output_handle, "fasta")
            logging.info(f"Shifted reference added to {output_file}")
        except Exception as e:
            logging.error(f"Error generating shifted reference: {e}")
            raise

    # Remove the sanitized copy of the input file
    os.remove(sanitized_input_file)
    logging.debug(f"Removed sanitized copy: {sanitized_input_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert a plasmid file (GenBank or xDNA) to a FASTA file and optionally generate a shifted reference"
    )
    parser.add_argument("-i", "--input_file", help="Input plasmid file", required=True)
    parser.add_argument("-o", "--output_file", help="Output FASTA file", required=True)
    parser.add_argument(
        "-t",
        "--file_type",
        choices=["genbank", "xdna", "snapgene"],
        help="Type of input file: 'genbank', 'xdna', or 'snapgene'",
        required=False,
    )
    parser.add_argument(
        "-sb",
        "--shift_bases",
        type=int,
        default=DEFAULT_SHIFT_BASES,
        help="Number of bases to shift in the shifted reference",
    )
    parser.add_argument(
        "-g",
        "--generate_shifted",
        action="store_true",
        help="Generate a shifted reference sequence",
    )
    parser.add_argument(
        "-w", "--overwrite", action="store_true", help="Overwrite existing output file"
    )
    parser.add_argument(
        "-pg",
        "--preserve_genbank",
        action="store_true",
        help="Preserve GenBank file in output directory for gene detection",
    )
    add_logging_args(parser)
    args = parser.parse_args()

    configure_logging_from_args(args)

    # Auto-detect format if not specified
    file_type = args.file_type
    if file_type is None:
        from .format_detection import detect_format

        file_type = detect_format(args.input_file)
        logging.info(f"Auto-detected format: {file_type}")

    convert_plasmidfile_to_fasta(
        args.input_file,
        args.output_file,
        file_type,
        args.shift_bases,
        args.generate_shifted,
        args.overwrite,
        args.preserve_genbank,
    )
