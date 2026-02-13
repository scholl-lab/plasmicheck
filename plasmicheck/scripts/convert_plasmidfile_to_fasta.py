from __future__ import annotations

import json
import logging
import os
import re
import shutil
from typing import Any

from Bio import SeqIO

from .utils import setup_logging  # Import setup_logging function

# Resolve the path to config.json in the parent directory of the current script
config_path: str = os.path.join(os.path.dirname(os.path.dirname(__file__)), "config.json")

# Load configuration from JSON file
with open(config_path) as config_file:
    config: dict[str, Any] = json.load(config_file)

DEFAULT_SHIFT_BASES: int = config["shift_bases"]
CONVERSION_CONFIG: dict[str, Any] = config["conversion"]


def sanitize_filename(filename: str) -> str:
    """Sanitize the filename by replacing any non-alphanumeric characters with underscores, including spaces."""
    return re.sub(r"[^a-zA-Z0-9_-]", "_", filename)


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
    records: Any = SeqIO.parse(input_file, "xdna")
    count: int = SeqIO.write(records, output_file, "genbank")
    logging.info(f"Converted {count} records from xDNA to GenBank")


def convert_plasmidfile_to_fasta(
    input_file: str,
    output_file: str,
    file_type: str,
    shift_bases: int = DEFAULT_SHIFT_BASES,
    generate_shifted: bool = False,
    overwrite: bool = False,
) -> None:
    """Convert a plasmid file to a FASTA file, optionally generating a shifted reference."""
    logging.info(f"Starting conversion of {input_file} to {output_file} as {file_type} file")

    # Create a sanitized copy of the input file
    sanitized_input_file: str = create_sanitized_copy(input_file)

    if os.path.exists(output_file) and not overwrite:
        logging.warning(f"File {output_file} already exists. Skipping conversion.")
        return

    # Determine the input file type and convert accordingly
    try:
        if file_type == "genbank":
            SeqIO.convert(sanitized_input_file, "genbank", output_file, "fasta")
        elif file_type == "xdna":
            temp_gb_file: str = os.path.splitext(output_file)[0] + ".gb"
            convert_xdna_to_gb(sanitized_input_file, temp_gb_file)
            SeqIO.convert(temp_gb_file, "genbank", output_file, "fasta")
            os.remove(temp_gb_file)
        else:
            raise ValueError(f"Unsupported file type: {file_type}")
        logging.info(f"Conversion to FASTA completed successfully for {output_file}")
    except Exception as e:
        logging.error(f"Error during conversion: {e}")
        raise

    if generate_shifted:
        logging.info(f"Generating shifted reference for {output_file}")
        try:
            # Read the converted FASTA file
            records: list[Any] = list(SeqIO.parse(output_file, "fasta"))

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
        choices=["genbank", "xdna"],
        help="Type of input file: 'genbank' or 'xdna'",
        required=True,
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
    parser.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser.add_argument("--log-file", help="Set the log output file", default=None)
    args = parser.parse_args()

    # Setup logging with the specified log level and file
    setup_logging(log_level=args.log_level.upper(), log_file=args.log_file)

    convert_plasmidfile_to_fasta(
        args.input_file,
        args.output_file,
        args.file_type,
        args.shift_bases,
        args.generate_shifted,
        args.overwrite,
    )
