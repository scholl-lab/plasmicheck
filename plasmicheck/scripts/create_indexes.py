from __future__ import annotations

import json
import logging
import os
import subprocess
from typing import Any

from .utils import setup_logging  # Import setup_logging function

# Load configuration from JSON file
with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "config.json")) as config_file:
    config: dict[str, Any] = json.load(config_file)

MINIMAP2_OPTIONS: list[str] = config["indexing"]["minimap2_options"]
SAMTOOLS_OPTIONS: list[str] = config["indexing"]["samtools_options"]
MINIMAP2_THREADS: int = config["alignment"]["minimap2_threads"]
SAMTOOLS_THREADS: int = config["alignment"]["samtools_threads"]


def create_indexes(fasta_file: str, overwrite: bool = False) -> None:
    base_name = os.path.splitext(fasta_file)[0]

    # Create Minimap2 index
    minimap2_index = f"{base_name}.mmi"
    samtools_index = f"{fasta_file}.fai"

    if not overwrite and (os.path.exists(minimap2_index) or os.path.exists(samtools_index)):
        logging.error("Index files already exist. Use the --overwrite flag to overwrite them.")
        raise FileExistsError(
            "Index files already exist. Use the --overwrite flag to overwrite them."
        )

    # Create Minimap2 index
    logging.info(f"Creating Minimap2 index for {fasta_file}")
    minimap2_command = [
        "minimap2",
        "-t",
        str(MINIMAP2_THREADS),
        *MINIMAP2_OPTIONS,
        minimap2_index,
        fasta_file,
    ]
    try:
        subprocess.run(minimap2_command, check=True)
        logging.info(f"Minimap2 index created: {minimap2_index}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to create Minimap2 index: {e}")
        raise

    # Create Samtools index
    logging.info(f"Creating Samtools index for {fasta_file}")
    samtools_command = ["samtools", "faidx", *SAMTOOLS_OPTIONS, fasta_file]
    try:
        subprocess.run(samtools_command, check=True)
        logging.info(f"Samtools index created: {samtools_index}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to create Samtools index: {e}")
        raise


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create Minimap2 and Samtools indexes for a FASTA file"
    )
    parser.add_argument("-f", "--fasta_file", help="FASTA file to index", required=True)
    parser.add_argument(
        "-w", "--overwrite", action="store_true", help="Overwrite existing index files"
    )
    parser.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser.add_argument("--log-file", help="Set the log output file", default=None)
    args = parser.parse_args()

    # Setup logging with the specified log level and file
    setup_logging(log_level=args.log_level.upper(), log_file=args.log_file)

    create_indexes(args.fasta_file, args.overwrite)
