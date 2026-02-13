from __future__ import annotations

import logging
import os
import subprocess

from filelock import FileLock

from plasmicheck.config import get_config

from .utils import add_logging_args, configure_logging_from_args

_cfg = get_config()
MINIMAP2_OPTIONS: list[str] = _cfg["indexing"]["minimap2_options"]
SAMTOOLS_OPTIONS: list[str] = _cfg["indexing"]["samtools_options"]
MINIMAP2_THREADS: int = _cfg["alignment"]["minimap2_threads"]
SAMTOOLS_THREADS: int = _cfg["alignment"]["samtools_threads"]


def create_indexes(fasta_file: str, overwrite: bool = False) -> None:
    base_name = os.path.splitext(fasta_file)[0]
    minimap2_index = f"{base_name}.mmi"
    samtools_index = f"{fasta_file}.fai"

    # Fast-path: if indexes already exist and we are not overwriting, avoid taking the lock.
    if not overwrite and os.path.exists(minimap2_index) and os.path.exists(samtools_index):
        logging.info(f"Index files already exist for {fasta_file}, skipping.")
        return

    lock_path = f"{base_name}.lock"
    lock = FileLock(lock_path, timeout=600)
    with lock:
        # Re-check inside lock — another process may have finished while we waited
        if not overwrite and os.path.exists(minimap2_index) and os.path.exists(samtools_index):
            logging.info(f"Index files already exist for {fasta_file}, skipping.")
            return

        # Warn about partial index state (one exists, the other doesn't)
        has_minimap2 = os.path.exists(minimap2_index)
        has_samtools = os.path.exists(samtools_index)
        if not overwrite and (has_minimap2 != has_samtools):
            present = minimap2_index if has_minimap2 else samtools_index
            missing = samtools_index if has_minimap2 else minimap2_index
            logging.warning(
                f"Partial index state: {present} exists but {missing} is missing. "
                "Recreating both indexes."
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
    add_logging_args(parser)
    args = parser.parse_args()

    configure_logging_from_args(args)

    create_indexes(args.fasta_file, args.overwrite)
