from __future__ import annotations

import logging
import os
import subprocess

from filelock import FileLock

from plasmicheck.config import get_config

_cfg = get_config()
MINIMAP2_OPTIONS: list[str] = _cfg["indexing"]["minimap2_options"]
SAMTOOLS_OPTIONS: list[str] = _cfg["indexing"]["samtools_options"]
MINIMAP2_THREADS: int = _cfg["alignment"]["minimap2_threads"]
SAMTOOLS_THREADS: int = _cfg["alignment"]["samtools_threads"]
THREADS_MINIMAP2: int = MINIMAP2_THREADS
THREADS_SAMTOOLS: int = SAMTOOLS_THREADS
LOCK_TIMEOUT: int = 600


def create_indexes(fasta_file: str, overwrite: bool = False, threads: int | None = None) -> None:
    base_name = os.path.splitext(fasta_file)[0]
    minimap2_index = f"{base_name}.mmi"
    samtools_index = f"{fasta_file}.fai"
    num_threads = threads if threads is not None else THREADS_MINIMAP2

    # Fast-path: if indexes already exist and we are not overwriting, avoid taking the lock.
    if not overwrite and os.path.exists(minimap2_index) and os.path.exists(samtools_index):
        logging.info(f"Index files already exist for {fasta_file}, skipping.")
        return

    lock_path = f"{base_name}.lock"
    lock = FileLock(lock_path, timeout=LOCK_TIMEOUT)
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

        # Create Minimap2 index atomically via temporary file
        logging.info(f"Creating Minimap2 index for {fasta_file}")
        tmp_mmi = f"{minimap2_index}.tmp.{os.getpid()}"
        minimap2_command = [
            "minimap2",
            "-t",
            str(num_threads),
            *MINIMAP2_OPTIONS,
            tmp_mmi,
            fasta_file,
        ]
        try:
            subprocess.run(minimap2_command, check=True)
            if os.path.exists(tmp_mmi):
                os.replace(tmp_mmi, minimap2_index)
            logging.info(f"Minimap2 index created: {minimap2_index}")
        except subprocess.CalledProcessError as e:
            if os.path.exists(tmp_mmi):
                os.remove(tmp_mmi)
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
