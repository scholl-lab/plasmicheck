from __future__ import annotations

import argparse
import hashlib
import logging
import os
import re
import subprocess
import sys
import tarfile
import time  # Added for retry delay
from datetime import datetime, timezone
from typing import Any

from plasmicheck.config import get_config

_cfg = get_config()
SUPPORTED_FORMATS: dict[str, Any] = _cfg["supported_formats"]


def calculate_md5(file_path: str) -> str:
    """Calculate and return the MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def write_md5sum(file_path: str, file_type: str, output_folder: str) -> None:
    """Write the MD5 checksum of a file to md5sum.txt, skipping duplicates."""
    md5sum: str = calculate_md5(file_path)
    md5sum_file: str = os.path.join(output_folder, "md5sum.txt")
    new_line: str = f"{file_type}\t{file_path}\t{md5sum}\n"

    existing_lines: set[str] = set()
    if os.path.exists(md5sum_file):
        with open(md5sum_file) as f:
            existing_lines = set(f.readlines())

    if new_line not in existing_lines:
        with open(md5sum_file, "a") as f:
            f.write(new_line)


def sanitize_filename(filename: str) -> str:
    """Sanitize a filename by replacing non-alphanumeric characters with underscores."""
    return re.sub(r"[^a-zA-Z0-9_-]", "_", filename)


def setup_logging(log_level: int = logging.INFO, log_file: str | None = None) -> None:
    """Setup logging configuration, avoid adding multiple handlers."""
    logger: logging.Logger = logging.getLogger()
    if not logger.handlers:  # Avoid adding multiple handlers
        handler: logging.StreamHandler[Any] = logging.StreamHandler(sys.stdout)
        formatter: logging.Formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        if log_file:
            file_handler: logging.FileHandler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        logger.setLevel(log_level)


def add_logging_args(parser: argparse.ArgumentParser) -> None:
    """Add standard --log-level and --log-file arguments to a parser."""
    parser.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser.add_argument("--log-file", help="Set the log output file", default=None)


def configure_logging_from_args(args: argparse.Namespace) -> None:
    """Configure logging from parsed CLI arguments."""
    log_level: int = getattr(logging, args.log_level.upper(), logging.INFO)
    setup_logging(log_level=log_level, log_file=args.log_file)


def run_command(command: str) -> subprocess.CompletedProcess[str]:
    """Run a command using subprocess with retry logic and log the output."""
    # Fetch retry settings from config or use default values
    retries: int = _cfg.get("retry_settings", {}).get("retries", 3)
    delay: int = _cfg.get("retry_settings", {}).get("delay", 5)

    for attempt in range(retries):
        logging.info(f"Running command: {command} (Attempt {attempt + 1}/{retries})")
        try:
            result: subprocess.CompletedProcess[str] = subprocess.run(
                command,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
            )
            logging.debug(result.stdout)
            if result.stderr:
                logging.warning(result.stderr)
            return result  # Exit after successful run
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}: {e.stderr}")
            if attempt < retries - 1:
                logging.info(f"Retrying command in {delay} seconds...")
                time.sleep(delay)
            else:
                logging.error(f"Command failed after {retries} attempts.")
                raise
    msg: str = f"Command failed after {retries} attempts: {command}"
    raise RuntimeError(msg)


def archive_output_folder(output_folder: str, archive_name: str | None = None) -> str:
    """Archive and compress the output folder into a .tar.gz file."""
    # Extract the folder structure and create a default archive name if not provided
    if not archive_name:
        folder_parts: list[str] = output_folder.strip(os.sep).split(os.sep)
        if len(folder_parts) >= 2:
            sample_name, plasmid_name = folder_parts[-2], folder_parts[-1]
        else:
            sample_name, plasmid_name = "sample", "plasmid"

        current_date: str = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d")
        archive_name = f"{sample_name}.{plasmid_name}.{current_date}.tar.gz"

    # Define the full path to the archive file
    archive_path: str = os.path.join(os.path.dirname(output_folder), archive_name)

    logging.info(f"Archiving output folder {output_folder} to {archive_path}")

    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(output_folder, arcname=os.path.basename(output_folder))

    logging.info(f"Output folder archived successfully: {archive_path}")
    return archive_path


def validate_file_existence(file_paths: list[str]) -> None:
    """Check if all files in the list exist."""
    missing_files: list[str] = [file for file in file_paths if not os.path.exists(file)]
    if missing_files:
        logging.error(f"The following files are missing: {', '.join(missing_files)}")
        raise FileNotFoundError(f"Missing files: {', '.join(missing_files)}")
    logging.debug("All files exist.")


def validate_file_format(file_paths: list[str], expected_formats: list[str] | str) -> None:
    """Check if files have the correct format based on their extensions."""
    if not isinstance(expected_formats, list):
        expected_formats = [expected_formats]  # Ensure it's a list for uniformity

    expected_extensions: list[str] = []
    for fmt in expected_formats:
        fmt_extensions: Any = SUPPORTED_FORMATS.get(fmt.lower())
        if not fmt_extensions:
            raise ValueError(f"Unsupported format: {fmt}")
        if isinstance(fmt_extensions, str):
            fmt_extensions = [fmt_extensions]
        expected_extensions.extend(fmt_extensions)  # Collect all valid extensions

    invalid_files: list[str] = [
        file for file in file_paths if not any(file.endswith(ext) for ext in expected_extensions)
    ]
    if invalid_files:
        logging.error(f"The following files have an incorrect format: {', '.join(invalid_files)}")
        raise ValueError(f"Incorrect format for files: {', '.join(invalid_files)}")

    logging.debug(f"All files have the correct format: {', '.join(expected_formats)}")


def quality_control(file_paths: list[str], expected_formats: list[str] | str | None = None) -> None:
    """Perform quality control checks on the input files."""
    # Step 1: Validate that files exist
    validate_file_existence(file_paths)

    # Step 2: Validate file format, if expected formats are provided
    if expected_formats:
        validate_file_format(file_paths, expected_formats)

    logging.debug("Quality control checks passed.")
