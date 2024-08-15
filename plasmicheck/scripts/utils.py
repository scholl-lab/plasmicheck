import hashlib
import os
import re
import logging
import sys
import subprocess
import tarfile
import json
import time  # Added for retry delay
from datetime import datetime

# Load configuration from JSON file
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

SUPPORTED_FORMATS = config['supported_formats']

def calculate_md5(file_path):
    """Calculate and return the MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def write_md5sum(file_path, file_type, output_folder):
    """Write the MD5 checksum of a file to md5sum.txt in the output folder."""
    md5sum = calculate_md5(file_path)
    md5sum_file = os.path.join(output_folder, "md5sum.txt")
    
    with open(md5sum_file, "a") as f:
        f.write(f"{file_type}\t{file_path}\t{md5sum}\n")

def sanitize_filename(filename):
    """Sanitize a filename by replacing non-alphanumeric characters with underscores."""
    return re.sub(r'[^a-zA-Z0-9_-]', '_', filename)

def setup_logging(log_level=logging.INFO, log_file=None):
    """Setup logging configuration, avoid adding multiple handlers."""
    logger = logging.getLogger()
    if not logger.handlers:  # Avoid adding multiple handlers
        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        logger.setLevel(log_level)

def run_command(command):
    """Run a command using subprocess with retry logic and log the output."""
    # Fetch retry settings from config or use default values
    retries = config.get('retry_settings', {}).get('retries', 3)
    delay = config.get('retry_settings', {}).get('delay', 5)
    
    for attempt in range(retries):
        logging.info(f"Running command: {command} (Attempt {attempt + 1}/{retries})")
        try:
            result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
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

def archive_output_folder(output_folder, archive_name=None):
    """Archive and compress the output folder into a .tar.gz file."""
    # Extract the folder structure and create a default archive name if not provided
    if not archive_name:
        folder_parts = output_folder.strip(os.sep).split(os.sep)
        if len(folder_parts) >= 2:
            sample_name, plasmid_name = folder_parts[-2], folder_parts[-1]
        else:
            sample_name, plasmid_name = "sample", "plasmid"
        
        current_date = datetime.now().strftime("%Y-%m-%d")
        archive_name = f"{sample_name}.{plasmid_name}.{current_date}.tar.gz"

    # Define the full path to the archive file
    archive_path = os.path.join(os.path.dirname(output_folder), archive_name)

    logging.info(f"Archiving output folder {output_folder} to {archive_path}")

    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(output_folder, arcname=os.path.basename(output_folder))

    logging.info(f"Output folder archived successfully: {archive_path}")
    return archive_path

def validate_file_existence(file_paths):
    """Check if all files in the list exist."""
    missing_files = [file for file in file_paths if not os.path.exists(file)]
    if missing_files:
        logging.error(f"The following files are missing: {', '.join(missing_files)}")
        raise FileNotFoundError(f"Missing files: {', '.join(missing_files)}")
    logging.debug("All files exist.")

def validate_file_format(file_paths, expected_formats):
    """Check if files have the correct format based on their extensions."""
    if not isinstance(expected_formats, list):
        expected_formats = [expected_formats]  # Ensure it's a list for uniformity
    
    expected_extensions = []
    for format in expected_formats:
        format_extensions = SUPPORTED_FORMATS.get(format.lower())
        if not format_extensions:
            raise ValueError(f"Unsupported format: {format}")
        if isinstance(format_extensions, str):
            format_extensions = [format_extensions]
        expected_extensions.extend(format_extensions)  # Collect all valid extensions

    invalid_files = [file for file in file_paths if not any(file.endswith(ext) for ext in expected_extensions)]
    if invalid_files:
        logging.error(f"The following files have an incorrect format: {', '.join(invalid_files)}")
        raise ValueError(f"Incorrect format for files: {', '.join(invalid_files)}")
    
    logging.debug(f"All files have the correct format: {', '.join(expected_formats)}")

def quality_control(file_paths, expected_formats=None):
    """Perform quality control checks on the input files."""
    # Step 1: Validate that files exist
    validate_file_existence(file_paths)
    
    # Step 2: Validate file format, if expected formats are provided
    if expected_formats:
        validate_file_format(file_paths, expected_formats)
    
    logging.debug("Quality control checks passed.")
