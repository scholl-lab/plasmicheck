import hashlib
import os
import re
import logging
import sys
import subprocess
import tarfile
from datetime import datetime

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
    """Run a command using subprocess and log the output."""
    logging.info(f"Running command: {command}")
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.debug(result.stdout)
        logging.error(result.stderr)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}: {e.stderr}")
        raise

def archive_output_folder(output_folder, archive_name=None):
    """
    Archive and compress the output folder into a .tar.gz file.

    Args:
        output_folder (str): The path to the output folder to be archived.
        archive_name (str, optional): The name of the archive file. If None, a name will be generated based on the
                                      folder structure and current date.
    """
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
