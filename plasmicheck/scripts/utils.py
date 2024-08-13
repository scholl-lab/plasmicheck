import hashlib
import os
import re

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
