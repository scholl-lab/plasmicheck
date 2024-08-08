import os
import shutil
from Bio import SeqIO
import json
import re

# Resolve the path to config.json in the parent directory of the current script
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')

# Load configuration from JSON file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

DEFAULT_SHIFT_BASES = config['shift_bases']
CONVERSION_CONFIG = config['conversion']

def sanitize_filename(filename):
    """Sanitize the filename by replacing any non-alphanumeric characters with underscores, including spaces."""
    return re.sub(r'[^a-zA-Z0-9_-]', '_', filename)

def create_sanitized_copy(input_file):
    """Create a sanitized copy of the input file."""
    sanitized_filename = sanitize_filename(os.path.basename(input_file))
    sanitized_filepath = os.path.join(os.path.dirname(input_file), sanitized_filename)
    shutil.copyfile(input_file, sanitized_filepath)
    return sanitized_filepath

def convert_xdna_to_gb(input_file, output_file):
    # Convert xDNA to GenBank using SeqIO
    records = SeqIO.parse(input_file, "xdna")
    count = SeqIO.write(records, output_file, "genbank")
    print(f"Converted {count} records from xDNA to GenBank")

def convert_plasmidfile_to_fasta(input_file, output_file, file_type, shift_bases=DEFAULT_SHIFT_BASES, generate_shifted=False, overwrite=False):
    # Create a sanitized copy of the input file
    sanitized_input_file = create_sanitized_copy(input_file)

    if os.path.exists(output_file) and not overwrite:
        print(f"File {output_file} already exists. Skipping conversion.")
        return

    # Determine the input file type and convert accordingly
    if file_type == 'genbank':
        SeqIO.convert(sanitized_input_file, 'genbank', output_file, 'fasta')
    elif file_type == 'xdna':
        temp_gb_file = os.path.splitext(output_file)[0] + '.gb'
        convert_xdna_to_gb(sanitized_input_file, temp_gb_file)
        SeqIO.convert(temp_gb_file, 'genbank', output_file, 'fasta')
        os.remove(temp_gb_file)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    if generate_shifted:
        # Read the converted FASTA file
        records = list(SeqIO.parse(output_file, 'fasta'))

        with open(output_file, "a") as output_handle:
            for record in records:
                # Generate the shifted sequence
                shifted_seq = record.seq[shift_bases:] + record.seq[:shift_bases]
                shifted_record = record[:]
                shifted_record.id = f"{record.id}_shifted"
                shifted_record.description = f"{record.description} (shifted by {shift_bases} bases)"
                shifted_record.seq = shifted_seq
                # Write the shifted records to the same file
                SeqIO.write(shifted_record, output_handle, "fasta")

    # Remove the sanitized copy of the input file
    os.remove(sanitized_input_file)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert a plasmid file (GenBank or xDNA) to a FASTA file and optionally generate a shifted reference")
    parser.add_argument("-i", "--input_file", help="Input plasmid file", required=True)
    parser.add_argument("-o", "--output_file", help="Output FASTA file", required=True)
    parser.add_argument("-t", "--file_type", choices=['genbank', 'xdna'], help="Type of input file: 'genbank' or 'xdna'", required=True)
    parser.add_argument("-sb", "--shift_bases", type=int, default=DEFAULT_SHIFT_BASES, help="Number of bases to shift in the shifted reference")
    parser.add_argument("-g", "--generate_shifted", action="store_true", help="Generate a shifted reference sequence")
    parser.add_argument("-w", "--overwrite", action="store_true", help="Overwrite existing output file")
    args = parser.parse_args()

    convert_plasmidfile_to_fasta(args.input_file, args.output_file, args.file_type, args.shift_bases, args.generate_shifted, args.overwrite)
