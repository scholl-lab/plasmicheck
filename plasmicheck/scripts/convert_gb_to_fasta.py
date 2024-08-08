import os
from Bio import SeqIO
import json

# Resolve the path to config.json in the parent directory of the current script
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')

# Load configuration from JSON file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

DEFAULT_SHIFT_BASES = config['shift_bases']

def convert_gb_to_fasta(input_file, output_file, shift_bases=DEFAULT_SHIFT_BASES, generate_shifted=False, overwrite=False):
    if os.path.exists(output_file) and not overwrite:
        print(f"File {output_file} already exists. Skipping conversion.")
        return

    # Convert GenBank to FASTA
    SeqIO.convert(input_file, 'genbank', output_file, 'fasta')

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

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert a GenBank file to a FASTA file and optionally generate a shifted reference")
    parser.add_argument("input_file", help="Input GenBank file")
    parser.add_argument("output_file", help="Output FASTA file")
    parser.add_argument("--shift_bases", type=int, default=DEFAULT_SHIFT_BASES, help="Number of bases to shift in the shifted reference")
    parser.add_argument("--generate_shifted", action="store_true", help="Generate a shifted reference sequence")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output file")
    args = parser.parse_args()

    convert_gb_to_fasta(args.input_file, args.output_file, args.shift_bases, args.generate_shifted, args.overwrite)
