import os
import subprocess
import json

# Load configuration from JSON file
with open(os.path.join(os.path.dirname(__file__), '..', 'config.json'), 'r') as config_file:
    config = json.load(config_file)

MINIMAP2_OPTIONS = config['indexing']['minimap2_options']
SAMTOOLS_OPTIONS = config['indexing']['samtools_options']
MINIMAP2_THREADS = config['alignment']['minimap2_threads']
SAMTOOLS_THREADS = config['alignment']['samtools_threads']

def create_indexes(fasta_file, overwrite=False):
    base_name = os.path.splitext(fasta_file)[0]
    
    # Create Minimap2 index
    minimap2_index = f"{base_name}.mmi"
    samtools_index = f"{fasta_file}.fai"

    if not overwrite:
        if os.path.exists(minimap2_index) or os.path.exists(samtools_index):
            raise FileExistsError("Index files already exist. Use the --overwrite flag to overwrite them.")

    # Create Minimap2 index
    minimap2_command = ["minimap2", "-t", str(MINIMAP2_THREADS)] + MINIMAP2_OPTIONS + [minimap2_index, fasta_file]
    subprocess.run(minimap2_command, check=True)

    # Create Samtools index (no threads argument for faidx)
    samtools_command = ["samtools", "faidx"] + SAMTOOLS_OPTIONS + [fasta_file]
    subprocess.run(samtools_command, check=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Create Minimap2 and Samtools indexes for a FASTA file")
    parser.add_argument("-f", "--fasta_file", help="FASTA file to index", required=True)
    parser.add_argument("-w", "--overwrite", action="store_true", help="Overwrite existing index files")
    args = parser.parse_args()

    create_indexes(args.fasta_file, args.overwrite)
