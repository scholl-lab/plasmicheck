import os
import subprocess

def create_indexes(fasta_file, overwrite=False):
    base_name = os.path.splitext(fasta_file)[0]
    
    # Create Minimap2 index
    minimap2_index = f"{base_name}.mmi"
    samtools_index = f"{fasta_file}.fai"

    if not overwrite:
        if os.path.exists(minimap2_index) or os.path.exists(samtools_index):
            raise FileExistsError("Index files already exist. Use the --overwrite flag to overwrite them.")

    # Create Minimap2 index
    subprocess.run(["minimap2", "-d", minimap2_index, fasta_file], check=True)

    # Create Samtools index
    subprocess.run(["samtools", "faidx", fasta_file], check=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Create Minimap2 and Samtools indexes for a FASTA file")
    parser.add_argument("fasta_file", help="FASTA file to index")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing index files")
    args = parser.parse_args()

    create_indexes(args.fasta_file, args.overwrite)
