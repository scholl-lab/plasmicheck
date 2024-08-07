import subprocess
import pysam
import os
from Bio import SeqIO

def find_fasta_file(base_name):
    extensions = ['.fasta', '.fa', '.fna', '.fsa', '.ffn']
    for ext in extensions:
        fasta_file = base_name + ext
        if os.path.isfile(fasta_file):
            return fasta_file
    raise FileNotFoundError(f"No FASTA file found for base name {base_name} with extensions {extensions}")

def spliced_alignment(human_index, plasmid_fasta, output_bam, padding=1000):
    # Perform spliced alignment
    subprocess.run(
        f"minimap2 -ax splice {human_index} {plasmid_fasta} | samtools view -h -F 4 - | samtools sort -o {output_bam}",
        shell=True,
        check=True
    )

    # Index the BAM file
    subprocess.run(f"samtools index {output_bam}", shell=True, check=True)

    # Extract spanned human reference regions with padding
    bamfile = pysam.AlignmentFile(output_bam, "rb")
    spanned_regions = []

    for read in bamfile.fetch():
        if read.reference_start is None or read.reference_end is None:
            continue
        start = max(0, read.reference_start - padding)
        end = read.reference_end + padding
        spanned_regions.append((read.reference_name, start, end))

    bamfile.close()

    # Deduplicate and sort the spanned regions
    spanned_regions = sorted(set(spanned_regions))

    return spanned_regions

def extract_human_reference(human_fasta, spanned_regions, output_fasta):
    with open(human_fasta, "r") as human_handle, open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(human_handle, "fasta"):
            for region in spanned_regions:
                if record.id == region[0]:
                    SeqIO.write(record[region[1]:region[2]], output_handle, "fasta")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Perform spliced alignment and extract human reference regions")
    parser.add_argument("output_fasta", help="Output FASTA file for the extracted human reference regions")
    parser.add_argument("human_index", help="Minimap2 index for the human reference genome")
    parser.add_argument("plasmid_fasta", help="FASTA file of the plasmid reference")
    parser.add_argument("output_bam", help="Output BAM file for the spliced alignment")
    parser.add_argument("--human_fasta", help="FASTA file of the human reference genome", default=None)
    parser.add_argument("--padding", type=int, default=1000, help="Padding to add to both sides of the spanned regions (default: 1000)")

    args = parser.parse_args()

    if args.human_fasta is None:
        base_name = os.path.splitext(args.human_index)[0]
        args.human_fasta = find_fasta_file(base_name)

    spanned_regions = spliced_alignment(args.human_index, args.plasmid_fasta, args.output_bam, args.padding)
    extract_human_reference(args.human_fasta, spanned_regions, args.output_fasta)
