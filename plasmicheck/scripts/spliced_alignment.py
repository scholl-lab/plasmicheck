import os
import pysam
from Bio import SeqIO
import json
import logging
from .utils import setup_logging, run_command  # Import setup_logging and run_command functions

# Load configuration from JSON file
with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json'), 'r') as config_file:
    config = json.load(config_file)

PADDING_DEFAULT = config['padding']
FASTA_EXTENSIONS = config['alignment']['fasta_extensions']
MINIMAP2_THREADS = config['alignment']['minimap2_threads']
SAMTOOLS_THREADS = config['alignment']['samtools_threads']

def find_fasta_file(base_name):
    logging.info(f"Finding FASTA file with base name: {base_name}")
    for ext in FASTA_EXTENSIONS:
        fasta_file = base_name + ext
        if os.path.isfile(fasta_file):
            logging.debug(f"Found FASTA file: {fasta_file}")
            return fasta_file
    logging.error(f"No FASTA file found for base name {base_name} with extensions {FASTA_EXTENSIONS}")
    raise FileNotFoundError(f"No FASTA file found for base name {base_name} with extensions {FASTA_EXTENSIONS}")

def spliced_alignment(human_index, plasmid_fasta, output_bam, padding=PADDING_DEFAULT):
    logging.info(f"Performing spliced alignment with {plasmid_fasta} to {human_index} and writing to {output_bam}")
    run_command(f"minimap2 -t {MINIMAP2_THREADS} -ax splice {human_index} {plasmid_fasta} | samtools view -@ {SAMTOOLS_THREADS} -h -F 4 - | samtools sort -@ {SAMTOOLS_THREADS} -o {output_bam}")

    logging.info(f"Indexing BAM file: {output_bam}")
    run_command(f"samtools index {output_bam}")

    logging.info(f"Extracting spanned human reference regions with padding: {padding}")
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
    logging.debug(f"Spanned regions: {spanned_regions}")

    return spanned_regions

def extract_human_reference(human_fasta, spanned_regions, output_fasta):
    logging.info(f"Extracting human reference regions from {human_fasta} to {output_fasta}")
    with open(human_fasta, "r") as human_handle, open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(human_handle, "fasta"):
            for region in spanned_regions:
                if record.id == region[0]:
                    logging.debug(f"Writing region {region} to {output_fasta}")
                    SeqIO.write(record[region[1]:region[2]], output_handle, "fasta")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Perform spliced alignment and extract human reference regions")
    parser.add_argument("-o", "--output_fasta", help="Output FASTA file for the extracted human reference regions", required=True)
    parser.add_argument("-i", "--human_index", help="Minimap2 index for the human reference genome", required=True)
    parser.add_argument("-p", "--plasmid_fasta", help="FASTA file of the plasmid reference", required=True)
    parser.add_argument("-b", "--output_bam", help="Output BAM file for the spliced alignment", required=True)
    parser.add_argument("-hf", "--human_fasta", help="FASTA file of the human reference genome", default=None)
    parser.add_argument("-d", "--padding", type=int, default=PADDING_DEFAULT, help=f"Padding to add to both sides of the spanned regions (default: {PADDING_DEFAULT})")
    parser.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser.add_argument("--log-file", help="Set the log output file", default=None)

    args = parser.parse_args()

    # Setup logging with the specified log level and file
    setup_logging(log_level=args.log_level.upper(), log_file=args.log_file)

    if args.human_fasta is None:
        base_name = os.path.splitext(args.human_index)[0]
        args.human_fasta = find_fasta_file(base_name)

    spanned_regions = spliced_alignment(args.human_index, args.plasmid_fasta, args.output_bam, args.padding)
    extract_human_reference(args.human_fasta, spanned_regions, args.output_fasta)
