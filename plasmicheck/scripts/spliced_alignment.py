import os
import pysam
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
    try:
        run_command(f"minimap2 -t {MINIMAP2_THREADS} -ax splice {human_index} {plasmid_fasta} | samtools view -@ {SAMTOOLS_THREADS} -h -F 4 - | samtools sort -@ {SAMTOOLS_THREADS} -o {output_bam}")
    except Exception as e:
        logging.error(f"Failed to perform spliced alignment: {e}")
        raise

    if not os.path.isfile(output_bam):
        logging.error(f"Spliced alignment output BAM file not found: {output_bam}")
        raise FileNotFoundError(f"Expected BAM file at {output_bam} but it was not found.")

    logging.info(f"Indexing BAM file: {output_bam}")
    try:
        run_command(f"samtools index {output_bam}")
    except Exception as e:
        logging.error(f"Failed to index BAM file: {e}")
        raise

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

    if not spanned_regions:
        logging.warning(f"No spanned regions were found in the BAM file: {output_bam}")
    else:
        logging.debug(f"Spanned regions: {spanned_regions}")

    # Deduplicate and sort the spanned regions
    spanned_regions = sorted(set(spanned_regions))
    return spanned_regions

def extract_human_reference(human_fasta, spanned_regions, output_fasta):
    logging.info(f"Extracting human reference regions from {human_fasta} to {output_fasta}")
    
    # Use pysam to open the FASTA file
    fasta = pysam.FastaFile(human_fasta)

    with open(output_fasta, "w") as output_handle:
        for region in spanned_regions:
            chrom = region[0]
            start = region[1]
            end = region[2]
            logging.debug(f"Fetching region {chrom}:{start}-{end}")
            sequence = fasta.fetch(chrom, start, end)
            
            # Write the fetched sequence to the output FASTA file
            output_handle.write(f">{chrom}:{start}-{end}\n{sequence}\n")
    
    fasta.close()
    logging.info(f"Human reference regions have been extracted to {output_fasta}")

def extract_plasmid_cDNA_positions(plasmid_fasta, bam_file, output_file):
    logging.info(f"Extracting cDNA start and end positions from {bam_file} and saving to {output_file}")
    
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    cDNA_start = None
    cDNA_end = None

    for read in bamfile.fetch():
        if read.is_unmapped:
            continue

        # Initialize position on the read
        pos_on_read = 0  # Position on the read (plasmid-derived)

        for operation, length in read.cigartuples:
            if operation in {4, 5}:  # Soft clip (S) or Hard clip (H)
                pos_on_read += length  # Skip over clipped regions without considering them as part of cDNA

            elif operation == 0:  # Match or mismatch (M)
                if cDNA_start is None:
                    # Set cDNA_start as the first aligned position on the read
                    cDNA_start = pos_on_read
                
                # Update cDNA_end to the last position in this match/mismatch
                cDNA_end = pos_on_read + length - 1
                pos_on_read += length

            elif operation == 1:  # Insertion in read (I)
                pos_on_read += length  # Move along the read, but do not count as aligned

            elif operation in {2, 3}:  # Deletion (D) or Skipped region from the reference (N)
                # Only move along the reference, no change to the read position
                continue

    bamfile.close()

    if cDNA_start is None or cDNA_end is None:
        logging.error("Failed to identify cDNA start and end positions.")
        raise ValueError("Could not determine cDNA start and end positions in the plasmid.")

    logging.debug(f"cDNA start: {cDNA_start}, cDNA end: {cDNA_end}")

    # Compute the INSERT_REGION dynamically
    insert_region = (cDNA_start, cDNA_end)
    
    with open(output_file, "w") as f:
        f.write(f"cDNA start position in plasmid: {cDNA_start}\n")
        f.write(f"cDNA end position in plasmid: {cDNA_end}\n")
        f.write(f"INSERT_REGION: {insert_region}\n")
    
    logging.info(f"cDNA positions and INSERT_REGION saved to {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Perform spliced alignment and extract human reference regions")
    parser.add_argument("-o", "--output_fasta", help="Output FASTA file for the extracted human reference regions", required=True)
    parser.add_argument("-i", "--human_index", help="Minimap2 index for the human reference genome", required=True)
    parser.add_argument("-p", "--plasmid_fasta", help="FASTA file of the plasmid reference", required=True)
    parser.add_argument("-b", "--output_bam", help="Output BAM file for the spliced alignment", required=True)
    parser.add_argument("-hf", "--human_fasta", help="FASTA file of the human reference genome", default=None)
    parser.add_argument("-d", "--padding", type=int, default=PADDING_DEFAULT, help=f"Padding to add to both sides of the spanned regions (default: {PADDING_DEFAULT})")
    parser.add_argument("--cDNA_output", help="Output file for cDNA start and end positions in the plasmid reference", default=None)
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

    if args.cDNA_output is None:
        args.cDNA_output = os.path.join(os.path.dirname(args.output_fasta), "cDNA_positions.txt")

    logging.info(f"Saving cDNA positions to: {args.cDNA_output}")
    extract_plasmid_cDNA_positions(args.plasmid_fasta, args.output_bam, args.cDNA_output)

    # Verify that the cDNA positions file was created successfully
    if not os.path.isfile(args.cDNA_output):
        logging.error(f"Failed to create cDNA positions file: {args.cDNA_output}")
        raise FileNotFoundError(f"Expected cDNA positions file at {args.cDNA_output} but it was not found.")
    else:
        logging.info(f"cDNA positions file successfully created: {args.cDNA_output}")
