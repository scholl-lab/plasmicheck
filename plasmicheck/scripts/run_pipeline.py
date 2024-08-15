import os
import subprocess
import re
import json
import sys
import logging
from .convert_plasmidfile_to_fasta import convert_plasmidfile_to_fasta
from .create_indexes import create_indexes
from .spliced_alignment import spliced_alignment, extract_human_reference, extract_plasmid_cDNA_positions
from .align_reads import align_reads
from .compare_alignments import compare_alignments
from .generate_report import main as generate_report, DEFAULT_THRESHOLD
from .utils import write_md5sum, sanitize_filename, setup_logging, archive_output_folder

# Resolve the path to config.json in the parent directory of the current script
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')

# Load configuration from JSON file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

def read_file_list(file_path):
    """Read a file containing newline-separated file paths and return a list of files."""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return [line.strip() for line in lines if line.strip()]

def get_file_list(file_or_list):
    """Determine if input is a single file or a list of files."""
    if os.path.isfile(file_or_list) and file_or_list.endswith('.txt'):
        return read_file_list(file_or_list)
    else:
        return [file_or_list]

def run_pipeline(human_fasta, plasmid_files, sequencing_files, output_folder, keep_intermediate=True, shift_bases=500, generate_shifted=False, overwrite=False, padding=1000, threshold=DEFAULT_THRESHOLD, md5_level="all", cDNA_output=None, archive_output=False):
    logging.info("Starting the pipeline...")
    # log the input parameters
    logging.info("Input parameters:")
    logging.info(f"Human reference FASTA: {human_fasta}")
    logging.info(f"Plasmid files: {plasmid_files}")
    logging.info(f"Sequencing files: {sequencing_files}")
    logging.info(f"Output folder: {output_folder}")
    
    plasmid_files = get_file_list(plasmid_files)
    sequencing_files = get_file_list(sequencing_files)

    for plasmid_file in plasmid_files:
        plasmid_file_type = 'genbank' if plasmid_file.endswith('.gb') or plasmid_file.endswith('.gbk') else 'xdna'
        logging.info(f"Processing plasmid file: {plasmid_file}")
        
        for sequencing_file in sequencing_files:
            logging.info(f"Processing sequencing file: {sequencing_file}")
            
            if sequencing_file.endswith('.bam'):
                sequencing_file_type = 'bam'
                fastq2 = None
            elif sequencing_file.endswith('.fastq'):
                sequencing_file_type = 'interleaved_fastq'
                fastq2 = None
            elif ',' in sequencing_file:
                sequencing_file_type = 'paired_fastq'
                sequencing_file, fastq2 = sequencing_file.split(',')
            else:
                logging.error("Unsupported sequencing file type. Must be .bam, .fastq, or paired FASTQ files separated by a comma.")
                raise ValueError("Unsupported sequencing file type. Must be .bam, .fastq, or paired FASTQ files separated by a comma.")

            bam_basename = sanitize_filename(os.path.splitext(os.path.basename(sequencing_file))[0])
            file_basename = sanitize_filename(os.path.splitext(os.path.basename(plasmid_file))[0])
            output_subfolder = os.path.join(output_folder, bam_basename, file_basename)

            if not os.path.exists(output_subfolder):
                logging.debug(f"Creating output subfolder: {output_subfolder}")
                os.makedirs(output_subfolder)

            # Step 1: Convert the plasmid file to a FASTA file or check if it exists
            plasmid_fasta = os.path.join(output_subfolder, os.path.splitext(os.path.basename(plasmid_file))[0] + ".fasta")
            if not os.path.exists(plasmid_fasta) or overwrite:
                logging.info(f"Converting {plasmid_file} to FASTA format.")
                logging.debug(f"Running: convert_plasmidfile_to_fasta({plasmid_file}, {plasmid_fasta}, {plasmid_file_type}, {shift_bases}, {generate_shifted}, {overwrite})")
                convert_plasmidfile_to_fasta(plasmid_file, plasmid_fasta, plasmid_file_type, shift_bases, generate_shifted, overwrite)
            if md5_level in ["all", "intermediate"]:
                write_md5sum(plasmid_fasta, "intermediate", output_subfolder)
            if md5_level in ["all"]:
                write_md5sum(plasmid_file, "input", output_subfolder)

            # Step 2: Generate indices for the human and plasmid FASTA files or check if they exist
            human_index = os.path.join(os.path.dirname(human_fasta), os.path.splitext(os.path.basename(human_fasta))[0] + ".mmi")
            plasmid_index = os.path.join(output_subfolder, os.path.splitext(os.path.basename(plasmid_fasta))[0] + ".mmi")
            if not os.path.exists(human_index) or overwrite:
                logging.info(f"Creating index for human reference: {human_fasta}")
                logging.debug(f"Running: create_indexes({human_fasta}, {overwrite})")
                create_indexes(human_fasta, overwrite)
            if not os.path.exists(plasmid_index) or overwrite:
                logging.info(f"Creating index for plasmid FASTA: {plasmid_fasta}")
                logging.debug(f"Running: create_indexes({plasmid_fasta}, {overwrite})")
                create_indexes(plasmid_fasta, overwrite)
            if md5_level in ["all", "intermediate"]:
                write_md5sum(plasmid_index, "intermediate", output_subfolder)
            if md5_level in ["all"]:
                write_md5sum(human_index, "intermediate", output_subfolder)

            # Step 3: Perform spliced alignment and extract human reference regions
            spliced_bam = os.path.join(output_subfolder, "spliced_alignment.bam")
            spliced_fasta = os.path.join(output_subfolder, "spliced_reference.fasta")
            logging.info(f"Performing spliced alignment for {sequencing_file} and {plasmid_file}...")
            logging.debug(f"Running: spliced_alignment({human_index}, {plasmid_fasta}, {spliced_bam}, {padding})")
            spanned_regions = spliced_alignment(human_index, plasmid_fasta, spliced_bam, padding)
            logging.info("Extracting human reference regions...")
            logging.debug(f"Running: extract_human_reference({human_fasta}, {spanned_regions}, {spliced_fasta})")
            extract_human_reference(human_fasta, spanned_regions, spliced_fasta)
            if md5_level in ["all", "intermediate"]:
                write_md5sum(spliced_bam, "intermediate", output_subfolder)
                write_md5sum(spliced_fasta, "intermediate", output_subfolder)

            # Generate index for the spliced reference
            spliced_index = os.path.join(output_subfolder, os.path.splitext(os.path.basename(spliced_fasta))[0] + ".mmi")
            if not os.path.exists(spliced_index) or overwrite:
                logging.info("Creating index for spliced reference...")
                logging.debug(f"Running: create_indexes({spliced_fasta}, {overwrite})")
                create_indexes(spliced_fasta, overwrite)
            if md5_level in ["all", "intermediate"]:
                write_md5sum(spliced_index, "intermediate", output_subfolder)

            # Step 4: Extract cDNA positions
            unique_cDNA_output = cDNA_output or os.path.join(output_subfolder, "cDNA_positions.txt")
            logging.info("Extracting cDNA positions from the spliced alignment...")
            logging.debug(f"Running: extract_plasmid_cDNA_positions({plasmid_fasta}, {spliced_bam}, {unique_cDNA_output})")
            extract_plasmid_cDNA_positions(plasmid_fasta, spliced_bam, unique_cDNA_output)
            if md5_level in ["all", "intermediate"]:
                write_md5sum(unique_cDNA_output, "intermediate", output_subfolder)

            # Step 5: Align reads to the plasmid and the spliced reference
            plasmid_bam = os.path.join(output_subfolder, "plasmid_alignment.bam")
            spliced_human_bam = os.path.join(output_subfolder, "spliced_human_alignment.bam")
            logging.info("Aligning reads to the plasmid and spliced human reference...")
            logging.debug(f"Running: align_reads({plasmid_index}, {sequencing_file}, {plasmid_bam}, 'plasmid', {fastq2})")
            align_reads(plasmid_index, sequencing_file, plasmid_bam, "plasmid", fastq2)
            logging.debug(f"Running: align_reads({spliced_index}, {sequencing_file}, {spliced_human_bam}, 'human', {fastq2})")
            align_reads(spliced_index, sequencing_file, spliced_human_bam, "human", fastq2)
            if md5_level in ["all", "intermediate"]:
                write_md5sum(plasmid_bam, "intermediate", output_subfolder)
                write_md5sum(spliced_human_bam, "intermediate", output_subfolder)

            # Step 6: Compare the two alignments
            comparison_output = os.path.join(output_subfolder, "comparison_result")
            logging.info("Comparing the two alignments...")
            logging.debug(f"Running: compare_alignments({plasmid_bam}, {spliced_human_bam}, {comparison_output})")
            compare_alignments(plasmid_bam, spliced_human_bam, comparison_output)

            # Step 7: Generate report
            logging.info("Generating report...")
            reads_assignment_file = f"{comparison_output}.reads_assignment.tsv"
            summary_file = f"{comparison_output}.summary.tsv"
            command_line = ' '.join(sys.argv)
            logging.debug(f"Running: generate_report({reads_assignment_file}, {summary_file}, {output_subfolder}, {threshold}, {command_line}, {human_fasta}, {plasmid_file}, {sequencing_file})")
            generate_report(reads_assignment_file, summary_file, output_subfolder, threshold, command_line, human_fasta, plasmid_file, sequencing_file)

            # Write MD5 checksums for the output files
            if md5_level in ["all", "intermediate", "output"]:
                write_md5sum(reads_assignment_file, "output", output_subfolder)
                write_md5sum(summary_file, "output", output_subfolder)

            # Step 8: Optionally delete intermediate files
            if not keep_intermediate:
                logging.info("Deleting intermediate files...")
                os.remove(plasmid_bam)
                os.remove(spliced_human_bam)
                os.remove(spliced_bam)

            # Step 9: Optionally archive the output folder
            if archive_output:
                archive_output_folder(output_subfolder)

    logging.info("Pipeline completed successfully.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the full pipeline to detect and quantify plasmid DNA contamination in sequencing data")
    parser.add_argument("-hf", "--human_fasta", help="Human reference FASTA file", required=True)
    parser.add_argument("-pf", "--plasmid_files", help="Plasmid files (single file or a file containing paths to multiple files)", required=True)
    parser.add_argument("-sf", "--sequencing_files", help="Sequencing files (single file or a file containing paths to multiple files)", required=True)
    parser.add_argument("-o", "--output_folder", help="Folder to write all outputs and intermediate files", required=True)
    parser.add_argument("-k", "--keep_intermediate", type=bool, default=True, help="Keep intermediate files (default: True)")
    parser.add_argument("-sb", "--shift_bases", type=int, default=500, help="Number of bases to shift in the shifted reference (default: 500)")
    parser.add_argument("-g", "--generate_shifted", action="store_true", help="Generate a shifted reference sequence")
    parser.add_argument("-w", "--overwrite", action="store_true", help="Overwrite existing output files")
    parser.add_argument("-d", "--padding", type=int, default=1000, help="Padding to add to both sides of the spanned regions (default: 1000)")
    parser.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")
    parser.add_argument("-md5", "--md5_level", type=str, choices=["all", "intermediate", "output"], default="intermediate", help="Level of MD5 checksum calculation (default: intermediate)")
    parser.add_argument("--cDNA_output", help="Output file for cDNA start and end positions in the plasmid reference", default=None)
    parser.add_argument("--archive_output", action="store_true", help="Archive and compress the output folder into a .tar.gz file")
    parser.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser.add_argument("--log-file", help="Set the log output file", default=None)

    args = parser.parse_args()

    setup_logging(log_level=getattr(logging, args.log_level.upper(), None), log_file=args.log_file)

    run_pipeline(args.human_fasta, args.plasmid_files, args.sequencing_files, args.output_folder, args.keep_intermediate, args.shift_bases, args.generate_shifted, args.overwrite, args.padding, args.threshold, args.md5_level, args.cDNA_output, args.archive_output)
