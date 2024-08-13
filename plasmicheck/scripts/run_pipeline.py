import os
import subprocess
import re
import json
import sys
from .convert_plasmidfile_to_fasta import convert_plasmidfile_to_fasta
from .create_indexes import create_indexes
from .spliced_alignment import spliced_alignment, extract_human_reference
from .align_reads import align_reads
from .compare_alignments import compare_alignments
from .generate_report import main as generate_report, DEFAULT_THRESHOLD

# Resolve the path to config.json in the parent directory of the current script
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')

# Load configuration from JSON file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

def sanitize_filename(filename):
    return re.sub(r'[^a-zA-Z0-9_-]', '_', filename)

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

def run_pipeline(human_fasta, plasmid_files, sequencing_files, output_folder, keep_intermediate=True, shift_bases=500, generate_shifted=False, overwrite=False, padding=1000, threshold=DEFAULT_THRESHOLD):
    plasmid_files = get_file_list(plasmid_files)
    sequencing_files = get_file_list(sequencing_files)

    for plasmid_file in plasmid_files:
        plasmid_file_type = 'genbank' if plasmid_file.endswith('.gb') or plasmid_file.endswith('.gbk') else 'xdna'
        for sequencing_file in sequencing_files:
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
                raise ValueError("Unsupported sequencing file type. Must be .bam, .fastq, or paired FASTQ files separated by a comma.")

            bam_basename = sanitize_filename(os.path.splitext(os.path.basename(sequencing_file))[0])
            file_basename = sanitize_filename(os.path.splitext(os.path.basename(plasmid_file))[0])
            output_subfolder = os.path.join(output_folder, bam_basename, file_basename)

            if not os.path.exists(output_subfolder):
                os.makedirs(output_subfolder)

            # Step 1: Convert the plasmid file to a FASTA file or check if it exists
            plasmid_fasta = os.path.join(output_subfolder, os.path.splitext(os.path.basename(plasmid_file))[0] + ".fasta")
            if not os.path.exists(plasmid_fasta) or overwrite:
                convert_plasmidfile_to_fasta(plasmid_file, plasmid_fasta, plasmid_file_type, shift_bases, generate_shifted, overwrite)

            # Step 2: Generate indices for the human and plasmid FASTA files or check if they exist
            human_index = os.path.join(os.path.dirname(human_fasta), os.path.splitext(os.path.basename(human_fasta))[0] + ".mmi")
            plasmid_index = os.path.join(output_subfolder, os.path.splitext(os.path.basename(plasmid_fasta))[0] + ".mmi")
            if not os.path.exists(human_index) or overwrite:
                create_indexes(human_fasta, overwrite)
            if not os.path.exists(plasmid_index) or overwrite:
                create_indexes(plasmid_fasta, overwrite)

            # Step 3: Perform spliced alignment and extract human reference regions
            spliced_bam = os.path.join(output_subfolder, "spliced_alignment.bam")
            spliced_fasta = os.path.join(output_subfolder, "spliced_reference.fasta")
            spanned_regions = spliced_alignment(human_index, plasmid_fasta, spliced_bam, padding)
            extract_human_reference(human_fasta, spanned_regions, spliced_fasta)

            # Generate index for the spliced reference
            spliced_index = os.path.join(output_subfolder, os.path.splitext(os.path.basename(spliced_fasta))[0] + ".mmi")
            if not os.path.exists(spliced_index) or overwrite:
                create_indexes(spliced_fasta, overwrite)

            # Step 4: Align reads to the plasmid and the spliced reference
            plasmid_bam = os.path.join(output_subfolder, "plasmid_alignment.bam")
            spliced_human_bam = os.path.join(output_subfolder, "spliced_human_alignment.bam")
            align_reads(plasmid_index, sequencing_file, plasmid_bam, "plasmid", fastq2)
            align_reads(spliced_index, sequencing_file, spliced_human_bam, "human", fastq2)

            # Step 5: Compare the two alignments
            comparison_output = os.path.join(output_subfolder, "comparison_result")
            compare_alignments(plasmid_bam, spliced_human_bam, comparison_output)

            # Step 6: Generate report
            reads_assignment_file = f"{comparison_output}.reads_assignment.tsv"
            summary_file = f"{comparison_output}.summary.tsv"
            command_line = ' '.join(sys.argv)
            generate_report(reads_assignment_file, summary_file, output_subfolder, threshold, command_line, human_fasta, plasmid_file, sequencing_file)

            # Step 7: Optionally delete intermediate files
            if not keep_intermediate:
                os.remove(plasmid_bam)
                os.remove(spliced_human_bam)
                os.remove(spliced_bam)

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

    args = parser.parse_args()

    run_pipeline(args.human_fasta, args.plasmid_files, args.sequencing_files, args.output_folder, args.keep_intermediate, args.shift_bases, args.generate_shifted, args.overwrite, args.padding, args.threshold)
