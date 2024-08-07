import os
import subprocess
import re
from .convert_gb_to_fasta import convert_gb_to_fasta
from .create_indexes import create_indexes
from .spliced_alignment import spliced_alignment, extract_human_reference
from .align_reads import align_reads
from .compare_alignments import compare_alignments
from .generate_report import main as generate_report, DEFAULT_THRESHOLD

def sanitize_filename(filename):
    return re.sub(r'[^a-zA-Z0-9_-]', '_', filename)

def run_pipeline(human_fasta, plasmid_gb, sequencing_file, output_folder, file_type, fastq2=None, keep_intermediate=True, shift_bases=500, generate_shifted=False, overwrite=False, padding=1000, threshold=DEFAULT_THRESHOLD):
    bam_basename = sanitize_filename(os.path.splitext(os.path.basename(sequencing_file))[0])
    gb_basename = sanitize_filename(os.path.splitext(os.path.basename(plasmid_gb))[0])
    output_subfolder = os.path.join(output_folder, bam_basename, gb_basename)

    if not os.path.exists(output_subfolder):
        os.makedirs(output_subfolder)

    # Step 1: Convert the GenBank plasmid file to a FASTA file or check if it exists
    plasmid_fasta = os.path.join(output_subfolder, os.path.splitext(os.path.basename(plasmid_gb))[0] + ".fasta")
    if not os.path.exists(plasmid_fasta) or overwrite:
        convert_gb_to_fasta(plasmid_gb, plasmid_fasta, shift_bases, generate_shifted, overwrite)

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
    align_reads(plasmid_index, sequencing_file, plasmid_bam, "plasmid", file_type, fastq2)
    align_reads(spliced_index, sequencing_file, spliced_human_bam, "human", file_type, fastq2)

    # Step 5: Compare the two alignments
    comparison_output = os.path.join(output_subfolder, "comparison_result")
    compare_alignments(plasmid_bam, spliced_human_bam, comparison_output)

    # Step 6: Generate report
    reads_assignment_file = f"{comparison_output}.reads_assignment.tsv"
    summary_file = f"{comparison_output}.summary.tsv"
    generate_report(reads_assignment_file, summary_file, output_subfolder, threshold, human_fasta, plasmid_gb, sequencing_file)

    # Step 7: Optionally delete intermediate files
    if not keep_intermediate:
        os.remove(plasmid_bam)
        os.remove(spliced_human_bam)
        os.remove(spliced_bam)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the full pipeline to detect and quantify plasmid DNA contamination in sequencing data")
    parser.add_argument("human_fasta", help="Human reference FASTA file")
    parser.add_argument("plasmid_gb", help="GenBank plasmid file")
    parser.add_argument("sequencing_file", help="Sequencing file (BAM, interleaved FASTQ, or first FASTQ file for paired FASTQ)")
    parser.add_argument("output_folder", help="Folder to write all outputs and intermediate files")
    parser.add_argument("file_type", help="Type of input file: 'bam', 'interleaved_fastq', or 'paired_fastq'")
    parser.add_argument("--fastq2", help="Second FASTQ file for paired FASTQ input", default=None)
    parser.add_argument("--keep_intermediate", action="store_true", help="Keep intermediate files (default: delete them)")
    parser.add_argument("--shift_bases", type=int, default=500, help="Number of bases to shift in the shifted reference (default: 500)")
    parser.add_argument("--generate_shifted", action="store_true", help="Generate a shifted reference sequence")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--padding", type=int, default=1000, help="Padding to add to both sides of the spanned regions (default: 1000)")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")

    args = parser.parse_args()

    run_pipeline(args.human_fasta, args.plasmid_gb, args.sequencing_file, args.output_folder, args.file_type, args.fastq2, args.keep_intermediate, args.shift_bases, args.generate_shifted, args.overwrite, args.padding, args.threshold)
