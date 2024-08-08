import argparse
import sys
import shutil
import json
import os

# Load configuration from JSON file
with open(os.path.join(os.path.dirname(__file__), 'config.json'), 'r') as config_file:
    config = json.load(config_file)

REQUIRED_TOOLS = config['required_tools']
REQUIRED_PYTHON_PACKAGES = config['required_python_packages']
DEFAULT_THRESHOLD = config['default_threshold']

def check_tools():
    missing_tools = []
    for tool in REQUIRED_TOOLS:
        if not shutil.which(tool):
            missing_tools.append(tool)
    return missing_tools

def check_python_packages():
    missing_packages = []
    for package in REQUIRED_PYTHON_PACKAGES:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    return missing_packages

def check_requirements():
    missing_tools = check_tools()
    missing_packages = check_python_packages()
    
    if missing_tools or missing_packages:
        print("The following requirements are not satisfied:")
        if missing_tools:
            print("Missing tools:", ", ".join(missing_tools))
        if missing_packages:
            print("Missing Python packages:", ", ".join(missing_packages))
        sys.exit(1)

def main(DEFAULT_THRESHOLD=DEFAULT_THRESHOLD):
    parser = argparse.ArgumentParser(description="PlasmiCheck: Detect and quantify plasmid DNA contamination in sequencing data")
    subparsers = parser.add_subparsers(dest="command")

    # Convert Command
    parser_convert = subparsers.add_parser("convert", help="Convert a GenBank file to a FASTA file and optionally generate a shifted reference")
    parser_convert.add_argument("input_file", help="Input GenBank file")
    parser_convert.add_argument("output_file", help="Output FASTA file")
    parser_convert.add_argument("--shift_bases", type=int, default=500, help="Number of bases to shift in the shifted reference (default: 500)")
    parser_convert.add_argument("--generate_shifted", action="store_true", help="Generate a shifted reference sequence")
    parser_convert.add_argument("--overwrite", action="store_true", help="Overwrite existing output file")

    # Index Command
    parser_index = subparsers.add_parser("index", help="Create Minimap2 and Samtools indexes for a FASTA file")
    parser_index.add_argument("fasta_file", help="FASTA file to index")
    parser_index.add_argument("--overwrite", action="store_true", help="Overwrite existing index files")

    # Align Command
    parser_align = subparsers.add_parser("align", help="Align reads to a reference and generate a BAI index")
    parser_align.add_argument("reference_index", help="Minimap2 index for the reference genome")
    parser_align.add_argument("input_file", help="Input file (BAM, interleaved FASTQ, or first FASTQ file for paired FASTQ)")
    parser_align.add_argument("output_bam", help="Output BAM file for alignment")
    parser_align.add_argument("alignment_type", help="Type of alignment: 'human' or 'plasmid'")
    parser_align.add_argument("file_type", help="Type of input file: 'bam', 'interleaved_fastq', or 'paired_fastq'")
    parser_align.add_argument("--fastq2", help="Second FASTQ file for paired FASTQ input", default=None)

    # Compare Command
    parser_compare = subparsers.add_parser("compare", help="Compare alignments and assign reads")
    parser_compare.add_argument("plasmid_bam", help="BAM file for plasmid alignment")
    parser_compare.add_argument("human_bam", help="BAM file for human alignment")
    parser_compare.add_argument("output_basename", help="Basename for output files")

    # Spliced Command
    parser_spliced = subparsers.add_parser("spliced", help="Perform spliced alignment and extract human reference regions")
    parser_spliced.add_argument("output_fasta", help="Output FASTA file for the extracted human reference regions")
    parser_spliced.add_argument("human_index", help="Minimap2 index for the human reference genome")
    parser_spliced.add_argument("plasmid_fasta", help="FASTA file of the plasmid reference")
    parser_spliced.add_argument("output_bam", help="Output BAM file for the spliced alignment")
    parser_spliced.add_argument("--human_fasta", help="FASTA file of the human reference genome", default=None)
    parser_spliced.add_argument("--padding", type=int, default=1000, help="Padding to add to both sides of the spanned regions (default: 1000)")

    # Pipeline Command
    parser_pipeline = subparsers.add_parser("pipeline", help="Run the full pipeline to detect and quantify plasmid DNA contamination in sequencing data")
    parser_pipeline.add_argument("human_fasta", help="Human reference FASTA file")
    parser_pipeline.add_argument("plasmid_gb", help="GenBank plasmid file")
    parser_pipeline.add_argument("sequencing_file", help="Sequencing file (BAM, interleaved FASTQ, or first FASTQ file for paired FASTQ)")
    parser_pipeline.add_argument("output_folder", help="Folder to write all outputs and intermediate files")
    parser_pipeline.add_argument("file_type", help="Type of input file: 'bam', 'interleaved_fastq', or 'paired_fastq'")
    parser_pipeline.add_argument("--fastq2", help="Second FASTQ file for paired FASTQ input", default=None)
    parser_pipeline.add_argument("--keep_intermediate", action="store_true", help="Keep intermediate files (default: delete them)")
    parser_pipeline.add_argument("--shift_bases", type=int, default=500, help="Number of bases to shift in the shifted reference (default: 500)")
    parser_pipeline.add_argument("--generate_shifted", action="store_true", help="Generate a shifted reference sequence")
    parser_pipeline.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    parser_pipeline.add_argument("--padding", type=int, default=1000, help="Padding to add to both sides of the spanned regions (default: 1000)")
    parser_pipeline.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")

    # Report Command
    parser_report = subparsers.add_parser("report", help="Generate a visualized HTML/PDF report from alignment comparison results")
    parser_report.add_argument("reads_assignment_file", help="Reads assignment file (reads_assignment.tsv)")
    parser_report.add_argument("summary_file", help="Summary file (summary.tsv)")
    parser_report.add_argument("output_folder", help="Folder to write the report and plots")
    parser_report.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")

    args = parser.parse_args()

    check_requirements()

    if args.command == "convert":
        from .scripts.convert_gb_to_fasta import convert_gb_to_fasta
        convert_gb_to_fasta(args.input_file, args.output_file, args.shift_bases, args.generate_shifted, args.overwrite)
    elif args.command == "index":
        from .scripts.create_indexes import create_indexes
        create_indexes(args.fasta_file, args.overwrite)
    elif args.command == "align":
        from .scripts.align_reads import align_reads
        align_reads(args.reference_index, args.input_file, args.output_bam, args.alignment_type, args.file_type, args.fastq2)
    elif args.command == "compare":
        from .scripts.compare_alignments import compare_alignments
        compare_alignments(args.plasmid_bam, args.human_bam, args.output_basename)
    elif args.command == "spliced":
        from .scripts.spliced_alignment import spliced_alignment, extract_human_reference, find_fasta_file
        if args.human_fasta is None:
            base_name = os.path.splitext(args.human_index)[0]
            args.human_fasta = find_fasta_file(base_name)
        spanned_regions = spliced_alignment(args.human_index, args.plasmid_fasta, args.output_bam, args.padding)
        extract_human_reference(args.human_fasta, spanned_regions, args.output_fasta)
    elif args.command == "pipeline":
        from .scripts.run_pipeline import run_pipeline, DEFAULT_THRESHOLD
        run_pipeline(args.human_fasta, args.plasmid_gb, args.sequencing_file, args.output_folder, args.file_type, args.fastq2, args.keep_intermediate, args.shift_bases, args.generate_shifted, args.overwrite, args.padding, args.threshold)
    elif args.command == "report":
        from .scripts.generate_report import main as generate_report
        generate_report(args.reads_assignment_file, args.summary_file, args.output_folder, args.threshold)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
