import argparse
import sys
import shutil
import json
import os
import logging

# Load configuration from JSON file
config_path = os.path.join(os.path.dirname(__file__), 'config.json')
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

REQUIRED_TOOLS = config['required_tools']
REQUIRED_PYTHON_PACKAGES = config['required_python_packages']
DEFAULT_THRESHOLD = config['default_threshold']
VERSION = config['version']

def setup_logging(log_level=logging.INFO, log_file=None):
    """Set up logging configuration."""
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=log_file,
        filemode='a'  # Append mode for logging to a file
    )
    if not log_file:  # If no log file, also log to stdout
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

def print_logo():
    logo_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'static', 'img', 'plasmicheck_ascii.txt')
    with open(logo_path, 'r') as f:
        print(f.read())

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
        logging.error("The following requirements are not satisfied:")
        if missing_tools:
            logging.error("Missing tools: %s", ", ".join(missing_tools))
        if missing_packages:
            logging.error("Missing Python packages: %s", ", ".join(missing_packages))
        sys.exit(1)

def main(DEFAULT_THRESHOLD=DEFAULT_THRESHOLD):
    parser = argparse.ArgumentParser(description="plasmicheck: Detect and quantify plasmid DNA contamination in sequencing data")
    parser.add_argument('--log-level', help="Set the logging level", default="INFO")
    parser.add_argument('--log-file', help="Set the log output file", default=None)
    
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {VERSION}')
    subparsers = parser.add_subparsers(dest="command")

    # Print the logo and version when --version is called
    if '-v' in sys.argv or '--version' in sys.argv:
        print_logo()
        print(f"Version: {VERSION}")
        sys.exit(0)

    # Convert Command
    parser_convert = subparsers.add_parser("convert", help="Convert a plasmid file to a FASTA file and optionally generate a shifted reference")
    parser_convert.add_argument("-i", "--input_file", help="Input plasmid file", required=True)
    parser_convert.add_argument("-o", "--output_file", help="Output FASTA file", required=True)
    parser_convert.add_argument("-t", "--file_type", choices=['genbank', 'xdna'], help="Type of input file: 'genbank' or 'xdna'", required=True)
    parser_convert.add_argument("-sb", "--shift_bases", type=int, default=500, help="Number of bases to shift in the shifted reference (default: 500)")
    parser_convert.add_argument("-g", "--generate_shifted", action="store_true", help="Generate a shifted reference sequence")
    parser_convert.add_argument("-w", "--overwrite", action="store_true", help="Overwrite existing output file")
    parser_convert.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser_convert.add_argument("--log-file", help="Set the log output file", default=None)

    # Index Command
    parser_index = subparsers.add_parser("index", help="Create Minimap2 and Samtools indexes for a FASTA file")
    parser_index.add_argument("-f", "--fasta_file", help="FASTA file to index", required=True)
    parser_index.add_argument("-w", "--overwrite", action="store_true", help="Overwrite existing index files")
    parser_index.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser_index.add_argument("--log-file", help="Set the log output file", default=None)

    # Align Command
    parser_align = subparsers.add_parser("align", help="Align reads to a reference and generate a BAI index")
    parser_align.add_argument("-r", "--reference_index", help="Minimap2 index for the reference genome", required=True)
    parser_align.add_argument("-i", "--input_files", nargs='+', help="Input file(s) (BAM, interleaved FASTQ, or paired FASTQ files)", required=True)
    parser_align.add_argument("-o", "--output_bam", help="Output BAM file for alignment", required=True)
    parser_align.add_argument("-a", "--alignment_type", help="Type of alignment: 'human' or 'plasmid'", required=True)
    parser_align.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser_align.add_argument("--log-file", help="Set the log output file", default=None)

    # Compare Command
    parser_compare = subparsers.add_parser("compare", help="Compare alignments and assign reads")
    parser_compare.add_argument("-p", "--plasmid_bam", help="BAM file for plasmid alignment", required=True)
    parser_compare.add_argument("-m", "--human_bam", help="BAM file for human alignment", required=True)
    parser_compare.add_argument("-o", "--output_basename", help="Basename for output files", required=True)
    parser_compare.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")
    parser_compare.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser_compare.add_argument("--log-file", help="Set the log output file", default=None)

    # Spliced Command
    parser_spliced = subparsers.add_parser("spliced", help="Perform spliced alignment and extract human reference regions")
    parser_spliced.add_argument("-o", "--output_fasta", help="Output FASTA file for the extracted human reference regions", required=True)
    parser_spliced.add_argument("-i", "--human_index", help="Minimap2 index for the human reference genome", required=True)
    parser_spliced.add_argument("-p", "--plasmid_fasta", help="FASTA file of the plasmid reference", required=True)
    parser_spliced.add_argument("-b", "--output_bam", help="Output BAM file for the spliced alignment", required=True)
    parser_spliced.add_argument("-hf", "--human_fasta", help="FASTA file of the human reference genome", default=None)
    parser_spliced.add_argument("-d", "--padding", type=int, default=1000, help="Padding to add to both sides of the spanned regions (default: 1000)")
    parser_spliced.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser_spliced.add_argument("--log-file", help="Set the log output file", default=None)

    # Pipeline Command
    parser_pipeline = subparsers.add_parser("pipeline", help="Run the full pipeline to detect and quantify plasmid DNA contamination in sequencing data")
    parser_pipeline.add_argument("-hf", "--human_fasta", help="Human reference FASTA file", required=True)
    parser_pipeline.add_argument("-pf", "--plasmid_files", help="Plasmid files (single file or a file containing paths to multiple files)", required=True)
    parser_pipeline.add_argument("-sf", "--sequencing_files", help="Sequencing files (single file or a file containing paths to multiple files)", required=True)
    parser_pipeline.add_argument("-o", "--output_folder", help="Folder to write all outputs and intermediate files", required=True)
    parser_pipeline.add_argument("-k", "--keep_intermediate", type=bool, default=True, help="Keep intermediate files (default: True)")
    parser_pipeline.add_argument("-sb", "--shift_bases", type=int, default=500, help="Number of bases to shift in the shifted reference (default: 500)")
    parser_pipeline.add_argument("-g", "--generate_shifted", action="store_true", help="Generate a shifted reference sequence")
    parser_pipeline.add_argument("-w", "--overwrite", action="store_true", help="Overwrite existing output files")
    parser_pipeline.add_argument("-d", "--padding", type=int, default=1000, help="Padding to add to both sides of the spanned regions (default: 1000)")
    parser_pipeline.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")
    parser_pipeline.add_argument("-md5", "--md5_level", type=str, choices=["all", "input", "intermediate", "output"], default="intermediate", help="Level of MD5 checksum calculation (default: all)")
    parser_pipeline.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser_pipeline.add_argument("--log-file", help="Set the log output file", default=None)

    # Report Command
    parser_report = subparsers.add_parser("report", help="Generate a visualized HTML/PDF report from alignment comparison results")
    parser_report.add_argument("-r", "--reads_assignment_file", help="Reads assignment file (reads_assignment.tsv)", required=True)
    parser_report.add_argument("-s", "--summary_file", help="Summary file (summary.tsv)", required=True)
    parser_report.add_argument("-o", "--output_folder", help="Folder to write the report and plots", required=True)
    parser_report.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")
    parser_report.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser_report.add_argument("--log-file", help="Set the log output file", default=None)

    # Summary Reports Command
    parser_summary_reports = subparsers.add_parser("summary_reports", help="Generate summary reports for multiple samples and plasmids.")
    parser_summary_reports.add_argument("-i", "--input_dir", help="Directory containing compare outputs", required=True)
    parser_summary_reports.add_argument("-o", "--output_dir", help="Directory to save the plots", required=True)
    parser_summary_reports.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")
    parser_summary_reports.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser_summary_reports.add_argument("--log-file", help="Set the log output file", default=None)

    args = parser.parse_args()

    # Setup logging based on command-line arguments
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    setup_logging(log_level=log_level, log_file=args.log_file)

    check_requirements()

    if args.command == "convert":
        from .scripts.convert_plasmidfile_to_fasta import convert_plasmidfile_to_fasta
        convert_plasmidfile_to_fasta(args.input_file, args.output_file, args.file_type, args.shift_bases, args.generate_shifted, args.overwrite)
    elif args.command == "index":
        from .scripts.create_indexes import create_indexes
        create_indexes(args.fasta_file, args.overwrite)
    elif args.command == "align":
        from .scripts.align_reads import align_reads
        align_reads(args.reference_index, args.input_files, args.output_bam, args.alignment_type)
    elif args.command == "compare":
        from .scripts.compare_alignments import compare_alignments
        compare_alignments(args.plasmid_bam, args.human_bam, args.output_basename, args.threshold)
    elif args.command == "spliced":
        from .scripts.spliced_alignment import spliced_alignment, extract_human_reference, find_fasta_file
        if args.human_fasta is None:
            base_name = os.path.splitext(args.human_index)[0]
            args.human_fasta = find_fasta_file(base_name)
        spanned_regions = spliced_alignment(args.human_index, args.plasmid_fasta, args.output_bam, args.padding)
        extract_human_reference(args.human_fasta, spanned_regions, args.output_fasta)
    elif args.command == "pipeline":
        from .scripts.run_pipeline import run_pipeline
        run_pipeline(args.human_fasta, args.plasmid_files, args.sequencing_files, args.output_folder, args.keep_intermediate, args.shift_bases, args.generate_shifted, args.overwrite, args.padding, args.threshold, args.md5_level)
    elif args.command == "report":
        from .scripts.generate_report import main as generate_report
        command_line = ' '.join(sys.argv)
        generate_report(args.reads_assignment_file, args.summary_file, args.output_folder, args.threshold, command_line)
    elif args.command == "summary_reports":
        from .scripts.generate_summary_reports import main as generate_summary_reports
        generate_summary_reports(args.input_dir, args.output_dir, args.threshold)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
