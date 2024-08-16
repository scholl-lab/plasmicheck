import pandas as pd
import plotly.express as px
import os
import json
import base64
import logging
from jinja2 import Environment, FileSystemLoader
from datetime import datetime

# Import the version from version.py
from plasmicheck.version import __version__ as VERSION

from .utils import setup_logging  # Import the setup_logging function

# Resolve the path to config.json in the parent directory of the current script
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')

# Load configuration from JSON file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

DEFAULT_THRESHOLD = config['default_threshold']
UNCLEAR_RANGE = config['unclear_range']
PLOT_SAMPLE_REPORT = config['plot_sample_report']
TEMPLATE_DIR = config['paths']['template_dir']
LOGO_PATH = config['paths']['logo_path']
PLOT_DIMENSIONS = PLOT_SAMPLE_REPORT.get('figsize', {'width': 800, 'height': 600})

# Setup logging for the libraries used in this script
logging.getLogger('jinja2').setLevel(logging.ERROR)
logging.getLogger('fontTools').setLevel(logging.ERROR)

def load_data(reads_assignment_file, summary_file):
    logging.info(f"Loading data from {reads_assignment_file} and {summary_file}")
    reads_df = pd.read_csv(reads_assignment_file, sep='\t')
    summary_df = pd.read_csv(summary_file, sep='\t')
    return reads_df, summary_df

def generate_plots(reads_df, output_folder):
    logging.info("Generating plots")
    counts = reads_df['AssignedTo'].value_counts().to_dict()

    # Ensure default values for width and height
    width = PLOT_SAMPLE_REPORT.get('figsize', {}).get('width', 600)  # Convert inches to pixels
    height = PLOT_SAMPLE_REPORT.get('figsize', {}).get('height', 400)  # Convert inches to pixels

    # Validate that width and height are integers
    if not isinstance(width, int):
        logging.warning(f"Invalid width value: {width}. Defaulting to 600.")
        width = 600
    if not isinstance(height, int):
        logging.warning(f"Invalid height value: {height}. Defaulting to 400.")
        height = 400

    # Create directory for plots if not exist
    plots_dir = os.path.join(output_folder, 'plots')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    # Box plot using Plotly
    boxplot_df = reads_df.copy()
    fig_box = px.box(
        boxplot_df, 
        x='AssignedTo', 
        y='PlasmidScore', 
        points='all',
        color='AssignedTo',
        title=f"{PLOT_SAMPLE_REPORT['title_box_plot']}<br>(Total Reads: {len(reads_df)})",
        labels={"PlasmidScore": PLOT_SAMPLE_REPORT['box_plot_y_label'], "AssignedTo": PLOT_SAMPLE_REPORT['box_plot_x_label']},
        width=width, 
        height=height
    )
    
    boxplot_filename_interactive = os.path.join(plots_dir, PLOT_SAMPLE_REPORT['output_box_plot_filename'].replace('.png', '.html'))
    fig_box.write_html(boxplot_filename_interactive)
    
    boxplot_filename_png = os.path.join(plots_dir, PLOT_SAMPLE_REPORT['output_box_plot_filename'])
    fig_box.write_image(boxplot_filename_png)

    # Scatter plot using Plotly
    fig_scatter = px.scatter(
        reads_df, 
        x='PlasmidScore', 
        y='HumanScore', 
        color='AssignedTo',
        title=f"{PLOT_SAMPLE_REPORT['title_scatter_plot']}<br>(Total Reads: {len(reads_df)})",
        labels={"PlasmidScore": PLOT_SAMPLE_REPORT['scatter_plot_x_label'], "HumanScore": PLOT_SAMPLE_REPORT['scatter_plot_y_label']},
        width=width, 
        height=height
    )
    
    scatter_filename_interactive = os.path.join(plots_dir, PLOT_SAMPLE_REPORT['output_scatter_plot_filename'].replace('.png', '.html'))
    fig_scatter.write_html(scatter_filename_interactive)

    scatter_filename_png = os.path.join(plots_dir, PLOT_SAMPLE_REPORT['output_scatter_plot_filename'])
    fig_scatter.write_image(scatter_filename_png)

    logging.info("Plots generated.")
    
    return boxplot_filename_interactive, boxplot_filename_png, scatter_filename_interactive, scatter_filename_png

def encode_image_to_base64(image_path):
    logging.info(f"Encoding image {image_path} to base64")
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded_string}"

def extract_verdict_from_summary(summary_df):
    verdict_row = summary_df[summary_df['Category'] == 'Verdict']
    if not verdict_row.empty:
        return verdict_row['Count'].values[0]
    else:
        return "Verdict not found in summary file"

def generate_report(summary_df, output_folder, verdict, ratio, threshold, unclear, command_line, human_fasta="None", plasmid_gb="None", sequencing_file="None", boxplot_filename_interactive=None, boxplot_filename_png=None, scatter_filename_interactive=None, scatter_filename_png=None):
    logging.info("Generating report")
    env = Environment(loader=FileSystemLoader(os.path.join('plasmicheck', TEMPLATE_DIR)))
    template = env.get_template('report_template.html')

    logo_base64 = encode_image_to_base64(os.path.join('plasmicheck', LOGO_PATH))

    verdict_color = "green" if "Sample is not contaminated with plasmid DNA" in verdict else "orange" if "Sample contamination status is unclear" in verdict else "red"

    # Read the interactive HTML plot content
    with open(boxplot_filename_interactive, 'r') as f:
        boxplot_html_content = f.read()

    with open(scatter_filename_interactive, 'r') as f:
        scatter_html_content = f.read()

    # Convert the static PNG files to Base64
    boxplot_png_base64 = encode_image_to_base64(boxplot_filename_png)
    scatter_png_base64 = encode_image_to_base64(scatter_filename_png)

    # Render the interactive HTML report
    html_content_interactive = template.render(
        summary_df=summary_df.to_html(classes='table table-striped'),
        box_plot=boxplot_html_content,
        scatter_plot=scatter_html_content,
        verdict=verdict,
        ratio=f"{ratio:.3f}",
        threshold=threshold,
        unclear_range=unclear,
        verdict_color=verdict_color,
        version=VERSION,
        logo_base64=logo_base64,
        human_fasta=human_fasta,
        plasmid_gb=plasmid_gb,
        sequencing_file=sequencing_file,
        run_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line,
        interactive=True
    )

    # Render the non-interactive HTML report
    html_content_non_interactive = template.render(
        summary_df=summary_df.to_html(classes='table table-striped'),
        box_plot=boxplot_png_base64,
        scatter_plot=scatter_png_base64,
        verdict=verdict,
        ratio=f"{ratio:.3f}",
        threshold=threshold,
        unclear_range=unclear,
        verdict_color=verdict_color,
        version=VERSION,
        logo_base64=logo_base64,
        human_fasta=human_fasta,
        plasmid_gb=plasmid_gb,
        sequencing_file=sequencing_file,
        run_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line,
        interactive=False
    )

    # Save interactive HTML report
    html_report_interactive = os.path.join(output_folder, 'report_interactive.html')
    with open(html_report_interactive, 'w') as f:
        f.write(html_content_interactive)

    # Save non-interactive HTML report
    html_report_non_interactive = os.path.join(output_folder, 'report_non_interactive.html')
    with open(html_report_non_interactive, 'w') as f:
        f.write(html_content_non_interactive)

def main(reads_assignment_file, summary_file, output_folder, threshold=DEFAULT_THRESHOLD, unclear=UNCLEAR_RANGE, human_fasta="None", plasmid_gb="None", sequencing_file="None", command_line=""):
    reads_df, summary_df = load_data(reads_assignment_file, summary_file)

    # Generate plots
    boxplot_filename_interactive, boxplot_filename_png, scatter_filename_interactive, scatter_filename_png = generate_plots(reads_df, output_folder)

    # Extract the verdict from the summary file
    verdict = extract_verdict_from_summary(summary_df)

    plasmid_count = int(summary_df[summary_df['Category'] == 'Plasmid']['Count'].values[0])
    human_count = int(summary_df[summary_df['Category'] == 'Human']['Count'].values[0])
    ratio = plasmid_count / human_count if human_count != 0 else float('inf')

    # Generate report
    generate_report(
        summary_df, 
        output_folder, 
        verdict, 
        ratio, 
        threshold, 
        UNCLEAR_RANGE, 
        command_line, 
        human_fasta, 
        plasmid_gb, 
        sequencing_file,
        boxplot_filename_interactive=boxplot_filename_interactive,
        boxplot_filename_png=boxplot_filename_png,
        scatter_filename_interactive=scatter_filename_interactive,
        scatter_filename_png=scatter_filename_png
    )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate a visualized HTML report from alignment comparison results")
    parser.add_argument("-r", "--reads_assignment_file", help="Reads assignment file (reads_assignment.tsv)", required=True)
    parser.add_argument("-s", "--summary_file", help="Summary file (summary.tsv)", required=True)
    parser.add_argument("-o", "--output_folder", help="Folder to write the report and plots", required=True)
    parser.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")
    parser.add_argument("-hf", "--human_fasta", default="None", help="Human reference FASTA file")
    parser.add_argument("-pg", "--plasmid_gb", default="None", help="GenBank plasmid file")
    parser.add_argument("-sf", "--sequencing_file", default="None", help="Sequencing file (BAM, interleaved FASTQ, or first FASTQ file for paired FASTQ)")
    parser.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser.add_argument("--log-file", help="Set the log output file", default=None)

    args = parser.parse_args()
    command_line = ' '.join(sys.argv)

    setup_logging(log_level=args.log_level.upper(), log_file=args.log_file)  # Setup logging with arguments

    main(
        args.reads_assignment_file, 
        args.summary_file, 
        args.output_folder, 
        args.threshold, 
        UNCLEAR_RANGE, 
        args.human_fasta, 
        args.plasmid_gb, 
        args.sequencing_file, 
        command_line
    )
