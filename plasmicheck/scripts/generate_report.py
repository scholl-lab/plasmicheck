import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from datetime import datetime
import json
import os
import base64

# Resolve the path to config.json in the parent directory of the current script
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')

# Load configuration from JSON file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

DEFAULT_THRESHOLD = config['default_threshold']
UNCLEAR_RANGE = config['unclear_range']
VERSION = config['version']
PLOT_SIZE = config['plot_size']
TEMPLATE_DIR = config['paths']['template_dir']
LOGO_PATH = config['paths']['logo_path']

def load_data(reads_assignment_file, summary_file):
    reads_df = pd.read_csv(reads_assignment_file, sep='\t')
    summary_df = pd.read_csv(summary_file, sep='\t')
    return reads_df, summary_df

def generate_plots(reads_df, output_folder):
    counts = reads_df['AssignedTo'].value_counts().to_dict()

    plt.figure(figsize=(PLOT_SIZE['width'], PLOT_SIZE['height']))
    sns.boxplot(data=reads_df, x='AssignedTo', y='PlasmidScore')
    plt.title(f"Box Plot of Plasmid Scores by Assignment\n(Total Reads: {len(reads_df)})")
    plt.xlabel("Assigned To")
    plt.ylabel("Plasmid Score")
    for category, count in counts.items():
        plt.text(category, 0, f'N={count}', ha='center', va='bottom')
    plt.savefig(f"{output_folder}/box_plot.png")
    plt.close()

    plt.figure(figsize=(PLOT_SIZE['width'], PLOT_SIZE['height']))
    sns.scatterplot(data=reads_df, x='PlasmidScore', y='HumanScore', hue='AssignedTo')
    plt.title(f"Scatter Plot of Plasmid vs. Human Scores by Assignment\n(Total Reads: {len(reads_df)})")
    plt.xlabel("Plasmid Score")
    plt.ylabel("Human Score")
    plt.legend(title='Assigned To', loc='best')
    plt.savefig(f"{output_folder}/scatter_plot.png")
    plt.close()

def encode_image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded_string}"

def generate_report(summary_df, output_folder, verdict, ratio, threshold, unclear_range, command_line, human_fasta="None", plasmid_gb="None", sequencing_file="None"):
    env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = env.get_template('report_template.html')

    logo_base64 = encode_image_to_base64(LOGO_PATH)

    verdict_color = "green" if "not contaminated" in verdict else "orange" if "unclear" in verdict else "red"

    html_content = template.render(
        summary_df=summary_df.to_html(classes='table table-striped'),
        box_plot="box_plot.png",
        scatter_plot="scatter_plot.png",
        verdict=verdict,
        ratio=f"{ratio:.3f}",
        threshold=threshold,
        unclear_range=unclear_range,
        verdict_color=verdict_color,
        version=VERSION,
        logo_base64=logo_base64,
        human_fasta=human_fasta,
        plasmid_gb=plasmid_gb,
        sequencing_file=sequencing_file,
        run_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line
    )

    html_report = f"{output_folder}/report.html"
    with open(html_report, 'w') as f:
        f.write(html_content)

    HTML(html_report).write_pdf(f"{output_folder}/report.pdf")

def main(reads_assignment_file, summary_file, output_folder, threshold=DEFAULT_THRESHOLD, unclear_range=UNCLEAR_RANGE, human_fasta="None", plasmid_gb="None", sequencing_file="None", command_line=""):
    reads_df, summary_df = load_data(reads_assignment_file, summary_file)
    generate_plots(reads_df, output_folder)

    plasmid_count = int(summary_df[summary_df['Category'] == 'Plasmid']['Count'].values[0])
    human_count = int(summary_df[summary_df['Category'] == 'Human']['Count'].values[0])
    ratio = plasmid_count / human_count if human_count != 0 else float('inf')
    
    verdict = (
        "Sample is contaminated with plasmid DNA" if ratio > threshold else
        "Sample is not contaminated with plasmid DNA"
    )

    generate_report(
        summary_df, 
        output_folder, 
        verdict, 
        ratio, 
        threshold, 
        unclear_range, 
        command_line, 
        human_fasta, 
        plasmid_gb, 
        sequencing_file
    )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate a visualized HTML/PDF report from alignment comparison results")
    parser.add_argument("-r", "--reads_assignment_file", help="Reads assignment file (reads_assignment.tsv)", required=True)
    parser.add_argument("-s", "--summary_file", help="Summary file (summary.tsv)", required=True)
    parser.add_argument("-o", "--output_folder", help="Folder to write the report and plots", required=True)
    parser.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")
    parser.add_argument("-hf", "--human_fasta", default="None", help="Human reference FASTA file")
    parser.add_argument("-pg", "--plasmid_gb", default="None", help="GenBank plasmid file")
    parser.add_argument("-sf", "--sequencing_file", default="None", help="Sequencing file (BAM, interleaved FASTQ, or first FASTQ file for paired FASTQ)")

    args = parser.parse_args()
    command_line = ' '.join(sys.argv)
    main(args.reads_assignment_file, args.summary_file, args.output_folder, args.threshold, UNCLEAR_RANGE, args.human_fasta, args.plasmid_gb, args.sequencing_file, command_line)
