import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from datetime import datetime

# Resolve the path to config.json in the parent directory of the current script
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')

# Load configuration from JSON file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

DEFAULT_THRESHOLD = config['default_threshold']
PLOT_CONFIG = config['plot_summary']
VERSION = config['version']

def find_tsv_files(input_dir, pattern):
    """Recursively find all tsv files matching the pattern in the input directory."""
    tsv_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(pattern):
                tsv_files.append(os.path.join(root, file))
    return tsv_files

def read_compare_outputs(input_dir):
    reads_assignment_files = find_tsv_files(input_dir, '.reads_assignment.tsv')
    summary_files = find_tsv_files(input_dir, '.summary.tsv')

    if not reads_assignment_files or not summary_files:
        raise ValueError("No valid TSV files found in the specified directory structure.")
    
    reads_df_list = []
    for file in reads_assignment_files:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(file)))
        plasmid_name = os.path.basename(os.path.dirname(file))
        df = pd.read_csv(file, sep='\t')
        df['Sample'] = sample_name
        df['Plasmid'] = plasmid_name
        reads_df_list.append(df)

    summary_df_list = []
    for file in summary_files:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(file)))
        plasmid_name = os.path.basename(os.path.dirname(file))
        df = pd.read_csv(file, sep='\t')
        df['Sample'] = sample_name
        df['Plasmid'] = plasmid_name
        df = df.rename(columns={"Count": "Value"})  # Rename Count to Value
        summary_df_list.append(df)
    
    reads_df = pd.concat(reads_df_list, ignore_index=True)
    summary_df = pd.concat(summary_df_list, ignore_index=True)
    
    return reads_df, summary_df

def create_plots(reads_df, summary_df, output_dir, threshold, plot_config):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Split summary dataframe into separate dataframes for different categories
    combined_df = summary_df[summary_df['Category'].isin(['Plasmid', 'Human', 'Tied'])].copy()
    verdict_df = summary_df[summary_df['Category'] == 'Verdict'].copy()
    ratio_df = summary_df[summary_df['Category'] == 'Ratio'].copy()

    # Ensure correct data types for combined and ratio dataframes
    combined_df.loc[:, 'Value'] = pd.to_numeric(combined_df['Value'], errors='coerce')
    ratio_df.loc[:, 'Value'] = pd.to_numeric(ratio_df['Value'], errors='coerce')

    # Create heatmap data for ratios
    ratio_data = ratio_df.pivot(index="Sample", columns="Plasmid", values="Value")

    # Apply threshold for coloring
    heatmap_data = ratio_data.apply(lambda x: x.map(lambda y: 1 if y > threshold else 0))

    plt.figure(figsize=(plot_config['figsize']['width'], plot_config['figsize']['height']))  # Adjust figure size
    cmap = sns.color_palette([plot_config['colors']['not_contaminated'], plot_config['colors']['contaminated']])  # Color-blind safe colors
    sns.heatmap(heatmap_data, annot=ratio_data, fmt=".2f", cmap=cmap, cbar=False, linewidths=plot_config['linewidths'], linecolor=plot_config['linecolor'])
    plt.title(plot_config['title'])
    plt.xticks(rotation=plot_config['xticks_rotation'], ha=plot_config['xticks_ha'])
    plt.yticks(rotation=plot_config['yticks_rotation'])
    plt.tight_layout()  # Adjust layout to fit labels
    plot_filename = os.path.join(output_dir, plot_config['output_filename'])
    plt.savefig(plot_filename)
    plt.close()
    
    return plot_filename, combined_df, verdict_df, ratio_df

def generate_report(combined_df, verdict_df, ratio_df, plot_filename, output_folder, threshold, command_line):
    # Load the template
    env = Environment(loader=FileSystemLoader('templates'))
    template = env.get_template('summary_template.html')

    # Render the template
    html_content = template.render(
        combined_df=combined_df.to_html(classes='table table-striped table-bordered', index=False),
        verdict_df=verdict_df.to_html(classes='table table-striped table-bordered', index=False),
        ratio_df=ratio_df.to_html(classes='table table-striped table-bordered', index=False),
        plot_image=os.path.basename(plot_filename),  # Use the basename to correctly link the plot
        threshold=threshold,
        version=VERSION,
        run_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line
    )

    # Save HTML report
    html_report = os.path.join(output_folder, 'summary_report.html')
    with open(html_report, 'w') as f:
        f.write(html_content)

    # Convert HTML to PDF
    HTML(html_report).write_pdf(os.path.join(output_folder, 'summary_report.pdf'))

def main(input_dir, output_dir, threshold=DEFAULT_THRESHOLD):
    reads_df, summary_df = read_compare_outputs(input_dir)
    plot_filename, combined_df, verdict_df, ratio_df = create_plots(reads_df, summary_df, output_dir, threshold, PLOT_CONFIG)

    generate_report(
        combined_df, 
        verdict_df, 
        ratio_df, 
        plot_filename, 
        output_dir, 
        threshold, 
        ' '.join(os.sys.argv)
    )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Merge and generate summary reports for multiple samples and plasmids.")
    parser.add_argument("-i", "--input_dir", help="Directory containing compare outputs", required=True)
    parser.add_argument("-o", "--output_dir", help="Directory to save the plots and reports", required=True)
    parser.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")

    args = parser.parse_args()
    main(args.input_dir, args.output_dir, args.threshold)
