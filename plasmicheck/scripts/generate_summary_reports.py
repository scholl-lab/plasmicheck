import os
import pandas as pd
import json
import logging
import base64
import plotly.express as px
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from datetime import datetime
import numpy as np

from .utils import setup_logging  # Import setup_logging function

# Resolve the path to config.json in the parent directory of the current script
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')

# Load configuration from JSON file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

DEFAULT_THRESHOLD = config['default_threshold']
UNCLEAR_RANGE = config['unclear_range']
PLOT_CONFIG = config['plot_summary']
VERSION = config['version']
TEMPLATE_DIR = config['paths']['template_dir']
LOGO_PATH = config['paths']['logo_path']
TABLE_SORTING = config['table_sorting']

# Replace magic numbers with config values
PLOT_DIMENSIONS = PLOT_CONFIG.get('plot_dimensions', {'width': 1200, 'height': 1200})
ROUND_DECIMALS = PLOT_CONFIG.get('round_decimals', 3)
LOG_OFFSET = PLOT_CONFIG.get('log_offset', 1e-9)
MARKER_STYLE = PLOT_CONFIG.get('marker_style', {'size': 8, 'opacity': 0.7})

# Setup logging for the libraries used in this script
logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('seaborn').setLevel(logging.ERROR)
logging.getLogger('jinja2').setLevel(logging.ERROR)
logging.getLogger('weasyprint').setLevel(logging.ERROR)
logging.getLogger('fontTools').setLevel(logging.ERROR)

def encode_image_to_base64(image_path):
    logging.info(f"Encoding image {image_path} to base64")
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded_string}"

def find_tsv_files(input_dir, pattern):
    """Recursively find all tsv files matching the pattern in the input directory."""
    tsv_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(pattern):
                tsv_files.append(os.path.join(root, file))
    return tsv_files

def read_compare_outputs(input_dir):
    logging.info(f"Reading compare outputs from {input_dir}")
    reads_assignment_files = find_tsv_files(input_dir, '.reads_assignment.tsv')
    summary_files = find_tsv_files(input_dir, '.summary.tsv')

    if not reads_assignment_files or not summary_files:
        logging.error("No valid TSV files found in the specified directory structure.")
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

def apply_sorting(df, sorting_config):
    """Apply sorting to a DataFrame based on the given sorting configuration."""
    if sorting_config:
        columns = sorting_config.get("columns", [])
        ascending = sorting_config.get("ascending", [True] * len(columns))
        df = df.sort_values(by=columns, ascending=ascending)
    return df

def save_tables_as_tsv_and_excel(combined_df, verdict_df, ratio_df, p_values_df, output_dir):
    """Save the combined, verdict, ratio, and p-value dataframes as TSV and Excel files."""
    data_dir = os.path.join(output_dir, 'data')
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # Apply sorting to each DataFrame
    combined_df = apply_sorting(combined_df, TABLE_SORTING.get('combined'))
    verdict_df = apply_sorting(verdict_df, TABLE_SORTING.get('verdict'))
    ratio_df = apply_sorting(ratio_df, TABLE_SORTING.get('ratio'))
    p_values_df = apply_sorting(p_values_df, TABLE_SORTING.get('p_value'))

    # Save TSV files
    combined_df.to_csv(os.path.join(data_dir, 'combined_table.tsv'), sep='\t', index=False)
    verdict_df.to_csv(os.path.join(data_dir, 'verdict_table.tsv'), sep='\t', index=False)
    ratio_df.to_csv(os.path.join(data_dir, 'ratio_table.tsv'), sep='\t', index=False)
    p_values_df.to_csv(os.path.join(data_dir, 'p_value_table.tsv'), sep='\t', index=False)

    # Save Excel file with multiple sheets
    excel_file = os.path.join(data_dir, 'summary_report_tables.xlsx')
    with pd.ExcelWriter(excel_file, engine='xlsxwriter') as writer:
        combined_df.to_excel(writer, sheet_name='Combined', index=False)
        verdict_df.to_excel(writer, sheet_name='Verdicts', index=False)
        ratio_df.to_excel(writer, sheet_name='Ratios', index=False)
        p_values_df.to_excel(writer, sheet_name='P-Values', index=False)
    
    logging.info("TSV and Excel files have been saved.")

def calculate_and_plot_variations(ratio_df, output_dir):
    """Calculate variations and generate boxplots for contamination ratios by plasmid."""
    logging.info("Calculating variations and generating boxplots.")
    
    # Convert the ratio_df to the appropriate format
    ratio_df['Value'] = pd.to_numeric(ratio_df['Value'], errors='coerce')
    boxplot_data = ratio_df.pivot(index="Sample", columns="Plasmid", values="Value")
    
    # Drop NaNs
    boxplot_data = boxplot_data.dropna(how='all', axis=1)

    # Calculate p-values and FDR correction
    p_values_list = []
    for plasmid in boxplot_data.columns:
        data = boxplot_data[plasmid].dropna()
        if len(data) > 1:
            # Compare each sample to the mean of others
            for i in range(len(data)):
                t_stat, p_value = stats.ttest_1samp(data.drop(data.index[i]), data.iloc[i])
                p_values_list.append({
                    'Sample': data.index[i],
                    'Plasmid': plasmid,
                    'p_value': p_value,
                })

    # Convert p-values list to DataFrame
    p_values_df = pd.DataFrame(p_values_list)

    # Apply FDR correction
    p_values_df['p_value_corrected'] = multipletests(p_values_df['p_value'], method='fdr_bh')[1]

    # Save p-values as a DataFrame in the data directory
    data_dir = os.path.join(output_dir, 'data')
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    
    p_value_table_filename = os.path.join(data_dir, 'p_value_table.tsv')
    p_values_df.to_csv(p_value_table_filename, sep='\t', index=False)

    # Interactive boxplot using Plotly
    boxplot_df = ratio_df.copy()
    boxplot_df['log_value'] = boxplot_df['Value'].apply(lambda x: np.log10(x + LOG_OFFSET))  # Avoid log(0) issues
    fig = px.box(
        boxplot_df, 
        x="Plasmid", 
        y="log_value", 
        points="all", 
        hover_data={'Sample': True, 'Value': True, 'log_value': False},  # Show original Value and hide log_value on hover
        title="Contamination Ratios by Plasmid (Log Scale)"
    )
    fig.update_layout(yaxis_title="Log of Value")
    
    # Adjust the marker style to be similar to your example
    fig.update_traces(marker=MARKER_STYLE)

    plots_dir = os.path.join(output_dir, 'plots')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    boxplot_filename = os.path.join(plots_dir, 'boxplot_contamination_ratios.html')
    fig.write_html(boxplot_filename)

    # Save non-interactive PNG plot
    boxplot_filename_png = os.path.join(plots_dir, 'boxplot_contamination_ratios.png')
    fig.write_image(boxplot_filename_png)
    
    logging.info("Boxplots and statistical calculations completed.")
    
    return boxplot_filename, boxplot_filename_png, p_value_table_filename, p_values_df

def create_plots(reads_df, summary_df, output_dir, threshold, unclear_range, plot_config, substring_to_remove=None):
    logging.info("Creating plots")
    plots_dir = os.path.join(output_dir, 'plots')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
    
    # Split summary dataframe into separate dataframes for different categories
    combined_df = summary_df[summary_df['Category'].isin(['Plasmid', 'Human', 'Tied'])].copy()
    verdict_df = summary_df[summary_df['Category'] == 'Verdict'].copy()
    ratio_df = summary_df[summary_df['Category'] == 'Ratio'].copy()

    # Ensure correct data types for combined and ratio dataframes
    combined_df.loc[:, 'Value'] = pd.to_numeric(combined_df['Value'], errors='coerce')
    ratio_df.loc[:, 'Value'] = pd.to_numeric(ratio_df['Value'], errors='coerce')

    # Remove the specified substring from sample names
    if substring_to_remove:
        ratio_df['Sample'] = ratio_df['Sample'].str.replace(substring_to_remove, '', regex=False)
        combined_df['Sample'] = combined_df['Sample'].str.replace(substring_to_remove, '', regex=False)
        verdict_df['Sample'] = verdict_df['Sample'].str.replace(substring_to_remove, '', regex=False)

    # Round the ratio data to 3 decimal places
    ratio_data = ratio_df.pivot(index="Sample", columns="Plasmid", values="Value").round(ROUND_DECIMALS)
    
    # Interactive heatmap using Plotly
    fig = px.imshow(
        ratio_data, 
        color_continuous_scale=px.colors.sequential.YlGnBu[::-1],  # Inverted color scale
        text_auto=True,
        width=PLOT_DIMENSIONS['width'],  # Set the width of the plot
        height=PLOT_DIMENSIONS['height']   # Set the height of the plot
    )
    fig.update_layout(
        title=plot_config['title'],
        xaxis_title="Plasmid",
        yaxis_title="Sample",
    )

    heatmap_filename_interactive = os.path.join(plots_dir, plot_config['output_filename'].replace('.png', '.html'))
    heatmap_filename_png = os.path.join(plots_dir, plot_config['output_filename'])

    # Save interactive HTML plot
    fig.write_html(heatmap_filename_interactive)

    # Save non-interactive PNG plot
    fig.write_image(heatmap_filename_png)

    return heatmap_filename_interactive, heatmap_filename_png, combined_df, verdict_df, ratio_df

def generate_report(combined_df, verdict_df, ratio_df, heatmap_filename_interactive, heatmap_filename_png, boxplot_filename_interactive, boxplot_filename_png, p_value_table_filename, output_folder, threshold, unclear_range, command_line):
    logging.info("Generating summary reports")

    # Apply sorting to each DataFrame
    combined_df = apply_sorting(combined_df, TABLE_SORTING.get('combined'))
    verdict_df = apply_sorting(verdict_df, TABLE_SORTING.get('verdict'))
    ratio_df = apply_sorting(ratio_df, TABLE_SORTING.get('ratio'))
    p_values_df = pd.read_csv(p_value_table_filename, sep='\t')  # Load the p-value table
    p_values_df = apply_sorting(p_values_df, TABLE_SORTING.get('p_value'))

    env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = env.get_template('summary_template.html')

    # Read the Plotly HTML files content for interactive report
    with open(heatmap_filename_interactive, 'r') as f:
        heatmap_html_content = f.read()

    with open(boxplot_filename_interactive, 'r') as f:
        boxplot_html_content = f.read()

    # Read the PNG files as base64 for non-interactive report
    heatmap_png_base64 = encode_image_to_base64(heatmap_filename_png)
    boxplot_png_base64 = encode_image_to_base64(boxplot_filename_png)

    # Convert sorted DataFrames to HTML
    combined_html = combined_df.to_html(classes='table table-striped table-bordered', index=False)
    verdict_html = verdict_df.to_html(classes='table table-striped table-bordered', index=False)
    ratio_html = ratio_df.to_html(classes='table table-striped table-bordered', index=False)
    p_value_html = p_values_df.to_html(classes='table table-striped table-bordered', index=False)

    # Encode the logo
    logo_base64 = encode_image_to_base64(LOGO_PATH)

    # Render interactive HTML report
    html_content_interactive = template.render(
        combined_df=combined_html,
        verdict_df=verdict_html,
        ratio_df=ratio_html,
        heatmap_content=heatmap_html_content,  # Plotly HTML content
        boxplot_content=boxplot_html_content,  # Plotly HTML content
        p_value_table=p_value_html,
        logo_base64=logo_base64,
        threshold=threshold,
        unclear_range=unclear_range,
        version=VERSION,
        run_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line,
        interactive=True  # Pass True for interactive report
    )

    # Render non-interactive HTML report
    html_content_non_interactive = template.render(
        combined_df=combined_html,
        verdict_df=verdict_html,
        ratio_df=ratio_html,
        heatmap_content=heatmap_png_base64,  # Base64-encoded PNG image
        boxplot_content=boxplot_png_base64,  # Base64-encoded PNG image
        p_value_table=p_value_html,
        logo_base64=logo_base64,
        threshold=threshold,
        unclear_range=unclear_range,
        version=VERSION,
        run_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line,
        interactive=False  # Pass False for non-interactive report
    )

    # Save interactive HTML report
    html_report_interactive = os.path.join(output_folder, 'summary_report_interactive.html')
    with open(html_report_interactive, 'w') as f:
        f.write(html_content_interactive)

    # Save non-interactive HTML report
    html_report_non_interactive = os.path.join(output_folder, 'summary_report_non_interactive.html')
    with open(html_report_non_interactive, 'w') as f:
        f.write(html_content_non_interactive)

    # Convert non-interactive HTML to PDF
    HTML(html_report_non_interactive).write_pdf(os.path.join(output_folder, 'summary_report.pdf'))

def main(input_dir, output_dir, threshold=DEFAULT_THRESHOLD, unclear_range=UNCLEAR_RANGE, substring_to_remove=None):
    reads_df, summary_df = read_compare_outputs(input_dir)
    heatmap_filename_interactive, heatmap_filename_png, combined_df, verdict_df, ratio_df = create_plots(
        reads_df, summary_df, output_dir, threshold, unclear_range, PLOT_CONFIG, substring_to_remove
    )

    # Save the tables as TSV and Excel
    boxplot_filename_interactive, boxplot_filename_png, p_value_table_filename, stats_results = calculate_and_plot_variations(ratio_df, output_dir)

    # Save all tables to both TSV and Excel
    save_tables_as_tsv_and_excel(combined_df, verdict_df, ratio_df, stats_results, output_dir)

    generate_report(
        combined_df, 
        verdict_df, 
        ratio_df, 
        heatmap_filename_interactive, 
        heatmap_filename_png,
        boxplot_filename_interactive, 
        boxplot_filename_png,
        p_value_table_filename,
        output_dir, 
        threshold, 
        unclear_range,
        ' '.join(os.sys.argv)
    )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Merge and generate summary reports for multiple samples and plasmids.")
    parser.add_argument("-i", "--input_dir", help="Directory containing compare outputs", required=True)
    parser.add_argument("-o", "--output_dir", help="Directory to save the plots and reports", required=True)
    parser.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")
    parser.add_argument("--substring_to_remove", help="Substring to remove from sample names", default=None)
    parser.add_argument("--log-level", help="Set the logging level", default="INFO")
    parser.add_argument("--log-file", help="Set the log output file", default=None)

    args = parser.parse_args()

    setup_logging(log_level=args.log_level.upper(), log_file=args.log_file)  # Setup logging with arguments

    main(args.input_dir, args.output_dir, args.threshold, substring_to_remove=args.substring_to_remove)
