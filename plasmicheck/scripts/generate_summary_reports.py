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

# Import the version from version.py
from plasmicheck.version import __version__ as VERSION

from .utils import setup_logging  # Import setup_logging function

# Load configuration from JSON file
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

# Constants from config
DEFAULT_THRESHOLD = config['default_threshold']
UNCLEAR_RANGE = config['unclear_range']
PLOT_CONFIG = config['plot_summary']
TEMPLATE_DIR = config['paths']['template_dir']
LOGO_PATH = config['paths']['logo_path']
TABLE_SORTING = config['table_sorting']
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


# Helper Functions
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


def read_compare_outputs(input_dir, substring_to_remove=None):
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
        
        # Remove the substring from the sample name if provided
        if substring_to_remove:
            sample_name = sample_name.replace(substring_to_remove, '')

        df = pd.read_csv(file, sep='\t', low_memory=False)  # Disable low_memory mode to avoid DtypeWarning
        df['Sample'] = sample_name
        df['Plasmid'] = plasmid_name
        reads_df_list.append(df)

    summary_df_list = []
    for file in summary_files:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(file)))
        plasmid_name = os.path.basename(os.path.dirname(file))
        
        # Remove the substring from the sample name if provided
        if substring_to_remove:
            sample_name = sample_name.replace(substring_to_remove, '')

        df = pd.read_csv(file, sep='\t', low_memory=False)  # Disable low_memory mode to avoid DtypeWarning
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


# Data Manipulation Functions
def calculate_variations(ratio_df):
    """Calculate p-values and apply FDR correction for contamination ratios by plasmid."""
    logging.info("Calculating variations for contamination ratios.")

    # Convert 'Value' to numeric, ensuring all entries are valid numbers
    ratio_df.loc[:, 'Value'] = pd.to_numeric(ratio_df['Value'], errors='coerce')
    boxplot_data = ratio_df.pivot(index="Sample", columns="Plasmid", values="Value")

    # Drop NaNs
    boxplot_data = boxplot_data.dropna(how='all', axis=1)

    p_values_list = []
    for plasmid in boxplot_data.columns:
        # Ensure the data is numeric and drop NaN values
        data = pd.to_numeric(boxplot_data[plasmid], errors='coerce').dropna()

        # Only run the t-test if there are enough data points
        if len(data) > 1:
            for i in range(len(data)):
                # Exclude the current sample and run the one-sample t-test
                try:
                    t_stat, p_value = stats.ttest_1samp(data.drop(data.index[i]), data.iloc[i])
                    p_values_list.append({
                        'Sample': data.index[i],
                        'Plasmid': plasmid,
                        'p_value': p_value,
                    })
                except Exception as e:
                    logging.error(f"Error during t-test for {plasmid}, sample {data.index[i]}: {e}")

    p_values_df = pd.DataFrame(p_values_list)
    if not p_values_df.empty:
        p_values_df['p_value_corrected'] = multipletests(p_values_df['p_value'], method='fdr_bh')[1]

    return boxplot_data, p_values_df

# Plotting Functions
def plot_boxplot(boxplot_data, output_dir):
    """Generate and save boxplots for contamination ratios by plasmid."""
    logging.info("Generating boxplots for contamination ratios.")

    boxplot_df = pd.melt(boxplot_data.reset_index(), id_vars=['Sample'], var_name='Plasmid', value_name='Value')
    boxplot_df['log_value'] = boxplot_df['Value'].apply(lambda x: np.log10(x + LOG_OFFSET))

    fig = px.box(
        boxplot_df, 
        x="Plasmid", 
        y="log_value", 
        points="all", 
        hover_data={'Sample': True, 'Value': True, 'log_value': False},
        title="Contamination Ratios by Plasmid (Log Scale)"
    )
    
    # Increase the height of the boxplot
    fig.update_layout(
        yaxis_title="Log of Value",
        height=800  # Adjust this value to make the boxplot taller
    )
    
    fig.update_traces(marker=MARKER_STYLE)

    plots_dir = os.path.join(output_dir, 'plots')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    boxplot_filename = os.path.join(plots_dir, 'boxplot_contamination_ratios.html')
    fig.write_html(boxplot_filename)

    boxplot_filename_png = os.path.join(plots_dir, 'boxplot_contamination_ratios.png')
    fig.write_image(boxplot_filename_png)
    
    logging.info("Boxplots generated.")
    
    return boxplot_filename, boxplot_filename_png

def plot_heatmap(ratio_df, output_dir, threshold, unclear_range, plot_config):
    """Generate and save heatmap for contamination ratios by plasmid without the color bar."""
    logging.info("Generating heatmap for contamination ratios.")

    # Use .loc[] to avoid SettingWithCopyWarning
    ratio_df.loc[:, 'Value'] = pd.to_numeric(ratio_df['Value'], errors='coerce')
    
    # Pivot the DataFrame for the heatmap
    ratio_data = ratio_df.pivot(index="Sample", columns="Plasmid", values="Value").round(ROUND_DECIMALS)

    # Define a custom color scale using Plotly's built-in scaling
    color_scale = [
        [0.0, plot_config['colors']['not_contaminated']],  # Start
        [unclear_range['lower_bound'] / ratio_data.max().max(), plot_config['colors']['unclear']],  # Lower bound
        [unclear_range['upper_bound'] / ratio_data.max().max(), plot_config['colors']['unclear']],  # Upper bound
        [1.0, plot_config['colors']['contaminated']]  # End
    ]

    fig = px.imshow(
        ratio_data,
        color_continuous_scale=color_scale,
        height=PLOT_DIMENSIONS['height'],
        aspect='auto'  # Adjust the aspect ratio to fit the plot dimensions
    )
    
    # Update layout for axis labels and title
    fig.update_layout(
        title=plot_config['title'],
        xaxis_title="Plasmid",
        yaxis_title="Sample",
        xaxis=dict(tickangle=plot_config['xticks_rotation'], tickfont=dict(size=9), tickmode='linear'),
        yaxis=dict(tickfont=dict(size=9))
    )
    
    # Disable the color bar
    fig.update_coloraxes(showscale=False)

    plots_dir = os.path.join(output_dir, 'plots')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    heatmap_filename_interactive = os.path.join(plots_dir, plot_config['output_filename'].replace('.png', '.html'))
    fig.write_html(heatmap_filename_interactive)

    heatmap_filename_png = os.path.join(plots_dir, plot_config['output_filename'])
    fig.write_image(heatmap_filename_png)

    logging.info("Heatmap generated.")
    
    return heatmap_filename_interactive, heatmap_filename_png

# Report Generation Functions
def generate_report(combined_df, verdict_df, ratio_df, heatmap_filename_interactive, heatmap_filename_png, boxplot_filename_interactive, boxplot_filename_png, p_value_table_filename, output_folder, threshold, unclear_range, command_line):
    logging.info("Generating summary reports")

    # Apply sorting to each DataFrame
    combined_df = apply_sorting(combined_df, TABLE_SORTING.get('combined'))
    verdict_df = apply_sorting(verdict_df, TABLE_SORTING.get('verdict'))
    ratio_df = apply_sorting(ratio_df, TABLE_SORTING.get('ratio'))
    p_values_df = pd.read_csv(p_value_table_filename, sep='\t')  
    p_values_df = apply_sorting(p_values_df, TABLE_SORTING.get('p_value'))

    env = Environment(loader=FileSystemLoader(os.path.join('plasmicheck', TEMPLATE_DIR)))
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
    logo_base64 = encode_image_to_base64(os.path.join('plasmicheck', LOGO_PATH))

    # Render interactive HTML report
    html_content_interactive = template.render(
        combined_df=combined_html,
        verdict_df=verdict_html,
        ratio_df=ratio_html,
        heatmap_content=heatmap_html_content,
        boxplot_content=boxplot_html_content,
        p_value_table=p_value_html,
        logo_base64=logo_base64,
        threshold=threshold,
        unclear_range=unclear_range,
        version=VERSION,
        run_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line,
        interactive=True  
    )

    # Render non-interactive HTML report
    html_content_non_interactive = template.render(
        combined_df=combined_html,
        verdict_df=verdict_html,
        ratio_df=ratio_html,
        heatmap_content=heatmap_png_base64,
        boxplot_content=boxplot_png_base64,
        p_value_table=p_value_html,
        logo_base64=logo_base64,
        threshold=threshold,
        unclear_range=unclear_range,
        version=VERSION,
        run_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line,
        interactive=False  
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
    reads_df, summary_df = read_compare_outputs(input_dir, substring_to_remove=substring_to_remove)

    # Calculate variations
    boxplot_data, p_values_df = calculate_variations(summary_df[summary_df['Category'] == 'Ratio'])

    # Generate plots
    boxplot_filename_interactive, boxplot_filename_png = plot_boxplot(boxplot_data, output_dir)
    heatmap_filename_interactive, heatmap_filename_png = plot_heatmap(summary_df[summary_df['Category'] == 'Ratio'], output_dir, threshold, unclear_range, PLOT_CONFIG)

    # Save tables
    save_tables_as_tsv_and_excel(summary_df[summary_df['Category'].isin(['Plasmid', 'Human', 'Tied'])], 
                                 summary_df[summary_df['Category'] == 'Verdict'], 
                                 summary_df[summary_df['Category'] == 'Ratio'], 
                                 p_values_df, output_dir)

    # Generate report
    generate_report(summary_df[summary_df['Category'].isin(['Plasmid', 'Human', 'Tied'])], 
                    summary_df[summary_df['Category'] == 'Verdict'], 
                    summary_df[summary_df['Category'] == 'Ratio'], 
                    heatmap_filename_interactive, heatmap_filename_png, 
                    boxplot_filename_interactive, boxplot_filename_png, 
                    os.path.join(output_dir, 'data', 'p_value_table.tsv'), 
                    output_dir, threshold, unclear_range, ' '.join(os.sys.argv))


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

    main(args.input_dir, args.output_dir, args.threshold, unclear_range=UNCLEAR_RANGE, substring_to_remove=args.substring_to_remove)
