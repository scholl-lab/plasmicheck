import subprocess
import json
import os

# Load configuration from JSON file
with open(os.path.join(os.path.dirname(__file__), '..', 'config.json'), 'r') as config_file:
    config = json.load(config_file)

MINIMAP2_THREADS = config['alignment']['minimap2_threads']
SAMTOOLS_THREADS = config['alignment']['samtools_threads']

def align_reads(reference_index, input_file, output_bam, alignment_type, fastq2=None):
    if alignment_type not in ['human', 'plasmid']:
        raise ValueError("alignment_type must be 'human' or 'plasmid'")

    if input_file.endswith('.bam'):
        command = f"samtools fasta {input_file} | minimap2 -t {MINIMAP2_THREADS} -ax sr {reference_index} - | samtools view -@ {SAMTOOLS_THREADS} -h -F 4 - | samtools sort -@ {SAMTOOLS_THREADS} -o {output_bam}"
    elif input_file.endswith('.fastq'):
        if fastq2:
            command = f"minimap2 -t {MINIMAP2_THREADS} -ax sr {reference_index} {input_file} {fastq2} | samtools view -@ {SAMTOOLS_THREADS} -h -F 4 - | samtools sort -@ {SAMTOOLS_THREADS} -o {output_bam}"
        else:
            command = f"minimap2 -t {MINIMAP2_THREADS} -ax sr {reference_index} {input_file} | samtools view -@ {SAMTOOLS_THREADS} -h -F 4 - | samtools sort -@ {SAMTOOLS_THREADS} -o {output_bam}"
    else:
        raise ValueError("Unsupported input file type. Must be .bam, .fastq, or paired FASTQ files.")

    # Run the alignment command
    subprocess.run(command, shell=True, check=True)

    # Generate the BAI index
    subprocess.run(f"samtools index {output_bam}", shell=True, check=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Align reads to a reference genome and generate a BAI index")
    parser.add_argument("-r", "--reference_index", help="Minimap2 index for the reference genome", required=True)
    parser.add_argument("-i", "--input_file", help="Input file (BAM, interleaved FASTQ, or first FASTQ file for paired FASTQ)", required=True)
    parser.add_argument("-o", "--output_bam", help="Output BAM file for alignment", required=True)
    parser.add_argument("-a", "--alignment_type", help="Type of alignment: 'human' or 'plasmid'", required=True)
    parser.add_argument("--fastq2", help="Second FASTQ file for paired FASTQ input", default=None)
    args = parser.parse_args()

    align_reads(args.reference_index, args.input_file, args.output_bam, args.alignment_type, args.fastq2)
