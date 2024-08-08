import pysam
import json
import os

# Resolve the path to config.json in the parent directory of the current script
config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.json')

# Load configuration from JSON file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

# Configuration variables for scoring
MATE_BONUS = config['scoring']['mate_bonus']
CLIPPING_PENALTY = config['scoring']['clipping_penalty']
MISMATCH_PENALTY = config['scoring']['mismatch_penalty']
DEFAULT_THRESHOLD = config['default_threshold']

def calculate_alignment_score(read):
    """
    Calculate a custom alignment score based on multiple criteria.
    """
    if read.is_unmapped:
        return 0

    score = read.mapping_quality

    # Penalize for clipping
    for operation, length in read.cigartuples:
        if operation in {4, 5}:  # 4: soft clipping, 5: hard clipping
            score -= length * CLIPPING_PENALTY

    # Penalize for mismatches
    mismatches = read.get_tag("NM") if read.has_tag("NM") else 0
    score -= mismatches * MISMATCH_PENALTY

    # Bonus if mate is mapped
    if not read.mate_is_unmapped:
        score += MATE_BONUS

    return score

def compare_alignments(plasmid_bam, human_bam, output_basename, threshold=DEFAULT_THRESHOLD):
    plasmid_samfile = pysam.AlignmentFile(plasmid_bam, "rb")
    human_samfile = pysam.AlignmentFile(human_bam, "rb")

    # Create dictionaries to store reads by query name
    plasmid_reads = {read.query_name: read for read in plasmid_samfile.fetch(until_eof=True)}
    human_reads = {read.query_name: read for read in human_samfile.fetch(until_eof=True)}

    assigned_counts = {"Plasmid": 0, "Human": 0, "Tied": 0}

    with open(f"{output_basename}.reads_assignment.tsv", 'w') as outfile:
        outfile.write("ReadID\tAssignedTo\tPlasmidScore\tHumanScore\n")
        for query_name, plasmid_read in plasmid_reads.items():
            plasmid_score = calculate_alignment_score(plasmid_read)
            human_read = human_reads.get(query_name)
            human_score = calculate_alignment_score(human_read) if human_read else 0

            if plasmid_score > human_score:
                assigned_to = "Plasmid"
                assigned_counts["Plasmid"] += 1
            elif human_score > plasmid_score:
                assigned_to = "Human"
                assigned_counts["Human"] += 1
            else:
                assigned_to = "Tied"
                assigned_counts["Tied"] += 1

            outfile.write(f"{query_name}\t{assigned_to}\t{plasmid_score}\t{human_score}\n")

    plasmid_count = assigned_counts["Plasmid"]
    human_count = assigned_counts["Human"]
    ratio = plasmid_count / human_count if human_count != 0 else float('inf')
    verdict = "Sample is contaminated with plasmid DNA" if ratio > threshold else "Sample is not contaminated with plasmid DNA"

    with open(f"{output_basename}.summary.tsv", 'w') as summary_file:
        summary_file.write("Category\tCount\n")
        for category, count in assigned_counts.items():
            summary_file.write(f"{category}\t{count}\n")
        summary_file.write(f"Verdict\t{verdict}\n")
        summary_file.write(f"Ratio\t{ratio}\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compare alignments and assign reads to plasmid or human reference")
    parser.add_argument("plasmid_bam", help="BAM file for plasmid alignment")
    parser.add_argument("human_bam", help="BAM file for human alignment")
    parser.add_argument("output_basename", help="Basename for output files")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})")
    args = parser.parse_args()

    compare_alignments(args.plasmid_bam, args.human_bam, args.output_basename, args.threshold)
