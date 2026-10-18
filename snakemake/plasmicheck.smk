import itertools
import os

from plasmicheck.scripts.utils import sanitize_filename

# ----------------------------------------------------------------------------------- #
# Load BAM and plasmid file paths
with open("data/bams.txt") as f:
    bam_files = [line.strip() for line in f if line.strip()]
with open("data/plasmids.txt") as f:
    plasmid_files = [line.strip() for line in f if line.strip()]

# Generate all combinations of BAM and plasmid files and store them in a dictionary
jobs = {}
for bam, plasmid in itertools.product(bam_files, plasmid_files):
    bam_basename = sanitize_filename(os.path.splitext(os.path.basename(bam))[0])
    plasmid_basename = sanitize_filename(os.path.splitext(os.path.basename(plasmid))[0])
    key = f"{bam_basename}__{plasmid_basename}"
    if key in jobs:
        raise ValueError(f"Duplicate job key: {key}")
    log_file = f"logs/{key}.log"
    jobs[key] = {
        "bam": bam,
        "plasmid": plasmid,
        "bam_basename": bam_basename,
        "plasmid_basename": plasmid_basename,
        "log": log_file,
    }

# ----------------------------------------------------------------------------------- #
# Helper functions
def get_mem_from_threads(wildcards, threads):
    """Calculate memory allocation based on the number of threads."""
    return threads * 2200  # 2.2 GB per thread

# ----------------------------------------------------------------------------------- #
# Define the rules
rule all:
    input:
        expand("logs/{bam_basename}__{plasmid_basename}.log",
               zip,
               bam_basename = [jobs[key]['bam_basename'] for key in jobs.keys()],
               plasmid_basename = [jobs[key]['plasmid_basename'] for key in jobs.keys()])

rule run_plasmicheck_pipeline:
    input:
        bam = lambda wildcards: f"{jobs[wildcards.bam_basename + '__' + wildcards.plasmid_basename]['bam']}",
        plasmid = lambda wildcards: f"{jobs[wildcards.bam_basename + '__' + wildcards.plasmid_basename]['plasmid']}",
    output:
        log = "logs/{bam_basename}__{plasmid_basename}.log",
    threads: 8
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = os.environ.get('TMPDIR', '/tmp')
    conda:
        "plasmicheck"
    shell:
        """
        plasmicheck pipeline -hf reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
            -pf {input.plasmid:q} -sf1 {input.bam:q} \
            -o output --threads {threads} --log-file {output.log:q} --archive_output
        """
