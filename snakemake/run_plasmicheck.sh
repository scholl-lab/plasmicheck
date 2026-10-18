#!/bin/bash
#
#SBATCH --job-name=sm_plasmicheck_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=1200M
#SBATCH --output=slurm_logs/%x-%j.log

# based on:
# https://hpc-docs.cubi.bihealth.org/best-practice/temp-files/#tmpdir-and-the-scheduler
# https://bihealth.github.io/bih-cluster/slurm/snakemake/#custom-logging-directory

set -euo pipefail

export TMPDIR="${HOME}/scratch/tmp"
mkdir -p "$TMPDIR"
export TMPDIR="$(mktemp -d -p "$TMPDIR")"
trap 'rm -rf "$TMPDIR"' EXIT

mkdir -p slurm_logs
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

CONFIGFILE="${1:-config.example.yaml}"
PROFILE="${SNAKEMAKE_PROFILE:-cubi-v1}"
JOBS="${SNAKEMAKE_JOBS:-100}"

CONFIG_ARGS=()
if [ -f "$CONFIGFILE" ]; then
    CONFIG_ARGS=(--configfile "$CONFIGFILE")
fi

date
srun snakemake -s plasmicheck.smk "${CONFIG_ARGS[@]}" --use-conda --profile="$PROFILE" -j"$JOBS"
date
