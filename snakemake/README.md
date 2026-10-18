# PlasmiCheck Snakemake Cohort Workflow

Config-driven Snakemake workflow (`snakemake/plasmicheck.smk`) for screening a cohort of coordinate-sorted BAMs (e.g., 97 BQSR exome BAMs) against a panel of plasmid constructs (e.g., 21 constructs) on SLURM without creating a `BAM × plasmid` job cross-product.

## Files

| File | Description |
| :--- | :--- |
| `plasmicheck.smk` | 3-rule cohort workflow (`prepare_constructs` → `scan_sample` → `summary_reports`) with legacy fallback to `data/bams.txt` and `data/plasmids.txt` |
| `config.example.yaml` | Documented example configuration for SLURM cohort execution |
| `config.schema.yaml` | JSON Schema (Draft 2020-12) for workflow configuration and per-rule resource bounds (`runtime <= 120` min) |
| `run_plasmicheck.sh` | SLURM submission wrapper script |

## 3-Rule DAG Architecture

1. **`prepare_constructs` (1 job):** Indexes `human_fasta` and populates the content-addressed `CacheManager` (`cache_dir`) with all plasmid construct indexes, spliced alignments, and spliced human reference indexes once.
2. **`scan_sample` (1 job per sample BAM):** Runs `plasmicheck pipeline -hf <human_fasta> -pf <plasmids_file> -sf1 <bam> --sample-name <sample> -o <output_dir> --threads <threads>` once per sample, extracting the BAM once and screening all constructs in a single pipeline run (<2 h walltime) before verifying all expected per-construct outputs.
3. **`summary_reports` (1 job):** Aggregates the explicit manifest of configured `(sample, construct)` comparison outputs via `plasmicheck summary_reports -i <manifest> -o <summary_dir>`.

## Usage

```bash
# Dry-run with example configuration
snakemake -s snakemake/plasmicheck.smk --configfile snakemake/config.example.yaml -n

# Run on a SLURM cluster with a Snakemake profile
snakemake -s snakemake/plasmicheck.smk \
  --configfile snakemake/config.example.yaml \
  --use-conda \
  --profile cubi-v1 \
  -j 100
```

## Configuration Reference (`config.schema.yaml`)

- **`samples`:** Sample sheet TSV (`sample\tbam` or `corrected_id` with `bam_folder`), mapping `{sample_id: bam_path}`, list of BAM paths, or newline-delimited BAM list file (fallback: `data/bams.txt`).
- **`plasmids_file` / `constructs`:** Path to newline-delimited construct list (fallback: `data/plasmids.txt`) or explicit list of construct files (`.gb`, `.xdna`, `.dna`, `.fasta`).
- **`human_fasta`:** Path to the human reference genome FASTA.
- **Optional flags forwarded to `plasmicheck pipeline`:** `regions_bed`, `include_unmapped` (default `true`), `threshold`, `padding`, `static_report`, `archive_output`, `overwrite`, `extra_args`.
- **`resources`:** Per-rule overrides (`prepare_constructs`, `scan_sample`, `summary_reports`) for `threads`, `mem_mb`, and `runtime` (minutes, max `120`).
