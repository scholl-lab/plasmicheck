# Snakemake cohort workflow for PlasmiCheck plasmid contamination detection (#129).
# Extracts each sample BAM once and aligns against all plasmid constructs in a single pass.

import csv
import os
import shlex
from typing import Any

from plasmicheck.scripts.utils import sanitize_filename

DEFAULT_HUMAN_FASTA = "reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
DEFAULT_BAMS_FILE = "data/bams.txt"
DEFAULT_PLASMIDS_FILE = "data/plasmids.txt"

DEFAULT_RESOURCES: dict[str, dict[str, int]] = {
    "prepare_constructs": {"threads": 4, "mem_mb": 16000, "runtime": 60},
    "scan_sample": {"threads": 4, "mem_mb": 8000, "runtime": 60},
    "summary_reports": {"threads": 2, "mem_mb": 4000, "runtime": 30},
}


def get_mem_from_threads(wildcards: Any, threads: int) -> int:
    """Calculate default memory allocation (MB) based on thread count."""
    return int(threads) * 2200


def validate_config(cfg: dict[str, Any]) -> None:
    """Validate workflow configuration dictionary and resource walltime limits."""
    if not isinstance(cfg, dict):
        raise ValueError("Configuration must be a dictionary.")

    samples_val = cfg.get("samples")
    if samples_val is not None:
        if isinstance(samples_val, (dict, list, str)):
            if not samples_val:
                raise ValueError("'samples' must not be empty.")
        else:
            raise ValueError("'samples' must be a dict, list, or file path string.")

    constructs_val = cfg.get("constructs")
    if constructs_val is not None:
        if isinstance(constructs_val, (list, str)):
            if not constructs_val:
                raise ValueError("'constructs' must not be empty.")
        else:
            raise ValueError("'constructs' must be a list of paths or a file path string.")

    threshold_val = cfg.get("threshold")
    if threshold_val is not None and (
        isinstance(threshold_val, bool)
        or not isinstance(threshold_val, (int, float))
        or threshold_val <= 0
    ):
        raise ValueError("'threshold' must be a positive number.")

    padding_val = cfg.get("padding")
    if padding_val is not None and (
        isinstance(padding_val, bool) or not isinstance(padding_val, int) or padding_val < 0
    ):
        raise ValueError("'padding' must be a non-negative integer.")

    resources_val = cfg.get("resources")
    if resources_val is not None:
        if not isinstance(resources_val, dict):
            raise ValueError("'resources' must be a dictionary.")
        for rule_name, rule_res in resources_val.items():
            if rule_name not in DEFAULT_RESOURCES:
                raise ValueError(f"Unknown rule in 'resources': {rule_name}")
            if not isinstance(rule_res, dict):
                raise ValueError(f"Resources for '{rule_name}' must be a dictionary.")
            for key, val in rule_res.items():
                if key not in ("threads", "mem_mb", "runtime"):
                    raise ValueError(f"Unknown resource key '{key}' for rule '{rule_name}'.")
                if isinstance(val, bool) or not isinstance(val, int) or val <= 0:
                    raise ValueError(f"Resource '{rule_name}.{key}' must be a positive integer.")
                if key == "runtime" and val > 120:
                    raise ValueError(
                        f"Resource '{rule_name}.runtime' ({val} min) exceeds the 120 min (<2 h) walltime limit."
                    )


def _read_lines(file_path: str) -> list[str]:
    """Read non-empty, non-comment lines from a text file."""
    with open(file_path, encoding="utf-8") as fh:
        return [
            line.strip()
            for line in fh
            if line.strip() and not line.strip().startswith("#")
        ]


def _bams_list_to_sample_map(bam_files: list[str]) -> dict[str, str]:
    """Convert a list of BAM file paths to {sanitized_sample_id: bam_path}."""
    if not bam_files:
        raise ValueError("No BAM files provided.")
    samples: dict[str, str] = {}
    for bam in bam_files:
        sample_id = sanitize_filename(os.path.splitext(os.path.basename(bam))[0])
        if sample_id in samples:
            raise ValueError(f"Duplicate sample key: {sample_id}")
        samples[sample_id] = bam
    return samples


def _parse_samples_tsv(tsv_path: str, cfg: dict[str, Any]) -> dict[str, str]:
    """Parse a TSV sample sheet ('sample\\tbam' or 'corrected_id' + bam_folder)."""
    lines = _read_lines(tsv_path)
    if not lines:
        raise ValueError(f"Sample sheet is empty: {tsv_path}")

    reader = list(csv.reader(lines, delimiter="\t"))
    first_row = [col.strip() for col in reader[0]]
    lower_first = [col.lower() for col in first_row]

    sample_cols = ("sample", "sample_id", "corrected_id")
    bam_cols = ("bam", "bam_path", "bam_file")
    has_header = any(c in lower_first for c in (*sample_cols, *bam_cols))

    samples: dict[str, str] = {}
    if has_header:
        sample_idx = next((lower_first.index(c) for c in sample_cols if c in lower_first), None)
        bam_idx = next((lower_first.index(c) for c in bam_cols if c in lower_first), None)
        if sample_idx is None:
            raise ValueError(
                f"Sample sheet {tsv_path} must contain a 'sample' or 'corrected_id' column."
            )
        bam_folder = cfg.get("bam_folder") or cfg.get("bam_dir")
        bam_ext = str(cfg.get("bam_ext", ".merged.dedup.bqsr.bam"))
        if bam_idx is None and not bam_folder:
            raise ValueError(
                f"Sample sheet {tsv_path} has no 'bam' column; set 'bam_folder' in config."
            )
        for row in reader[1:]:
            if not row or len(row) <= sample_idx:
                continue
            sample_id = row[sample_idx].strip()
            if not sample_id:
                continue
            if bam_idx is not None:
                if len(row) <= bam_idx or not row[bam_idx].strip():
                    raise ValueError(f"Missing BAM path for sample '{sample_id}' in {tsv_path}")
                bam_path = row[bam_idx].strip()
            else:
                bam_path = os.path.join(str(bam_folder), f"{sample_id}{bam_ext}")
            if sample_id in samples and samples[sample_id] != bam_path:
                raise ValueError(f"Duplicate sample key: {sample_id}")
            samples[sample_id] = bam_path
    else:
        for row in reader:
            if len(row) < 2:
                raise ValueError(
                    f"Expected 'sample\\tbam' columns in {tsv_path}, got: {row!r}"
                )
            sample_id, bam_path = row[0].strip(), row[1].strip()
            if sample_id in samples:
                raise ValueError(f"Duplicate sample key: {sample_id}")
            samples[sample_id] = bam_path

    if not samples:
        raise ValueError(f"No samples resolved from {tsv_path}")
    return samples


def resolve_samples(cfg: dict[str, Any]) -> dict[str, str]:
    """Resolve {sample_id: bam_path} from config dict, TSV sheet, or BAM list file."""
    samples_src = cfg.get("samples")
    if samples_src is None:
        samples_src = cfg.get("bams_file", DEFAULT_BAMS_FILE)

    if isinstance(samples_src, dict):
        if not samples_src:
            raise ValueError("No samples configured.")
        return {str(k).strip(): str(v).strip() for k, v in samples_src.items()}

    if isinstance(samples_src, list):
        return _bams_list_to_sample_map([str(x).strip() for x in samples_src if str(x).strip()])

    if isinstance(samples_src, str):
        if samples_src.endswith(".tsv"):
            return _parse_samples_tsv(samples_src, cfg)
        lines = _read_lines(samples_src)
        if lines and "\t" in lines[0]:
            return _parse_samples_tsv(samples_src, cfg)
        return _bams_list_to_sample_map(lines)

    raise ValueError(f"Unsupported 'samples' configuration: {samples_src!r}")


def resolve_plasmids(cfg: dict[str, Any], output_dir: str) -> tuple[str, list[str]]:
    """Resolve (plasmids_file_path, construct_files_list) from config or legacy fallback."""
    plasmids_file = cfg.get("plasmids_file")
    constructs = cfg.get("constructs")

    if plasmids_file:
        pf = str(plasmids_file)
        files = _read_lines(pf) if (os.path.isfile(pf) and pf.endswith(".txt")) else [pf]
        return pf, files

    if constructs is not None:
        if isinstance(constructs, str):
            pf = constructs
            files = _read_lines(pf) if (os.path.isfile(pf) and pf.endswith(".txt")) else [pf]
            return pf, files
        if isinstance(constructs, list):
            files = [str(x).strip() for x in constructs if str(x).strip()]
            if not files:
                raise ValueError("'constructs' list must not be empty.")
            os.makedirs(output_dir, exist_ok=True)
            pf = os.path.join(output_dir, ".plasmids.txt")
            with open(pf, "w", encoding="utf-8") as fh:
                fh.write("\n".join(files) + "\n")
            return pf, files

    pf = DEFAULT_PLASMIDS_FILE
    files = _read_lines(pf) if os.path.isfile(pf) else []
    return pf, files


def get_rule_resource(cfg: dict[str, Any], rule_name: str, key: str) -> int:
    """Return per-rule resource value ('threads', 'mem_mb', 'runtime') with config override."""
    defaults = DEFAULT_RESOURCES[rule_name]
    rule_cfg = cfg.get("resources", {}).get(rule_name, {})
    return int(rule_cfg.get(key, defaults[key]))


def build_pipeline_flags(cfg: dict[str, Any]) -> str:
    """Build optional CLI flags for 'plasmicheck pipeline' from config."""
    flags: list[str] = []

    if cfg.get("cache_dir"):
        flags.append(f"--cache-dir {shlex.quote(str(cfg['cache_dir']))}")

    regions_bed = cfg.get("regions_bed")
    if regions_bed:
        flags.append(f"--regions {shlex.quote(str(regions_bed))}")
        if cfg.get("include_unmapped", True):
            flags.append("--include-unmapped")
    elif cfg.get("include_unmapped") is True and "regions_bed" not in cfg:
        flags.append("--include-unmapped")

    if cfg.get("target_bed"):
        flags.append(f"--target-bed {shlex.quote(str(cfg['target_bed']))}")
    if cfg.get("gtf"):
        flags.append(f"--gtf {shlex.quote(str(cfg['gtf']))}")
    if cfg.get("vcf"):
        flags.append(f"--vcf {shlex.quote(str(cfg['vcf']))}")
    if cfg.get("threshold") is not None:
        flags.append(f"--threshold {shlex.quote(str(cfg['threshold']))}")
    if cfg.get("padding") is not None:
        flags.append(f"--padding {shlex.quote(str(cfg['padding']))}")
    if cfg.get("static_report"):
        flags.append("--static-report")
    if cfg.get("archive_output"):
        flags.append("--archive_output")
    if cfg.get("overwrite"):
        flags.append("--overwrite")
    if cfg.get("extra_args"):
        flags.append(str(cfg["extra_args"]).strip())

    return " ".join(flags)


def build_summary_flags(cfg: dict[str, Any]) -> str:
    """Build optional CLI flags for 'plasmicheck summary_reports' from config."""
    flags: list[str] = []
    if cfg.get("threshold") is not None:
        flags.append(f"--threshold {shlex.quote(str(cfg['threshold']))}")
    if cfg.get("static_report"):
        flags.append("--static-report")
    if cfg.get("substring_to_remove"):
        flags.append(f"--substring_to_remove {shlex.quote(str(cfg['substring_to_remove']))}")
    return " ".join(flags)


# ----------------------------------------------------------------------------------- #
# Resolve configuration and inputs
# ----------------------------------------------------------------------------------- #

config = globals().get("config", {})
validate_config(config)

HUMAN_FASTA = str(config.get("human_fasta", DEFAULT_HUMAN_FASTA))
OUTPUT_DIR = str(config.get("output_dir", "output"))
SUMMARY_DIR = str(config.get("summary_dir", os.path.join(OUTPUT_DIR, "summary")))
LOG_DIR = str(config.get("log_dir", "logs"))
CACHE_DIR = str(config.get("cache_dir", os.path.expanduser("~/.plasmicheck/cache")))
CONDA_ENV = str(config.get("conda_env", "plasmicheck"))

SAMPLE_BAMS = resolve_samples(config)
SAMPLES = list(SAMPLE_BAMS.keys())
PLASMIDS_FILE, CONSTRUCT_FILES = resolve_plasmids(config, OUTPUT_DIR)
bam_files = list(SAMPLE_BAMS.values())
plasmid_files = CONSTRUCT_FILES

PIPELINE_FLAGS = build_pipeline_flags(config)
SUMMARY_FLAGS = build_summary_flags(config)


# ----------------------------------------------------------------------------------- #
# Rules (3-rule cohort DAG: prepare_constructs -> scan_sample -> summary_reports)
# ----------------------------------------------------------------------------------- #

rule all:
    """Default cohort target: per-sample contamination scans and cohort summary reports."""
    input:
        done=expand(f"{LOG_DIR}/{{sample}}.done", sample=SAMPLES),
        summary_done=f"{LOG_DIR}/summary_reports.done",
        verdict_tsv=f"{SUMMARY_DIR}/data/verdict_table.tsv",
        ratio_tsv=f"{SUMMARY_DIR}/data/ratio_table.tsv",
        report_html=f"{SUMMARY_DIR}/summary_report_interactive.html",


rule prepare_constructs:
    """Build human reference index and stage plasmid construct references once."""
    input:
        human_fasta=HUMAN_FASTA,
        plasmids_file=PLASMIDS_FILE,
        constructs=CONSTRUCT_FILES,
    output:
        done=f"{LOG_DIR}/prepare_constructs.done",
    params:
        cache_dir=CACHE_DIR,
    threads: get_rule_resource(config, "prepare_constructs", "threads")
    resources:
        mem_mb=get_rule_resource(config, "prepare_constructs", "mem_mb"),
        runtime=get_rule_resource(config, "prepare_constructs", "runtime"),
        tmpdir=os.environ.get("TMPDIR", "/tmp"),
    conda:
        CONDA_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output.done:q})" {params.cache_dir:q}
        plasmicheck index -f {input.human_fasta:q}
        if [ -n "{input.constructs}" ]; then
            for plasmid in {input.constructs:q}; do
                name="$(basename "${{plasmid%.*}}")"
                plasmicheck convert -i "$plasmid" -o "{params.cache_dir:q}/${{name}}.fasta" -w
                plasmicheck index -f "{params.cache_dir:q}/${{name}}.fasta" -w
            done
        elif [ -f {input.plasmids_file:q} ]; then
            while IFS= read -r plasmid || [ -n "$plasmid" ]; do
                case "$plasmid" in ""|\#*) continue ;; esac
                name="$(basename "${{plasmid%.*}}")"
                plasmicheck convert -i "$plasmid" -o "{params.cache_dir:q}/${{name}}.fasta" -w
                plasmicheck index -f "{params.cache_dir:q}/${{name}}.fasta" -w
            done < {input.plasmids_file:q}
        fi
        touch {output.done:q}
        """


rule scan_sample:
    """Screen a single sample BAM against all plasmid constructs in one pass."""
    input:
        bam=lambda wildcards: SAMPLE_BAMS[wildcards.sample],
        human_fasta=HUMAN_FASTA,
        plasmids_file=PLASMIDS_FILE,
        prepared=f"{LOG_DIR}/prepare_constructs.done",
    output:
        done=f"{LOG_DIR}/{{sample}}.done",
        log=f"{LOG_DIR}/{{sample}}.log",
    params:
        output_dir=OUTPUT_DIR,
        extra_flags=PIPELINE_FLAGS,
    threads: get_rule_resource(config, "scan_sample", "threads")
    resources:
        mem_mb=get_rule_resource(config, "scan_sample", "mem_mb"),
        runtime=get_rule_resource(config, "scan_sample", "runtime"),
        tmpdir=os.environ.get("TMPDIR", "/tmp"),
    conda:
        CONDA_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output.done:q})" {params.output_dir:q}
        plasmicheck pipeline \
            -hf {input.human_fasta:q} \
            -pf {input.plasmids_file:q} \
            -sf1 {input.bam:q} \
            -o {params.output_dir:q} \
            --threads {threads} \
            --log-file {output.log:q} \
            {params.extra_flags}
        touch {output.done:q}
        """


rule summary_reports:
    """Aggregate per-sample plasmid scan results into cohort summary tables and reports."""
    input:
        done=expand(f"{LOG_DIR}/{{sample}}.done", sample=SAMPLES),
    output:
        done=f"{LOG_DIR}/summary_reports.done",
        verdict_tsv=f"{SUMMARY_DIR}/data/verdict_table.tsv",
        ratio_tsv=f"{SUMMARY_DIR}/data/ratio_table.tsv",
        report_html=f"{SUMMARY_DIR}/summary_report_interactive.html",
    params:
        input_dir=OUTPUT_DIR,
        summary_dir=SUMMARY_DIR,
        extra_flags=SUMMARY_FLAGS,
    threads: get_rule_resource(config, "summary_reports", "threads")
    resources:
        mem_mb=get_rule_resource(config, "summary_reports", "mem_mb"),
        runtime=get_rule_resource(config, "summary_reports", "runtime"),
        tmpdir=os.environ.get("TMPDIR", "/tmp"),
    conda:
        CONDA_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output.done:q})" {params.summary_dir:q}
        plasmicheck summary_reports \
            -i {params.input_dir:q} \
            -o {params.summary_dir:q} \
            {params.extra_flags}
        touch {output.done:q}
        """
