# Snakemake cohort workflow for PlasmiCheck plasmid contamination detection (#129).
# Extracts each sample BAM once and screens all plasmid constructs in a single pipeline run.

import csv
import os
import re
import shlex
from pathlib import Path
from typing import Any

try:
    from plasmicheck.output_paths import OutputPaths
    from plasmicheck.scripts.utils import sanitize_filename
except ImportError:  # pragma: no cover — cluster head-node without plasmicheck installed

    def sanitize_filename(name: str) -> str:
        return re.sub(r"[^a-zA-Z0-9_-]", "_", name)

    class OutputPaths:  # type: ignore[no-redef]
        def __init__(self, output_root: Path, sample_name: str, plasmid_name: str) -> None:
            self.comparison_prefix = str(
                Path(output_root) / sample_name / plasmid_name / "data" / "comparison_result"
            )

DEFAULT_HUMAN_FASTA = "reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
DEFAULT_BAMS_FILE = "data/bams.txt"
DEFAULT_PLASMIDS_FILE = "data/plasmids.txt"

DEFAULT_RESOURCES: dict[str, dict[str, int]] = {
    "prepare_constructs": {"threads": 4, "mem_mb": 16000, "runtime": 30},
    "scan_sample": {"threads": 4, "mem_mb": 4000, "runtime": 20},
    "summary_reports": {"threads": 2, "mem_mb": 4000, "runtime": 30},
}


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

    shift_bases_val = cfg.get("shift_bases")
    if shift_bases_val is not None and (
        isinstance(shift_bases_val, bool)
        or not isinstance(shift_bases_val, int)
        or shift_bases_val < 0
    ):
        raise ValueError("'shift_bases' must be a non-negative integer.")

    generate_shifted_val = cfg.get("generate_shifted")
    if generate_shifted_val is not None and not isinstance(generate_shifted_val, bool):
        raise ValueError("'generate_shifted' must be a boolean.")

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


def _validate_sample_map(raw_samples: dict[str, str]) -> dict[str, str]:
    """Sanitize and validate uniqueness of sample IDs and BAM paths."""
    if not raw_samples:
        raise ValueError("No samples configured.")
    validated: dict[str, str] = {}
    for raw_id, raw_bam in raw_samples.items():
        sample_id = sanitize_filename(str(raw_id).strip())
        bam_path = str(raw_bam).strip()
        if not sample_id:
            raise ValueError(f"Invalid empty sample key: {raw_id!r}")
        if not bam_path:
            raise ValueError(f"Missing BAM path for sample '{sample_id}'")
        if sample_id in validated:
            raise ValueError(f"Duplicate sample key: {sample_id}")
        validated[sample_id] = bam_path
    return validated


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
    return _validate_sample_map(samples)


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
            raw_id = row[sample_idx].strip()
            if not raw_id:
                continue
            sample_id = sanitize_filename(raw_id)
            if bam_idx is not None:
                if len(row) <= bam_idx or not row[bam_idx].strip():
                    raise ValueError(f"Missing BAM path for sample '{sample_id}' in {tsv_path}")
                bam_path = row[bam_idx].strip()
            else:
                bam_path = os.path.join(str(bam_folder), f"{raw_id}{bam_ext}")
            if sample_id in samples and samples[sample_id] != bam_path:
                raise ValueError(f"Duplicate sample key: {sample_id}")
            samples[sample_id] = bam_path
    else:
        for row in reader:
            if len(row) < 2:
                raise ValueError(
                    f"Expected 'sample\\tbam' columns in {tsv_path}, got: {row!r}"
                )
            sample_id = sanitize_filename(row[0].strip())
            bam_path = row[1].strip()
            if sample_id in samples:
                raise ValueError(f"Duplicate sample key: {sample_id}")
            samples[sample_id] = bam_path

    if not samples:
        raise ValueError(f"No samples resolved from {tsv_path}")
    return _validate_sample_map(samples)


def resolve_samples(cfg: dict[str, Any]) -> dict[str, str]:
    """Resolve {sample_id: bam_path} from config dict, TSV sheet, or BAM list file."""
    samples_src = cfg.get("samples")
    if samples_src is None:
        samples_src = cfg.get("bams_file", DEFAULT_BAMS_FILE)

    if isinstance(samples_src, dict):
        if not samples_src:
            raise ValueError("No samples configured.")
        sanitized: dict[str, str] = {}
        for k, v in samples_src.items():
            sid = sanitize_filename(str(k).strip())
            if not sid:
                raise ValueError(f"Invalid empty sample key: {k!r}")
            if sid in sanitized:
                raise ValueError(f"Duplicate sample key: {sid}")
            sanitized[sid] = str(v).strip()
        return _validate_sample_map(sanitized)

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


def _construct_names(constructs: list[str]) -> list[str]:
    """Return unique sanitized construct names in order."""
    names: list[str] = []
    seen: set[str] = set()
    for p in constructs:
        name = sanitize_filename(os.path.splitext(os.path.basename(p))[0])
        if name and name not in seen:
            seen.add(name)
            names.append(name)
    return names


def expected_sample_outputs(output_dir: str, sample: str, constructs: list[str]) -> list[str]:
    """Return expected (.summary.tsv, .reads_assignment.tsv) paths for a sample across all constructs."""
    sample_id = sanitize_filename(sample)
    expected: list[str] = []
    for plasmid_name in _construct_names(constructs):
        paths = OutputPaths(
            output_root=Path(output_dir),
            sample_name=sample_id,
            plasmid_name=plasmid_name,
        )
        expected.append(f"{paths.comparison_prefix}.summary.tsv")
        expected.append(f"{paths.comparison_prefix}.reads_assignment.tsv")
    return expected


def verify_sample_outputs(output_dir: str, sample: str, constructs: list[str]) -> None:
    """Verify that every expected construct output for a sample exists and is non-empty."""
    expected = expected_sample_outputs(output_dir, sample, constructs)
    if not expected:
        raise RuntimeError(f"No constructs configured for sample '{sample}'.")
    missing = [p for p in expected if not os.path.isfile(p) or os.path.getsize(p) == 0]
    if missing:
        raise RuntimeError(
            f"Incomplete construct scan for sample '{sample}'; missing outputs: {', '.join(missing)}"
        )


def expected_cohort_summary_files(
    output_dir: str, samples: list[str], constructs: list[str]
) -> list[str]:
    """Return explicit .summary.tsv paths for the configured cohort (samples x constructs)."""
    summary_files: list[str] = []
    construct_names = _construct_names(constructs)
    for sample in samples:
        sample_id = sanitize_filename(sample)
        for plasmid_name in construct_names:
            paths = OutputPaths(
                output_root=Path(output_dir),
                sample_name=sample_id,
                plasmid_name=plasmid_name,
            )
            summary_files.append(f"{paths.comparison_prefix}.summary.tsv")
    return summary_files


def write_summary_manifest(
    manifest_path: str, output_dir: str, samples: list[str], constructs: list[str]
) -> list[str]:
    """Write explicit manifest of configured .summary.tsv paths and return the list."""
    files = expected_cohort_summary_files(output_dir, samples, constructs)
    os.makedirs(os.path.dirname(os.path.abspath(manifest_path)), exist_ok=True)
    with open(manifest_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(files) + ("\n" if files else ""))
    return files


def get_rule_resource(cfg: dict[str, Any], rule_name: str, key: str) -> int:
    """Return per-rule resource value ('threads', 'mem_mb', 'runtime') with config override."""
    defaults = DEFAULT_RESOURCES[rule_name]
    rule_cfg = cfg.get("resources", {}).get(rule_name, {})
    return int(rule_cfg.get(key, defaults[key]))


def resolve_gtf(cfg: dict[str, Any]) -> str | None:
    """Resolve optional GTF/MANE annotation path from config or aliases."""
    gtf_val = cfg.get("gtf_file") or cfg.get("gtf") or cfg.get("mane_gtf")
    return str(gtf_val) if gtf_val else None


def build_prepare_kwargs(cfg: dict[str, Any]) -> dict[str, Any]:
    """Build keyword arguments for prepare_references matching scan_sample pipeline flags."""
    kwargs: dict[str, Any] = {
        "cache_dir": str(cfg.get("cache_dir") or os.path.expanduser("~/.plasmicheck/cache")),
        "padding": int(cfg["padding"]) if cfg.get("padding") is not None else 1000,
        "shift_bases": int(cfg["shift_bases"]) if cfg.get("shift_bases") is not None else 500,
        "generate_shifted": bool(cfg.get("generate_shifted", False)),
        "overwrite": bool(cfg.get("overwrite", False)),
        "gtf_file": resolve_gtf(cfg),
    }
    extra = cfg.get("extra_args")
    if extra:
        tokens: list[str] = []
        for raw in shlex.split(str(extra)):
            tokens.extend(raw.split("=", 1) if raw.startswith("-") and "=" in raw else [raw])
        idx = 0
        while idx < len(tokens):
            tok = tokens[idx]
            if tok == "--gtf" and idx + 1 < len(tokens):
                kwargs["gtf_file"] = tokens[idx + 1]
                idx += 2
            elif tok in ("-d", "--padding") and idx + 1 < len(tokens):
                kwargs["padding"] = int(tokens[idx + 1])
                idx += 2
            elif tok in ("-sb", "--shift_bases") and idx + 1 < len(tokens):
                kwargs["shift_bases"] = int(tokens[idx + 1])
                idx += 2
            elif tok == "--cache-dir" and idx + 1 < len(tokens):
                kwargs["cache_dir"] = tokens[idx + 1]
                idx += 2
            elif tok in ("-g", "--generate_shifted"):
                kwargs["generate_shifted"] = True
                idx += 1
            elif tok in ("-w", "--overwrite"):
                kwargs["overwrite"] = True
                idx += 1
            else:
                idx += 1
    return kwargs


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
        for key in ("target_bed_max_intervals", "target_bed_seed"):
            if cfg.get(key) is not None:
                flags.append(f"--{key.replace('_', '-')} {int(cfg[key])}")
    gtf_val = resolve_gtf(cfg)
    if gtf_val:
        flags.append(f"--gtf {shlex.quote(gtf_val)}")
    vcf_val = cfg.get("vcf_file") or cfg.get("vcf")
    if vcf_val:
        flags.append(f"--vcf {shlex.quote(str(vcf_val))}")
    if cfg.get("threshold") is not None:
        flags.append(f"--threshold {shlex.quote(str(cfg['threshold']))}")
    if cfg.get("padding") is not None:
        flags.append(f"--padding {shlex.quote(str(cfg['padding']))}")
    if cfg.get("shift_bases") is not None:
        flags.append(f"--shift_bases {shlex.quote(str(cfg['shift_bases']))}")
    if cfg.get("generate_shifted"):
        flags.append("--generate_shifted")
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
    vcf_val = cfg.get("vcf_file") or cfg.get("vcf")
    if vcf_val:
        flags.append(f"--vcf {shlex.quote(str(vcf_val))}")
    return " ".join(flags)


# ----------------------------------------------------------------------------------- #
# Resolve configuration and inputs
# ----------------------------------------------------------------------------------- #

config = globals().get("config", {})
validate_config(config)

PREPARE_KWARGS = build_prepare_kwargs(config)
HUMAN_FASTA = str(config.get("human_fasta", DEFAULT_HUMAN_FASTA))
OUTPUT_DIR = str(config.get("output_dir", "output"))
SUMMARY_DIR = str(config.get("summary_dir", os.path.join(OUTPUT_DIR, "summary")))
LOG_DIR = str(config.get("log_dir", "logs"))
CACHE_DIR = str(PREPARE_KWARGS["cache_dir"])
CONDA_ENV = str(config.get("conda_env", "plasmicheck"))
GTF_FILE = PREPARE_KWARGS["gtf_file"]

SAMPLE_BAMS = resolve_samples(config)
SAMPLES = list(SAMPLE_BAMS.keys())

# Restrict {sample} to configured IDs so it cannot match prepare_constructs/summary_reports markers.
SAMPLE_WILDCARD_RE = "|".join(re.escape(s) for s in SAMPLES) or "(?!)"
PLASMIDS_FILE, CONSTRUCT_FILES = resolve_plasmids(config, OUTPUT_DIR)

PIPELINE_FLAGS = build_pipeline_flags(config)
SUMMARY_FLAGS = build_summary_flags(config)
SUMMARY_MANIFEST = os.path.join(OUTPUT_DIR, ".summary_manifest.txt")


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
    """Build human reference index and populate CacheManager for all plasmid constructs once."""
    input:
        human_fasta=HUMAN_FASTA,
        plasmids_file=PLASMIDS_FILE,
        constructs=CONSTRUCT_FILES,
    output:
        done=f"{LOG_DIR}/prepare_constructs.done",
    params:
        cache_dir=CACHE_DIR,
        padding=PREPARE_KWARGS["padding"],
        shift_bases=PREPARE_KWARGS["shift_bases"],
        generate_shifted="1" if PREPARE_KWARGS["generate_shifted"] else "0",
        overwrite="1" if PREPARE_KWARGS["overwrite"] else "0",
        gtf_file=GTF_FILE or "",
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
        python3 -c 'import sys; from plasmicheck.scripts.run_pipeline import prepare_references; prepare_references(sys.argv[1], sys.argv[2], cache_dir=sys.argv[3], padding=int(sys.argv[4]), overwrite=sys.argv[5] == "1", threads=int(sys.argv[6]), gtf_file=sys.argv[7] or None, shift_bases=int(sys.argv[8]), generate_shifted=sys.argv[9] == "1")' \
            {input.human_fasta:q} \
            {input.plasmids_file:q} \
            {params.cache_dir:q} \
            {params.padding} \
            {params.overwrite} \
            {threads} \
            {params.gtf_file:q} \
            {params.shift_bases} \
            {params.generate_shifted}
        touch {output.done:q}
        """


rule scan_sample:
    """Screen a single sample BAM against all plasmid constructs in one pipeline run."""
    wildcard_constraints:
        sample=SAMPLE_WILDCARD_RE,
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
        expected_outputs=lambda wildcards: expected_sample_outputs(
            OUTPUT_DIR, wildcards.sample, CONSTRUCT_FILES
        ),
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
        rm -f {params.expected_outputs:q}
        plasmicheck pipeline \
            -hf {input.human_fasta:q} \
            -pf {input.plasmids_file:q} \
            -sf1 {input.bam:q} \
            --sample-name {wildcards.sample:q} \
            -o {params.output_dir:q} \
            --threads {threads} \
            --log-file {output.log:q} \
            {params.extra_flags}
        for f in {params.expected_outputs:q}; do
            [ -s "$f" ] || {{ echo "Missing or empty construct output: $f" >&2; exit 1; }}
        done
        touch {output.done:q}
        """


rule summary_reports:
    """Aggregate configured per-sample plasmid scan results into cohort summary tables and reports."""
    input:
        done=expand(f"{LOG_DIR}/{{sample}}.done", sample=SAMPLES),
        plasmids_file=PLASMIDS_FILE,
    output:
        done=f"{LOG_DIR}/summary_reports.done",
        verdict_tsv=f"{SUMMARY_DIR}/data/verdict_table.tsv",
        ratio_tsv=f"{SUMMARY_DIR}/data/ratio_table.tsv",
        report_html=f"{SUMMARY_DIR}/summary_report_interactive.html",
    params:
        summary_dir=SUMMARY_DIR,
        manifest=SUMMARY_MANIFEST,
        summary_files=lambda wildcards: expected_cohort_summary_files(
            OUTPUT_DIR, SAMPLES, CONSTRUCT_FILES
        ),
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
        mkdir -p "$(dirname {output.done:q})" "$(dirname {params.manifest:q})" {params.summary_dir:q}
        printf '%s\n' {params.summary_files:q} > {params.manifest:q}
        plasmicheck summary_reports \
            -i {params.manifest:q} \
            -o {params.summary_dir:q} \
            {params.extra_flags}
        touch {output.done:q}
        """
