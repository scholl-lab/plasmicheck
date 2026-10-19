from __future__ import annotations

import contextlib
import hashlib
import logging
import os
import secrets
import shutil
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from collections.abc import Generator
from dataclasses import dataclass, field
from pathlib import Path
from shlex import quote
from typing import Any

import pysam

from ..cache_manager import CacheManager, _resolve_md5
from ..config import get_config
from ..output_paths import OutputPaths
from ..thread_config import allocate_threads, detect_cpu_count
from . import utils as _utils
from .align_reads import align_reads, bam_to_fastq_cmd, has_bam_index
from .compare_alignments import (
    TARGET_BED_MAX_INTERVALS,
    TARGET_BED_SEED,
    compare_alignments,
    compute_bed_median_depth,
)
from .convert_plasmidfile_to_fasta import convert_plasmidfile_to_fasta
from .create_indexes import create_indexes
from .generate_report import DEFAULT_THRESHOLD
from .generate_report import main as generate_report
from .spliced_alignment import (
    extract_human_reference,
    extract_plasmid_cdna_positions,
    spliced_alignment,
)
from .utils import (
    archive_output_folder,
    calculate_md5,
    quality_control,
    run_command,
    sanitize_filename,
    write_md5sum,
)

_orig_run_command = _utils.run_command


def write_igv_session(paths: OutputPaths) -> None:
    """Write IGV session XML with relative paths to plasmid reference and BAM."""
    paths.igv_session.parent.mkdir(parents=True, exist_ok=True)
    session = ET.Element(
        "Session", genome=paths.plasmid_fasta.relative_to(paths.base).as_posix(), version="8"
    )
    resources = ET.SubElement(session, "Resources")
    ET.SubElement(resources, "Resource", path=paths.plasmid_bam.relative_to(paths.base).as_posix())
    ET.ElementTree(session).write(paths.igv_session, encoding="UTF-8", xml_declaration=True)


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SequencingInput:
    """A single or paired-end sequencing input."""

    file1: str
    file2: str | None = None

    @property
    def is_paired(self) -> bool:
        return self.file2 is not None


@dataclass
class PipelineStep:
    """A single step in the pipeline execution plan."""

    name: str
    description: str
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)
    skip: bool = False
    skip_reason: str = ""


@dataclass
class PipelinePlan:
    """Full execution plan for a pipeline run."""

    human_fasta: str
    plasmid_files: list[str]
    sequencing_inputs: list[SequencingInput]
    output_folder: str
    overwrite: bool
    combinations: list[tuple[str, SequencingInput]] = field(default_factory=list)
    combination_steps: dict[str, list[PipelineStep]] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    built_indexes: set[str] = field(default_factory=set)

    @property
    def total_steps(self) -> int:
        return sum(len(steps) for steps in self.combination_steps.values())

    @property
    def skipped_steps(self) -> int:
        return sum(1 for steps in self.combination_steps.values() for s in steps if s.skip)


@dataclass
class CombinationResult:
    """Result of processing a single combination (plasmid x sequencing file)."""

    combo_label: str
    success: bool
    duration: float
    error: Exception | None = None


# ---------------------------------------------------------------------------
# File list helpers
# ---------------------------------------------------------------------------


def read_file_list(file_path: str) -> list[str]:
    """Read a file containing newline-separated file paths and return a list of files."""
    with open(file_path) as file:
        lines: list[str] = file.readlines()
    return [line.strip() for line in lines if line.strip()]


def get_file_list(file_or_list: str) -> list[str]:
    """Determine if input is a single file or a list of files."""
    if os.path.isfile(file_or_list) and file_or_list.endswith(".txt"):
        return read_file_list(file_or_list)
    else:
        return [file_or_list]


def resolve_sequencing_inputs(
    sequencing_files_r1: str,
    sequencing_files_r2: str | None = None,
) -> list[SequencingInput]:
    """Normalize CLI sequencing file arguments into a list of SequencingInput.

    Args:
        sequencing_files_r1: Forward (R1) FASTQ/BAM file or file list (.txt).
        sequencing_files_r2: Reverse (R2) FASTQ file or file list (.txt) for paired-end.
    """
    r1_list = get_file_list(sequencing_files_r1)
    r2_list = get_file_list(sequencing_files_r2) if sequencing_files_r2 else []
    if sequencing_files_r2 and len(r1_list) != len(r2_list):
        raise ValueError(
            f"-sf1 has {len(r1_list)} files but -sf2 has {len(r2_list)} files; "
            "counts must match for paired-end."
        )
    if r2_list:
        return [SequencingInput(r1, r2) for r1, r2 in zip(r1_list, r2_list, strict=True)]
    return [SequencingInput(r1) for r1 in r1_list]


# ---------------------------------------------------------------------------
# Plan builder
# ---------------------------------------------------------------------------


def build_plan(
    human_fasta: str,
    plasmid_files: str,
    output_folder: str,
    overwrite: bool,
    sequencing_files_r1: str,
    sequencing_files_r2: str | None = None,
    cdna_output: str | None = None,
    sample_name: str | None = None,
) -> PipelinePlan:
    """Build an execution plan without running anything."""
    plasmid_file_list = get_file_list(plasmid_files)
    seq_inputs = resolve_sequencing_inputs(
        sequencing_files_r1=sequencing_files_r1,
        sequencing_files_r2=sequencing_files_r2,
    )
    if sample_name is not None:
        if not sanitize_filename(sample_name.strip()):
            raise ValueError("--sample-name must not be empty.")
        if len(seq_inputs) > 1:
            raise ValueError("--sample-name can only be used with a single sequencing input.")

    plan = PipelinePlan(
        human_fasta=human_fasta,
        plasmid_files=plasmid_file_list,
        sequencing_inputs=seq_inputs,
        output_folder=output_folder,
        overwrite=overwrite,
    )

    if overwrite:
        plan.warnings.append("Overwrite enabled — existing files will be replaced.")

    combo_idx = 0
    for seq_input in seq_inputs:
        for plasmid_file in plasmid_file_list:
            combo_idx += 1
            sequencing_file = seq_input.file1
            bam_basename = (
                sanitize_filename(sample_name.strip())
                if sample_name
                else sanitize_filename(os.path.splitext(os.path.basename(sequencing_file))[0])
            )
            file_basename = sanitize_filename(os.path.splitext(os.path.basename(plasmid_file))[0])
            combo_label = f"{os.path.basename(plasmid_file)} x {os.path.basename(sequencing_file)}"
            plan.combinations.append((combo_label, seq_input))

            # Create OutputPaths instance for this combination
            paths = OutputPaths(
                output_root=Path(output_folder),
                sample_name=bam_basename,
                plasmid_name=file_basename,
            )

            # Compute output paths using OutputPaths dataclass
            output_sub = str(paths.base)
            plasmid_fasta = str(paths.plasmid_fasta)
            human_index = os.path.join(
                os.path.dirname(human_fasta),
                os.path.splitext(os.path.basename(human_fasta))[0] + ".mmi",
            )
            plasmid_index = str(paths.plasmid_index)
            spliced_bam = str(paths.spliced_bam)
            spliced_fasta = str(paths.spliced_fasta)
            unique_cdna = cdna_output or str(paths.cdna_positions)
            plasmid_bam = str(paths.plasmid_bam)
            spliced_human_bam = str(paths.spliced_human_bam)
            comparison_output = paths.comparison_prefix

            steps: list[PipelineStep] = []

            # Step 1: Convert
            convert_skip = not overwrite and os.path.exists(plasmid_fasta)
            steps.append(
                PipelineStep(
                    name="convert",
                    description=f"Convert {os.path.basename(plasmid_file)} to FASTA",
                    inputs=[plasmid_file],
                    outputs=[plasmid_fasta],
                    skip=convert_skip,
                    skip_reason="Output FASTA already exists" if convert_skip else "",
                )
            )

            # Step 2: Index human + plasmid
            human_idx_skip = not overwrite and os.path.exists(human_index)
            steps.append(
                PipelineStep(
                    name="index",
                    description="Create indexes (human + plasmid references)",
                    inputs=[human_fasta, plasmid_fasta],
                    outputs=[human_index, plasmid_index],
                    skip=human_idx_skip and os.path.exists(plasmid_index),
                    skip_reason="Indexes already exist" if human_idx_skip else "",
                )
            )

            # Step 3: Spliced alignment
            steps.append(
                PipelineStep(
                    name="spliced",
                    description="Spliced alignment + extract human reference",
                    inputs=[human_index, plasmid_fasta],
                    outputs=[spliced_bam, spliced_fasta],
                )
            )

            # Step 4: cDNA positions
            steps.append(
                PipelineStep(
                    name="cdna",
                    description="Extract cDNA positions",
                    inputs=[plasmid_fasta, spliced_bam],
                    outputs=[unique_cdna],
                )
            )

            # Step 5: Align reads
            align_inputs = [sequencing_file]
            if seq_input.file2:
                align_inputs.append(seq_input.file2)
            steps.append(
                PipelineStep(
                    name="align",
                    description="Align reads to plasmid + human reference",
                    inputs=align_inputs,
                    outputs=[plasmid_bam, spliced_human_bam],
                )
            )

            # Step 6: Compare
            steps.append(
                PipelineStep(
                    name="compare",
                    description="Compare alignments and assign reads",
                    inputs=[plasmid_bam, spliced_human_bam],
                    outputs=[
                        f"{comparison_output}.reads_assignment.tsv",
                        f"{comparison_output}.summary.tsv",
                    ],
                )
            )

            # Step 7: Report
            steps.append(
                PipelineStep(
                    name="report",
                    description=f"Generate report in {os.path.basename(output_sub)}/",
                    inputs=[
                        f"{comparison_output}.reads_assignment.tsv",
                        f"{comparison_output}.summary.tsv",
                    ],
                    outputs=[],
                )
            )

            plan.combination_steps[combo_label] = steps

    return plan


def print_plan(plan: PipelinePlan) -> None:
    """Print the execution plan to stdout."""
    print("PlasmiCheck Pipeline Dry-Run")
    print("=" * 40)
    print()

    # Environment
    print("Environment:")
    minimap2 = shutil.which("minimap2")
    samtools = shutil.which("samtools")
    print(f"  [{'OK' if minimap2 else 'MISSING'}] minimap2")
    print(f"  [{'OK' if samtools else 'MISSING'}] samtools")
    human_exists = os.path.exists(plan.human_fasta)
    human_size = os.path.getsize(plan.human_fasta) if human_exists else 0
    print(
        f"  [{'OK' if human_exists else 'MISSING'}] Human reference: "
        f"{plan.human_fasta} ({human_size:,} bytes)"
    )
    print(f"  Output directory: {plan.output_folder}")
    print()

    # Inputs
    print("Inputs:")
    print(
        f"  Plasmids ({len(plan.plasmid_files)}): "
        + ", ".join(os.path.basename(p) for p in plan.plasmid_files)
    )
    print(
        f"  Sequencing inputs ({len(plan.sequencing_inputs)}): "
        + ", ".join(os.path.basename(s.file1) for s in plan.sequencing_inputs)
    )
    print()

    # Warnings
    for w in plan.warnings:
        print(f"  WARNING: {w}")
    if plan.warnings:
        print()

    # Steps
    n_combos = len(plan.combinations)
    human_index_path = os.path.join(
        os.path.dirname(plan.human_fasta),
        os.path.splitext(os.path.basename(plan.human_fasta))[0] + ".mmi",
    )
    print(f"Execution Plan ({n_combos} combinations, {plan.total_steps} steps):")
    print(f"  Human index: {human_index_path} (built once, shared across all combinations)")
    print()
    for idx, (combo_label, _seq_input) in enumerate(plan.combinations, 1):
        print(f"  [{idx}/{n_combos}] {combo_label}")
        for step_idx, step in enumerate(plan.combination_steps[combo_label], 1):
            skip_tag = f"  [SKIP: {step.skip_reason}]" if step.skip else ""
            print(f"    {step_idx}. {step.description}{skip_tag}")
        print()

    active = plan.total_steps - plan.skipped_steps
    print(f"Summary: {plan.total_steps} steps ({plan.skipped_steps} skipped, {active} to run)")


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Progress tracking
# ---------------------------------------------------------------------------


class _ProgressTracker:
    """Thin wrapper around Rich Progress (or no-op) that stores the task ID."""

    def __init__(self, progress: Any, task_id: int) -> None:
        self._progress = progress
        self._task_id = task_id

    def advance(self) -> None:
        self._progress.advance(self._task_id)

    @property
    def console(self) -> Any:
        return self._progress.console


class _NoOpConsole:
    @staticmethod
    def print(*_args: Any, **_kwargs: Any) -> None:
        pass


class _NoOpProgress:
    """Minimal stand-in when progress display is disabled."""

    console = _NoOpConsole()

    def advance(self, task_id: int = 0) -> None:
        pass


@contextlib.contextmanager
def _progress_context(enabled: bool, total: int) -> Generator[_ProgressTracker, None, None]:
    """Yield a progress tracker (Rich) or a no-op."""
    if not enabled:
        yield _ProgressTracker(_NoOpProgress(), 0)
        return

    from rich.progress import (
        BarColumn,
        Progress,
        SpinnerColumn,
        TextColumn,
        TimeElapsedColumn,
    )

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
    ) as progress:
        task_id = progress.add_task("Pipeline", total=total)
        yield _ProgressTracker(progress, task_id)


# ---------------------------------------------------------------------------
# Pipeline execution
# ---------------------------------------------------------------------------


def _extract_spliced_fasta_regions(prepared_plasmids: dict[str, dict[str, Any]]) -> list[str]:
    """Extract deduplicated 'chrom:start-end' region strings from prepared spliced FASTA headers."""
    regions: list[str] = []
    for info in prepared_plasmids.values():
        spliced_fasta = info.get("spliced_fasta")
        if not spliced_fasta or not os.path.exists(spliced_fasta):
            continue
        with open(spliced_fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    if ":" in header and "-" in header and header not in regions:
                        regions.append(header)
    return regions


def _prepare_plasmid_in_dir(
    plasmid_file: str,
    p_dir: Path,
    file_basename: str,
    human_fasta: str,
    human_index_path: str,
    *,
    cache: CacheManager | None,
    get_md5: Any,
    shift_bases: int = 500,
    generate_shifted: bool = False,
    overwrite: bool = False,
    padding: int = 1000,
    mm2_threads: int = 2,
    sam_threads: int = 2,
    gtf_file: str | None = None,
    gtf_md5: str = "",
) -> str:
    """Convert, index, and spliced-align a single plasmid into p_dir with CacheManager backing."""
    p_dir.mkdir(parents=True, exist_ok=True)
    from .format_detection import detect_format

    plasmid_file_type = detect_format(plasmid_file)
    logging.info(f"Processing plasmid file: {plasmid_file}")

    # Step 1: Convert the plasmid file to a FASTA file
    staged_fasta = p_dir / f"{file_basename}.fasta"
    if not os.path.exists(staged_fasta) or overwrite:
        logging.info(f"Converting {plasmid_file} to FASTA format.")
        convert_plasmidfile_to_fasta(
            plasmid_file,
            str(staged_fasta),
            plasmid_file_type,
            shift_bases,
            generate_shifted,
            overwrite,
            preserve_genbank=True,
        )

    # Step 2: Generate indices
    staged_index = p_dir / f"{file_basename}.mmi"
    staged_fai = p_dir / f"{file_basename}.fasta.fai"
    cached_index = None
    if cache and not overwrite:
        try:
            cached_index = cache.get_index(str(staged_fasta))
        except (OSError, PermissionError) as e:
            logging.warning(f"Cache read error: {e}")

    if cached_index:
        shutil.copy2(str(cached_index["mmi"]), str(staged_index))
        shutil.copy2(str(cached_index["fai"]), str(staged_fai))
    else:
        create_indexes(str(staged_fasta), overwrite, threads=mm2_threads)
        if cache:
            try:
                cache.put_index(str(staged_fasta), str(staged_index), str(staged_fai))
            except (OSError, PermissionError) as e:
                logging.warning(f"Cache write error: {e}")

    # Step 3: Spliced alignment
    staged_spliced_bam = p_dir / "spliced_alignment.bam"
    staged_spliced_fasta = p_dir / "spliced_reference.fasta"
    staged_cdna = p_dir / "cDNA_positions.txt"

    cached_spliced = None
    plasmid_md5 = ""
    human_md5 = ""
    if cache:
        try:
            plasmid_md5 = get_md5(str(staged_fasta))
            human_md5 = get_md5(human_fasta)
            if not overwrite:
                cached_spliced = cache.get_spliced(
                    plasmid_md5, human_md5, padding=padding, gtf_md5=gtf_md5
                )
        except (OSError, PermissionError) as e:
            logging.warning(f"Cache read error: {e}")

    if cached_spliced:
        shutil.copy2(str(cached_spliced["bam"]), str(staged_spliced_bam))
        shutil.copy2(str(cached_spliced["fasta"]), str(staged_spliced_fasta))
        shutil.copy2(str(cached_spliced["cdna_positions"]), str(staged_cdna))
        if os.path.exists(str(staged_spliced_bam)):
            pysam.index(str(staged_spliced_bam))
    else:
        logging.info(f"Performing spliced alignment for {plasmid_file}...")
        spanned_regions = spliced_alignment(
            human_index_path,
            str(staged_fasta),
            str(staged_spliced_bam),
            padding,
            threads_minimap2=mm2_threads,
            threads_samtools=sam_threads,
        )
        logging.info("Extracting human reference regions...")
        extract_human_reference(human_fasta, spanned_regions, str(staged_spliced_fasta))

        logging.info("Extracting cDNA positions from spliced alignment...")
        extract_plasmid_cdna_positions(
            str(staged_fasta),
            str(staged_spliced_bam),
            str(staged_cdna),
            padding=padding,
            gtf_file=gtf_file,
        )

        if cache:
            try:
                cache.put_spliced(
                    plasmid_md5,
                    human_md5,
                    str(staged_spliced_bam),
                    str(staged_spliced_fasta),
                    str(staged_cdna),
                    padding=padding,
                    gtf_md5=gtf_md5,
                )
            except (OSError, PermissionError) as e:
                logging.warning(f"Cache write error: {e}")

    # Step 4: Index spliced FASTA (with cache support)
    staged_spliced_index = p_dir / "spliced_reference.mmi"
    staged_spliced_fai = p_dir / "spliced_reference.fasta.fai"
    cached_spliced_idx = None
    if cache and not overwrite:
        try:
            cached_spliced_idx = cache.get_index(str(staged_spliced_fasta))
        except (OSError, PermissionError) as e:
            logging.warning(f"Cache read error: {e}")

    if cached_spliced_idx:
        shutil.copy2(str(cached_spliced_idx["mmi"]), str(staged_spliced_index))
        shutil.copy2(str(cached_spliced_idx["fai"]), str(staged_spliced_fai))
    else:
        create_indexes(str(staged_spliced_fasta), overwrite=True, threads=mm2_threads)
        if cache:
            try:
                cache.put_index(
                    str(staged_spliced_fasta),
                    str(staged_spliced_index),
                    str(staged_spliced_fai),
                )
            except (OSError, PermissionError) as e:
                logging.warning(f"Cache write error: {e}")

    return str(staged_spliced_fasta)


def prepare_references(
    human_fasta: str,
    plasmid_files: str | list[str],
    *,
    cache_dir: str | None = None,
    shift_bases: int = 500,
    generate_shifted: bool = False,
    overwrite: bool = False,
    padding: int = 1000,
    threads: int | None = None,
    gtf_file: str | None = None,
) -> None:
    """Build human reference index and populate CacheManager for all plasmid constructs."""
    plasmid_file_list = (
        get_file_list(plasmid_files)
        if isinstance(plasmid_files, str)
        else [str(p).strip() for p in plasmid_files if str(p).strip()]
    )
    if not plasmid_file_list:
        raise ValueError("No plasmid files provided for reference preparation.")

    total_threads = threads if threads is not None else detect_cpu_count()[0]
    mm2_threads, sam_threads = allocate_threads(total_threads)

    quality_control([human_fasta], expected_formats=["fasta"])
    quality_control(plasmid_file_list, expected_formats=["genbank", "xdna", "snapgene", "fasta"])

    cache_path = Path(cache_dir) if cache_dir else Path.home() / ".plasmicheck" / "cache"
    cache = CacheManager(cache_path)

    human_index_path = os.path.join(
        os.path.dirname(human_fasta),
        os.path.splitext(os.path.basename(human_fasta))[0] + ".mmi",
    )
    create_indexes(human_fasta, overwrite, threads=mm2_threads)

    md5_cache: dict[str, str] = {}

    def _get_md5(path: str) -> str:
        if path not in md5_cache:
            if os.path.exists(path):
                md5_cache[path] = calculate_md5(path)
            else:
                md5_cache[path] = hashlib.md5(path.encode()).hexdigest()
        return md5_cache[path]

    gtf_md5 = _resolve_md5(gtf_file) if gtf_file else ""

    with tempfile.TemporaryDirectory(prefix="plasmicheck_prep_") as tmp_dir:
        prep_dir = Path(tmp_dir)
        prepared: set[str] = set()
        for plasmid_file in plasmid_file_list:
            if plasmid_file in prepared:
                continue
            file_basename = sanitize_filename(os.path.splitext(os.path.basename(plasmid_file))[0])
            _prepare_plasmid_in_dir(
                plasmid_file,
                prep_dir / file_basename,
                file_basename,
                human_fasta,
                human_index_path,
                cache=cache,
                get_md5=_get_md5,
                shift_bases=shift_bases,
                generate_shifted=generate_shifted,
                overwrite=overwrite,
                padding=padding,
                mm2_threads=mm2_threads,
                sam_threads=sam_threads,
                gtf_file=gtf_file,
                gtf_md5=gtf_md5,
            )
            prepared.add(plasmid_file)


def run_pipeline(
    human_fasta: str,
    plasmid_files: str,
    output_folder: str = "output",
    keep_intermediate: bool = True,
    shift_bases: int = 500,
    generate_shifted: bool = False,
    overwrite: bool = False,
    padding: int = 1000,
    threshold: float = DEFAULT_THRESHOLD,
    md5_level: str = "all",
    cdna_output: str | None = None,
    archive_output: bool = False,
    sequencing_files_r1: str | None = None,
    sequencing_files_r2: str | None = None,
    dry_run: bool = False,
    progress: bool = False,
    static_report: bool = False,
    plotly_mode: str = "directory",
    plot_backend: str = "matplotlib",
    threads: int | None = None,
    cache_dir: str | None = None,
    no_cache: bool = False,
    anonymize: bool = False,
    regions: str | list[str] | None = None,
    include_unmapped: bool = False,
    gtf_file: str | None = None,
    vcf_file: str | None = None,
    target_bed: str | None = None,
    dp_sample_median: float | None = None,
    dp_locus: float | None = None,
    r_cohort: float = 0.60,
    regions_bed: str | None = None,
    target_bed_max_intervals: int = TARGET_BED_MAX_INTERVALS,
    target_bed_seed: int = TARGET_BED_SEED,
    sample_name: str | None = None,
) -> None:
    if sequencing_files_r1 is None:
        raise ValueError("-sf1 (sequencing_files_r1) is required.")
    if regions is None and regions_bed is not None:
        if not os.path.isfile(regions_bed):
            raise FileNotFoundError(f"Regions BED file not found: {regions_bed}")
        regions = regions_bed

    anonymize_key = secrets.token_bytes(16) if anonymize else None
    if anonymize and keep_intermediate:
        logging.warning(
            "--anonymize deletes intermediate BAMs: they hold original read names and sequences."
        )
        keep_intermediate = False

    md5_cache: dict[str, str] = {}

    def _get_md5(path: str) -> str:
        if path not in md5_cache:
            if os.path.exists(path):
                md5_cache[path] = calculate_md5(path)
            else:
                md5_cache[path] = hashlib.md5(path.encode()).hexdigest()
        return md5_cache[path]

    # Build the plan regardless of dry-run
    plan = build_plan(
        human_fasta=human_fasta,
        plasmid_files=plasmid_files,
        output_folder=output_folder,
        overwrite=overwrite,
        sequencing_files_r1=sequencing_files_r1,
        sequencing_files_r2=sequencing_files_r2,
        cdna_output=cdna_output,
        sample_name=sample_name,
    )

    # Initialize cache
    cache: CacheManager | None = None
    if not no_cache:
        try:
            cache_path = Path(cache_dir) if cache_dir else Path.home() / ".plasmicheck" / "cache"
            cache = CacheManager(cache_path)
            logging.info(f"Cache initialized at {cache_path}")
        except (OSError, PermissionError) as e:
            logging.warning(f"Cache directory inaccessible: {e}. Disabling cache for this run.")
            cache = None

    if dry_run:
        print_plan(plan)
        return

    logging.info("Starting the pipeline...")
    logging.info("Input parameters:")
    logging.info(f"Human reference FASTA: {human_fasta}")
    logging.info(f"Plasmid files: {plasmid_files}")
    logging.info(f"Output folder: {output_folder}")
    logging.info(f"Sequencing inputs: {len(plan.sequencing_inputs)} entries")

    # Thread detection and allocation
    if threads is not None:
        total_threads = threads
        thread_source = f"CLI --threads={threads}"
    else:
        total_threads, thread_source = detect_cpu_count()

    mm2_threads, sam_threads = allocate_threads(total_threads)
    sort_memory = get_config().get("alignment", {}).get("samtools_sort_memory", "2G")

    logging.info(f"Using {total_threads} threads ({thread_source})")
    logging.info(f"  minimap2: {mm2_threads} threads, samtools: {sam_threads} threads")
    logging.info(f"  samtools sort memory: {sort_memory} per thread")

    # Perform quality control checks
    logging.info("Performing quality control checks on input files...")
    quality_control([human_fasta], expected_formats=["fasta"])
    quality_control(plan.plasmid_files, expected_formats=["genbank", "xdna", "snapgene", "fasta"])

    fastq_extensions = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    for seq_input in plan.sequencing_inputs:
        if seq_input.file1.endswith(".bam"):
            quality_control([seq_input.file1], expected_formats=["bam"])
        elif seq_input.file2:
            quality_control([seq_input.file1, seq_input.file2], expected_formats=["fastq"])
        elif seq_input.file1.endswith(fastq_extensions):
            quality_control([seq_input.file1], expected_formats=["fastq"])
        else:
            raise ValueError(
                f"Unsupported sequencing file type: {seq_input.file1}. "
                "Must be .bam, .fastq, .fastq.gz, or paired FASTQ files."
            )

    # Upfront human reference indexing (shared across all combinations)
    human_index_path = os.path.join(
        os.path.dirname(human_fasta),
        os.path.splitext(os.path.basename(human_fasta))[0] + ".mmi",
    )
    logging.info("Building human reference index (shared across all combinations)...")
    create_indexes(human_fasta, overwrite, threads=mm2_threads)
    plan.built_indexes.add(human_index_path)
    logging.info(f"Built human reference index: {human_index_path}")

    # Track results for batch resilience
    results: list[CombinationResult] = []

    with (
        contextlib.ExitStack() as cleanup_stack,
        _progress_context(progress, plan.total_steps) as prog,
    ):
        prep_dir = Path(
            cleanup_stack.enter_context(tempfile.TemporaryDirectory(prefix="plasmicheck_prep_"))
        )
        gtf_md5 = _resolve_md5(gtf_file) if gtf_file else ""
        prepared_plasmids: dict[str, dict[str, Any]] = {}
        prep_errors: dict[str, Exception] = {}

        # Hoist plasmid preparation before the sample loop so all plasmid spliced
        # references are ready before extracting region-restricted BAM reads.
        for plasmid_file in plan.plasmid_files:
            if plasmid_file in prepared_plasmids or plasmid_file in prep_errors:
                continue
            file_basename = sanitize_filename(os.path.splitext(os.path.basename(plasmid_file))[0])
            p_dir = prep_dir / file_basename
            try:
                staged_spliced_fasta = _prepare_plasmid_in_dir(
                    plasmid_file,
                    p_dir,
                    file_basename,
                    human_fasta,
                    human_index_path,
                    cache=cache,
                    get_md5=_get_md5,
                    shift_bases=shift_bases,
                    generate_shifted=generate_shifted,
                    overwrite=overwrite,
                    padding=padding,
                    mm2_threads=mm2_threads,
                    sam_threads=sam_threads,
                    gtf_file=gtf_file,
                    gtf_md5=gtf_md5,
                )
                prepared_plasmids[plasmid_file] = {
                    "p_dir": p_dir,
                    "file_basename": file_basename,
                    "spliced_fasta": staged_spliced_fasta,
                }
            except Exception as e:
                prep_errors[plasmid_file] = e

        auto_spliced_regions = _extract_spliced_fasta_regions(prepared_plasmids)

        for seq_input in plan.sequencing_inputs:
            sequencing_file = seq_input.file1
            fastq2 = seq_input.file2
            bam_basename = (
                sanitize_filename(sample_name.strip())
                if sample_name
                else sanitize_filename(os.path.splitext(os.path.basename(sequencing_file))[0])
            )
            sample_temp_fastq: str | None = None
            eff_regions: str | list[str] | None = None
            try:
                if not seq_input.is_paired and sequencing_file.endswith(".bam"):
                    if regions is not None:
                        eff_regions = regions
                    elif has_bam_index(sequencing_file) and auto_spliced_regions:
                        eff_regions = auto_spliced_regions

                    if os.path.exists(sequencing_file) and (
                        bool(eff_regions)
                        or (len(plan.plasmid_files) > 1 and os.path.getsize(sequencing_file) > 64)
                    ):
                        sample_temp_fastq = str(
                            prep_dir
                            / f".tmp_{sanitize_filename(os.path.basename(sequencing_file))}_{os.getpid()}.fastq"
                        )
                        logging.info(
                            f"Pre-collating paired reads from BAM once for {sequencing_file}..."
                        )
                        collate_cmd = (
                            f"{bam_to_fastq_cmd(sequencing_file, sam_threads, eff_regions, include_unmapped=include_unmapped)} "
                            f"> {quote(sample_temp_fastq)}"
                        )
                        _cmd_runner = (
                            run_command
                            if run_command is not _orig_run_command
                            else _utils.run_command
                        )
                        _cmd_runner(collate_cmd)

                sample_median_dp = dp_sample_median
                if (
                    sample_median_dp is None
                    and target_bed is not None
                    and sequencing_file.endswith(".bam")
                ):
                    sample_median_dp = compute_bed_median_depth(
                        sequencing_file, target_bed, target_bed_max_intervals, target_bed_seed
                    )

                align_input = sample_temp_fastq if sample_temp_fastq else sequencing_file
                align_regions = eff_regions if not sample_temp_fastq else None

                for plasmid_file in plan.plasmid_files:
                    combo_label = (
                        f"{os.path.basename(plasmid_file)} x {os.path.basename(sequencing_file)}"
                    )
                    start_time = time.time()

                    try:
                        logging.info(f"Processing sequencing file: {sequencing_file}")

                        file_basename = sanitize_filename(
                            os.path.splitext(os.path.basename(plasmid_file))[0]
                        )

                        # Create OutputPaths instance for this combination
                        paths = OutputPaths(
                            output_root=Path(output_folder),
                            sample_name=bam_basename,
                            plasmid_name=file_basename,
                        )
                        paths.ensure_directories()
                        output_subfolder = str(paths.base)

                        if plasmid_file in prep_errors:
                            raise prep_errors[plasmid_file]

                        p_dir = prep_dir / file_basename

                        # Copy staged artifacts into OutputPaths
                        plasmid_fasta = str(paths.plasmid_fasta)
                        plasmid_index = str(paths.plasmid_index)
                        spliced_bam = str(paths.spliced_bam)
                        spliced_fasta = str(paths.spliced_fasta)
                        spliced_index = str(paths.spliced_index)
                        unique_cdna_output = cdna_output or str(paths.cdna_positions)

                        def _copy_if_exists(src: Path | str, dst: Path | str) -> None:
                            if os.path.exists(src):
                                shutil.copy2(str(src), str(dst))

                        _copy_if_exists(p_dir / f"{file_basename}.fasta", plasmid_fasta)
                        _copy_if_exists(
                            p_dir / f"{file_basename}.fasta.fai", f"{plasmid_fasta}.fai"
                        )
                        _copy_if_exists(p_dir / f"{file_basename}.mmi", plasmid_index)
                        _copy_if_exists(p_dir / f"{file_basename}.gb", paths.plasmid_genbank)

                        _copy_if_exists(p_dir / "spliced_alignment.bam", spliced_bam)
                        _copy_if_exists(p_dir / "spliced_alignment.bam.bai", f"{spliced_bam}.bai")
                        _copy_if_exists(p_dir / "spliced_reference.fasta", spliced_fasta)
                        _copy_if_exists(
                            p_dir / "spliced_reference.fasta.fai", f"{spliced_fasta}.fai"
                        )
                        _copy_if_exists(p_dir / "spliced_reference.mmi", spliced_index)
                        _copy_if_exists(p_dir / "cDNA_positions.txt", unique_cdna_output)

                        if md5_level in ["all", "intermediate"]:
                            if os.path.exists(plasmid_fasta):
                                write_md5sum(plasmid_fasta, "intermediate", output_subfolder)
                            if os.path.exists(plasmid_index):
                                write_md5sum(plasmid_index, "intermediate", output_subfolder)
                            if os.path.exists(human_index_path) and md5_level in ["all"]:
                                write_md5sum(human_index_path, "intermediate", output_subfolder)
                            if os.path.exists(spliced_bam):
                                write_md5sum(spliced_bam, "intermediate", output_subfolder)
                            if os.path.exists(spliced_fasta):
                                write_md5sum(spliced_fasta, "intermediate", output_subfolder)
                            if os.path.exists(spliced_index):
                                write_md5sum(spliced_index, "intermediate", output_subfolder)
                            if os.path.exists(unique_cdna_output):
                                write_md5sum(unique_cdna_output, "intermediate", output_subfolder)
                        if md5_level in ["all", "input"] and os.path.exists(plasmid_file):
                            write_md5sum(plasmid_file, "input", output_subfolder)

                        prog.advance()  # Step 1: convert
                        prog.advance()  # Step 2: index
                        prog.advance()  # Step 3: spliced
                        prog.advance()  # Step 4: cDNA

                        # Step 5: Align reads
                        plasmid_bam = str(paths.plasmid_bam)
                        spliced_human_bam = str(paths.spliced_human_bam)
                        logging.info("Aligning reads to plasmid and spliced human reference...")
                        if (
                            total_threads >= 4
                            and os.path.exists(plasmid_index)
                            and os.path.exists(spliced_index)
                        ):
                            from concurrent.futures import ThreadPoolExecutor

                            half_mm2 = max(1, mm2_threads // 2)
                            half_sam = max(1, sam_threads // 2)
                            with ThreadPoolExecutor(max_workers=2) as pool:
                                f1 = pool.submit(
                                    align_reads,
                                    plasmid_index,
                                    align_input,
                                    plasmid_bam,
                                    "plasmid",
                                    fastq2,
                                    minimap2_threads=half_mm2,
                                    samtools_threads=half_sam,
                                    samtools_sort_memory=sort_memory,
                                    regions=align_regions,
                                    include_unmapped=include_unmapped,
                                )
                                f2 = pool.submit(
                                    align_reads,
                                    spliced_index,
                                    align_input,
                                    spliced_human_bam,
                                    "human",
                                    fastq2,
                                    minimap2_threads=half_mm2,
                                    samtools_threads=half_sam,
                                    samtools_sort_memory=sort_memory,
                                    regions=align_regions,
                                    include_unmapped=include_unmapped,
                                )
                                f1.result()
                                f2.result()
                        else:
                            align_reads(
                                plasmid_index,
                                align_input,
                                plasmid_bam,
                                "plasmid",
                                fastq2,
                                minimap2_threads=mm2_threads,
                                samtools_threads=sam_threads,
                                samtools_sort_memory=sort_memory,
                                regions=align_regions,
                                include_unmapped=include_unmapped,
                            )
                            align_reads(
                                spliced_index,
                                align_input,
                                spliced_human_bam,
                                "human",
                                fastq2,
                                minimap2_threads=mm2_threads,
                                samtools_threads=sam_threads,
                                samtools_sort_memory=sort_memory,
                                regions=align_regions,
                                include_unmapped=include_unmapped,
                            )
                        if md5_level in ["all", "intermediate"]:
                            if os.path.exists(plasmid_bam):
                                write_md5sum(plasmid_bam, "intermediate", output_subfolder)
                            if os.path.exists(spliced_human_bam):
                                write_md5sum(spliced_human_bam, "intermediate", output_subfolder)
                        prog.advance()

                        # Step 6: Compare alignments
                        comparison_output = paths.comparison_prefix
                        logging.info("Comparing the two alignments...")

                        from .format_detection import detect_format

                        p_fmt = detect_format(plasmid_file)
                        genbank_file_path: str | None = None
                        if p_fmt == "genbank":
                            genbank_file_path = plasmid_file  # Original .gb file
                        elif p_fmt in ["xdna", "snapgene"]:
                            candidate_gb = str(paths.plasmid_genbank)
                            if os.path.exists(candidate_gb):
                                genbank_file_path = candidate_gb

                        compare_alignments(
                            plasmid_bam,
                            spliced_human_bam,
                            comparison_output,
                            threshold=threshold,
                            genbank_file=genbank_file_path,
                            cdna_positions_file=unique_cdna_output,
                            threads=sam_threads,
                            anonymize=anonymize,
                            anonymize_key=anonymize_key,
                            vcf_file=vcf_file,
                            target_bed=target_bed,
                            target_bed_max_intervals=target_bed_max_intervals,
                            target_bed_seed=target_bed_seed,
                            dp_sample_median=sample_median_dp,
                            dp_locus=dp_locus,
                            r_cohort=r_cohort,
                            sample_bam=sequencing_file
                            if sequencing_file.endswith(".bam")
                            else None,
                        )
                        prog.advance()

                        # IGV session XML generation (Issue #49) before report so report can link it
                        if keep_intermediate and paths.plasmid_bam.exists():
                            write_igv_session(paths)

                        # Step 7: Generate report
                        logging.info("Generating report...")
                        reads_assignment_file = f"{comparison_output}.reads_assignment.tsv"
                        summary_file = f"{comparison_output}.summary.tsv"

                        def _shown(p: str) -> str:
                            return (
                                os.path.basename(p.rstrip(os.sep))
                                if anonymize and os.sep in p
                                else p
                            )

                        command_line = " ".join(_shown(a) for a in sys.argv)
                        report_extra_kwargs: dict[str, Any] = (
                            {"sample_name": bam_basename} if sample_name else {}
                        )
                        generate_report(
                            reads_assignment_file,
                            summary_file,
                            output_subfolder,
                            threshold=threshold,
                            human_fasta=_shown(human_fasta),
                            plasmid_gb=_shown(plasmid_file),
                            sequencing_file=_shown(sequencing_file),
                            command_line=command_line,
                            static_report=static_report,
                            plotly_mode=plotly_mode,
                            plot_backend=plot_backend,
                            reports_dir=str(paths.reports_dir),
                            anonymize=anonymize,
                            **report_extra_kwargs,
                        )

                        # Write MD5 checksums for the output files
                        if md5_level in ["all", "intermediate", "output"]:
                            if os.path.exists(reads_assignment_file):
                                write_md5sum(reads_assignment_file, "output", output_subfolder)
                            if os.path.exists(summary_file):
                                write_md5sum(summary_file, "output", output_subfolder)
                        prog.advance()

                        # Step 8: Optionally delete intermediate files
                        if not keep_intermediate:
                            logging.info("Deleting intermediate files...")
                            intermediate_files = [
                                plasmid_bam,
                                spliced_human_bam,
                                spliced_bam,
                                f"{spliced_bam}.bai",
                                f"{spliced_fasta}.fai",
                                f"{plasmid_bam}.bai",
                                f"{spliced_human_bam}.bai",
                            ]
                            for f in intermediate_files:
                                if os.path.exists(f):
                                    with contextlib.suppress(OSError):
                                        os.remove(f)

                        # Step 9: Optionally archive the output folder
                        if archive_output:
                            archive_output_folder(output_subfolder)

                        # Record successful completion
                        duration = time.time() - start_time
                        results.append(
                            CombinationResult(
                                combo_label=combo_label, success=True, duration=duration, error=None
                            )
                        )
                        logging.info(f"{combo_label}: {duration:.1f}s")

                    except Exception as e:
                        # Record failure and continue to next combination
                        duration = time.time() - start_time
                        results.append(
                            CombinationResult(
                                combo_label=combo_label, success=False, duration=duration, error=e
                            )
                        )
                        logging.error(f"{combo_label} FAILED after {duration:.1f}s: {e}")
            finally:
                if sample_temp_fastq:
                    with contextlib.suppress(OSError):
                        if os.path.exists(sample_temp_fastq):
                            os.remove(sample_temp_fastq)

    # Log batch summary
    successful = [r for r in results if r.success]
    failed = [r for r in results if not r.success]

    logging.info(f"Pipeline complete: {len(successful)} succeeded, {len(failed)} failed")

    if failed:
        logging.warning("Failed combinations:")
        for result in failed:
            logging.warning(f"  - {result.combo_label}: {result.error}")

    # If ALL combinations failed, raise an error
    if len(failed) == len(results) and results:
        error_summary = "; ".join(f"{r.combo_label}: {r.error}" for r in failed)
        raise RuntimeError(f"All {len(failed)} combinations failed. Errors: {error_summary}")

    if successful:
        logging.info("Pipeline completed successfully for successful combinations.")
