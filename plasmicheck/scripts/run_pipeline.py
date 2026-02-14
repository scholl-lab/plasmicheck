from __future__ import annotations

import contextlib
import logging
import os
import shutil
import sys
from collections.abc import Generator
from dataclasses import dataclass, field
from typing import Any

from .align_reads import align_reads
from .compare_alignments import compare_alignments
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
    add_logging_args,
    archive_output_folder,
    configure_logging_from_args,
    quality_control,
    sanitize_filename,
    write_md5sum,
)

# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SequencingInput:
    """A single or paired-end sequencing input."""

    file1: str
    file2: str | None = None


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

    @property
    def total_steps(self) -> int:
        return sum(len(steps) for steps in self.combination_steps.values())

    @property
    def skipped_steps(self) -> int:
        return sum(1 for steps in self.combination_steps.values() for s in steps if s.skip)


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
) -> PipelinePlan:
    """Build an execution plan without running anything."""
    plasmid_file_list = get_file_list(plasmid_files)
    seq_inputs = resolve_sequencing_inputs(
        sequencing_files_r1=sequencing_files_r1,
        sequencing_files_r2=sequencing_files_r2,
    )

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
    for plasmid_file in plasmid_file_list:
        for seq_input in seq_inputs:
            combo_idx += 1
            sequencing_file = seq_input.file1
            bam_basename = sanitize_filename(os.path.splitext(os.path.basename(sequencing_file))[0])
            file_basename = sanitize_filename(os.path.splitext(os.path.basename(plasmid_file))[0])
            output_sub = os.path.join(output_folder, bam_basename, file_basename)
            combo_label = f"{os.path.basename(plasmid_file)} x {os.path.basename(sequencing_file)}"
            plan.combinations.append((combo_label, seq_input))

            # Compute output paths
            plasmid_fasta = os.path.join(
                output_sub, os.path.splitext(os.path.basename(plasmid_file))[0] + ".fasta"
            )
            human_index = os.path.join(
                os.path.dirname(human_fasta),
                os.path.splitext(os.path.basename(human_fasta))[0] + ".mmi",
            )
            plasmid_index = os.path.join(
                output_sub, os.path.splitext(os.path.basename(plasmid_fasta))[0] + ".mmi"
            )
            spliced_bam = os.path.join(output_sub, "spliced_alignment.bam")
            spliced_fasta = os.path.join(output_sub, "spliced_reference.fasta")
            unique_cdna = cdna_output or os.path.join(output_sub, "cDNA_positions.txt")
            plasmid_bam = os.path.join(output_sub, "plasmid_alignment.bam")
            spliced_human_bam = os.path.join(output_sub, "spliced_human_alignment.bam")
            comparison_output = os.path.join(output_sub, "comparison_result")

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
    print(f"Execution Plan ({n_combos} combinations, {plan.total_steps} steps):")
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
) -> None:
    if sequencing_files_r1 is None:
        raise ValueError("-sf1 (sequencing_files_r1) is required.")

    # Build the plan regardless of dry-run
    plan = build_plan(
        human_fasta=human_fasta,
        plasmid_files=plasmid_files,
        output_folder=output_folder,
        overwrite=overwrite,
        sequencing_files_r1=sequencing_files_r1,
        sequencing_files_r2=sequencing_files_r2,
        cdna_output=cdna_output,
    )

    if dry_run:
        print_plan(plan)
        return

    logging.info("Starting the pipeline...")
    logging.info("Input parameters:")
    logging.info(f"Human reference FASTA: {human_fasta}")
    logging.info(f"Plasmid files: {plasmid_files}")
    logging.info(f"Output folder: {output_folder}")
    logging.info(f"Sequencing inputs: {len(plan.sequencing_inputs)} entries")

    # Perform quality control checks
    logging.info("Performing quality control checks on input files...")
    quality_control([human_fasta], expected_formats=["fasta"])
    quality_control(plan.plasmid_files, expected_formats=["genbank", "xdna"])

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

    with _progress_context(progress, plan.total_steps) as prog:
        for plasmid_file in plan.plasmid_files:
            plasmid_file_type: str = (
                "genbank"
                if plasmid_file.endswith(".gb") or plasmid_file.endswith(".gbk")
                else "xdna"
            )
            logging.info(f"Processing plasmid file: {plasmid_file}")

            for seq_input in plan.sequencing_inputs:
                sequencing_file = seq_input.file1
                fastq2 = seq_input.file2
                logging.info(f"Processing sequencing file: {sequencing_file}")

                bam_basename = sanitize_filename(
                    os.path.splitext(os.path.basename(sequencing_file))[0]
                )
                file_basename = sanitize_filename(
                    os.path.splitext(os.path.basename(plasmid_file))[0]
                )
                output_subfolder = os.path.join(output_folder, bam_basename, file_basename)

                os.makedirs(output_subfolder, exist_ok=True)

                # Step 1: Convert the plasmid file to a FASTA file
                plasmid_fasta = os.path.join(
                    output_subfolder,
                    os.path.splitext(os.path.basename(plasmid_file))[0] + ".fasta",
                )
                if not os.path.exists(plasmid_fasta) or overwrite:
                    logging.info(f"Converting {plasmid_file} to FASTA format.")
                    convert_plasmidfile_to_fasta(
                        plasmid_file,
                        plasmid_fasta,
                        plasmid_file_type,
                        shift_bases,
                        generate_shifted,
                        overwrite,
                    )
                if md5_level in ["all", "intermediate"]:
                    write_md5sum(plasmid_fasta, "intermediate", output_subfolder)
                if md5_level in ["all", "input"]:
                    write_md5sum(plasmid_file, "input", output_subfolder)
                prog.advance()

                # Step 2: Generate indices
                human_index = os.path.join(
                    os.path.dirname(human_fasta),
                    os.path.splitext(os.path.basename(human_fasta))[0] + ".mmi",
                )
                plasmid_index = os.path.join(
                    output_subfolder,
                    os.path.splitext(os.path.basename(plasmid_fasta))[0] + ".mmi",
                )
                create_indexes(human_fasta, overwrite)
                create_indexes(plasmid_fasta, overwrite)
                if md5_level in ["all", "intermediate"]:
                    write_md5sum(plasmid_index, "intermediate", output_subfolder)
                if md5_level in ["all"]:
                    write_md5sum(human_index, "intermediate", output_subfolder)
                prog.advance()

                # Step 3: Spliced alignment
                spliced_bam = os.path.join(output_subfolder, "spliced_alignment.bam")
                spliced_fasta = os.path.join(output_subfolder, "spliced_reference.fasta")
                logging.info(
                    f"Performing spliced alignment for {sequencing_file} and {plasmid_file}..."
                )
                spanned_regions = spliced_alignment(
                    human_index, plasmid_fasta, spliced_bam, padding
                )
                logging.info("Extracting human reference regions...")
                extract_human_reference(human_fasta, spanned_regions, spliced_fasta)
                if md5_level in ["all", "intermediate"]:
                    write_md5sum(spliced_bam, "intermediate", output_subfolder)
                    write_md5sum(spliced_fasta, "intermediate", output_subfolder)

                spliced_index = os.path.join(
                    output_subfolder,
                    os.path.splitext(os.path.basename(spliced_fasta))[0] + ".mmi",
                )
                create_indexes(spliced_fasta, overwrite)
                if md5_level in ["all", "intermediate"]:
                    write_md5sum(spliced_index, "intermediate", output_subfolder)
                prog.advance()

                # Step 4: cDNA positions
                unique_cdna_output = cdna_output or os.path.join(
                    output_subfolder, "cDNA_positions.txt"
                )
                logging.info("Extracting cDNA positions from the spliced alignment...")
                extract_plasmid_cdna_positions(plasmid_fasta, spliced_bam, unique_cdna_output)
                if md5_level in ["all", "intermediate"]:
                    write_md5sum(unique_cdna_output, "intermediate", output_subfolder)
                prog.advance()

                # Step 5: Align reads
                plasmid_bam = os.path.join(output_subfolder, "plasmid_alignment.bam")
                spliced_human_bam = os.path.join(output_subfolder, "spliced_human_alignment.bam")
                logging.info("Aligning reads to the plasmid and spliced human reference...")
                align_reads(plasmid_index, sequencing_file, plasmid_bam, "plasmid", fastq2)
                align_reads(spliced_index, sequencing_file, spliced_human_bam, "human", fastq2)
                if md5_level in ["all", "intermediate"]:
                    write_md5sum(plasmid_bam, "intermediate", output_subfolder)
                    write_md5sum(spliced_human_bam, "intermediate", output_subfolder)
                prog.advance()

                # Step 6: Compare alignments
                comparison_output = os.path.join(output_subfolder, "comparison_result")
                logging.info("Comparing the two alignments...")
                compare_alignments(plasmid_bam, spliced_human_bam, comparison_output)
                prog.advance()

                # Step 7: Generate report
                logging.info("Generating report...")
                reads_assignment_file = f"{comparison_output}.reads_assignment.tsv"
                summary_file = f"{comparison_output}.summary.tsv"
                command_line = " ".join(sys.argv)
                generate_report(
                    reads_assignment_file,
                    summary_file,
                    output_subfolder,
                    threshold=threshold,
                    human_fasta=human_fasta,
                    plasmid_gb=plasmid_file,
                    sequencing_file=sequencing_file,
                    command_line=command_line,
                    static_report=static_report,
                    plotly_mode=plotly_mode,
                    output_root=output_folder,
                )

                # Write MD5 checksums for the output files
                if md5_level in ["all", "intermediate", "output"]:
                    write_md5sum(reads_assignment_file, "output", output_subfolder)
                    write_md5sum(summary_file, "output", output_subfolder)
                prog.advance()

                # Step 8: Optionally delete intermediate files
                if not keep_intermediate:
                    logging.info("Deleting intermediate files...")
                    os.remove(plasmid_bam)
                    os.remove(spliced_human_bam)
                    os.remove(spliced_bam)

                # Step 9: Optionally archive the output folder
                if archive_output:
                    archive_output_folder(output_subfolder)

    logging.info("Pipeline completed successfully.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the full pipeline to detect and quantify plasmid DNA contamination in sequencing data"
    )
    parser.add_argument("-hf", "--human_fasta", help="Human reference FASTA file", required=True)
    parser.add_argument(
        "-pf",
        "--plasmid_files",
        help="Plasmid files (single file or a file containing paths to multiple files)",
        required=True,
    )
    parser.add_argument(
        "-sf1",
        "--sequencing_files_r1",
        help="Forward (R1) FASTQ/BAM file or file list (.txt)",
        required=True,
    )
    parser.add_argument(
        "-sf2",
        "--sequencing_files_r2",
        help="Reverse (R2) FASTQ files or file list (.txt) for paired-end",
        default=None,
    )
    parser.add_argument(
        "-o",
        "--output_folder",
        help="Folder to write all outputs and intermediate files",
        required=True,
    )
    parser.add_argument(
        "-k",
        "--keep_intermediate",
        type=bool,
        default=True,
        help="Keep intermediate files (default: True)",
    )
    parser.add_argument(
        "-sb",
        "--shift_bases",
        type=int,
        default=500,
        help="Number of bases to shift in the shifted reference (default: 500)",
    )
    parser.add_argument(
        "-g",
        "--generate_shifted",
        action="store_true",
        help="Generate a shifted reference sequence",
    )
    parser.add_argument(
        "-w", "--overwrite", action="store_true", help="Overwrite existing output files"
    )
    parser.add_argument(
        "-d",
        "--padding",
        type=int,
        default=1000,
        help="Padding to add to both sides of the spanned regions (default: 1000)",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help=f"Threshold for contamination verdict (default: {DEFAULT_THRESHOLD})",
    )
    parser.add_argument(
        "-md5",
        "--md5_level",
        type=str,
        choices=["all", "intermediate", "output"],
        default="intermediate",
        help="Level of MD5 checksum calculation (default: intermediate)",
    )
    parser.add_argument(
        "--cDNA_output",
        help="Output file for cDNA start and end positions in the plasmid reference",
        default=None,
    )
    parser.add_argument(
        "--archive_output",
        action="store_true",
        help="Archive and compress the output folder into a .tar.gz file",
    )
    parser.add_argument(
        "-n",
        "--dry_run",
        action="store_true",
        help="Show execution plan without running anything",
    )
    parser.add_argument(
        "--no_progress",
        action="store_true",
        help="Disable progress bar",
    )
    parser.add_argument(
        "--static-report",
        action="store_true",
        help="Generate static PNG reports alongside interactive HTML (opt-in, slower)",
    )
    parser.add_argument(
        "--plotly-mode",
        choices=["cdn", "directory", "embedded"],
        default="directory",
        help="Plotly.js inclusion mode for interactive reports (default: directory)",
    )
    add_logging_args(parser)

    args = parser.parse_args()

    configure_logging_from_args(args)

    run_pipeline(
        args.human_fasta,
        args.plasmid_files,
        args.output_folder,
        args.keep_intermediate,
        args.shift_bases,
        args.generate_shifted,
        args.overwrite,
        args.padding,
        args.threshold,
        args.md5_level,
        args.cDNA_output,
        args.archive_output,
        sequencing_files_r1=args.sequencing_files_r1,
        sequencing_files_r2=args.sequencing_files_r2,
        dry_run=args.dry_run,
        progress=sys.stderr.isatty() and not args.no_progress,
        static_report=args.static_report,
        plotly_mode=args.plotly_mode,
    )
