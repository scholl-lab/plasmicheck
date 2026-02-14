#!/usr/bin/env python3
"""Benchmark PlasmiCheck pipeline per-step performance.

This script times each pipeline step separately to identify performance bottlenecks
and validate optimization efforts in Phases 5-7 of v0.32.0.

Usage:
    python scripts/benchmark.py                     # 3 iterations on synthetic dataset
    python scripts/benchmark.py -n 5                # 5 iterations
    python scripts/benchmark.py -o results.md       # Custom output file
    python scripts/benchmark.py --dataset-dir PATH  # Custom dataset
    python scripts/benchmark.py --steps convert,align  # Benchmark specific steps
"""

import argparse
import logging
import os
import platform
import shutil
import statistics
import sys
import tempfile
import time
from collections.abc import Generator
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

# Add parent directory to path for imports
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from plasmicheck.scripts.align_reads import align_reads  # noqa: E402
from plasmicheck.scripts.compare_alignments import compare_alignments  # noqa: E402
from plasmicheck.scripts.convert_plasmidfile_to_fasta import (  # noqa: E402
    convert_plasmidfile_to_fasta,
)
from plasmicheck.scripts.create_indexes import create_indexes  # noqa: E402
from plasmicheck.scripts.generate_report import main as generate_report  # noqa: E402
from plasmicheck.scripts.spliced_alignment import (  # noqa: E402
    extract_human_reference,
    extract_plasmid_cdna_positions,
    spliced_alignment,
)
from plasmicheck.scripts.utils import sanitize_filename  # noqa: E402
from plasmicheck.version import __version__  # noqa: E402


@contextmanager
def timed_step(name: str, timings: dict[str, float]) -> Generator[None, None, None]:
    """Context manager to time a pipeline step."""
    start = time.perf_counter()
    try:
        yield
    finally:
        timings[name] = time.perf_counter() - start


def check_external_tools() -> None:
    """Check that required external tools are available."""
    missing = []
    for tool in ["minimap2", "samtools"]:
        if shutil.which(tool) is None:
            missing.append(tool)

    if missing:
        print(f"ERROR: Required tools not found on PATH: {', '.join(missing)}")
        print("\nPlease install the following tools:")
        print("  - minimap2: https://github.com/lh3/minimap2")
        print("  - samtools: http://www.htslib.org/download/")
        sys.exit(1)


def count_reads(fastq_path: str) -> int:
    """Count reads in a FASTQ file."""
    with open(fastq_path) as f:
        return sum(1 for line in f) // 4


def run_pipeline_step(
    step_name: str,
    human_fasta: str,
    plasmid_gb: str,
    fastq_r1: str,
    fastq_r2: str | None,
    output_dir: str,
    timings: dict[str, float],
    enabled_steps: set[str] | None = None,
) -> None:
    """Run all pipeline steps with per-step timing.

    Args:
        step_name: Name of this timing run (for logging)
        human_fasta: Path to human reference FASTA
        plasmid_gb: Path to plasmid GenBank file
        fastq_r1: Path to R1 FASTQ
        fastq_r2: Path to R2 FASTQ (or None for single-end)
        output_dir: Output directory for this iteration
        timings: Dictionary to store timing results
        enabled_steps: Set of steps to time (None = all steps)
    """
    # Suppress pipeline logging during benchmark
    logging.getLogger().setLevel(logging.WARNING)

    # Compute paths matching run_pipeline.py logic
    bam_basename = sanitize_filename(os.path.splitext(os.path.basename(fastq_r1))[0])
    file_basename = sanitize_filename(os.path.splitext(os.path.basename(plasmid_gb))[0])
    output_subfolder = os.path.join(output_dir, bam_basename, file_basename)
    os.makedirs(output_subfolder, exist_ok=True)

    plasmid_fasta = os.path.join(
        output_subfolder, os.path.splitext(os.path.basename(plasmid_gb))[0] + ".fasta"
    )
    human_index = os.path.join(
        os.path.dirname(human_fasta), os.path.splitext(os.path.basename(human_fasta))[0] + ".mmi"
    )
    plasmid_index = os.path.join(
        output_subfolder, os.path.splitext(os.path.basename(plasmid_fasta))[0] + ".mmi"
    )
    spliced_bam = os.path.join(output_subfolder, "spliced_alignment.bam")
    spliced_fasta = os.path.join(output_subfolder, "spliced_reference.fasta")
    spliced_index = os.path.join(
        output_subfolder, os.path.splitext(os.path.basename(spliced_fasta))[0] + ".mmi"
    )
    cdna_positions = os.path.join(output_subfolder, "cDNA_positions.txt")
    plasmid_bam = os.path.join(output_subfolder, "plasmid_alignment.bam")
    spliced_human_bam = os.path.join(output_subfolder, "spliced_human_alignment.bam")
    comparison_output = os.path.join(output_subfolder, "comparison_result")

    def should_time(step: str) -> bool:
        return enabled_steps is None or step in enabled_steps

    # Step 1: Convert plasmid file to FASTA
    if should_time("convert"):
        with timed_step("convert", timings):
            convert_plasmidfile_to_fasta(
                plasmid_gb,
                plasmid_fasta,
                "genbank",
                shift_bases=500,
                generate_shifted=False,
                overwrite=True,
            )
    else:
        # Still run but don't time
        convert_plasmidfile_to_fasta(
            plasmid_gb, plasmid_fasta, "genbank", shift_bases=500, generate_shifted=False, overwrite=True
        )

    # Step 2: Create indexes
    if should_time("index"):
        with timed_step("index", timings):
            create_indexes(human_fasta, overwrite=True)
            create_indexes(plasmid_fasta, overwrite=True)
    else:
        create_indexes(human_fasta, overwrite=True)
        create_indexes(plasmid_fasta, overwrite=True)

    # Step 3: Spliced alignment (includes extract_human_reference, create_indexes, extract_cdna_positions)
    if should_time("spliced"):
        with timed_step("spliced", timings):
            spanned_regions = spliced_alignment(human_index, plasmid_fasta, spliced_bam, padding=1000)
            extract_human_reference(human_fasta, spanned_regions, spliced_fasta)
            create_indexes(spliced_fasta, overwrite=True)
            extract_plasmid_cdna_positions(plasmid_fasta, spliced_bam, cdna_positions)
    else:
        spanned_regions = spliced_alignment(human_index, plasmid_fasta, spliced_bam, padding=1000)
        extract_human_reference(human_fasta, spanned_regions, spliced_fasta)
        create_indexes(spliced_fasta, overwrite=True)
        extract_plasmid_cdna_positions(plasmid_fasta, spliced_bam, cdna_positions)

    # Step 4: Align reads to plasmid and spliced human reference
    if should_time("align"):
        with timed_step("align", timings):
            align_reads(plasmid_index, fastq_r1, plasmid_bam, "plasmid", fastq2=fastq_r2)
            align_reads(spliced_index, fastq_r1, spliced_human_bam, "human", fastq2=fastq_r2)
    else:
        align_reads(plasmid_index, fastq_r1, plasmid_bam, "plasmid", fastq2=fastq_r2)
        align_reads(spliced_index, fastq_r1, spliced_human_bam, "human", fastq2=fastq_r2)

    # Step 5: Compare alignments
    if should_time("compare"):
        with timed_step("compare", timings):
            compare_alignments(plasmid_bam, spliced_human_bam, comparison_output, threshold=0.8)
    else:
        compare_alignments(plasmid_bam, spliced_human_bam, comparison_output, threshold=0.8)

    # Step 6: Generate report
    if should_time("report"):
        reads_assignment_file = f"{comparison_output}.reads_assignment.tsv"
        summary_file = f"{comparison_output}.summary.tsv"
        with timed_step("report", timings):
            generate_report(
                reads_assignment_file,
                summary_file,
                output_subfolder,
                threshold=0.8,
                human_fasta=human_fasta,
                plasmid_gb=plasmid_gb,
                sequencing_file=fastq_r1,
                command_line="benchmark",
            )
    else:
        reads_assignment_file = f"{comparison_output}.reads_assignment.tsv"
        summary_file = f"{comparison_output}.summary.tsv"
        generate_report(
            reads_assignment_file,
            summary_file,
            output_subfolder,
            threshold=0.8,
            human_fasta=human_fasta,
            plasmid_gb=plasmid_gb,
            sequencing_file=fastq_r1,
            command_line="benchmark",
        )


def run_benchmark(
    dataset_dir: str, iterations: int, enabled_steps: set[str] | None = None
) -> dict[str, list[float]]:
    """Run benchmark iterations and collect timing data.

    Args:
        dataset_dir: Path to dataset directory
        iterations: Number of iterations to run
        enabled_steps: Set of steps to time (None = all steps)

    Returns:
        Dictionary mapping step name to list of timings
    """
    # Locate dataset files
    human_fasta = os.path.join(dataset_dir, "human_ref.fasta")
    plasmid_gb = os.path.join(dataset_dir, "plasmid.gb")
    fastq_r1 = os.path.join(dataset_dir, "contaminated_R1.fastq")
    fastq_r2_path = os.path.join(dataset_dir, "contaminated_R2.fastq")

    # Validate dataset
    for path in [human_fasta, plasmid_gb, fastq_r1]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Dataset file not found: {path}")

    # Check if R2 exists (optional for single-end)
    fastq_r2: str | None = fastq_r2_path if os.path.exists(fastq_r2_path) else None

    # Count reads for reporting
    read_count = count_reads(fastq_r1)

    print(f"Dataset: {dataset_dir}")
    print(f"Read count: {read_count}")
    print(f"Iterations: {iterations} (+ 1 warm-up)")
    if enabled_steps:
        print(f"Steps: {', '.join(sorted(enabled_steps))}")
    print()

    # Collect timings per step
    all_timings: dict[str, list[float]] = {
        "convert": [],
        "index": [],
        "spliced": [],
        "align": [],
        "compare": [],
        "report": [],
        "total": [],
    }

    # Warm-up iteration (not counted)
    print("Running warm-up iteration...", end=" ", flush=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        warmup_timings: dict[str, float] = {}
        run_pipeline_step(
            "warmup",
            human_fasta,
            plasmid_gb,
            fastq_r1,
            fastq_r2,
            tmpdir,
            warmup_timings,
            enabled_steps,
        )
    print("done")

    # Timed iterations
    for i in range(1, iterations + 1):
        print(f"Running iteration {i}/{iterations}...", end=" ", flush=True)
        with tempfile.TemporaryDirectory() as tmpdir:
            iteration_timings: dict[str, float] = {}
            iter_start = time.perf_counter()

            run_pipeline_step(
                f"iteration-{i}",
                human_fasta,
                plasmid_gb,
                fastq_r1,
                fastq_r2,
                tmpdir,
                iteration_timings,
                enabled_steps,
            )

            iter_total = time.perf_counter() - iter_start
            iteration_timings["total"] = iter_total

            # Record timings
            for step, timing in iteration_timings.items():
                all_timings[step].append(timing)

        print(f"done ({iter_total:.2f}s)")

    print()
    return all_timings


def format_benchmark_report(
    timings: dict[str, list[float]], dataset_dir: str, iterations: int, enabled_steps: set[str] | None
) -> str:
    """Format benchmark results as Markdown.

    Args:
        timings: Dictionary of step timings
        dataset_dir: Path to dataset directory
        iterations: Number of iterations run
        enabled_steps: Set of steps that were timed

    Returns:
        Markdown formatted report
    """
    # Count reads for report
    fastq_r1 = os.path.join(dataset_dir, "contaminated_R1.fastq")
    read_count = count_reads(fastq_r1)

    # Compute statistics
    def stats(values: list[float]) -> tuple[float, float, float, float]:
        if not values:
            return (0.0, 0.0, 0.0, 0.0)
        mean_val = statistics.mean(values)
        std_val = statistics.stdev(values) if len(values) > 1 else 0.0
        min_val = min(values)
        max_val = max(values)
        return (mean_val, std_val, min_val, max_val)

    # Build report
    lines = []
    lines.append("# PlasmiCheck Benchmark Results")
    lines.append("")
    lines.append(f"**Version:** {__version__}")
    lines.append(f"**Date:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"**Dataset:** {os.path.basename(dataset_dir)} ({read_count} reads)")
    lines.append(f"**Iterations:** {iterations} (+ 1 warm-up)")
    lines.append(f"**System:** {platform.system()} {platform.machine()}")
    if enabled_steps:
        lines.append(f"**Steps:** {', '.join(sorted(enabled_steps))}")
    lines.append("")
    lines.append("## Per-Step Timing")
    lines.append("")
    lines.append("| Step | Mean (s) | Std Dev (s) | Min (s) | Max (s) | % of Total |")
    lines.append("|------|----------|-------------|---------|---------|------------|")

    # Calculate total for percentage computation
    total_mean = statistics.mean(timings["total"]) if timings["total"] else 0.0

    # Report each step
    step_order = ["convert", "index", "spliced", "align", "compare", "report"]
    for step in step_order:
        mean_val, std_val, min_val, max_val = stats(timings[step])
        pct = (mean_val / total_mean * 100) if total_mean > 0 else 0.0

        # Mark steps that weren't timed
        if enabled_steps and step not in enabled_steps:
            lines.append(
                f"| {step} | — | — | — | — | — |"
            )
        else:
            lines.append(
                f"| {step} | {mean_val:.3f} | {std_val:.3f} | {min_val:.3f} | {max_val:.3f} | {pct:.1f}% |"
            )

    # Total row
    total_mean, total_std, total_min, total_max = stats(timings["total"])
    lines.append(
        f"| **Total** | **{total_mean:.3f}** | **{total_std:.3f}** | **{total_min:.3f}** | **{total_max:.3f}** | **100%** |"
    )

    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- Warm-up iteration excluded from statistics")
    lines.append("- Index generation includes both human and plasmid references")
    lines.append("- Spliced step includes alignment, extraction, indexing, and cDNA position extraction")
    lines.append("- Align step includes both plasmid and spliced human alignments")
    lines.append("- Timings include subprocess overhead (minimap2, samtools)")
    if enabled_steps:
        lines.append("- Steps marked with — were not timed (but still executed as dependencies)")
    lines.append("")

    return "\n".join(lines)


def main() -> None:
    """Main entry point for benchmark script."""
    parser = argparse.ArgumentParser(
        description="Benchmark PlasmiCheck pipeline per-step performance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                           # 3 iterations on synthetic dataset
  %(prog)s -n 5                      # 5 iterations
  %(prog)s -o results.md             # Custom output file
  %(prog)s --dataset-dir PATH        # Custom dataset
  %(prog)s --steps convert,align     # Benchmark specific steps
        """,
    )

    default_dataset = str(REPO_ROOT / "tests" / "data" / "synthetic")

    parser.add_argument(
        "-n",
        "--iterations",
        type=int,
        default=3,
        help="Number of iterations to run (default: 3)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="BENCHMARK.md",
        help="Output file path (default: BENCHMARK.md)",
    )
    parser.add_argument(
        "--dataset-dir",
        default=default_dataset,
        help=f"Dataset directory (default: {default_dataset})",
    )
    parser.add_argument(
        "--steps",
        help="Comma-separated list of steps to time (default: all). Valid: convert,index,spliced,align,compare,report",
    )

    args = parser.parse_args()

    # Check external tools
    check_external_tools()

    # Parse step filter
    enabled_steps: set[str] | None = None
    if args.steps:
        valid_steps = {"convert", "index", "spliced", "align", "compare", "report"}
        enabled_steps = {s.strip() for s in args.steps.split(",")}
        invalid = enabled_steps - valid_steps
        if invalid:
            print(f"ERROR: Invalid steps: {', '.join(invalid)}")
            print(f"Valid steps: {', '.join(sorted(valid_steps))}")
            sys.exit(1)

    # Run benchmark
    try:
        timings = run_benchmark(args.dataset_dir, args.iterations, enabled_steps)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    # Generate and write report
    report = format_benchmark_report(timings, args.dataset_dir, args.iterations, enabled_steps)

    with open(args.output, "w") as f:
        f.write(report)

    print(f"Benchmark results written to: {args.output}")


if __name__ == "__main__":
    main()
