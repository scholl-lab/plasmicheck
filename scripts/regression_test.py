#!/usr/bin/env python3
"""Regression test script for PlasmiCheck pipeline.

This standalone script validates that pipeline outputs (contamination ratios,
read assignments, verdicts) remain unchanged after code modifications.

Usage:
    python scripts/regression_test.py --generate    # Create baseline from v0.31.0
    python scripts/regression_test.py               # Compare against baseline
    python scripts/regression_test.py --baseline-dir PATH  # Custom cache location

Note: First run may be slower due to index generation.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


def get_repo_root() -> Path:
    """Get repository root (parent of scripts directory)."""
    return Path(__file__).resolve().parent.parent


def get_synthetic_data_dir() -> Path:
    """Get path to synthetic test data directory."""
    return get_repo_root() / "tests" / "data" / "synthetic"


def check_external_tools() -> tuple[bool, str]:
    """Check if minimap2 and samtools are available.

    Returns:
        Tuple of (success, error_message)
    """
    minimap2_available = shutil.which("minimap2") is not None
    samtools_available = shutil.which("samtools") is not None

    if not minimap2_available and not samtools_available:
        return False, "minimap2 and samtools not found on PATH"
    elif not minimap2_available:
        return False, "minimap2 not found on PATH"
    elif not samtools_available:
        return False, "samtools not found on PATH"

    return True, ""


def run_pipeline_scenario(
    scenario_name: str,
    r1_filename: str,
    r2_filename: str,
    output_dir: Path,
) -> dict[str, str]:
    """Run pipeline for a single scenario and extract TSV outputs.

    Args:
        scenario_name: Name of scenario ("contaminated" or "not_contaminated")
        r1_filename: Read 1 FASTQ filename
        r2_filename: Read 2 FASTQ filename
        output_dir: Temporary output directory

    Returns:
        Dictionary with keys "summary" and "reads_assignment" containing TSV text
    """
    from plasmicheck.scripts.run_pipeline import run_pipeline

    synthetic_data = get_synthetic_data_dir()

    human_fasta = str(synthetic_data / "human_ref.fasta")
    plasmid_gb = str(synthetic_data / "plasmid.gb")
    r1 = str(synthetic_data / r1_filename)
    r2 = str(synthetic_data / r2_filename)

    # Suppress pipeline logging during runs
    logging.getLogger("plasmicheck").setLevel(logging.WARNING)

    run_pipeline(
        human_fasta=human_fasta,
        plasmid_files=plasmid_gb,
        output_folder=str(output_dir),
        keep_intermediate=True,
        overwrite=True,
        threshold=0.8,
        sequencing_files_r1=r1,
        sequencing_files_r2=r2,
    )

    # Extract output files
    summary_files = list(output_dir.rglob("*.summary.tsv"))
    reads_assignment_files = list(output_dir.rglob("*.reads_assignment.tsv"))

    if not summary_files:
        raise RuntimeError(f"No summary.tsv found for {scenario_name}")
    if not reads_assignment_files:
        raise RuntimeError(f"No reads_assignment.tsv found for {scenario_name}")

    results = {
        "summary": summary_files[0].read_text(),
        "reads_assignment": reads_assignment_files[0].read_text(),
    }

    return results


def generate_baseline(baseline_dir: Path) -> None:
    """Generate baseline from current pipeline output.

    Args:
        baseline_dir: Directory to store baseline.json
    """
    from plasmicheck.version import __version__

    print(f"Generating baseline with PlasmiCheck v{__version__}...")

    baseline_dir.mkdir(parents=True, exist_ok=True)

    baseline_data: dict[str, Any] = {
        "version": __version__,
        "scenarios": {},
    }

    # Run both scenarios
    scenarios = [
        ("contaminated", "contaminated_R1.fastq", "contaminated_R2.fastq"),
        ("not_contaminated", "not_contaminated_R1.fastq", "not_contaminated_R2.fastq"),
    ]

    for scenario_name, r1_file, r2_file in scenarios:
        print(f"  Running scenario: {scenario_name}...")
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir)
            results = run_pipeline_scenario(scenario_name, r1_file, r2_file, output_dir)
            baseline_data["scenarios"][scenario_name] = results

    # Atomic write using tempfile + os.replace
    fd, temp_path = tempfile.mkstemp(dir=baseline_dir, suffix=".json")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(baseline_data, f, indent=2)

        baseline_file = baseline_dir / "baseline.json"
        os.replace(temp_path, baseline_file)
        print(f"\nBaseline saved to: {baseline_file}")
        print(f"Version: {__version__}")
    except Exception:
        # Clean up temp file on error
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        raise


def compare_summary_tsv(baseline_text: str, current_text: str) -> list[str]:
    """Compare summary.tsv files.

    Args:
        baseline_text: Baseline TSV content
        current_text: Current TSV content

    Returns:
        List of error messages (empty if match)
    """
    errors: list[str] = []

    # Parse with pandas
    from io import StringIO
    baseline_df = pd.read_csv(StringIO(baseline_text), sep="\t", dtype=str)
    current_df = pd.read_csv(StringIO(current_text), sep="\t", dtype=str)

    if len(baseline_df) != len(current_df):
        errors.append(f"Row count mismatch: baseline={len(baseline_df)}, current={len(current_df)}")
        return errors

    # Compare each row
    for baseline_row, current_row in zip(baseline_df.itertuples(), current_df.itertuples(), strict=True):
        category = baseline_row.Category

        if str(category).lower() == "verdict":
            # Exact string match (case-insensitive)
            baseline_val = str(baseline_row.Count).strip().lower()
            current_val = str(current_row.Count).strip().lower()
            if baseline_val != current_val:
                errors.append(f"Verdict mismatch: baseline='{baseline_val}', current='{current_val}'")

        elif category in ["Plasmid", "Human", "Tied"]:
            # Exact integer match
            baseline_int = int(baseline_row.Count)  # type: ignore[arg-type]
            current_int = int(current_row.Count)  # type: ignore[arg-type]
            if baseline_int != current_int:
                errors.append(f"{category} count mismatch: baseline={baseline_int}, current={current_int}")

        elif category in ["Ratio", "CoverageOutsideINSERT"]:
            # Float comparison with tolerance
            baseline_float = float(baseline_row.Count)  # type: ignore[arg-type]
            current_float = float(current_row.Count)  # type: ignore[arg-type]
            if abs(baseline_float - current_float) > 0.001:
                errors.append(
                    f"{category} divergence: baseline={baseline_float:.6f}, "
                    f"current={current_float:.6f}, diff={abs(baseline_float - current_float):.6f}"
                )

        elif category == "MismatchesNearINSERT":
            # String match (it's a dict representation)
            baseline_val = str(baseline_row.Count).strip()
            current_val = str(current_row.Count).strip()
            if baseline_val != current_val:
                errors.append(f"MismatchesNearINSERT mismatch: baseline='{baseline_val}', current='{current_val}'")

    return errors


def compare_reads_assignment_tsv(baseline_text: str, current_text: str) -> list[str]:
    """Compare reads_assignment.tsv files.

    Args:
        baseline_text: Baseline TSV content
        current_text: Current TSV content

    Returns:
        List of error messages (empty if match)
    """
    errors: list[str] = []

    # Parse with pandas
    from io import StringIO
    baseline_df = pd.read_csv(StringIO(baseline_text), sep="\t")
    current_df = pd.read_csv(StringIO(current_text), sep="\t")

    if len(baseline_df) != len(current_df):
        errors.append(f"Row count mismatch: baseline={len(baseline_df)}, current={len(current_df)}")
        return errors

    # Compare AssignedTo column (exact match)
    baseline_assigned = baseline_df["AssignedTo"].tolist()
    current_assigned = current_df["AssignedTo"].tolist()

    mismatches = sum(1 for b, c in zip(baseline_assigned, current_assigned, strict=True) if b != c)
    if mismatches > 0:
        errors.append(f"AssignedTo mismatch: {mismatches}/{len(baseline_df)} reads differ")

    # Compare PlasmidScore and HumanScore (float tolerance)
    for col in ["PlasmidScore", "HumanScore"]:
        if col in baseline_df.columns and col in current_df.columns:
            baseline_scores = baseline_df[col].astype(float)
            current_scores = current_df[col].astype(float)

            max_diff = (baseline_scores - current_scores).abs().max()
            if max_diff > 0.001:
                errors.append(f"{col} max difference: {max_diff:.6f} (tolerance: 0.001)")

    return errors


def compare_against_baseline(baseline_dir: Path) -> bool:
    """Compare current pipeline output against baseline.

    Args:
        baseline_dir: Directory containing baseline.json

    Returns:
        True if all tests pass, False otherwise
    """
    from plasmicheck.version import __version__

    baseline_file = baseline_dir / "baseline.json"

    if not baseline_file.exists():
        print(f"ERROR: No baseline found at {baseline_file}")
        print("Run with --generate first to create a baseline.")
        return False

    # Load baseline
    with open(baseline_file) as f:
        baseline_data = json.load(f)

    baseline_version = baseline_data.get("version", "unknown")
    if baseline_version != __version__:
        print(f"WARNING: Baseline was generated with v{baseline_version}, current version is v{__version__}")
        print()

    print(f"Comparing against baseline (v{baseline_version})...")

    all_pass = True

    # Run both scenarios
    scenarios = [
        ("contaminated", "contaminated_R1.fastq", "contaminated_R2.fastq"),
        ("not_contaminated", "not_contaminated_R1.fastq", "not_contaminated_R2.fastq"),
    ]

    for scenario_name, r1_file, r2_file in scenarios:
        print(f"\n{'='*70}")
        print(f"Scenario: {scenario_name}")
        print('='*70)

        # Run pipeline
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir)
            current_results = run_pipeline_scenario(scenario_name, r1_file, r2_file, output_dir)

        baseline_results = baseline_data["scenarios"][scenario_name]

        # Compare summary.tsv
        print("Comparing summary.tsv...")
        summary_errors = compare_summary_tsv(
            baseline_results["summary"],
            current_results["summary"]
        )

        if summary_errors:
            print("  FAIL - Summary differences detected:")
            for error in summary_errors:
                print(f"    - {error}")
            all_pass = False
        else:
            print("  PASS - Summary matches baseline")

        # Compare reads_assignment.tsv
        print("Comparing reads_assignment.tsv...")
        reads_errors = compare_reads_assignment_tsv(
            baseline_results["reads_assignment"],
            current_results["reads_assignment"]
        )

        if reads_errors:
            print("  FAIL - Reads assignment differences detected:")
            for error in reads_errors:
                print(f"    - {error}")
            all_pass = False
        else:
            print("  PASS - Reads assignment matches baseline")

    print(f"\n{'='*70}")
    if all_pass:
        print("RESULT: All tests PASSED")
        print('='*70)
    else:
        print("RESULT: Some tests FAILED")
        print('='*70)

    return all_pass


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Regression test for PlasmiCheck pipeline",
        epilog="Note: First run may be slower due to index generation."
    )
    parser.add_argument(
        "--generate",
        action="store_true",
        help="Generate baseline from current pipeline output"
    )
    parser.add_argument(
        "--baseline-dir",
        type=Path,
        default=get_repo_root() / "scripts" / ".regression_cache",
        help="Directory to store/load baseline (default: scripts/.regression_cache/)"
    )

    args = parser.parse_args()

    # Check for external tools
    tools_available, error_msg = check_external_tools()
    if not tools_available:
        print(f"ERROR: {error_msg}")
        print("\nRegression tests require minimap2 and samtools to be installed and available on PATH.")
        print("Please install these tools and try again.")
        return 1

    if args.generate:
        try:
            generate_baseline(args.baseline_dir)
            return 0
        except Exception as e:
            print(f"\nERROR: Failed to generate baseline: {e}")
            return 1
    else:
        try:
            success = compare_against_baseline(args.baseline_dir)
            return 0 if success else 1
        except Exception as e:
            print(f"\nERROR: Failed to compare against baseline: {e}")
            import traceback
            traceback.print_exc()
            return 1


if __name__ == "__main__":
    sys.exit(main())
