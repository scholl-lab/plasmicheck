"""Shared test fixtures for PlasmiCheck test suite."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Mock pysam when it is not installable (e.g. Windows — no wheels on PyPI).
# Source modules like compare_alignments.py and spliced_alignment.py have
# `import pysam` at the top level, which would fail at test-collection time.
# All unit tests already use mock read objects, so a MagicMock placeholder
# is sufficient for the modules to import.
# ---------------------------------------------------------------------------
try:
    import pysam  # noqa: F401
except ImportError:
    sys.modules["pysam"] = MagicMock()


@pytest.fixture
def sample_config() -> dict[str, Any]:
    """Load the project config.json."""
    config_path = Path(__file__).resolve().parent.parent / "plasmicheck" / "config.json"
    with open(config_path) as f:
        return json.load(f)


@pytest.fixture
def tmp_output_dir(tmp_path: Path) -> Path:
    """Provide a temporary output directory."""
    out = tmp_path / "output"
    out.mkdir()
    return out


@pytest.fixture
def sample_reads_assignment_df() -> pd.DataFrame:
    """DataFrame matching reads_assignment.tsv schema."""
    return pd.DataFrame(
        {
            "ReadID": ["read1", "read2", "read3", "read4", "read5"],
            "AssignedTo": ["Plasmid", "Human", "Tied", "Plasmid", "Human"],
            "PlasmidScore": [70, 30, 50, 65, 20],
            "HumanScore": [30, 70, 50, 25, 60],
            "PlasmidCIGAR": ["100M", "50M50S", "75M25S", "100M", "30M70S"],
            "HumanCIGAR": ["50M50S", "100M", "75M25S", "30M70S", "100M"],
            "PlasmidMapQ": [60, 20, 40, 55, 10],
            "HumanMapQ": [20, 60, 40, 15, 55],
        }
    )


@pytest.fixture
def sample_summary_df() -> pd.DataFrame:
    """DataFrame matching summary.tsv schema."""
    return pd.DataFrame(
        {
            "Category": [
                "Plasmid",
                "Human",
                "Tied",
                "Verdict",
                "Ratio",
                "CoverageOutsideINSERT",
                "MismatchesNearINSERT",
            ],
            "Count": [
                100,
                50,
                10,
                "Sample is contaminated with plasmid DNA",
                2.0,
                0.1234,
                "{'with_mismatches_or_clipping': 5, 'without_mismatches_or_clipping': 10}",
            ],
        }
    )


@pytest.fixture
def mock_pysam_read() -> Any:
    """Factory fixture for creating mock pysam aligned reads."""

    def _make_read(
        *,
        is_unmapped: bool = False,
        mapping_quality: int = 60,
        cigartuples: list[tuple[int, int]] | None = None,
        mate_is_unmapped: bool = False,
        nm_tag: int | None = 0,
        query_name: str = "test_read",
        cigarstring: str | None = "100M",
        reference_start: int = 0,
        reference_end: int = 100,
        query_length: int = 100,
    ) -> MagicMock:
        read = MagicMock()
        read.is_unmapped = is_unmapped
        read.mapping_quality = mapping_quality
        read.cigartuples = cigartuples if cigartuples is not None else [(0, 100)]
        read.mate_is_unmapped = mate_is_unmapped
        read.query_name = query_name
        read.cigarstring = cigarstring
        read.reference_start = reference_start
        read.reference_end = reference_end
        read.query_length = query_length

        if nm_tag is not None:
            read.has_tag.return_value = True
            read.get_tag.return_value = nm_tag
        else:
            read.has_tag.return_value = False
            read.get_tag.return_value = 0

        return read

    return _make_read


@pytest.fixture
def cdna_positions_file(tmp_path: Path) -> Path:
    """Write a temp cDNA_positions.txt and return its path."""
    f = tmp_path / "cDNA_positions.txt"
    f.write_text(
        "cDNA start position in plasmid: 100\n"
        "cDNA end position in plasmid: 500\n"
        "INSERT_REGION: (100, 500)\n"
    )
    return f


@pytest.fixture
def real_data_dir() -> Path | None:
    """Path to tests/data/real/ -- returns None if missing (skip tests that need it)."""
    p = Path(__file__).resolve().parent / "data" / "real"
    if p.is_dir():
        return p
    return None
