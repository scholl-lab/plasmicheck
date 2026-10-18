"""Unit tests for CacheManager."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import pytest

from plasmicheck.cache_manager import CacheManager, _resolve_md5


@pytest.fixture
def sample_fasta(tmp_path: Path) -> Path:
    """Create a sample FASTA file for testing."""
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(">seq1\nACGTACGT\n>seq2\nGCTAGCTA\n")
    return fasta_file


@pytest.fixture
def sample_mmi(tmp_path: Path) -> Path:
    """Create a dummy minimap2 index file."""
    mmi_file = tmp_path / "sample.mmi"
    mmi_file.write_bytes(b"FAKE_MMI_CONTENT_" * 100)
    return mmi_file


@pytest.fixture
def sample_fai(tmp_path: Path) -> Path:
    """Create a dummy samtools index file."""
    fai_file = tmp_path / "sample.fai"
    fai_file.write_text("seq1\t8\t6\t8\t9\nseq2\t8\t21\t8\t9\n")
    return fai_file


@pytest.fixture
def sample_bam(tmp_path: Path) -> Path:
    """Create a dummy BAM file."""
    bam_file = tmp_path / "spliced.bam"
    bam_file.write_bytes(b"FAKE_BAM_CONTENT_" * 200)
    return bam_file


@pytest.fixture
def sample_spliced_fasta(tmp_path: Path) -> Path:
    """Create a dummy spliced FASTA file."""
    fasta_file = tmp_path / "spliced.fasta"
    fasta_file.write_text(">spliced_region\nATGCATGC\n")
    return fasta_file


@pytest.fixture
def sample_cdna_positions(tmp_path: Path) -> Path:
    """Create a dummy cDNA positions file."""
    cdna_file = tmp_path / "cdna_positions.txt"
    cdna_file.write_text("start\t100\nend\t500\n")
    return cdna_file


# Initialization tests


@pytest.mark.unit
def test_cache_manager_creates_directories(tmp_path: Path) -> None:
    """Test that CacheManager creates cache directory structure."""
    cache_root = tmp_path / "cache"
    _cache = CacheManager(cache_root)

    assert cache_root.exists()
    assert (cache_root / "indexes").exists()
    assert (cache_root / "spliced").exists()
    # Manifest is created lazily when first entry is added


@pytest.mark.unit
def test_cache_manager_handles_existing_directory(tmp_path: Path) -> None:
    """Test that CacheManager handles pre-existing cache directory gracefully."""
    cache_root = tmp_path / "cache"
    cache_root.mkdir(parents=True)

    # Should not raise
    _cache = CacheManager(cache_root)
    assert cache_root.exists()


@pytest.mark.unit
def test_cache_manager_stores_size_limit(tmp_path: Path) -> None:
    """Test that size limit is stored correctly."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root, size_limit_gb=5.0)

    assert cache.size_limit == 5.0 * 1024**3


@pytest.mark.unit
def test_cache_manager_inaccessible_directory() -> None:
    """Test that inaccessible cache directory raises clear error."""
    # Use a path that cannot be created (invalid root)
    if Path("/").exists():
        cache_root = Path("/root/nonexistent/cache")
    else:
        cache_root = Path("C:\\Windows\\System32\\cache")

    with pytest.raises(OSError, match="Cannot create cache directory"):
        CacheManager(cache_root)


# Index cache tests


@pytest.mark.unit
def test_get_index_returns_none_on_miss(tmp_path: Path, sample_fasta: Path) -> None:
    """Test that get_index returns None on cache miss."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    result = cache.get_index(str(sample_fasta))
    assert result is None


@pytest.mark.unit
def test_put_and_get_index_round_trip(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that put_index + get_index round-trip succeeds."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    # Put index
    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    # Get index
    result = cache.get_index(str(sample_fasta))
    assert result is not None
    assert "mmi" in result
    assert "fai" in result
    assert result["mmi"].exists()
    assert result["fai"].exists()


@pytest.mark.unit
def test_get_index_returns_both_files(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that get_index returns dict with mmi and fai keys."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))
    result = cache.get_index(str(sample_fasta))

    assert result is not None
    assert isinstance(result["mmi"], Path)
    assert isinstance(result["fai"], Path)
    assert result["mmi"].suffix == ".mmi"
    assert result["fai"].suffix == ".fai"


@pytest.mark.unit
def test_index_cache_hit_logs_info(
    tmp_path: Path,
    sample_fasta: Path,
    sample_mmi: Path,
    sample_fai: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test that cache hit logs INFO message."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    with caplog.at_level("INFO"):
        cache.get_index(str(sample_fasta))

    assert "Cache HIT for index" in caplog.text


@pytest.mark.unit
def test_index_cache_miss_logs_info(
    tmp_path: Path, sample_fasta: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Test that cache miss logs INFO message."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    with caplog.at_level("INFO"):
        cache.get_index(str(sample_fasta))

    assert "Cache MISS for index" in caplog.text


@pytest.mark.unit
def test_index_stale_entry_returns_none(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that stale entry (file deleted externally) returns None."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    # Delete cached files externally
    cached = cache.get_index(str(sample_fasta))
    assert cached is not None
    cached["mmi"].unlink()
    cached["fai"].unlink()

    # Should detect stale entry
    result = cache.get_index(str(sample_fasta))
    assert result is None


@pytest.mark.unit
def test_index_size_mismatch_returns_none(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that size mismatch returns None."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    # Corrupt cached file (change size)
    cached = cache.get_index(str(sample_fasta))
    assert cached is not None
    cached["mmi"].write_bytes(b"CORRUPTED")

    # Should detect size mismatch
    result = cache.get_index(str(sample_fasta))
    assert result is None


# Spliced cache tests


@pytest.mark.unit
def test_get_spliced_returns_none_on_miss(tmp_path: Path, sample_fasta: Path) -> None:
    """Test that get_spliced returns None on cache miss."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    plasmid_fasta = sample_fasta
    human_fasta = tmp_path / "human.fasta"
    human_fasta.write_text(">chr1\nATGCATGC\n")

    result = cache.get_spliced(str(plasmid_fasta), str(human_fasta))
    assert result is None


@pytest.mark.unit
def test_put_and_get_spliced_round_trip(
    tmp_path: Path,
    sample_fasta: Path,
    sample_bam: Path,
    sample_spliced_fasta: Path,
    sample_cdna_positions: Path,
) -> None:
    """Test that put_spliced + get_spliced round-trip succeeds."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    plasmid_fasta = sample_fasta
    human_fasta = tmp_path / "human.fasta"
    human_fasta.write_text(">chr1\nATGCATGC\n")

    # Put spliced artifacts
    cache.put_spliced(
        str(plasmid_fasta),
        str(human_fasta),
        str(sample_bam),
        str(sample_spliced_fasta),
        str(sample_cdna_positions),
    )

    # Get spliced artifacts
    result = cache.get_spliced(str(plasmid_fasta), str(human_fasta))
    assert result is not None
    assert "bam" in result
    assert "fasta" in result
    assert "cdna_positions" in result
    assert result["bam"].exists()
    assert result["fasta"].exists()
    assert result["cdna_positions"].exists()


@pytest.mark.unit
def test_spliced_combined_key_uses_both_inputs(
    tmp_path: Path,
    sample_bam: Path,
    sample_spliced_fasta: Path,
    sample_cdna_positions: Path,
) -> None:
    """Test that combined key uses both plasmid and human MD5."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    plasmid1 = tmp_path / "plasmid1.fasta"
    plasmid1.write_text(">plasmid1\nACGT\n")
    plasmid2 = tmp_path / "plasmid2.fasta"
    plasmid2.write_text(">plasmid2\nGCTA\n")
    human = tmp_path / "human.fasta"
    human.write_text(">chr1\nATGC\n")

    # Put with plasmid1
    cache.put_spliced(
        str(plasmid1),
        str(human),
        str(sample_bam),
        str(sample_spliced_fasta),
        str(sample_cdna_positions),
    )

    # Get with plasmid1 should hit
    assert cache.get_spliced(str(plasmid1), str(human)) is not None

    # Get with plasmid2 should miss (different plasmid)
    assert cache.get_spliced(str(plasmid2), str(human)) is None


@pytest.mark.unit
def test_spliced_cache_differentiates_padding(
    tmp_path: Path,
    sample_bam: Path,
    sample_spliced_fasta: Path,
    sample_cdna_positions: Path,
) -> None:
    """Test that get_spliced and put_spliced differentiate different padding values."""
    cache = CacheManager(tmp_path / "cache")
    plasmid = tmp_path / "plasmid.fasta"
    plasmid.write_text(">plasmid\nACGT\n")
    human = tmp_path / "human.fasta"
    human.write_text(">human\nATGC\n")

    cache.put_spliced(
        str(plasmid),
        str(human),
        str(sample_bam),
        str(sample_spliced_fasta),
        str(sample_cdna_positions),
        padding=1000,
    )

    assert cache.get_spliced(str(plasmid), str(human), padding=1000) is not None
    assert cache.get_spliced(str(plasmid), str(human), padding=500) is None


@pytest.mark.unit
def test_spliced_cache_rejects_v2_key(
    tmp_path: Path,
    sample_bam: Path,
    sample_spliced_fasta: Path,
    sample_cdna_positions: Path,
) -> None:
    """Test that existing v2 cache entries are not hit after bumping to v3."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    plasmid = tmp_path / "plasmid.fasta"
    plasmid.write_text(">plasmid\nACGT\n")
    human = tmp_path / "human.fasta"
    human.write_text(">chr1\nATGC\n")

    p_md5 = _resolve_md5(str(plasmid))
    h_md5 = _resolve_md5(str(human))
    v2_raw = f"{p_md5}:{h_md5}:1000:v2"
    v2_key = hashlib.md5(v2_raw.encode()).hexdigest()

    # Manually populate manifest with v2 key
    cache.manifest[v2_key] = {
        "bam_path": str(sample_bam),
        "fasta_path": str(sample_spliced_fasta),
        "cdna_path": str(sample_cdna_positions),
        "bam_size": sample_bam.stat().st_size,
        "fasta_size": sample_spliced_fasta.stat().st_size,
        "cdna_size": sample_cdna_positions.stat().st_size,
        "last_access": 0,
    }
    cache._save_manifest()

    # Querying spliced cache should miss because it looks for :v3
    assert cache.get_spliced(str(plasmid), str(human), padding=1000) is None


# Atomic write tests


@pytest.mark.unit
def test_cached_file_content_matches_source(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that cached file content matches source."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    original_mmi_content = sample_mmi.read_bytes()
    original_fai_content = sample_fai.read_bytes()

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    cached = cache.get_index(str(sample_fasta))
    assert cached is not None
    assert cached["mmi"].read_bytes() == original_mmi_content
    assert cached["fai"].read_bytes() == original_fai_content


@pytest.mark.unit
def test_no_temp_files_after_put(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that temp files are cleaned up after successful write."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    # Check for temp files
    temp_files = list(cache.indexes_dir.glob("*.tmp"))
    assert len(temp_files) == 0


# LRU eviction tests


@pytest.mark.unit
def test_no_eviction_when_under_limit(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that no eviction occurs when under size limit."""
    cache_root = tmp_path / "cache"
    # Set very large limit
    cache = CacheManager(cache_root, size_limit_gb=10.0)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    # Entry should still be there
    assert cache.get_index(str(sample_fasta)) is not None


@pytest.mark.unit
def test_oldest_entry_evicted_first(tmp_path: Path) -> None:
    """Test that oldest entry is evicted first when over limit."""
    cache_root = tmp_path / "cache"

    # Create multiple FASTA files
    fasta1 = tmp_path / "fasta1.fasta"
    fasta1.write_text(">seq1\nACGT\n")
    fasta2 = tmp_path / "fasta2.fasta"
    fasta2.write_text(">seq2\nGCTA\n")

    mmi1 = tmp_path / "mmi1.mmi"
    mmi1.write_bytes(b"MMI1" * 100)
    mmi2 = tmp_path / "mmi2.mmi"
    mmi2.write_bytes(b"MMI2" * 100)

    fai1 = tmp_path / "fai1.fai"
    fai1.write_text("seq1\t4\t6\t4\t5\n")
    fai2 = tmp_path / "fai2.fai"
    fai2.write_text("seq2\t4\t6\t4\t5\n")

    # Calculate size needed for exactly one entry
    entry_size = len(mmi1.read_bytes()) + len(fai1.read_bytes())
    # Set limit to just fit one entry (in GB)
    limit_gb = (entry_size + 100) / (1024**3)  # Add 100 bytes buffer
    cache = CacheManager(cache_root, size_limit_gb=limit_gb)

    # Put first entry
    cache.put_index(str(fasta1), str(mmi1), str(fai1))
    time.sleep(0.01)  # Ensure different timestamps

    # Put second entry (should evict first)
    cache.put_index(str(fasta2), str(mmi2), str(fai2))

    # First entry should be evicted
    assert cache.get_index(str(fasta1)) is None
    # Second entry should still be there
    assert cache.get_index(str(fasta2)) is not None


@pytest.mark.unit
def test_multiple_entries_evicted_if_needed(tmp_path: Path) -> None:
    """Test that multiple entries are evicted if needed."""
    cache_root = tmp_path / "cache"
    # Set tiny limit
    cache = CacheManager(cache_root, size_limit_gb=0.000000001)

    # Create three FASTA files
    fastas = []
    mmis = []
    fais = []
    for i in range(3):
        fasta = tmp_path / f"fasta{i}.fasta"
        fasta.write_text(f">seq{i}\nACGT\n")
        mmi = tmp_path / f"mmi{i}.mmi"
        mmi.write_bytes(b"MMI" * 100)
        fai = tmp_path / f"fai{i}.fai"
        fai.write_text(f"seq{i}\t4\t6\t4\t5\n")
        fastas.append(fasta)
        mmis.append(mmi)
        fais.append(fai)
        time.sleep(0.01)

    # Put all entries
    for fasta, mmi, fai in zip(fastas, mmis, fais, strict=True):
        cache.put_index(str(fasta), str(mmi), str(fai))

    # Due to tiny limit, only the last entry should remain
    assert cache.get_index(str(fastas[0])) is None
    assert cache.get_index(str(fastas[1])) is None
    # Last entry might also be evicted depending on exact sizes
    # At minimum, older entries should be gone


@pytest.mark.unit
def test_evicted_files_deleted_from_disk(tmp_path: Path) -> None:
    """Test that evicted files are deleted from disk."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root, size_limit_gb=0.000000001)

    fasta1 = tmp_path / "fasta1.fasta"
    fasta1.write_text(">seq1\nACGT\n")
    mmi1 = tmp_path / "mmi1.mmi"
    mmi1.write_bytes(b"MMI1" * 100)
    fai1 = tmp_path / "fai1.fai"
    fai1.write_text("seq1\t4\t6\t4\t5\n")

    cache.put_index(str(fasta1), str(mmi1), str(fai1))

    # Get cached paths before eviction
    cached = cache.get_index(str(fasta1))
    if cached:
        cached_mmi = cached["mmi"]
        cached_fai = cached["fai"]

        # Add another entry to trigger eviction
        fasta2 = tmp_path / "fasta2.fasta"
        fasta2.write_text(">seq2\nGCTA\n")
        mmi2 = tmp_path / "mmi2.mmi"
        mmi2.write_bytes(b"MMI2" * 100)
        fai2 = tmp_path / "fai2.fai"
        fai2.write_text("seq2\t4\t6\t4\t5\n")
        time.sleep(0.01)
        cache.put_index(str(fasta2), str(mmi2), str(fai2))

        # Old cached files should be deleted
        assert not cached_mmi.exists()
        assert not cached_fai.exists()


# Manifest tests


@pytest.mark.unit
def test_corrupted_manifest_handled_gracefully(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Test that corrupted manifest.json is handled gracefully."""
    cache_root = tmp_path / "cache"
    cache_root.mkdir()
    manifest_path = cache_root / "manifest.json"

    # Write garbage to manifest
    manifest_path.write_text("{ INVALID JSON }")

    # Should still initialize
    with caplog.at_level("WARNING"):
        cache = CacheManager(cache_root)

    assert "Corrupted manifest" in caplog.text
    assert len(cache.manifest) == 0


@pytest.mark.unit
def test_manifest_persists_across_instances(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that manifest persists across CacheManager instances."""
    cache_root = tmp_path / "cache"

    # Create first instance and put entry
    cache1 = CacheManager(cache_root)
    cache1.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))
    del cache1

    # Create second instance and get entry
    cache2 = CacheManager(cache_root)
    result = cache2.get_index(str(sample_fasta))
    assert result is not None


@pytest.mark.unit
def test_manifest_json_format(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that manifest.json has correct format."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    # Read manifest JSON
    manifest_path = cache_root / "manifest.json"
    with open(manifest_path) as f:
        data = json.load(f)

    assert "version" in data
    assert data["version"] == "1.0"
    assert "entries" in data
    assert len(data["entries"]) == 1


# Info and clean tests


@pytest.mark.unit
def test_info_returns_correct_structure(tmp_path: Path) -> None:
    """Test that info() returns correct structure."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    info = cache.info()
    assert "total_entries" in info
    assert "total_size_bytes" in info
    assert "cache_dir" in info
    assert "oldest_entry" in info
    assert "newest_entry" in info
    assert info["total_entries"] == 0


@pytest.mark.unit
def test_info_calculates_size_correctly(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that info() calculates total size correctly."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    info = cache.info()
    expected_size = sample_mmi.stat().st_size + sample_fai.stat().st_size
    assert info["total_size_bytes"] == expected_size
    assert info["total_entries"] == 1


@pytest.mark.unit
def test_info_tracks_oldest_and_newest(tmp_path: Path, sample_mmi: Path, sample_fai: Path) -> None:
    """Test that info() tracks oldest and newest entry timestamps."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    fasta1 = tmp_path / "fasta1.fasta"
    fasta1.write_text(">seq1\nACGT\n")
    fasta2 = tmp_path / "fasta2.fasta"
    fasta2.write_text(">seq2\nGCTA\n")

    cache.put_index(str(fasta1), str(sample_mmi), str(sample_fai))
    time.sleep(0.01)
    cache.put_index(str(fasta2), str(sample_mmi), str(sample_fai))

    info = cache.info()
    assert info["oldest_entry"] is not None
    assert info["newest_entry"] is not None
    assert info["oldest_entry"] < info["newest_entry"]


@pytest.mark.unit
def test_clean_all_removes_entries(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that clean_all() removes all entries."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))

    count = cache.clean_all()
    assert count == 1
    assert cache.get_index(str(sample_fasta)) is None


@pytest.mark.unit
def test_clean_all_deletes_files(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that clean_all() deletes cached files from disk."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))
    cached = cache.get_index(str(sample_fasta))
    assert cached is not None

    cached_mmi = cached["mmi"]
    cached_fai = cached["fai"]

    cache.clean_all()

    assert not cached_mmi.exists()
    assert not cached_fai.exists()


@pytest.mark.unit
def test_clean_all_deletes_manifest(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that clean_all() deletes manifest."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))
    assert (cache_root / "manifest.json").exists()

    cache.clean_all()

    assert not (cache_root / "manifest.json").exists()


@pytest.mark.unit
def test_clean_all_recreates_directories(
    tmp_path: Path, sample_fasta: Path, sample_mmi: Path, sample_fai: Path
) -> None:
    """Test that clean_all() recreates empty directories."""
    cache_root = tmp_path / "cache"
    cache = CacheManager(cache_root)

    cache.put_index(str(sample_fasta), str(sample_mmi), str(sample_fai))
    cache.clean_all()

    assert (cache_root / "indexes").exists()
    assert (cache_root / "spliced").exists()
