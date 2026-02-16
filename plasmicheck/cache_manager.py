"""Content-addressed cache for pipeline artifacts.

This module provides CacheManager for caching minimap2 indexes and spliced alignment
artifacts. Uses MD5-based content addressing, JSON manifest for metadata tracking,
atomic file writes, and LRU eviction when size limits are exceeded.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import tempfile
import time
from collections import OrderedDict
from pathlib import Path
from typing import Any

from plasmicheck.scripts.utils import calculate_md5

logger = logging.getLogger(__name__)


class CacheManager:
    """Content-addressed cache for pipeline artifacts.

    Stores minimap2 indexes (.mmi + .fai) and spliced alignment artifacts (BAM, FASTA,
    cDNA positions) with LRU eviction when size limit exceeded.

    Args:
        cache_root: Root directory for cache storage
        size_limit_gb: Maximum cache size in gigabytes (default: 10.0)

    Raises:
        OSError: If cache directory cannot be created
        PermissionError: If cache directory is not writable
    """

    def __init__(self, cache_root: Path, size_limit_gb: float = 10.0) -> None:
        """Initialize cache manager and create directory structure."""
        self.cache_root = Path(cache_root)
        self.indexes_dir = self.cache_root / "indexes"
        self.spliced_dir = self.cache_root / "spliced"
        self.manifest_path = self.cache_root / "manifest.json"
        self.size_limit = size_limit_gb * 1024**3  # Convert GB to bytes

        try:
            # Create cache directories
            self.indexes_dir.mkdir(parents=True, exist_ok=True)
            self.spliced_dir.mkdir(parents=True, exist_ok=True)
        except (OSError, PermissionError) as e:
            raise OSError(f"Cannot create cache directory {self.cache_root}: {e}") from e

        # Load or initialize manifest
        self.manifest: OrderedDict[str, dict[str, Any]] = self._load_manifest()

    def _load_manifest(self) -> OrderedDict[str, dict[str, Any]]:
        """Load manifest from disk or create new.

        Returns:
            OrderedDict of cache entries, empty if manifest doesn't exist or is corrupted
        """
        if self.manifest_path.exists():
            try:
                with open(self.manifest_path) as f:
                    data = json.load(f)
                    return OrderedDict(data.get("entries", {}))
            except (json.JSONDecodeError, OSError) as e:
                logger.warning(f"Corrupted manifest at {self.manifest_path}, creating new: {e}")
        return OrderedDict()

    def _save_manifest(self) -> None:
        """Atomically save manifest to disk.

        Uses tempfile + os.replace for atomic write to prevent corruption.
        """
        manifest_data = {"version": "1.0", "entries": dict(self.manifest)}

        # Create temp file in same directory for atomic replace
        with tempfile.NamedTemporaryFile(
            mode="w", dir=str(self.cache_root), delete=False, suffix=".tmp"
        ) as temp:
            json.dump(manifest_data, temp, indent=2)
            temp_name = temp.name

        os.replace(temp_name, str(self.manifest_path))

    def get_index(self, input_fasta: str) -> dict[str, Path] | None:
        """Get cached index for input FASTA.

        Args:
            input_fasta: Path to input FASTA file

        Returns:
            Dict with "mmi" and "fai" keys pointing to cached files, or None on miss
        """
        cache_key = calculate_md5(input_fasta)

        if cache_key in self.manifest:
            entry = self.manifest[cache_key]
            mmi_path = Path(entry["mmi_path"])
            fai_path = Path(entry["fai_path"])

            # Verify both files exist and sizes match
            if (
                mmi_path.exists()
                and fai_path.exists()
                and mmi_path.stat().st_size == entry["mmi_size"]
                and fai_path.stat().st_size == entry["fai_size"]
            ):
                # Cache hit: update access time and move to end (LRU)
                self.manifest.move_to_end(cache_key)
                self.manifest[cache_key]["last_access"] = time.time()
                self._save_manifest()
                logger.info(f"Cache HIT for index: {input_fasta}")
                return {"mmi": mmi_path, "fai": fai_path}
            else:
                # Stale entry: files deleted or size mismatch
                logger.warning(
                    f"Stale cache entry for {input_fasta} (files missing or size mismatch)"
                )
                del self.manifest[cache_key]
                self._save_manifest()

        logger.info(f"Cache MISS for index: {input_fasta}")
        return None

    def put_index(self, input_fasta: str, mmi_path: str, fai_path: str) -> None:
        """Store index in cache.

        Args:
            input_fasta: Path to input FASTA file (used for cache key)
            mmi_path: Path to minimap2 index file (.mmi)
            fai_path: Path to samtools index file (.fai)
        """
        cache_key = calculate_md5(input_fasta)
        cached_mmi = self.indexes_dir / f"{cache_key}.mmi"
        cached_fai = self.indexes_dir / f"{cache_key}.fai"

        # Copy both files to cache atomically
        self._atomic_copy(mmi_path, cached_mmi)
        self._atomic_copy(fai_path, cached_fai)

        # Update manifest with both file sizes
        mmi_size = cached_mmi.stat().st_size
        fai_size = cached_fai.stat().st_size
        self.manifest[cache_key] = {
            "mmi_path": str(cached_mmi),
            "fai_path": str(cached_fai),
            "mmi_size": mmi_size,
            "fai_size": fai_size,
            "created": time.time(),
            "last_access": time.time(),
            "type": "index",
        }
        self.manifest.move_to_end(cache_key)
        self._save_manifest()

        # Evict LRU if over size limit
        self._evict_if_needed()

        logger.info(f"Cached index for {input_fasta} ({cache_key[:8]}...)")

    def get_spliced(self, plasmid_fasta: str, human_fasta: str) -> dict[str, Path] | None:
        """Get cached spliced alignment artifacts.

        Args:
            plasmid_fasta: Path to plasmid FASTA file
            human_fasta: Path to human reference FASTA file

        Returns:
            Dict with "bam", "fasta", "cdna_positions" keys pointing to cached files,
            or None on miss
        """
        # Compute combined cache key
        plasmid_md5 = calculate_md5(plasmid_fasta)
        human_md5 = calculate_md5(human_fasta)
        cache_key = hashlib.md5(f"{plasmid_md5}:{human_md5}".encode()).hexdigest()

        if cache_key in self.manifest:
            entry = self.manifest[cache_key]
            bam_path = Path(entry["bam_path"])
            fasta_path = Path(entry["fasta_path"])
            cdna_path = Path(entry["cdna_path"])

            # Verify all files exist and sizes match
            if (
                bam_path.exists()
                and fasta_path.exists()
                and cdna_path.exists()
                and bam_path.stat().st_size == entry["bam_size"]
                and fasta_path.stat().st_size == entry["fasta_size"]
                and cdna_path.stat().st_size == entry["cdna_size"]
            ):
                # Cache hit: update access time and move to end (LRU)
                self.manifest.move_to_end(cache_key)
                self.manifest[cache_key]["last_access"] = time.time()
                self._save_manifest()
                logger.info(f"Cache HIT for spliced alignment: {plasmid_fasta} + {human_fasta}")
                return {"bam": bam_path, "fasta": fasta_path, "cdna_positions": cdna_path}
            else:
                # Stale entry: files deleted or size mismatch
                logger.warning(
                    "Stale cache entry for spliced alignment (files missing or size mismatch)"
                )
                del self.manifest[cache_key]
                self._save_manifest()

        logger.info(f"Cache MISS for spliced alignment: {plasmid_fasta} + {human_fasta}")
        return None

    def put_spliced(
        self,
        plasmid_fasta: str,
        human_fasta: str,
        bam_path: str,
        fasta_path: str,
        cdna_path: str,
    ) -> None:
        """Store spliced alignment artifacts in cache.

        Args:
            plasmid_fasta: Path to plasmid FASTA file (for cache key)
            human_fasta: Path to human reference FASTA file (for cache key)
            bam_path: Path to spliced alignment BAM file
            fasta_path: Path to spliced reference FASTA file
            cdna_path: Path to cDNA positions file
        """
        # Compute combined cache key
        plasmid_md5 = calculate_md5(plasmid_fasta)
        human_md5 = calculate_md5(human_fasta)
        cache_key = hashlib.md5(f"{plasmid_md5}:{human_md5}".encode()).hexdigest()

        cached_bam = self.spliced_dir / f"{cache_key}.bam"
        cached_fasta = self.spliced_dir / f"{cache_key}.fasta"
        cached_cdna = self.spliced_dir / f"{cache_key}.txt"

        # Copy all three files to cache atomically
        self._atomic_copy(bam_path, cached_bam)
        self._atomic_copy(fasta_path, cached_fasta)
        self._atomic_copy(cdna_path, cached_cdna)

        # Update manifest with all file sizes
        bam_size = cached_bam.stat().st_size
        fasta_size = cached_fasta.stat().st_size
        cdna_size = cached_cdna.stat().st_size

        self.manifest[cache_key] = {
            "bam_path": str(cached_bam),
            "fasta_path": str(cached_fasta),
            "cdna_path": str(cached_cdna),
            "bam_size": bam_size,
            "fasta_size": fasta_size,
            "cdna_size": cdna_size,
            "created": time.time(),
            "last_access": time.time(),
            "type": "spliced",
        }
        self.manifest.move_to_end(cache_key)
        self._save_manifest()

        # Evict LRU if over size limit
        self._evict_if_needed()

        logger.info(
            f"Cached spliced alignment for {plasmid_fasta} + {human_fasta} ({cache_key[:8]}...)"
        )

    def _atomic_copy(self, source: str, target: Path) -> None:
        """Copy file to target atomically.

        Args:
            source: Source file path
            target: Target file path

        Uses tempfile in same directory + os.replace for atomic write.
        """
        try:
            # Create temp file in same directory (same filesystem)
            with tempfile.NamedTemporaryFile(
                dir=str(target.parent), delete=False, suffix=".tmp"
            ) as temp:
                with open(source, "rb") as src:
                    shutil.copyfileobj(src, temp)
                temp_name = temp.name

            # Atomic replace
            os.replace(temp_name, str(target))
        except Exception:
            # Clean up temp file if copy failed
            if "temp_name" in locals() and os.path.exists(temp_name):
                os.unlink(temp_name)
            raise

    def _evict_if_needed(self) -> None:
        """Evict LRU entries if total size exceeds limit.

        Removes oldest entries first until total size is under limit.
        """
        # Calculate total size of all cached files
        total_size = 0
        for entry in self.manifest.values():
            if entry["type"] == "index":
                total_size += entry["mmi_size"] + entry["fai_size"]
            elif entry["type"] == "spliced":
                total_size += entry["bam_size"] + entry["fasta_size"] + entry["cdna_size"]

        # Evict oldest entries until under limit
        while total_size > self.size_limit and self.manifest:
            # Pop oldest (first) entry
            oldest_key, oldest_entry = self.manifest.popitem(last=False)

            # Delete cached files
            if oldest_entry["type"] == "index":
                mmi_path = Path(oldest_entry["mmi_path"])
                fai_path = Path(oldest_entry["fai_path"])
                mmi_path.unlink(missing_ok=True)
                fai_path.unlink(missing_ok=True)
                evicted_size = oldest_entry["mmi_size"] + oldest_entry["fai_size"]
            elif oldest_entry["type"] == "spliced":
                bam_path = Path(oldest_entry["bam_path"])
                fasta_path = Path(oldest_entry["fasta_path"])
                cdna_path = Path(oldest_entry["cdna_path"])
                bam_path.unlink(missing_ok=True)
                fasta_path.unlink(missing_ok=True)
                cdna_path.unlink(missing_ok=True)
                evicted_size = (
                    oldest_entry["bam_size"]
                    + oldest_entry["fasta_size"]
                    + oldest_entry["cdna_size"]
                )
            else:
                evicted_size = 0

            total_size -= evicted_size
            logger.warning(f"Cache evicted LRU entry: {oldest_key[:8]}...")

        # Save manifest after eviction
        if total_size <= self.size_limit:
            self._save_manifest()

    def info(self) -> dict[str, Any]:
        """Get cache information.

        Returns:
            Dict with total_entries, total_size_bytes, cache_dir, oldest_entry, newest_entry
        """
        total_size = 0
        oldest_time = None
        newest_time = None

        for entry in self.manifest.values():
            # Calculate total size
            if entry["type"] == "index":
                total_size += entry["mmi_size"] + entry["fai_size"]
            elif entry["type"] == "spliced":
                total_size += entry["bam_size"] + entry["fasta_size"] + entry["cdna_size"]

            # Track oldest and newest
            created = entry["created"]
            if oldest_time is None or created < oldest_time:
                oldest_time = created
            if newest_time is None or created > newest_time:
                newest_time = created

        return {
            "total_entries": len(self.manifest),
            "total_size_bytes": total_size,
            "cache_dir": str(self.cache_root),
            "oldest_entry": oldest_time,
            "newest_entry": newest_time,
        }

    def clean_all(self) -> int:
        """Delete all cached files and manifest.

        Returns:
            Number of entries removed
        """
        entry_count = len(self.manifest)

        # Delete all cached files
        for entry in self.manifest.values():
            if entry["type"] == "index":
                Path(entry["mmi_path"]).unlink(missing_ok=True)
                Path(entry["fai_path"]).unlink(missing_ok=True)
            elif entry["type"] == "spliced":
                Path(entry["bam_path"]).unlink(missing_ok=True)
                Path(entry["fasta_path"]).unlink(missing_ok=True)
                Path(entry["cdna_path"]).unlink(missing_ok=True)

        # Delete manifest
        self.manifest_path.unlink(missing_ok=True)

        # Clear in-memory manifest
        self.manifest.clear()

        # Recreate empty directories
        self.indexes_dir.mkdir(parents=True, exist_ok=True)
        self.spliced_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Cleaned all cache entries ({entry_count} entries removed)")
        return entry_count
