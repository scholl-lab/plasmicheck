# Phase 14: Cache Infrastructure - Research

**Researched:** 2026-02-15
**Domain:** Content-addressed file caching for bioinformatics pipeline artifacts (minimap2 indexes, spliced alignments)
**Confidence:** HIGH

## Summary

Phase 14 implements a content-addressed cache system to eliminate redundant computation in PlasmiCheck's pipeline. The cache stores minimap2 indexes and spliced alignment artifacts using MD5-based content addressing, supporting both within-run deduplication and cross-run persistence. The implementation follows established Python patterns: chunked MD5 hashing for large reference files, atomic file writes via tempfile + os.replace, LRU eviction when size limits are exceeded, and JSON manifest for metadata tracking.

**Primary recommendation:** Build a lightweight custom CacheManager class (not a third-party library like diskcache) that uses MD5 content-addressing, stores cache entries in `~/.plasmicheck/cache/{indexes,spliced}/`, tracks metadata in a JSON manifest, and integrates seamlessly with the existing pipeline orchestration in `run_pipeline.py`. The existing `built_indexes` set pattern already demonstrates in-memory deduplication; extend this to persistent disk cache with the same lookup semantics.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| hashlib | stdlib | MD5 computation | Already used in codebase (`utils.py:calculate_md5`); chunked reading pattern established |
| json | stdlib | Manifest storage | Simple, human-readable, sufficient for metadata tracking |
| os | stdlib | Atomic file operations | `os.replace()` provides cross-platform atomic rename |
| tempfile | stdlib | Atomic write pattern | `NamedTemporaryFile` in same directory ensures same-filesystem atomic replace |
| pathlib | stdlib | Path manipulation | Already used in `OutputPaths` dataclass (Phase 13) |
| collections.OrderedDict | stdlib | LRU tracking | O(1) access and eviction via `move_to_end()` method |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| time | stdlib | Last-access timestamps | Record cache hit times for LRU eviction |
| shutil | stdlib | Directory size calculation | Enforce cache size cap (10 GB default) |
| logging | stdlib | Cache hit/miss events | INFO-level events as specified in CONTEXT.md |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Custom CacheManager | diskcache library | diskcache adds external dependency, overkill for simple content-addressed cache; PlasmiCheck needs <200 LOC custom solution |
| JSON manifest | SQLite database | SQLite more complex, unnecessary overhead for ~10-100 cache entries; JSON is human-readable for debugging |
| MD5 hashing | SHA256 hashing | MD5 sufficient for content-addressing non-adversarial data; already used in codebase; faster for multi-GB references |

**Installation:**
No new dependencies required — all libraries are Python stdlib.

## Architecture Patterns

### Recommended Project Structure
```
plasmicheck/
├── cache_manager.py         # New: CacheManager class
├── scripts/
│   ├── run_pipeline.py      # Modified: integrate cache_manager
│   ├── create_indexes.py    # Modified: check cache before indexing
│   ├── spliced_alignment.py # Modified: check cache before alignment
│   └── utils.py             # Already has calculate_md5, no changes needed
└── cli.py                   # Modified: add `cache` subcommand

~/.plasmicheck/cache/        # User cache directory
├── manifest.json            # Cache metadata
├── indexes/
│   └── <md5_hash>.mmi       # Cached minimap2 indexes
└── spliced/
    ├── <combined_md5>.bam   # Cached spliced alignment BAM
    ├── <combined_md5>.fasta # Cached spliced reference FASTA
    └── <combined_md5>.txt   # Cached cDNA positions
```

### Pattern 1: Content-Addressed Cache Key Generation
**What:** Compute MD5 hash of input file contents to generate cache key
**When to use:** Before indexing (plasmid/human reference) or spliced alignment
**Example:**
```python
# Source: Existing pattern in utils.py (lines 21-27)
from plasmicheck.scripts.utils import calculate_md5

# Index cache: single file hash
plasmid_md5 = calculate_md5(plasmid_fasta)
cache_key = plasmid_md5  # e.g., "a3f2c8b1d9e4..."

# Spliced alignment cache: combined hash
plasmid_md5 = calculate_md5(plasmid_fasta)
human_md5 = calculate_md5(human_fasta)
cache_key = hashlib.md5(f"{plasmid_md5}:{human_md5}".encode()).hexdigest()
```

**Rationale:** MD5 is already used throughout PlasmiCheck for checksumming (`write_md5sum`). Consistency matters more than algorithm choice. The chunked reading pattern (4096-byte chunks) in `calculate_md5` handles multi-GB reference files efficiently without memory spikes.

### Pattern 2: Atomic Cache Write
**What:** Write to temp file in same directory, then atomic rename
**When to use:** When populating cache with newly created index or spliced alignment
**Example:**
```python
# Source: Python stdlib best practices
# https://python-atomicwrites.readthedocs.io/
# https://docs.python.org/3/library/tempfile.html
import os
import tempfile
from pathlib import Path

def cache_file_atomically(source_path: str, cache_path: Path) -> None:
    """Copy file to cache with atomic write."""
    cache_dir = cache_path.parent
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Create temp file in same directory (same filesystem)
    with tempfile.NamedTemporaryFile(
        dir=str(cache_dir),
        delete=False,
        suffix=".tmp"
    ) as temp:
        temp_path = temp.name
        # Copy source to temp
        with open(source_path, 'rb') as src:
            while chunk := src.read(4096 * 1024):  # 4 MB chunks
                temp.write(chunk)

    # Atomic replace
    os.replace(temp_path, str(cache_path))
```

**Rationale:** `os.replace()` is atomic on POSIX (Linux/macOS) and uses MoveFileEx on Windows (best-effort atomic). Creating temp file in same directory ensures same filesystem, which is required for atomic rename. This prevents corrupted cache entries if write is interrupted.

### Pattern 3: LRU Eviction with OrderedDict
**What:** Track cache entries in insertion order, evict oldest when size limit exceeded
**When to use:** After adding new cache entry, check total size and evict if needed
**Example:**
```python
# Source: LRU cache patterns with OrderedDict
# https://www.geeksforgeeks.org/python/lru-cache-in-python-using-ordereddict/
from collections import OrderedDict
from pathlib import Path

class CacheManifest:
    """Tracks cache entries with LRU eviction."""

    def __init__(self, manifest_path: Path, size_limit_gb: float = 10.0):
        self.manifest_path = manifest_path
        self.size_limit = size_limit_gb * 1024**3  # bytes
        self.entries: OrderedDict[str, dict] = OrderedDict()
        self._load()

    def access(self, cache_key: str) -> None:
        """Mark entry as recently used (move to end)."""
        if cache_key in self.entries:
            self.entries.move_to_end(cache_key)
            self.entries[cache_key]["last_access"] = time.time()

    def evict_lru(self, cache_dir: Path) -> None:
        """Remove least recently used entries until under size limit."""
        total_size = sum(e["size"] for e in self.entries.values())

        while total_size > self.size_limit and self.entries:
            # Pop from beginning (oldest)
            oldest_key, oldest_entry = self.entries.popitem(last=False)

            # Delete cache file(s)
            for file_path in oldest_entry["files"]:
                Path(file_path).unlink(missing_ok=True)

            total_size -= oldest_entry["size"]
            logging.warning(f"Cache evicted LRU entry: {oldest_key}")
```

**Rationale:** OrderedDict maintains insertion order and `move_to_end()` is O(1), making LRU tracking efficient. The manifest persists between runs, so LRU order survives process restarts. When a cache hit occurs, move that entry to the end (most recent).

### Pattern 4: Cache Integration in Pipeline
**What:** Check cache before expensive operations (indexing, spliced alignment)
**When to use:** In `create_indexes.py` and `spliced_alignment.py` before running minimap2
**Example:**
```python
# Source: Integration pattern for run_pipeline.py
from plasmicheck.cache_manager import CacheManager

def run_pipeline(..., cache_dir: str | None = None, no_cache: bool = False):
    # Initialize cache manager
    if no_cache:
        cache = None
    else:
        cache_path = Path(cache_dir) if cache_dir else Path.home() / ".plasmicheck" / "cache"
        cache = CacheManager(cache_path)

    # In indexing step
    plasmid_index = str(paths.plasmid_index)
    if cache and (cached_index := cache.get_index(plasmid_fasta)):
        logging.info(f"Cache HIT: reusing index for {plasmid_fasta}")
        shutil.copy(cached_index, plasmid_index)
    else:
        logging.info(f"Cache MISS: building index for {plasmid_fasta}")
        create_indexes(plasmid_fasta, overwrite)
        if cache:
            cache.put_index(plasmid_fasta, plasmid_index)
```

**Rationale:** Minimal changes to existing pipeline flow. Cache is checked before expensive operations; if cache miss, proceed with normal operation and populate cache. The `built_indexes` set pattern (lines 74, 507-511 of `run_pipeline.py`) already demonstrates in-memory deduplication; cache extends this to persistent storage.

### Anti-Patterns to Avoid
- **Don't validate cache entries on every read:** Trust the manifest (size check only). Re-hashing multi-GB files on every read defeats caching performance. Validate only on detection of corruption (size mismatch).
- **Don't use file locking for concurrent writes:** Multiple processes writing same cache key is rare (same input file). Atomic rename handles race condition: last writer wins, both writes succeed. File locking adds complexity and deadlock risk.
- **Don't cache output-specific files:** Cache only input-dependent artifacts (indexes, spliced alignments). Don't cache per-sample outputs (comparison TSVs, reports) which depend on sequencing data.
- **Don't hardcode cache paths:** Always allow override via `--cache-dir` flag. Users may need custom cache locations (network storage, different disk).

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Atomic file writes | Manual rename with error handling | `tempfile.NamedTemporaryFile` + `os.replace()` | Stdlib pattern handles same-filesystem requirement, Windows compatibility, cleanup on error |
| MD5 hashing large files | Load entire file into memory | Chunked reading with `hashlib.md5()` | Already implemented in `utils.py:calculate_md5`; handles multi-GB references efficiently |
| Directory traversal and size | Recursive os.walk with manual sum | `shutil.disk_usage()` for limits, manual walk for eviction | Stdlib provides disk usage; need manual walk to identify specific files for eviction |
| LRU eviction logic | Custom linked list implementation | `collections.OrderedDict` with `move_to_end()` | O(1) operations, maintained by stdlib, less error-prone than manual linked list |

**Key insight:** Python stdlib provides all necessary primitives for cache implementation. Third-party libraries like `diskcache` are overkill for PlasmiCheck's simple content-addressed cache (estimated <200 LOC for `CacheManager`). Custom implementation gives full control over cache key format (combined MD5 for spliced alignments), directory structure, and integration points.

## Common Pitfalls

### Pitfall 1: Cross-Filesystem Rename Failure
**What goes wrong:** `os.replace()` fails with OSError if temp file and target are on different filesystems
**Why it happens:** Default `tempfile.NamedTemporaryFile()` uses system temp directory (often `/tmp`), which may be on different filesystem than `~/.plasmicheck/cache/`
**How to avoid:** Always create temp file in same directory as target:
```python
# WRONG: may be different filesystem
temp = tempfile.NamedTemporaryFile(delete=False)

# CORRECT: guaranteed same filesystem
temp = tempfile.NamedTemporaryFile(dir=str(cache_dir), delete=False)
```
**Warning signs:** OSError with message "Invalid cross-device link" (Linux) or "The system cannot move the file to a different disk drive" (Windows)

### Pitfall 2: Cache Key Collision for Combined Hashes
**What goes wrong:** Two different input combinations produce same cache key
**Why it happens:** Naive concatenation without separator: `md5(plasmid_md5 + human_md5)` could collide if boundaries shift
**How to avoid:** Use delimiter in combined hash:
```python
# WRONG: collision risk
cache_key = hashlib.md5(f"{plasmid_md5}{human_md5}".encode()).hexdigest()

# CORRECT: delimiter prevents collision
cache_key = hashlib.md5(f"{plasmid_md5}:{human_md5}".encode()).hexdigest()
```
**Warning signs:** Cache hits when inputs changed; mysterious incorrect results from cached spliced alignments

### Pitfall 3: Manifest Corruption on Concurrent Writes
**What goes wrong:** Two pipeline runs update manifest.json simultaneously, corrupting JSON structure
**Why it happens:** Non-atomic read-modify-write of manifest file
**How to avoid:** Use atomic write pattern for manifest updates:
```python
# WRONG: non-atomic update
manifest = json.load(open("manifest.json"))
manifest[key] = value
json.dump(manifest, open("manifest.json", "w"))

# CORRECT: atomic update
with tempfile.NamedTemporaryFile(mode='w', dir=cache_dir, delete=False) as temp:
    json.dump(manifest, temp, indent=2)
os.replace(temp.name, manifest_path)
```
**Warning signs:** JSON decode errors on manifest load; corrupted manifest with partial writes

### Pitfall 4: Inaccessible Cache Directory Stops Pipeline
**What goes wrong:** Pipeline fails with exception when cache directory is unwritable (permissions, disk full)
**Why it happens:** Not catching OSError/PermissionError on cache operations
**How to avoid:** Graceful degradation — log warning and disable caching:
```python
try:
    cache = CacheManager(cache_path)
except (OSError, PermissionError) as e:
    logging.warning(f"Cache directory inaccessible: {e}. Disabling cache for this run.")
    cache = None
```
**Warning signs:** Pipeline crashes on cache initialization; permission denied errors in logs

### Pitfall 5: Stale Cache After minimap2 Version Upgrade
**What goes wrong:** Cached indexes built with minimap2 v2.24 used after upgrade to v2.28, causing alignment errors
**Why it happens:** Index format may change between minimap2 versions; manifest doesn't track tool versions
**How to avoid:** Include tool version in cache key OR manifest metadata:
```python
# Option 1: Version in cache key (recommended)
minimap2_version = subprocess.run(["minimap2", "--version"],
                                  capture_output=True, text=True).stdout.strip()
cache_key = hashlib.md5(f"{file_md5}:{minimap2_version}".encode()).hexdigest()

# Option 2: Version in manifest, validate on load
manifest_entry = {
    "cache_key": file_md5,
    "minimap2_version": minimap2_version,
    "created": time.time()
}
```
**Warning signs:** Alignment errors after minimap2 upgrade; minimap2 complaints about index format
**Decision:** Defer to planning phase — version tracking adds complexity. Consider scope for v0.34.0.

## Code Examples

Verified patterns from official sources:

### Chunked MD5 Hashing (Already Implemented)
```python
# Source: plasmicheck/scripts/utils.py lines 21-27
# Pattern verified: https://gist.github.com/aunyks/042c2798383f016939c40aa1be4f4aaf
import hashlib

def calculate_md5(file_path: str) -> str:
    """Calculate and return the MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
```

### Atomic File Write Pattern
```python
# Source: https://python-atomicwrites.readthedocs.io/
# https://code.activestate.com/recipes/579097-safely-and-atomically-write-to-a-file/
import os
import tempfile
from pathlib import Path

def atomic_write(content: bytes, target_path: Path) -> None:
    """Write content to file atomically."""
    target_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        dir=str(target_path.parent),
        delete=False
    ) as temp:
        temp.write(content)
        temp_name = temp.name

    os.replace(temp_name, str(target_path))
```

### Cache Manager Skeleton
```python
# Source: Design based on research findings
# LRU pattern: https://www.geeksforgeeks.org/python/lru-cache-in-python-using-ordereddict/
# Manifest pattern: dbt artifacts https://docs.getdbt.com/reference/artifacts/manifest-json

from __future__ import annotations
import json
import logging
import os
import shutil
import tempfile
import time
from collections import OrderedDict
from pathlib import Path
from plasmicheck.scripts.utils import calculate_md5

class CacheManager:
    """Content-addressed cache for pipeline artifacts."""

    def __init__(self, cache_root: Path, size_limit_gb: float = 10.0):
        self.cache_root = Path(cache_root)
        self.indexes_dir = self.cache_root / "indexes"
        self.spliced_dir = self.cache_root / "spliced"
        self.manifest_path = self.cache_root / "manifest.json"
        self.size_limit = size_limit_gb * 1024**3

        # Create cache directories
        self.indexes_dir.mkdir(parents=True, exist_ok=True)
        self.spliced_dir.mkdir(parents=True, exist_ok=True)

        # Load or initialize manifest
        self.manifest = self._load_manifest()

    def _load_manifest(self) -> OrderedDict:
        """Load manifest from disk or create new."""
        if self.manifest_path.exists():
            try:
                with open(self.manifest_path) as f:
                    data = json.load(f)
                    return OrderedDict(data.get("entries", {}))
            except (json.JSONDecodeError, OSError) as e:
                logging.warning(f"Corrupted manifest, creating new: {e}")
        return OrderedDict()

    def _save_manifest(self) -> None:
        """Atomically save manifest to disk."""
        manifest_data = {
            "version": "1.0",
            "entries": dict(self.manifest)
        }
        with tempfile.NamedTemporaryFile(
            mode='w',
            dir=str(self.cache_root),
            delete=False
        ) as temp:
            json.dump(manifest_data, temp, indent=2)
            temp_name = temp.name
        os.replace(temp_name, str(self.manifest_path))

    def get_index(self, input_fasta: str) -> Path | None:
        """Get cached index for input FASTA, or None if miss."""
        cache_key = calculate_md5(input_fasta)

        if cache_key in self.manifest:
            cached_path = Path(self.manifest[cache_key]["path"])
            if cached_path.exists():
                # Cache hit: update access time and move to end
                self.manifest.move_to_end(cache_key)
                self.manifest[cache_key]["last_access"] = time.time()
                self._save_manifest()
                logging.info(f"Cache HIT: {cache_key[:8]}... ({input_fasta})")
                return cached_path
            else:
                # Stale entry: file deleted outside cache manager
                del self.manifest[cache_key]
                self._save_manifest()

        logging.info(f"Cache MISS: {cache_key[:8]}... ({input_fasta})")
        return None

    def put_index(self, input_fasta: str, index_path: str) -> None:
        """Store index in cache."""
        cache_key = calculate_md5(input_fasta)
        cached_path = self.indexes_dir / f"{cache_key}.mmi"

        # Copy to cache atomically
        self._atomic_copy(index_path, cached_path)

        # Update manifest
        self.manifest[cache_key] = {
            "path": str(cached_path),
            "size": cached_path.stat().st_size,
            "created": time.time(),
            "last_access": time.time(),
            "type": "index"
        }
        self.manifest.move_to_end(cache_key)
        self._save_manifest()

        # Evict LRU if over size limit
        self._evict_if_needed()

    def _atomic_copy(self, source: str, target: Path) -> None:
        """Copy file to target atomically."""
        with tempfile.NamedTemporaryFile(
            dir=str(target.parent),
            delete=False
        ) as temp:
            with open(source, 'rb') as src:
                shutil.copyfileobj(src, temp)
            temp_name = temp.name
        os.replace(temp_name, str(target))

    def _evict_if_needed(self) -> None:
        """Evict LRU entries if total size exceeds limit."""
        total_size = sum(e["size"] for e in self.manifest.values())

        while total_size > self.size_limit and self.manifest:
            # Pop oldest (first) entry
            oldest_key, oldest_entry = self.manifest.popitem(last=False)
            cached_path = Path(oldest_entry["path"])

            if cached_path.exists():
                cached_path.unlink()
                total_size -= oldest_entry["size"]
                logging.warning(f"Cache evicted LRU entry: {oldest_key[:8]}...")

        self._save_manifest()
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Rebuild indexes every run | Cache and reuse indexes | v0.34.0 (this phase) | Dramatic runtime reduction for batch runs |
| In-memory dedup only (`built_indexes` set) | Persistent cache across runs | v0.34.0 (this phase) | Cross-run speedup, not just within-run |
| MD5 checksumming for verification only | MD5 as content-address cache key | v0.34.0 (this phase) | Leverage existing hashing for cache lookup |
| Manual cache management (if any) | Automatic LRU eviction | v0.34.0 (this phase) | Prevents unbounded cache growth |

**Deprecated/outdated:**
- None — this is a new feature, not replacing existing functionality

**Current best practices:**
- Python 3.11+ has `hashlib.file_digest()` for efficient file hashing, but PlasmiCheck targets Python 3.9+ (based on `from __future__ import annotations` usage). Stick with existing `calculate_md5` pattern.
- `os.replace()` preferred over `os.rename()` since Python 3.3+ for cross-platform atomic operations

## Open Questions

Things that couldn't be fully resolved:

1. **Tool version tracking in cache keys**
   - What we know: minimap2 index format may change between versions; using old index with new version could cause errors
   - What's unclear: How frequently does minimap2 index format change? Is version-aware caching worth the complexity in v0.34.0?
   - Recommendation: Start without version tracking. If index format incompatibility reported by users, add version tracking in future release. Include warning in cache info output about potential version mismatches.

2. **Cache warming strategy**
   - What we know: First run after cache clear will be slow (cache population)
   - What's unclear: Should PlasmiCheck provide `plasmicheck cache warm` command to pre-populate cache?
   - Recommendation: Defer to future release. Cache warming requires predicting which references will be used. Auto-population on first use is simpler.

3. **Concurrent pipeline runs sharing cache**
   - What we know: Atomic writes handle concurrent cache population (last writer wins, both succeed)
   - What's unclear: Race condition where Process A checks cache (miss), Process B populates cache, Process A redundantly rebuilds
   - Recommendation: Accept redundant work on rare race condition. File locking adds complexity and deadlock risk. Both processes succeed, one does wasted work. Log both cache puts, second will overwrite first atomically.

4. **XDG Base Directory Specification compliance**
   - What we know: Linux XDG standard says cache should use `$XDG_CACHE_HOME` (defaults to `~/.cache`)
   - What's unclear: Should default be `~/.cache/plasmicheck/` (XDG-compliant) or `~/.plasmicheck/cache/` (app-specific)?
   - Recommendation: Use `~/.plasmicheck/cache/` for consistency with potential future config files (`~/.plasmicheck/config.toml`). Users can override with `--cache-dir` if XDG compliance required. Document in help text.

## Sources

### Primary (HIGH confidence)
- Python stdlib documentation:
  - hashlib: https://docs.python.org/3/library/hashlib.html
  - tempfile: https://docs.python.org/3/library/tempfile.html
  - os.replace(): https://zetcode.com/python/os-replace/
  - collections.OrderedDict: https://docs.python.org/3/library/collections.html#collections.OrderedDict

- PlasmiCheck codebase (verified by reading):
  - `plasmicheck/scripts/utils.py` — `calculate_md5()` implementation (lines 21-27)
  - `plasmicheck/scripts/run_pipeline.py` — `built_indexes` pattern (lines 74, 507-511)
  - `plasmicheck/output_paths.py` — OutputPaths dataclass (Phase 13)

### Secondary (MEDIUM confidence)
- Python atomic writes patterns:
  - [Atomicwrites documentation](https://python-atomicwrites.readthedocs.io/)
  - [Safely and atomically write to a file](https://code.activestate.com/recipes/579097-safely-and-atomically-write-to-a-file/)

- LRU cache implementation patterns:
  - [LRU Cache in Python using OrderedDict - GeeksforGeeks](https://www.geeksforgeeks.org/python/lru-cache-in-python-using-ordereddict/)
  - [Implementing a Least Recently Used (LRU) Cache in Python - llego.dev](https://llego.dev/posts/implement-lru-cache-python/)

- File hashing best practices:
  - [Hash a large file in Python - GitHub Gist](https://gist.github.com/aunyks/042c2798383f016939c40aa1be4f4aaf)
  - [Python MD5 Hash: A Comprehensive Guide - CodeRivers](https://coderivers.org/blog/python-md5-hash/)

- Bioinformatics pipeline caching:
  - [minimap2 GitHub repository](https://github.com/lh3/minimap2) — index reuse patterns
  - [pbmm2 documentation](https://github.com/PacificBiosciences/pbmm2) — index pre-computation strategy

### Tertiary (LOW confidence)
- Diskcache library (not recommended for this use case):
  - [DiskCache Tutorial](https://grantjenks.com/docs/diskcache/tutorial.html)
  - [diskcache PyPI](https://pypi.org/project/diskcache/)
  - Reason for LOW confidence: Not applicable to PlasmiCheck's simple content-addressed needs; custom implementation preferred

- XDG Base Directory:
  - [XDG Base Directory Specification](https://xdgbasedirectoryspecification.com/)
  - [PyXDG documentation](https://pyxdg.readthedocs.io/en/latest/basedirectory.html)
  - Reason for LOW confidence: Compliance optional for PlasmiCheck; `--cache-dir` override sufficient

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all stdlib, no version dependencies, patterns verified in codebase
- Architecture: HIGH — atomic writes, LRU eviction, content-addressing are well-established patterns
- Pitfalls: HIGH — based on verified stdlib behavior and common Python file I/O issues

**Research date:** 2026-02-15
**Valid until:** 90 days (2026-05-16) — cache implementation patterns are stable; minimap2 index format changes infrequent
