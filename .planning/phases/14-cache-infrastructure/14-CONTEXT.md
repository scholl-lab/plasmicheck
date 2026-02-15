# Phase 14: Cache Infrastructure - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Content-addressed caching for minimap2 indexes and spliced alignments, enabling reuse within and across pipeline runs. Covers cache storage, within-run deduplication, cross-run persistence, CLI management commands, and graceful degradation. Does not include caching of alignment results or report outputs.

</domain>

<decisions>
## Implementation Decisions

### Cache key design
- Content-addressed using MD5 hash of file contents (consistent with existing `write_md5sum` usage in codebase)
- Index cache key: MD5 hash of the source plasmid/reference file
- Spliced alignment cache key: combined MD5 hashes of both plasmid file and human reference file
- Within-run deduplication uses the same cache mechanism — first build populates the cache, subsequent samples find the cache hit. No separate in-memory tracking needed.

### Storage & eviction policy
- Default cache directory: `~/.plasmicheck/cache/`
- Subdirectories: `indexes/` and `spliced/` (at minimum)
- JSON manifest file tracks entries, sizes, and last-access timestamps
- No automatic TTL-based eviction
- Size cap: 10 GB default, configurable
- Eviction strategy: LRU (least recently used) when size limit is exceeded

### CLI surface & user control
- Caching enabled by default; disable with `--no-cache` flag
- Override cache directory with `--cache-dir <path>` flag
- Cache hit/miss events logged at INFO level (always visible in normal output)
- `plasmicheck cache info` — summary only: total entries, total size, cache directory path, oldest/newest entry
- `plasmicheck cache clean --all` — clears entire cache. No selective filters in v0.34.0.

### Failure & edge cases
- Corrupted entries (incomplete file, wrong size): delete and rebuild silently, log a warning
- Concurrent writes: write to temp file, atomic rename into cache. No file locking.
- Inaccessible cache directory (permissions, disk full): log warning, disable caching for this run, pipeline continues normally
- Cache validation on read: trust the manifest (size check only). No re-hashing on every read.

### Claude's Discretion
- Internal cache directory structure beyond `indexes/` and `spliced/`
- Manifest file format details (exact JSON schema)
- Hash computation strategy (full file vs chunked reads for large references)
- Exact temp file naming convention for atomic writes
- Config.json integration for cache settings

</decisions>

<specifics>
## Specific Ideas

- MD5 is already used in the codebase (`write_md5sum` in utils.py) — stay consistent
- Cache key for spliced alignments must include both plasmid AND human reference hashes, since the spliced alignment depends on both inputs
- Within-run dedup is "free" — just a consequence of using the persistent cache. No special mechanism.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 14-cache-infrastructure*
*Context gathered: 2026-02-15*
