# Phase 13: Output Restructuring - Research

**Researched:** 2026-02-15
**Domain:** Python path management, dataclasses, file organization patterns
**Confidence:** HIGH

## Summary

This phase restructures pipeline outputs from a flat directory structure into logical subfolders (alignments/, references/, data/, reports/, plots/) while maintaining backward compatibility for summary report readers. The standard approach uses frozen dataclasses with pathlib.Path for type-safe, centralized path construction, and depth-agnostic file discovery via os.walk() for robust summary report reading.

Current state: PlasmiCheck creates outputs in `output/<sample>/<plasmid>/` with all files flat in the same directory. Path construction is scattered across pipeline scripts using string concatenation with `os.path.join()`. The summary report reader (`generate_summary_reports.py`) uses `os.walk()` to find `.reads_assignment.tsv` and `.summary.tsv` files, extracting sample/plasmid names from parent directory structure (lines 74-75, 90-91).

**Primary recommendation:** Use frozen dataclass with pathlib.Path properties for all path construction; maintain os.walk() for summary reader but make it robust to arbitrary nesting depth; introduce OutputPaths as single source of truth for all output file locations.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| dataclasses | stdlib (3.10+) | Type-safe configuration objects | Built-in, zero dependencies, excellent mypy support |
| pathlib | stdlib (3.10+) | Object-oriented path handling | Modern Python standard, cross-platform, cleaner than os.path |
| os.walk | stdlib | Recursive directory traversal | Standard library, depth-agnostic, lazy iteration |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| typing | stdlib | Type hints for Path | Already used in PlasmiCheck for strict mypy |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| pathlib | os.path strings | pathlib 3-5x slower in microbenchmarks but negligible in I/O-bound code; pathlib preferred for readability |
| os.walk | Path.rglob() | rglob() can be 40-100x slower on deeply nested directories (1000+ levels); os.walk safer for arbitrary structures |
| frozen dataclass | NamedTuple | Dataclass more extensible, supports methods and properties; NamedTuple only for lightweight immutable records |

**Installation:**
```bash
# No installation needed - all stdlib
```

## Architecture Patterns

### Recommended Project Structure
```
output/
├── <sample_name>/
│   └── <plasmid_name>/
│       ├── alignments/           # BAM files
│       │   ├── plasmid_alignment.bam
│       │   ├── spliced_human_alignment.bam
│       │   └── spliced_alignment.bam
│       ├── references/           # Reference sequences and indexes
│       │   ├── <plasmid>.fasta
│       │   ├── <plasmid>.mmi
│       │   ├── spliced_reference.fasta
│       │   └── spliced_reference.mmi
│       ├── data/                 # TSV outputs
│       │   ├── comparison_result.reads_assignment.tsv
│       │   ├── comparison_result.summary.tsv
│       │   └── cDNA_positions.txt
│       ├── reports/              # HTML reports
│       │   ├── report_interactive.html
│       │   └── report_non_interactive.html
│       └── plots/                # Plot files (if static)
│           └── *.png
```

### Pattern 1: Frozen Dataclass for Path Management
**What:** Centralized, immutable path configuration using dataclass with computed properties
**When to use:** When multiple modules need consistent path construction
**Example:**
```python
# Source: Dataclasses best practices (frozen, slots, kw_only)
from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True, slots=True)
class OutputPaths:
    """Centralized path construction for pipeline outputs.

    Args:
        output_root: Base output directory
        sample_name: Sample identifier
        plasmid_name: Plasmid identifier
    """
    output_root: Path
    sample_name: str
    plasmid_name: str

    @property
    def base(self) -> Path:
        """Sample/plasmid output directory."""
        return self.output_root / self.sample_name / self.plasmid_name

    @property
    def alignments(self) -> Path:
        """Alignment files directory."""
        return self.base / "alignments"

    @property
    def plasmid_bam(self) -> Path:
        """Plasmid alignment BAM file."""
        return self.alignments / "plasmid_alignment.bam"

    def ensure_dirs(self) -> None:
        """Create all output subdirectories."""
        for subdir in [self.alignments, self.references, self.data,
                       self.reports, self.plots]:
            subdir.mkdir(parents=True, exist_ok=True)
```

**Why frozen:** Makes dataclass immutable and hashable; prevents accidental modification during pipeline execution. Approximately 2.4x slower to instantiate but negligible in I/O-bound pipeline context.

**Why slots:** Reduces memory usage and speeds up attribute access (Python 3.10+).

**Why properties:** Computed paths ensure consistency; changing directory structure only requires updating properties, not all call sites.

### Pattern 2: Depth-Agnostic File Discovery
**What:** Robust file finding that works regardless of directory nesting depth
**When to use:** When reading files from user-provided or variable-depth directories
**Example:**
```python
# Source: os.walk() best practices for large trees
def find_tsv_files(input_dir: str, pattern: str) -> list[str]:
    """Recursively find all tsv files matching the pattern."""
    tsv_files: list[str] = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(pattern):
                tsv_files.append(os.path.join(root, file))
    return tsv_files
```

**Key insight:** This pattern already exists in `generate_summary_reports.py:49-56` and is correct. The parent directory parsing (lines 74-75, 90-91) assumes exactly 2 levels above the file, which is the constraint to relax.

### Pattern 3: Parent Directory Name Extraction (Depth-Agnostic)
**What:** Extract sample/plasmid names from file path without assuming depth
**When to use:** After directory restructuring when files may be at different depths
**Example:**
```python
# Source: Path manipulation patterns with pathlib
from pathlib import Path

def extract_sample_plasmid(file_path: str | Path) -> tuple[str, str]:
    """Extract sample and plasmid names from comparison file path.

    Expected path: .../sample_name/plasmid_name/.../file.tsv
    Works with any subdirectory depth under plasmid_name/.
    """
    path = Path(file_path)
    # Find the directory containing the file
    # Could be data/ subfolder or flat
    current = path.parent

    # Walk up to find structure matching sample/plasmid pattern
    # Look for pattern: parent.parent has multiple children (samples)
    #                  parent has this plasmid
    plasmid_name = current.name if current.name != "data" else current.parent.name
    # Continue walking up if we were in data/ subfolder
    sample_dir = current.parent if current.name != "data" else current.parent.parent
    sample_name = sample_dir.name

    return sample_name, plasmid_name
```

**Better approach:** Store metadata in the TSV files themselves or use consistent naming pattern.

### Anti-Patterns to Avoid
- **Hardcoded path strings:** Scattered `os.path.join(output, "plasmid_alignment.bam")` across modules makes refactoring error-prone
- **Depth assumptions:** `os.path.dirname(os.path.dirname(file))` breaks when directory structure changes
- **String paths in new code:** Use pathlib.Path for new code; only convert to str at system boundaries (subprocess calls)
- **Mutable path containers:** Path configuration objects should be immutable to prevent accidental modification during pipeline

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Path joining | String concatenation | pathlib.Path / operator | Cross-platform, handles edge cases (trailing slashes, absolute paths) |
| Directory creation | if not exists + mkdir | Path.mkdir(parents=True, exist_ok=True) | Idempotent, race-condition safe, handles parent creation |
| Recursive file finding | Manual recursion | os.walk() or Path.rglob() | Handles symlinks, permissions errors, filters efficiently |
| Path normalization | Manual string cleanup | Path.resolve() | Handles .., ., symlinks, relative-to-absolute conversion |
| File existence checks | try/except or os.path.exists | Path.exists() or Path.is_file() | More readable, distinguishes files vs directories |

**Key insight:** Python's stdlib path handling is battle-tested for edge cases (Windows paths, symlinks, permissions, Unicode, race conditions). Custom implementations miss these cases.

## Common Pitfalls

### Pitfall 1: Frozen Dataclass __post_init__ Assignment
**What goes wrong:** Trying to set computed fields in `__post_init__` on frozen dataclass raises FrozenInstanceError
**Why it happens:** frozen=True generates __setattr__ that prevents all attribute assignment after __init__
**How to avoid:** Use properties for computed values, or call `object.__setattr__(self, 'field', value)` in __post_init__
**Warning signs:** FrozenInstanceError during dataclass instantiation

**Example:**
```python
# WRONG
@dataclass(frozen=True)
class Paths:
    root: Path
    def __post_init__(self):
        self.subdir = self.root / "sub"  # FrozenInstanceError!

# RIGHT - use property
@dataclass(frozen=True)
class Paths:
    root: Path
    @property
    def subdir(self) -> Path:
        return self.root / "sub"

# ALSO RIGHT - use object.__setattr__ if truly needed
@dataclass(frozen=True)
class Paths:
    root: Path
    subdir: Path = field(init=False)
    def __post_init__(self):
        object.__setattr__(self, 'subdir', self.root / "sub")
```

### Pitfall 2: pathlib.Path.rglob() Performance on Deep Nesting
**What goes wrong:** Path.rglob() becomes 40-100x slower than glob.glob() on deeply nested directories (1000+ levels)
**Why it happens:** rglob() implementation has O(n²) behavior with depth due to internal recursion
**How to avoid:** Use os.walk() for arbitrary-depth recursive search; reserve rglob() for shallow (<100 levels) known structures
**Warning signs:** Summary report generation taking minutes on normal directory structures

**Example:**
```python
# SLOW on deep trees
def find_files_slow(root: Path) -> list[Path]:
    return list(root.rglob("*.tsv"))  # Can take 90+ seconds on 1000-level depth

# FAST - os.walk is optimized
def find_files_fast(root: Path) -> list[Path]:
    files = []
    for dirpath, _, filenames in os.walk(root):
        files.extend(Path(dirpath) / f for f in filenames if f.endswith(".tsv"))
    return files
```

### Pitfall 3: Mutable Fields in Frozen Dataclass
**What goes wrong:** frozen=True makes dataclass immutable, but fields containing mutable objects (list, dict) can still be modified
**Why it happens:** Immutability is shallow - prevents reassignment but not mutation of field contents
**How to avoid:** Use immutable types for fields (tuple, frozenset, str, int, Path); validate in __post_init__ if needed
**Warning signs:** "Immutable" dataclass state changes during execution

**Example:**
```python
# WRONG - mutable field in frozen dataclass
@dataclass(frozen=True)
class Config:
    paths: list[Path]  # Can still call paths.append()!

# RIGHT - use immutable tuple
@dataclass(frozen=True)
class Config:
    paths: tuple[Path, ...]  # Cannot modify after creation
```

### Pitfall 4: Backward Compatibility Breaking
**What goes wrong:** Changing output structure breaks existing summary reports, user scripts, or downstream tools
**Why it happens:** External tools may hardcode path assumptions (e.g., `output/sample/plasmid/comparison_result.summary.tsv`)
**How to avoid:**
- Keep file discovery depth-agnostic (already using os.walk - good)
- Document migration path for users
- Consider symlinks or deprecation warnings for old paths
**Warning signs:** Summary reports fail to find files after pipeline run

### Pitfall 5: os.path.join() vs pathlib.Path() Mixing
**What goes wrong:** Mixing string paths (os.path) with Path objects causes type errors or unexpected behavior
**Why it happens:** Path / str works, but str + Path doesn't; subprocess.run() needs str not Path
**How to avoid:**
- Convert at boundaries: Path objects internally, str() when calling subprocess
- Be consistent within a module
- Use Union[str, Path] types if accepting both
**Warning signs:** TypeError: unsupported operand type(s) for /: 'str' and 'Path'

**Example:**
```python
# WRONG - mixing types
base: str = "output"
subdir: Path = Path("alignments")
full = base / subdir  # TypeError!

# RIGHT - consistent types
base = Path("output")
subdir = Path("alignments")
full = base / subdir

# RIGHT - convert at boundary
subprocess.run(["minimap2", "-o", str(output_path)])
```

## Code Examples

Verified patterns from official sources:

### Creating OutputPaths Dataclass
```python
# Source: Python dataclasses documentation + pathlib best practices
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True, slots=True)
class OutputPaths:
    """Centralized output path construction for PlasmiCheck pipeline.

    All pipeline components should construct paths via this class.
    Provides type safety and single source of truth for output structure.

    Args:
        output_root: Base output directory (e.g., "output")
        sample_name: Sanitized sample identifier
        plasmid_name: Sanitized plasmid identifier
    """
    output_root: Path
    sample_name: str
    plasmid_name: str

    @property
    def base(self) -> Path:
        """Sample/plasmid base directory."""
        return self.output_root / self.sample_name / self.plasmid_name

    @property
    def alignments_dir(self) -> Path:
        return self.base / "alignments"

    @property
    def references_dir(self) -> Path:
        return self.base / "references"

    @property
    def data_dir(self) -> Path:
        return self.base / "data"

    @property
    def reports_dir(self) -> Path:
        return self.base / "reports"

    @property
    def plots_dir(self) -> Path:
        return self.base / "plots"

    # Alignment files
    @property
    def plasmid_bam(self) -> Path:
        return self.alignments_dir / "plasmid_alignment.bam"

    @property
    def spliced_human_bam(self) -> Path:
        return self.alignments_dir / "spliced_human_alignment.bam"

    @property
    def spliced_bam(self) -> Path:
        return self.alignments_dir / "spliced_alignment.bam"

    # Reference files
    @property
    def plasmid_fasta(self) -> Path:
        # Use plasmid_name stem for consistent naming
        return self.references_dir / f"{self.plasmid_name}.fasta"

    @property
    def plasmid_index(self) -> Path:
        return self.references_dir / f"{self.plasmid_name}.mmi"

    @property
    def spliced_fasta(self) -> Path:
        return self.references_dir / "spliced_reference.fasta"

    @property
    def spliced_index(self) -> Path:
        return self.references_dir / "spliced_reference.mmi"

    # Data files
    @property
    def comparison_reads(self) -> Path:
        return self.data_dir / "comparison_result.reads_assignment.tsv"

    @property
    def comparison_summary(self) -> Path:
        return self.data_dir / "comparison_result.summary.tsv"

    @property
    def cdna_positions(self) -> Path:
        return self.data_dir / "cDNA_positions.txt"

    # Reports
    @property
    def report_interactive(self) -> Path:
        return self.reports_dir / "report_interactive.html"

    @property
    def report_static(self) -> Path:
        return self.reports_dir / "report_non_interactive.html"

    def ensure_directories(self) -> None:
        """Create all output subdirectories."""
        for directory in [
            self.alignments_dir,
            self.references_dir,
            self.data_dir,
            self.reports_dir,
            self.plots_dir,
        ]:
            directory.mkdir(parents=True, exist_ok=True)
```

### Using OutputPaths in Pipeline
```python
# Source: Dataclass pattern usage in run_pipeline.py context
from pathlib import Path

def run_combination(
    plasmid_file: str,
    sequencing_file: str,
    output_folder: str,
) -> None:
    # Create paths object
    paths = OutputPaths(
        output_root=Path(output_folder),
        sample_name=sanitize_filename(Path(sequencing_file).stem),
        plasmid_name=sanitize_filename(Path(plasmid_file).stem),
    )

    # Ensure directories exist
    paths.ensure_directories()

    # Use paths throughout pipeline - no hardcoded strings
    convert_plasmidfile_to_fasta(
        plasmid_file,
        str(paths.plasmid_fasta),  # Convert Path to str at boundary
        # ...
    )

    align_reads(
        str(paths.plasmid_index),
        sequencing_file,
        str(paths.plasmid_bam),
        # ...
    )

    compare_alignments(
        str(paths.plasmid_bam),
        str(paths.spliced_human_bam),
        str(paths.data_dir / "comparison_result"),  # Prefix only, no .tsv
        # ...
    )
```

### Depth-Agnostic Summary Reader
```python
# Source: os.walk patterns + current generate_summary_reports.py
from pathlib import Path

def read_compare_outputs(
    input_dir: str,
    substring_to_remove: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read comparison outputs from any directory structure.

    Finds .reads_assignment.tsv and .summary.tsv files at any depth,
    extracts sample/plasmid from parent directories.

    Works with both flat and subfolder structures:
    - output/sample/plasmid/file.tsv (old)
    - output/sample/plasmid/data/file.tsv (new)
    """
    import pandas as pd

    reads_files = find_tsv_files(input_dir, ".reads_assignment.tsv")
    summary_files = find_tsv_files(input_dir, ".summary.tsv")

    if not reads_files or not summary_files:
        raise ValueError("No valid TSV files found")

    reads_dfs: list[pd.DataFrame] = []
    for file_path in reads_files:
        path = Path(file_path)

        # Extract sample/plasmid from path structure
        # New structure: .../sample/plasmid/data/file.tsv
        # Old structure: .../sample/plasmid/file.tsv
        # Strategy: plasmid is always 2 or 3 levels up
        if path.parent.name == "data":
            plasmid_name = path.parent.parent.name
            sample_name = path.parent.parent.parent.name
        else:
            plasmid_name = path.parent.name
            sample_name = path.parent.parent.name

        if substring_to_remove:
            sample_name = sample_name.replace(substring_to_remove, "")

        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        df["Sample"] = sample_name
        df["Plasmid"] = plasmid_name
        reads_dfs.append(df)

    # Same logic for summary files...

    return pd.concat(reads_dfs, ignore_index=True), summary_df
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| os.path string manipulation | pathlib.Path objects | Python 3.4+ (2014) | Cleaner code, cross-platform, but 3-5x slower (negligible in I/O) |
| Scattered path construction | Centralized dataclass | Modern Python (2020+) | Type safety, single source of truth, easier refactoring |
| @dataclass defaults | @dataclass(frozen=True, slots=True) | Python 3.10+ (2021) | Immutability, memory efficiency, hashability |
| glob.glob(recursive=True) | os.walk() for deep trees | Ongoing (rglob perf issue) | 40-100x faster on deeply nested directories |

**Deprecated/outdated:**
- Manual path concatenation: Use pathlib / operator
- os.path.* functions in new code: Use pathlib.Path methods (except os.walk for recursion)
- Mutable config objects: Use frozen dataclasses for configuration

## Open Questions

Things that couldn't be fully resolved:

1. **GenBank file preservation path**
   - What we know: convert_plasmidfile_to_fasta preserves .gb files in output_subfolder for resistance gene detection
   - What's unclear: Should .gb files go in references/ or data/ subfolder?
   - Recommendation: Place in references/ alongside .fasta; they're reference data not analysis output

2. **MD5 checksum file location**
   - What we know: utils.write_md5sum creates md5sum.txt in output_subfolder
   - What's unclear: Should md5sum.txt be in data/ or at base level?
   - Recommendation: Keep at base level (plasmid/sample dir) - it checksums files across all subdirs

3. **Archive behavior with subfolders**
   - What we know: archive_output_folder creates tar.gz of output_subfolder
   - What's unclear: Does it handle subdirectories correctly? Test needed.
   - Recommendation: Verify with integration test; archive should work unchanged

4. **Human reference index location**
   - What we know: Human index is created in same directory as human_fasta (reference/ top-level)
   - What's unclear: Should this remain separate from per-sample outputs?
   - Recommendation: Keep unchanged - human reference is shared across all samples, not sample-specific

## Sources

### Primary (HIGH confidence)
- [Python pathlib documentation](https://docs.python.org/3/library/pathlib.html) - Official pathlib API
- [Python dataclasses documentation](https://docs.python.org/3/library/dataclasses.html) - Official dataclass parameters
- [Python os.walk documentation](https://docs.python.org/3/library/os.html#os.walk) - Official os.walk API
- [Python pathlib: The Complete Guide for 2026](https://devtoolbox.dedyn.io/blog/python-pathlib-complete-guide) - Modern pathlib patterns
- [Python Dataclasses: The Complete Guide for 2026](https://devtoolbox.dedyn.io/blog/python-dataclasses-guide) - Modern dataclass patterns

### Secondary (MEDIUM confidence)
- [Pathlib vs. OS.Path: Choosing the Right Tool](https://www.oreateai.com/blog/pathlib-vs-ospath-choosing-the-right-tool-for-file-path-management-in-python/d20b79071092d10c7894d8ce93715e72) - Performance comparison
- [Dataclass patterns](https://grueter.dev/blog/dataclass-patterns/) - Best practices including frozen, slots
- [Path.rglob performance issues](https://github.com/python/cpython/issues/102613) - Known rglob() slowness on deep nesting
- [Choose Your Path: Python's pathlib vs os.path](https://medium.com/@barila/choose-your-path-pythons-pathlib-vs-os-path-4de0b1e752dd) - When to use each
- [Statically enforcing frozen data classes](https://rednafi.com/python/statically-enforcing-frozen-dataclasses/) - Frozen dataclass patterns

### Tertiary (LOW confidence)
- [Path Glob v.s. OS Walk](https://medium.com/@avq5ac1/path-glob-v-s-os-walk-4278d6add6ac) - Community comparison (verify performance claims)
- [7 Python Best Practices](https://python.plainenglish.io/7-python-best-practices-so-powerful-i-refactored-my-entire-codebase-in-a-weekend-9bb66fc2e77e) - General patterns (verify dataclass recommendations)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All stdlib, well-documented, used in existing codebase
- Architecture: HIGH - Patterns verified from official docs and existing code
- Pitfalls: HIGH - Documented in official sources and bug trackers

**Research date:** 2026-02-15
**Valid until:** 2026-03-15 (30 days - stable stdlib features)
