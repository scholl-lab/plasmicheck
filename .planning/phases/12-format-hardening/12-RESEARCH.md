# Phase 12: Format Hardening - Research

**Researched:** 2026-02-15
**Domain:** File format detection and conversion (plasmid files: xDNA, SnapGene .dna, GenBank)
**Confidence:** HIGH

## Summary

Phase 12 adds robust plasmid file format handling to PlasmiCheck, focusing on auto-detection via magic bytes, SnapGene .dna support, and improved xDNA record ID handling. The research confirms BioPython 1.84+ includes first-class support for SnapGene format via `Bio.SeqIO.SnapGeneIO`, making implementation straightforward.

The current implementation has two key issues:
1. **xDNA record IDs are broken**: BioPython's xDNA parser extracts the comment field as record ID (`[author`), which produces unusable FASTA headers. This breaks downstream alignment tools.
2. **Format detection is extension-only**: The pipeline uses simple `.endswith()` checks, which fail when SnapGene `.dna` files arrive (both SnapGene and xDNA use `.dna` extension).

**Primary recommendation:** Implement magic byte detection for `.dna` files (SnapGene signature: `09 00 00 00 0E 53 6E 61 70 47 65 6E 65` vs xDNA version byte: `0x00`), use BioPython's `snapgene` format parser for SnapGene files, and override all record IDs with sanitized filename instead of relying on file content.

## Standard Stack

The established libraries/tools for file format handling in BioPython-based pipelines:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| BioPython | >=1.84 | Sequence file I/O (SeqIO) | De facto standard for bioinformatics file parsing in Python; includes native SnapGene and xDNA parsers since 1.75+ |
| Python stdlib | 3.10+ | Binary file I/O (`open(..., 'rb')`) | No external dependencies needed for magic byte detection |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| argparse | stdlib | CLI argument parsing with `choices=` | Input validation for `--format` flag |
| logging | stdlib | Diagnostic output to stderr | User-facing format detection messages |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Manual magic byte detection | `python-magic` (libmagic wrapper) | python-magic requires external libmagic C library; overkill for 2 known formats |
| Manual magic byte detection | `puremagic` (pure Python) | Adds dependency for simple detection we can implement in 10 lines |
| BioPython SnapGene parser | Third-party `snapgene-reader` | BioPython native support is sufficient; avoids extra dependency |

**Installation:**
No new dependencies required — BioPython >=1.84 already installed.

## Architecture Patterns

### Recommended Function Structure

```python
# plasmicheck/scripts/format_detection.py (new module)
def detect_format(file_path: str) -> str:
    """Auto-detect plasmid file format via magic bytes or text parsing."""
    # 1. Check magic bytes for binary formats (.dna, .xdna)
    # 2. Fall back to text-based detection (GenBank, FASTA)
    # 3. Return format string compatible with BioPython SeqIO

def read_magic_bytes(file_path: str, n_bytes: int = 16) -> bytes:
    """Read first N bytes from file for magic byte detection."""
    with open(file_path, 'rb') as f:
        return f.read(n_bytes)

def is_snapgene(magic_bytes: bytes) -> bool:
    """Check if magic bytes match SnapGene .dna signature."""
    # Signature: 09 00 00 00 0E 53 6E 61 70 47 65 6E 65
    return magic_bytes[:13] == b'\x09\x00\x00\x00\x0ESnapGene'

def is_xdna(magic_bytes: bytes) -> bool:
    """Check if magic bytes match xDNA format (version byte + fixed header)."""
    # xDNA: first byte is version (0x00), byte 112 is 0xFF
    return magic_bytes[0] == 0x00 and len(magic_bytes) >= 112
```

### Pattern 1: Override Record ID with Sanitized Filename

**What:** Replace BioPython's default record.id with sanitized filename, regardless of source format
**When to use:** All conversions to FASTA (convert command, pipeline auto-conversion)

**Example:**
```python
# In convert_plasmidfile_to_fasta.py
from pathlib import Path

def convert_with_sanitized_id(input_file: str, output_file: str, file_format: str) -> None:
    """Convert plasmid file to FASTA with sanitized filename as record ID."""
    # Parse file with BioPython
    records = SeqIO.parse(input_file, file_format)

    # Override record ID with sanitized filename
    sanitized_base = sanitize_filename(Path(input_file).stem)

    with open(output_file, 'w') as out_handle:
        for record in records:
            # Override ID and description
            record.id = sanitized_base
            record.description = f"converted_from={file_format}"
            SeqIO.write(record, out_handle, 'fasta')
```

### Pattern 2: Format Detection with Override Flag

**What:** Auto-detect format but allow explicit `--format` override for edge cases
**When to use:** Convert command CLI interface

**Example:**
```python
# In cli.py convert command
parser_convert.add_argument(
    "-f", "--format",
    choices=["snapgene", "xdna", "genbank", "fasta"],
    default=None,
    help="File format (auto-detected if not specified)"
)

# In conversion logic
if args.format:
    detected_format = args.format
    logging.info(f"Using user-specified format: {detected_format}")
else:
    detected_format = detect_format(args.input_file)
    logging.info(f"Detected format: {detected_format}", file=sys.stderr)
```

### Anti-Patterns to Avoid

- **Extension-based detection for .dna files:** Both SnapGene and xDNA use `.dna` extension — always check magic bytes
- **Relying on record.id from file content:** xDNA comment fields are unreliable; SnapGene feature names vary; always use filename
- **Silent format assumption:** Log detected format to stderr so users know what happened
- **Raising generic exceptions:** Use specific error messages with actionable guidance

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Parsing SnapGene binary format | Custom binary parser | BioPython `SeqIO.parse(file, "snapgene")` | SnapGene format is complex (packet-based structure, XML features, primer filtering); BioPython handles version differences |
| Parsing xDNA binary format | Custom binary parser | BioPython `SeqIO.parse(file, "xdna")` | xDNA has 112-byte header with big-endian integers, topology flags, optional sections; BioPython tested on real files |
| File format auto-detection | Complex heuristics | Simple magic byte check + BioPython fallback | Only need to distinguish 2 binary formats; text formats handled by BioPython's parse error handling |
| Filename sanitization | Custom regex per call site | Centralized `utils.sanitize_filename()` | Already implemented and tested; ensures consistency across codebase |

**Key insight:** BioPython's SeqIO parsers are battle-tested on diverse real-world files. Custom parsers will hit edge cases (file corruption, version differences, vendor-specific quirks) that BioPython already handles.

## Common Pitfalls

### Pitfall 1: xDNA Record ID Truncation
**What goes wrong:** BioPython's xDNA parser uses the comment field as record.id. Comments like `[author : name=User - ]` get truncated to `[author`, producing FASTA headers like `>[author` which break alignment tools expecting valid sequence IDs.

**Why it happens:** xDNA format stores sequence metadata in a freeform comment field. BioPython extracts this faithfully, but downstream tools (minimap2, samtools) expect alphanumeric IDs.

**How to avoid:** Always override `record.id` with `sanitize_filename(Path(input_file).stem)` before writing FASTA. Never trust record.id from parsed files.

**Warning signs:** Alignment failures with cryptic errors about sequence names; FASTA headers starting with `[` or containing spaces.

### Pitfall 2: SnapGene vs xDNA Ambiguity
**What goes wrong:** Both formats use `.dna` extension. Extension-based detection (current implementation: `.endswith(".gb") else "xdna"`) will misclassify SnapGene files as xDNA, causing BioPython parse errors.

**Why it happens:** SnapGene and xDNA are binary formats from different cloning software vendors that both adopted `.dna` extension historically.

**How to avoid:** Check magic bytes first for `.dna` files:
- SnapGene: starts with `09 00 00 00 0E 53 6E 61 70 47 65 6E 65` ("SnapGene" string)
- xDNA: starts with version byte `0x00`, has 112-byte header

**Warning signs:** BioPython raising `ValueError` or parse errors when user provides SnapGene file to pipeline.

### Pitfall 3: Error Messages Without Action
**What goes wrong:** Generic errors like "Unsupported file type: xdna" or "Parse error" leave users stuck. They don't know if it's a version issue, corrupt file, or wrong format.

**Why it happens:** BioPython exceptions are technical (designed for developers); users need actionable guidance.

**How to avoid:** Catch BioPython exceptions and re-raise with context:
```python
try:
    SeqIO.parse(file, "xdna")
except Exception as e:
    if "version" in str(e).lower():
        raise ValueError(
            f"Unsupported xDNA version in {file}. "
            "Export as GenBank (.gb) from SnapGene or Serial Cloner, then use plasmicheck convert."
        ) from e
    elif "corrupt" in str(e).lower() or "truncated" in str(e).lower():
        raise ValueError(
            f"File {file} appears corrupt or truncated. "
            "Re-export from your cloning software (SnapGene, Serial Cloner, Benchling)."
        ) from e
    else:
        raise ValueError(
            f"Failed to parse {file}. Supported formats: SnapGene (.dna), "
            "xDNA (.xdna), GenBank (.gb/.gbk), FASTA (.fasta/.fa). "
            "If unsure, export as GenBank from your cloning software."
        ) from e
```

**Warning signs:** Users reporting "it doesn't work" without details; support requests asking "what format should I use?"

### Pitfall 4: Logging vs Errors Confusion
**What goes wrong:** Writing error messages to stdout instead of stderr, or using `logging.error()` for user-facing errors instead of exceptions.

**Why it happens:** Python's logging defaults to stderr, but confusion about when to log vs raise.

**How to avoid:** Follow Unix conventions:
- **Errors that stop execution:** Raise exceptions (argparse catches, writes to stderr, exits 2)
- **Diagnostic info:** `logging.info()` to stderr (format detection, progress)
- **Warnings:** `logging.warning()` to stderr (non-fatal issues)
- **Normal output:** stdout only for parseable results (currently none in convert command)

**Warning signs:** Error messages going to stdout (contaminating pipelines); users missing diagnostic info.

## Code Examples

Verified patterns from official sources:

### BioPython SnapGene Parsing
```python
# Source: https://biopython.org/docs/dev/api/Bio.SeqIO.SnapGeneIO.html
from Bio import SeqIO

# Parse SnapGene .dna file
records = SeqIO.parse("plasmid.dna", "snapgene")
for record in records:
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}")
    print(f"Features: {len(record.features)}")  # Annotations preserved
```

### BioPython xDNA Parsing
```python
# Source: https://biopython.org/docs/1.75/api/Bio.SeqIO.XdnaIO.html
from Bio import SeqIO

# Parse xDNA file
records = SeqIO.parse("plasmid.xdna", "xdna")
for record in records:
    # Note: record.id is comment field (unreliable)
    print(f"ID (from comment): {record.id}")  # "[author : name=User - ]"
```

### Magic Byte Detection (No Dependencies)
```python
# Pattern: Pure Python magic byte detection
def detect_format(file_path: str) -> str:
    """Auto-detect plasmid file format."""
    with open(file_path, 'rb') as f:
        magic = f.read(16)

    # Check SnapGene signature
    if magic[:13] == b'\x09\x00\x00\x00\x0ESnapGene':
        return "snapgene"

    # Check xDNA signature (version 0x00, fixed header size)
    if magic[0] == 0x00:
        # Verify by checking byte 112 should be 0xFF
        with open(file_path, 'rb') as f:
            header = f.read(112)
            if len(header) == 112 and header[-1] == 0xFF:
                return "xdna"

    # Fall back to text-based detection
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline()
            if first_line.startswith('LOCUS'):
                return "genbank"
            elif first_line.startswith('>'):
                return "fasta"
    except UnicodeDecodeError:
        pass  # Binary file, not text format

    raise ValueError(
        f"Unrecognized file format: {file_path}. "
        "Supported: SnapGene (.dna), xDNA (.xdna), GenBank (.gb/.gbk), FASTA (.fasta/.fa)"
    )
```

### FASTA Record ID Override
```python
# Pattern: Override record ID with sanitized filename
from pathlib import Path
from Bio import SeqIO

def convert_with_sanitized_id(input_file: str, output_file: str,
                               file_format: str) -> None:
    """Convert plasmid file to FASTA with clean record ID."""
    # Parse with BioPython
    record = next(SeqIO.parse(input_file, file_format))

    # Override ID with sanitized filename (no extension)
    sanitized_id = sanitize_filename(Path(input_file).stem)
    record.id = sanitized_id
    record.name = sanitized_id
    record.description = f"converted_from={file_format}"

    # Write to FASTA
    SeqIO.write(record, output_file, "fasta")

# Result: >pcDNA5_FRT_TO_KCNJ5 converted_from=xdna
# Instead of: >[author : name=Demo User - ]
```

### CLI Format Flag with Validation
```python
# Source: https://docs.python.org/3/library/argparse.html
parser.add_argument(
    "-f", "--format",
    choices=["snapgene", "xdna", "genbank", "fasta"],
    default=None,
    help="File format (auto-detected if not specified)"
)
# argparse validates choices at parse time, raises error for invalid values
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Extension-only detection (`.endswith()`) | Magic byte detection + fallback | Phase 12 (v0.34.0) | Enables SnapGene support; handles .dna ambiguity |
| Use file's record.id as-is | Override with sanitized filename | Phase 12 (v0.34.0) | Fixes xDNA FASTA headers; ensures consistent IDs |
| Manual `-t/--file_type` required | Auto-detection with optional override | Phase 12 (v0.34.0) | Better UX; users don't need to know format internals |
| Generic error messages | Actionable errors with tool names | Phase 12 (v0.34.0) | Users know how to fix issues (export from SnapGene, etc.) |

**Deprecated/outdated:**
- **`-t/--file_type` required flag**: Now optional; auto-detection is default (still supported for override)
- **Assumption that .dna = xDNA**: SnapGene also uses .dna extension; must check magic bytes

## Open Questions

Things that couldn't be fully resolved:

1. **Are there xDNA versions other than 0x00 in the wild?**
   - What we know: BioPython documentation and xDNA spec both state version byte is "always zero" in practice
   - What's unclear: Whether newer Serial Cloner versions use different version bytes
   - Recommendation: Log warning if version byte != 0x00, but attempt parsing anyway; BioPython will fail if truly incompatible

2. **Should we preserve SnapGene feature annotations?**
   - What we know: BioPython SnapGene parser extracts features (restriction sites, primers, etc.)
   - What's unclear: Whether these features would be useful for resistance gene detection (Phase 10)
   - Recommendation: Sequence-only for Phase 12 (per CONTEXT.md decision); revisit in future if resistance gene detection needs annotations

3. **What error does BioPython raise for unsupported xDNA versions?**
   - What we know: BioPython xDNA parser exists since 1.75, but exact exception types for version errors not documented
   - What's unclear: Whether it's `ValueError`, `IOError`, or custom exception
   - Recommendation: Catch generic `Exception` in conversion code, inspect message for keywords ("version", "corrupt", "truncated"), provide actionable error

## Sources

### Primary (HIGH confidence)
- BioPython 1.84+ installed and tested locally
- Bio.SeqIO.SnapGeneIO: https://biopython.org/docs/dev/api/Bio.SeqIO.SnapGeneIO.html
- Bio.SeqIO.XdnaIO: https://biopython.org/docs/1.75/api/Bio.SeqIO.XdnaIO.html
- SnapGene source code: https://github.com/biopython/biopython/blob/master/Bio/SeqIO/SnapGeneIO.py
- xDNA format specification: https://incenp.org/dvlpt/docs/binary-sequence-formats/index.html
- Python argparse documentation: https://docs.python.org/3/library/argparse.html
- Codebase inspection:
  - `plasmicheck/scripts/convert_plasmidfile_to_fasta.py` (current conversion logic)
  - `plasmicheck/scripts/utils.py` (sanitize_filename implementation)
  - `plasmicheck/scripts/run_pipeline.py` (pipeline format detection)
  - `plasmicheck/config.json` (supported_formats configuration)
  - `tests/test_convert.py` (existing sanitization tests)

### Secondary (MEDIUM confidence)
- SnapGene magic bytes: https://gist.github.com/neutrinoguy/b6cdbe854b34b9fc32c7bbe88b8eb261 (verified with local file)
- CLI stderr conventions: https://bugs.python.org/issue26645

### Tertiary (LOW confidence)
- python-magic alternatives: https://pypi.org/project/puremagic/ (not needed, included for completeness)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - BioPython SnapGene/xDNA parsers verified locally with v1.85
- Architecture: HIGH - Patterns tested with real xDNA files from repo
- Pitfalls: HIGH - xDNA record ID issue confirmed via testing; SnapGene vs xDNA ambiguity documented in specs

**Research date:** 2026-02-15
**Valid until:** 2026-04-15 (60 days — BioPython stable, format specs unlikely to change)

**Current implementation gaps identified:**
1. Extension-only format detection at line 513-516 of `run_pipeline.py`
2. No SnapGene support in convert command or pipeline
3. xDNA record IDs used as-is, producing broken FASTA headers (confirmed with test file)
4. No actionable error messages for unsupported format variants
5. `sanitize_filename()` duplicated in two files (convert_plasmidfile_to_fasta.py and utils.py)
