# Phase 12: Format Hardening - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Reliable conversion of diverse plasmid file formats (xDNA, SnapGene .dna) without manual intervention. Users can pass any supported binary plasmid format to `plasmicheck convert` or directly to the pipeline via `-pf`, and get correct FASTA output with clean record IDs. Unsupported formats produce actionable error messages.

</domain>

<decisions>
## Implementation Decisions

### Error messaging & guidance
- Actionable errors: include the specific problem, what the user can do, and name specific tools (SnapGene, Serial Cloner, Benchling)
- Unsupported xDNA version: "Unsupported xDNA version X. Export as GenBank (.gb) from SnapGene or Serial Cloner, then use plasmicheck convert."
- Unrecognized format: list all supported formats inline in the error message
- Corrupt/truncated files: suggest re-exporting from cloning software
- Error output destination: research CLI best practices during research phase and decide accordingly

### FASTA record ID convention
- Always use sanitized filename (without extension) as record ID, not the original record ID from the file
- Sanitization: alphanumeric + hyphens + underscores only [a-zA-Z0-9_-]
- Reuse existing `sanitize_filename` function from utils.py for consistency
- Include source format in FASTA description line: `>pcDNA5_FRT_TO_KCNJ5 converted_from=xdna`

### SnapGene support scope
- Full first-class support: accept .dna (SnapGene) as input to convert command
- Parse via BioPython's snapgene parser, produce FASTA
- Pipeline (-pf flag) accepts .dna files directly — auto-converts internally
- Sequence only — discard SnapGene feature annotations during conversion
- Auto-detect SnapGene vs xDNA by magic bytes, with --format override flag

### Format detection behavior
- Magic bytes checked first for binary files (.dna, .xdna)
- If binary detection fails, fall back to text-based parsing (try GenBank, then FASTA)
- --format flag accepts fixed choices: {snapgene, xdna, genbank, fasta} — validated at parse time
- Log detected format to stderr: "Detected format: SnapGene (.dna)" for user awareness
- Corrupt files get actionable error with re-export suggestion

### Claude's Discretion
- Error output destination (stderr only vs stderr + log file) — research best practices
- Exact magic byte detection implementation details
- How to handle edge cases in text-based fallback parsing
- Logging verbosity level conventions (align with existing pipeline patterns)

</decisions>

<specifics>
## Specific Ideas

- Error messages should name specific tools by name (SnapGene, Serial Cloner, Benchling) — users in this domain know these tools
- Follow the "actionable error" pattern consistently: problem + what to do + tool suggestion
- FASTA description line should include provenance: `converted_from=xdna` or `converted_from=snapgene`

</specifics>

<deferred>
## Deferred Ideas

- SnapGene feature annotation preservation for resistance gene detection — could be valuable but out of scope for format hardening (potential future enhancement)

</deferred>

---

*Phase: 12-format-hardening*
*Context gathered: 2026-02-15*
