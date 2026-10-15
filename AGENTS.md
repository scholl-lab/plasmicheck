# Working on PlasmiCheck

Detects plasmid contamination from sequencing data. Keep changes focused and
preserve scientific results, CLI compatibility, and report/TSV schemas.

## Commands

- Setup: `uv sync --locked` (Python 3.12 from `.python-version`).
- Fast loop: `uv run --locked pytest tests/test_<area>.py -q` or `make test-fast`.
- Before finishing: `make ci-check` (Ruff lint/format, strict mypy, coverage, build).
- Pipeline changes: also `make test-integration`; needs minimap2 and samtools >=1.19.
- Workflow changes: `make lint-actions` (Docker); alignment pipelines require Bash.
- Fix formatting: `make format`. Install hooks: `uv run --locked pre-commit install`.
- Update dependencies deliberately with `uv lock`; commit `pyproject.toml` and `uv.lock` together.

## Code map

- `plasmicheck/cli.py`: argument parsing and command dispatch.
- `plasmicheck/scripts/run_pipeline.py`: orchestration and output lifecycle.
- `plasmicheck/scripts/compare_alignments.py`: scoring, streaming merge, verdicts.
- `plasmicheck/scripts/`: alignment, conversion, coverage, reports, and shared utilities.
- `plasmicheck/scripts/plotting/`, `templates/`, `static/`: plots and report assets.
- `plasmicheck/{config,resources,output_paths,thread_config,cache_manager}.py`: shared infrastructure.
- `tests/`: portable tests; `test_integration.py` uses committed synthetic data.
- `scripts/`: optional benchmarks/regression tools; `snakemake/`: workflow integration.
- `.planning/`: historical milestone records, not current implementation instructions.

## Change rules

- Read callers and relevant tests before editing; reuse existing helpers.
- Keep committed documentation self-contained; omit external project comparisons.
- Keep Python 3.10 compatibility; default development and strict typing use 3.12.
- Test changed behavior with small synthetic inputs. Do not require private data.
- Maintain the 70% branch-aware coverage floor in `pyproject.toml`; do not hide
  untested modules or add blanket lint/type suppressions to pass checks.
- Read-name sorting and merge comparison must use the same lexical order.
- Propagate thresholds/options consistently from CLI to comparison and reports.
- Keep cache locks, atomic writes, cleanup, and validation when moving code.
- Never commit patient data, local references, generated reports, or credentials.
- Prefer independent parallel tasks with explicit file ownership when useful.
- Finish with changes, checks actually run, and remaining risks. Do not claim
  hosted CI or publication ran when only local checks were performed.

See [README.md](README.md#development) for setup and release details.
