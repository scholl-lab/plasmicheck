.PHONY: help install install-dev hooks clean lint lint-actions format format-check typecheck test test-fast test-integration build smoke-package check ci-check fix

UV ?= uv
RUN = $(UV) run --locked
JOBS ?= 4

help: ## Show available targets
	@awk 'BEGIN {FS = ":.*## "} /^[a-zA-Z_-]+:.*## / {printf "  %-20s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install locked runtime dependencies
	$(UV) sync --locked --no-dev

install-dev: ## Install locked development environment
	$(UV) sync --locked

hooks: ## Install pre-commit and pre-push git hooks
	$(RUN) pre-commit install --hook-type pre-commit --hook-type pre-push

lint: ## Check all Python files with the locked Ruff version
	$(RUN) ruff check .

lint-actions: ## Validate GitHub workflows (requires Docker)
	docker run --rm -v "$(CURDIR):/repo:ro" -w /repo rhysd/actionlint:1.7.12@sha256:b1934ee5f1c509618f2508e6eb47ee0d3520686341fec936f3b79331f9315667

format: ## Apply safe lint fixes and format all Python files
	$(RUN) ruff check --fix .
	$(RUN) ruff format .

format-check: ## Check formatting without modifying files
	$(RUN) ruff format --check .

typecheck: ## Strictly type-check package, scripts, and tests
	$(RUN) mypy

check: ## Run fast static checks (lint, format-check, typecheck) in parallel
	@$(MAKE) -j $(JOBS) --output-sync=target lint format-check typecheck

test: ## Run portable tests with the 70% branch-aware coverage floor
	$(RUN) pytest -m "not slow and not integration" --cov=plasmicheck --cov-report=term-missing:skip-covered --cov-report=xml

test-fast: ## Run portable tests without coverage for quick iteration
	$(RUN) pytest -m "not slow and not integration" -q

test-integration: ## Run synthetic pipeline tests; requires minimap2 and samtools >=1.19
	@command -v minimap2 >/dev/null || { echo "Install minimap2 to run integration tests"; exit 1; }
	@command -v samtools >/dev/null || { echo "Install samtools >=1.19 to run integration tests"; exit 1; }
	$(RUN) pytest -m integration

build: ## Build the source distribution and wheel
	rm -rf dist
	$(UV) build --no-sources

smoke-package: build ## Smoke-test built wheel and sdist in an isolated temp venv outside checkout
	@tmp="$$(mktemp -d)"; trap 'rm -rf "$$tmp"' EXIT; \
	$(UV) export --locked --no-dev --no-emit-project --format requirements-txt > "$$tmp/reqs.txt"; \
	for artifact in "$(CURDIR)"/dist/*.whl "$(CURDIR)"/dist/*.tar.gz; do \
		$(UV) venv --quiet --clear "$$tmp/venv"; \
		$(UV) pip sync --quiet --python "$$tmp/venv/bin/python" "$$tmp/reqs.txt"; \
		$(UV) pip install --quiet --python "$$tmp/venv/bin/python" --no-deps "$$artifact"; \
		(cd "$$tmp" && "$$tmp/venv/bin/plasmicheck" --version >/dev/null && "$$tmp/venv/bin/python" -c 'import sys; from pathlib import Path; import plasmicheck; from plasmicheck.resources import get_resource_path; assert Path(plasmicheck.__file__).is_relative_to(Path(sys.prefix)); [ None for r in ("config.json", "py.typed", "templates/report_template.html", "templates/summary_template.html", "static/img/plasmicheck_ascii.txt", "static/img/plasmicheck_logo.svg", "static/img/plasmicheck_logo_200px.png", "static/img/plasmicheck_logo_600px.png") if not get_resource_path(r).is_file() and (_ for _ in ()).throw(AssertionError(r)) ]'); \
	done

ci-check: check test smoke-package ## Run full CI gates locally (including lint-actions & integration when tools exist)
	@if command -v docker >/dev/null 2>&1; then $(MAKE) lint-actions; fi
	@if command -v minimap2 >/dev/null 2>&1 && command -v samtools >/dev/null 2>&1; then $(MAKE) test-integration; fi

fix: format lint ## Format and lint

clean: ## Remove project build outputs and tool caches
	rm -rf build dist .pytest_cache .mypy_cache .ruff_cache htmlcov .coverage coverage.xml
