.PHONY: help install install-dev clean lint lint-actions format format-check typecheck test test-fast test-integration build ci-check fix

UV ?= uv
RUN = $(UV) run --locked

help: ## Show available targets
	@awk 'BEGIN {FS = ":.*## "} /^[a-zA-Z_-]+:.*## / {printf "  %-20s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install locked runtime dependencies
	$(UV) sync --locked --no-dev

install-dev: ## Install locked development environment
	$(UV) sync --locked

lint: ## Check all Python files with the locked Ruff version
	$(RUN) ruff check .

lint-actions: ## Validate GitHub workflows (requires Docker)
	docker run --rm -v "$(CURDIR):/repo:ro" -w /repo rhysd/actionlint:1.7.12@sha256:b1934ee5f1c509618f2508e6eb47ee0d3520686341fec936f3b79331f9315667

format: ## Apply safe lint fixes and format all Python files
	$(RUN) ruff check --fix .
	$(RUN) ruff format .

format-check: ## Check formatting without modifying files
	$(RUN) ruff format --check .

typecheck: ## Strictly type-check the complete package
	$(RUN) mypy plasmicheck

test: ## Run portable tests with the 70% branch-aware coverage floor
	$(RUN) pytest -m "not slow and not integration" --cov=plasmicheck --cov-report=term-missing:skip-covered --cov-report=xml

test-fast: ## Run portable tests without coverage for quick iteration
	$(RUN) pytest -m "not slow and not integration" -q

test-integration: ## Run synthetic pipeline tests; requires minimap2 and samtools >=1.19
	@command -v minimap2 >/dev/null || { echo "Install minimap2 to run integration tests"; exit 1; }
	@command -v samtools >/dev/null || { echo "Install samtools >=1.19 to run integration tests"; exit 1; }
	$(RUN) pytest -m integration

build: ## Build the source distribution and wheel
	$(UV) build --no-sources

ci-check: lint format-check typecheck test build ## Run the portable CI checks

fix: format lint ## Format and lint

clean: ## Remove project build outputs and tool caches
	rm -rf build dist .pytest_cache .mypy_cache .ruff_cache htmlcov .coverage coverage.xml
