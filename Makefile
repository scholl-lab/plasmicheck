.PHONY: help install install-dev clean lint format format-check typecheck test test-fast ci-check fix

CYAN  := \033[0;36m
GREEN := \033[0;32m
YELLOW := \033[0;33m
NC    := \033[0m

PYTHON := python

help: ## Show this help message
	@echo "$(CYAN)Available targets:$(NC)"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  $(GREEN)%-20s$(NC) %s\n", $$1, $$2}'

##@ Installation

install: ## Install the package
	uv pip install -e .

install-dev: ## Install with dev dependencies
	uv pip install -e ".[dev]"

##@ Cleaning

clean: ## Clean build artifacts and caches
	@echo "$(YELLOW)Cleaning...$(NC)"
	rm -rf build/ dist/ *.egg-info .pytest_cache .mypy_cache .ruff_cache htmlcov/ .coverage coverage.xml
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	@echo "$(GREEN)Cleaned$(NC)"

##@ Code Quality

lint: ## Run ruff linter
	@echo "$(CYAN)Linting...$(NC)"
	$(PYTHON) -m ruff check plasmicheck/ tests/
	@echo "$(GREEN)Lint passed$(NC)"

format: ## Auto-format code with ruff
	@echo "$(CYAN)Formatting...$(NC)"
	$(PYTHON) -m ruff format plasmicheck/ tests/
	$(PYTHON) -m ruff check --fix plasmicheck/ tests/
	@echo "$(GREEN)Formatted$(NC)"

format-check: ## Check formatting without modifying
	@echo "$(CYAN)Checking format...$(NC)"
	$(PYTHON) -m ruff format --check --diff plasmicheck/ tests/
	@echo "$(GREEN)Format OK$(NC)"

typecheck: ## Run mypy type checker
	@echo "$(CYAN)Type checking...$(NC)"
	$(PYTHON) -m mypy plasmicheck/
	@echo "$(GREEN)Type check passed$(NC)"

##@ Testing

test: ## Run all tests with coverage
	@echo "$(CYAN)Testing...$(NC)"
	$(PYTHON) -m pytest tests/ --verbose --cov=plasmicheck --cov-report=term-missing --cov-report=xml
	@echo "$(GREEN)Tests passed$(NC)"

test-fast: ## Run fast tests only (no integration/slow)
	$(PYTHON) -m pytest -m "not slow and not integration" tests/ --tb=short -q

##@ CI Verification

ci-check: ## Run all CI checks locally (lint + format + typecheck + test)
	@echo "$(CYAN)Running CI checks...$(NC)"
	@echo ""
	@echo "$(CYAN)[1/4] Linting...$(NC)"
	@$(MAKE) lint
	@echo ""
	@echo "$(CYAN)[2/4] Format check...$(NC)"
	@$(MAKE) format-check
	@echo ""
	@echo "$(CYAN)[3/4] Type check...$(NC)"
	@$(MAKE) typecheck
	@echo ""
	@echo "$(CYAN)[4/4] Tests...$(NC)"
	@$(MAKE) test-fast
	@echo ""
	@echo "$(GREEN)ALL CI CHECKS PASSED$(NC)"

##@ Quick Commands

fix: format lint ## Auto-fix then lint
	@echo "$(GREEN)Fixed$(NC)"
