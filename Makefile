SHELL := /usr/bin/env bash

RUN ?= manifests/examples/sample.yaml
EXTRA ?=
PYTHON ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)

.PHONY: help run status validate validate-manifest test clean clean-cache clean-all

help:
	@echo "Targets:"
	@echo "  run                - run the pipeline (set RUN=manifests/<sample>.yaml)"
	@echo "  status             - print cache/log/stage status locations for a manifest"
	@echo "  validate           - run syntax checks"
	@echo "  validate-manifest  - validate a manifest (set RUN=manifests/<sample>.yaml)"
	@echo "  test               - run unit tests"
	@echo "  clean              - remove repo-facing results only"
	@echo "  clean-cache        - remove Linux-side cache/work state for this project"
	@echo "  clean-all          - remove both results and cache"

run:
	@test -n "$(RUN)" || (echo "Set RUN=manifests/<sample>.yaml"; exit 1)
	@bash run.sh "$(RUN)" $(EXTRA)

status:
	@test -n "$(RUN)" || (echo "Set RUN=manifests/<sample>.yaml"; exit 1)
	@PYTHONPATH=. $(PYTHON) -m workflow.core.manifest status --manifest "$(RUN)"

validate:
	@$(PYTHON) -m py_compile workflow/*.py workflow/core/*.py workflow/stages/*.py workflow/reporting/*.py tests/*.py
	@bash -n run.sh
	@echo "Validation checks passed."

validate-manifest:
	@test -n "$(RUN)" || (echo "Set RUN=manifests/<sample>.yaml"; exit 1)
	@PYTHONPATH=. $(PYTHON) -m workflow.core.manifest validate --manifest "$(RUN)" --schema manifests/schema/run_manifest.schema.json --output /tmp/genomics_manifest.validate.json
	@echo "Manifest is valid."

test:
	@$(PYTHON) -m unittest discover -s tests -p "test_*.py"

clean:
	@rm -rf output/results

clean-cache:
	@CACHE_ROOT="$$(PYTHONPATH=. $(PYTHON) -m workflow.core.manifest print-project-cache-root --manifest "$(RUN)")" && rm -rf "$$CACHE_ROOT"

clean-all: clean clean-cache
	@find workflow tests -type d -name "__pycache__" -prune -exec rm -rf {} +
