SHELL := /usr/bin/env bash

.PHONY: help run validate validate-manifest test clean

help:
	@echo "Targets:"
	@echo "  run                - run the pipeline (requires RUN=manifests/<sample>.yaml)"
	@echo "  validate           - run syntax checks"
	@echo "  validate-manifest  - validate a manifest (requires RUN=manifests/<sample>.yaml)"
	@echo "  test               - run unit tests"
	@echo "  clean              - remove generated output/"

run:
	@test -n "$(RUN)" || (echo "Set RUN=manifests/<sample>.yaml"; exit 1)
	@bash run.sh "$(RUN)"

validate:
	@python3 -m py_compile workflow/core/*.py workflow/monogenic/*.py workflow/prs/*.py tests/*.py
	@bash -n run.sh workflow/core/*.sh workflow/prs/*.sh workflow/monogenic/annotate_with_vep.sh
	@echo "Validation checks passed."

validate-manifest:
	@test -n "$(RUN)" || (echo "Set RUN=manifests/<sample>.yaml"; exit 1)
	@PYTHONPATH=. python3 -m workflow.core.validate_manifest --manifest "$(RUN)" --schema manifests/schema/run_manifest.schema.json --output /tmp/genomics_manifest.validate.json
	@echo "Manifest is valid."

test:
	@python3 -m unittest discover -s tests -p "test_*.py"

clean:
	@rm -rf output
	@find workflow tests -type d -name "__pycache__" -prune -exec rm -rf {} +
