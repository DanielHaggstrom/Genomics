SHELL := /usr/bin/env bash

.PHONY: help bootstrap validate prepublish clean monogenic prs pipeline

help:
	@echo "Targets:"
	@echo "  bootstrap  - create local config files from examples (if missing)"
	@echo "  validate   - run syntax checks for Python and shell scripts"
	@echo "  prepublish - scan tracked files for obvious secret/path leaks"
	@echo "  clean      - delete generated pipeline runs/artifacts"
	@echo "  monogenic  - run monogenic stage through workflow/run_pipeline.sh"
	@echo "  prs        - run PRS stage through workflow/run_pipeline.sh"
	@echo "  pipeline   - run full pipeline through workflow/run_pipeline.sh"

bootstrap:
	@test -f config/pipeline.env || cp config/pipeline.env.example config/pipeline.env
	@test -f config/secrets.env || cp config/secrets.env.example config/secrets.env
	@test -f .env || cp .env.example .env
	@test -f config/prs_samplesheet.csv || cp config/prs_samplesheet.example.csv config/prs_samplesheet.csv
	@test -f config/prs_pgs_ids.txt || cp config/prs_pgs_ids.example.txt config/prs_pgs_ids.txt
	@test -f config/prs_update_sex.tsv || cp config/prs_update_sex.example.tsv config/prs_update_sex.tsv
	@echo "Bootstrap complete. Edit local files before running."

validate:
	@python3 -m py_compile workflow/monogenic/*.py workflow/utils/*.py
	@bash -n workflow/prs/*.sh workflow/monogenic/annotate_with_vep.sh workflow/run_pipeline.sh workflow/utils/prepublish_check.sh workflow/utils/cleanup_runs.sh
	@echo "Validation checks passed."

prepublish:
	@bash workflow/utils/prepublish_check.sh

clean:
	@bash workflow/utils/cleanup_runs.sh

monogenic:
	@bash workflow/run_pipeline.sh --stage monogenic

prs:
	@bash workflow/run_pipeline.sh --stage prs

pipeline:
	@bash workflow/run_pipeline.sh --stage all
