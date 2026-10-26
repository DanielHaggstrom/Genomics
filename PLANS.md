# Current plan: pruning and resilient PRS execution

The value proposition is useful PRS calculation from one WGS sample, with costly
completed work reusable after interruption. Prefer deletion and explicit stage
boundaries. No full pipeline execution is authorized. Protect all original-data
files and all files of uncertain provenance; never clean their directory trees.

## Implementation

- [x] Replace the outer cache-disabled Nextflow wrapper with one local Python driver.
  Keep Nextflow only for the pinned upstream scorer's actual resumable tasks.
- [x] Default to PRS only; retain optional monogenic triage and PRS-only imputation.
  Consolidate three run examples into one. Remove remote reference probes,
  automatic reference downloading, redundant downloader fallbacks/cache scavenging,
  obsolete helper tests, and duplicate large VCF/report-manifest publication.
- [x] Handle Ctrl-C/SIGTERM, terminate active process groups/direct Docker containers,
  cancel queued tasks, stream bounded-memory logs, and lock a cache against two writers.
- [x] Write metadata atomically; gate generated VCF reuse on successful completion
  receipts with output/index fingerprints. Keep completed chromosome gVCFs shared
  across branches/score selections. Continue partial HTTP downloads.
- [x] Preserve the scientific fixes: quality-filtered explicit scaffold alleles,
  conservative missingness/indels, direct-call precedence over imputation, input-aware
  conversion caches, upstream resource limits, and honest score reporting.
- [x] Keep originals read-only; create missing indexes on cache copies and reject
  source/output overlap. Keep Windows/WSL helper migration deferred as requested.
- [x] Finish focused validation and independent diff review: 85 synthetic/unit tests
  pass, syntax checks and the generic example manifest validate, and the diff has
  no whitespace errors. No real genomic tools or pipeline stages were launched.

## Validation and acceptance

Use `make test`, `make validate`, and example manifest validation. Tests exercise
synthetic genotype arithmetic, stage wiring with mocked external tools, source
protection, report completeness, failed/incomplete outcomes, stop/restart/locking,
streaming logs, and real short Python subprocess cancellation. No personal genomics
pipeline, Docker, GATK, Beagle, PLINK, or Nextflow runs are authorized in this pass.

A passing suite is not sample validation. The remaining acceptance gate, after
explicit authorization, is real-tool integration on a bounded synthetic/truth-known
region, followed by one manageable additive autosomal score. Check genotype
concordance, missing high-weight variants, row-level allele matching and independent
score arithmetic. Keep published score weights intact; do not lower matching gates
to manufacture success. Validate ancestry calibration before presenting percentiles.

Resume boundaries are completed tasks/chromosomes, not the middle of a caller
algorithm. Unclean power loss/SIGKILL can leave processes running and stale status;
HTTP continuation depends on server support. Legacy outputs without receipts may
be recomputed once; old cache trees are preserved, not automatically deleted.
