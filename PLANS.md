# PLANS.md

This file is the shared execution plan for work that is too large for a quick patch.

It is intentionally narrow. It is not a product wishlist.

## North Star

Deliver one sane, reproducible research-use pipeline:

`CRAM/CRAI -> validated reference + staged Linux-native inputs -> monogenic triage + PRS scoring -> one report`

## What Success Looks Like

- one command from a local sample manifest to final results
- no legacy execution path
- no private data tracked in git
- no repo-root runtime sprawl
- stage-boundary reuse works and is easy to reason about
- performance is good enough for overnight personal-genome runs on a normal workstation

## Current Priorities

### 1. Stabilize the live end-to-end run

Goal:
- finish a full CRAM-first run on the real sample
- confirm that monogenic and PRS both complete from the same command

Definition of done:
- final report exists
- stage manifests are `ok`
- no manual rescue commands were required after launch

### 2. Improve throughput without widening scope

Goal:
- reduce wall-clock time while keeping the same canonical pipeline

Path of highest sanity:
- keep WSL-native staging as the default
- make worker selection less conservative, especially for monogenic
- keep shard reuse intact across reruns
- improve progress visibility before adding new features

Concrete next steps:
- expose clearer stage progress in logs or manifests
- define explicit operator commands for tuning and resuming runs
- investigate safe live tuning during a run; if true live changes are not reliable, support a documented graceful stop-and-resume path that preserves completed shard work

Operator-facing target behavior:
- support automatic worker defaults derived from `runtime.max_cpus`, `runtime.max_memory`, current shard plan, and observed free disk headroom
- keep manual override knobs narrow and explicit:
  - `analysis.monogenic.workers`
  - `analysis.prs.workers`
  - `analysis.monogenic.pairhmm_threads`
  - `analysis.prs.pairhmm_threads`
- add clear operational commands to the README once implemented, for example:
  - `bash run.sh manifests/local/<sample>.yaml`
  - `bash run.sh manifests/local/<sample>.yaml -resume`
  - `nextflow log <run-name>`
  - `docker stats`
- investigate whether safe live tuning is possible with supported tool/runtime commands; if not, document the sane path as:
  - let current shards finish
  - stop the run cleanly
  - adjust worker settings
  - relaunch with cache reuse

Explicit non-goals for this phase:
- no new workflow engine
- no Kubernetes/cloud detour
- no new UI

Status:
- automatic worker planning is now implemented for monogenic and PRS dense genotyping
- narrow per-stage override knobs are now part of the manifest schema
- live tuning during an active run is still future work

### 3. Harden cache and import behavior

Goal:
- make restarts and imports predictable

Concrete next steps:
- confirm completed shard reuse after interruption
- tighten cache-completeness checks where needed
- ensure imported artifacts are clearly marked in manifests and reports

### 4. Rescue PRS without widening scope too early

Goal:
- improve PRS compatibility using the existing WGS data before considering larger workflow changes

Path of highest sanity:
- first improve harmonization of scorefiles against the dense PRS genotype set
- only consider imputation if direct harmonized WGS scoring still cannot produce defensible overlap

Concrete next steps:
- add support for `hm_inferOtherAllele` when `other_allele` is absent in PGS Catalog scorefiles
- improve handling of multiallelic and otherwise awkward scorefile rows so more harmonized loci survive matching
- evaluate a GATK non-variant-sites PRS export path for direct WGS compatibility:
  - targeted `HaplotypeCaller` on PRS loci
  - gVCF-style output where appropriate
  - `GenotypeGVCFs --include-non-variant-sites true`
  - final PRS-facing VCF normalization without `<NON_REF>` leftovers
- compare that direct non-variant-site export path against the current dense PRS VCF path and keep only the saner one
- rerun PRS from cache after harmonization changes and compare match rates score-by-score
- if overlap is still insufficient, evaluate imputation as an explicit scope decision rather than a silent pipeline default

Imputation scope note:
- imputation is not "inventing a second genome"; it is a standard statistical inference step that uses linkage disequilibrium plus a reference panel to estimate genotypes at common variants in a form that many PRS pipelines expect
- for this repo, imputation should be treated as an optional compatibility layer for PRS, not as a required baseline for all analysis
- do not add imputation until we have measured the best direct-from-WGS harmonized result first
- if added later, keep it narrow:
  - PRS-only
  - clearly marked in manifests and reports
  - reusable as a stage-boundary artifact
  - documented with the exact reference panel, tool, and quality thresholds used

Direct-WGS compatibility note:
- the pgsc_calc discussion at `PGScatalog/pgsc_calc#123` suggests a standard route for including reference calls at PRS loci:
  - call targeted loci from BAM/CRAM
  - genotype with non-variant sites included
  - emit a PRS-facing VCF that contains both variant and homozygous-reference positions
- for this repo, keep that idea single-sample and cacheable; do not add multisample combine steps unless project scope changes

### 5. Clarify monogenic vs PRS execution modes

Goal:
- acknowledge that monogenic triage and PRS are different analysis products while keeping one sane repo surface

Path of highest sanity:
- keep one shared codebase and one shared CRAM/reference preparation path
- support three explicit execution modes:
  - monogenic only
  - PRS only
  - both
- avoid splitting into two competing repos or two unrelated pipeline implementations

Concrete next steps:
- make the supported execution modes explicit in README and manifests
- ensure stage toggles are first-class and easy to reason about
- keep shared preparation (`manifest -> context -> reference validation -> staging`) common to both analysis types
- make outputs clearly separated by analysis type even when launched together
- keep the unified report for `both`, but allow analysis-specific outputs to stand on their own when only one branch is enabled

Design guardrails:
- do not reintroduce legacy alternative execution paths
- do not duplicate monogenic and PRS infrastructure that should remain shared
- if user-facing commands diverge later, they should still dispatch into the same underlying pipeline surface

### 6. Keep reporting useful and minimal

Goal:
- one report that is readable, cautious, and actually helpful

Concrete next steps:
- keep one primary Markdown artifact
- preserve technical TSV/VCF artifacts as machine-readable companions
- improve wording only when it increases clarity, not length

### 7. Streamline code without forcing clever abstractions

Goal:
- reduce accidental complexity and dead code without hiding real genomics complexity

Current audit summary:
- the repo is no longer carrying two pipelines
- most remaining live code is concentrated where it should be:
  - `workflow/stages/monogenic.py`
  - `workflow/stages/prs.py`
  - `workflow/core/reference.py`
- the best remaining simplification opportunities are local refactors, not another architectural rewrite

Ranked quick wins:
- remove dead helpers that are no longer used:
  - `workflow/stages/monogenic.py::parse_summary_file`
  - `workflow/core/runtime.py::resolve_docker_bin`
- trim obviously duplicated small stage helpers where a single shared helper is clearer than two copies
- keep deleting vestigial code whenever a refactor makes a helper unnecessary

Ranked worthwhile medium refactors:
- unify annotated-VCF / CSQ walking in the monogenic stage
  - today both priority extraction and ClinVar triage parse the same annotated VCF shape separately
  - extract one shared iterator or row-normalization helper and keep the selection logic separate
- split `workflow/core/reference.py::cmd_prepare_context` into a few clear helper steps
  - input + disk/executable checks
  - CRAM header / shard planning
  - reference bundle resolution
  - final context assembly
- lightly deduplicate stage orchestration patterns shared by monogenic and PRS
  - disabled-stage short circuit
  - cache-hit reuse
  - standard failure payload shape
  - import provenance summary
  - do this only if the shared helper is simpler than the duplicated code

Areas to leave alone unless pain increases:
- `workflow/reporting/report.py`
  - manual Markdown is acceptable here; do not add a templating layer just to save lines
- the split between `workflow/core`, `workflow/stages`, and `workflow/reporting`
  - this is currently a sane boundary, not over-modularization
- tracked example manifests
  - the separate `both`, `monogenic-only`, and `prs-only` examples improve clarity more than they cost in lines
- `runtime.py` as the boring utility bucket
  - keep it consolidated unless a helper group clearly wants its own file

Guardrails for this cleanup phase:
- prefer fewer, clearer functions over more files
- do not refactor just to reduce LoC on a graph
- do not add abstractions that make stage behavior harder to follow
- only keep refactors that make the pipeline easier to trust, test, or resume

## Things We Should Not Add Unless Scope Changes

- clinical interpretation framework beyond research triage
- SV/CNV/trio workflows
- extra configuration surfaces
- multiple competing report styles
- general lab-pipeline tooling unrelated to this repo's core path

## Shared Decisions

- `AGENTS.md` is tracked and repo-wide.
- `PLANS.md` is tracked and shared.
- private notes, manifests, presets, results, and raw data stay untracked.
- the repo should favor fewer modules and fewer docs, as long as readability stays good.

## Status Log

- Done: repo reduced to one canonical CRAM-first pipeline surface.
- Done: reference probing, validation, and WSL-native staging are implemented.
- Done: PRS now builds a dense genotype VCF from CRAM before scoring.
- Done: monogenic now completes end-to-end from the CRAM-first path on the real sample.
- Done: PRS compatibility export now builds a sex-aware PLINK pfile from the dense PRS VCF.
- Done: PRS union-site building now supports `hm_inferOtherAllele` fallback and keeps valid multiallelic score rows instead of discarding them outright.
- Done: PRS dense-genotype cache lineage now changes when harmonization logic changes, so improved targeting is not masked by stale cached artifacts.
- Done: execution mode is now explicit and derived as `monogenic-only`, `prs-only`, or `both`, with tracked example manifests and report output that adapts to the selected mode.
- Done: resource-aware worker planning is now implemented for monogenic and PRS dense-genotype stages, with optional narrow manifest overrides for workers and `pairhmm_threads`.
- In progress: PRS direct-from-WGS harmonization; scoring still fails the configured overlap threshold in the last measured run.
- Next likely improvement: expose clearer stage progress and a documented stop-adjust-resume path before deciding whether any true live tuning support is worth adding.
