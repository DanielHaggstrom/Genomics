# PLANS.md

This file is the shared execution plan for work that is too large for a quick patch.

It is intentionally narrow. It is not a product wishlist.

## North Star

Deliver one sane, reproducible research-use pipeline:

`CRAM/CRAI -> validated reference + staged Linux-native inputs -> monogenic triage + PRS scaffold generation + optional PRS imputation + PRS scoring -> one report`

## What Success Looks Like

- one command from a local sample manifest to final results
- no legacy execution path
- no private data tracked in git
- no repo-root runtime sprawl
- stage-boundary reuse works and is easy to reason about
- performance is good enough for overnight personal-genome runs on a normal workstation
- PRS is calculated from a target-genome representation that matches how `pgsc_calc` is actually intended to be used

## Current Priorities

### 1. Rebuild PRS around the standard scaffold-first path

Goal:
- stop treating the PRS-targeted dense VCF as the canonical PRS substrate
- build a broad common-variant scaffold from the CRAM, then score from that scaffold directly or from an imputed scaffold

Why this is now the right path:
- the pipeline already proved that the underlying WGS data are present and usable
- the current blocker is target-genome representation, not missing DNA
- `pgsc_calc` documentation says the workflow works best with imputed genomes and standardized target data
- a real WGS prototype that clears the overlap bar uses a broad HGDP+1kGP-style scaffold, not a PRS-target-locus VCF

Path of highest sanity:
- keep the repo local-first
- keep monogenic unchanged
- make the canonical PRS path:
  - `CRAM -> PRS scaffold VCF -> pgsc_calc`
  - optional: `CRAM -> PRS scaffold VCF -> imputed PRS scaffold VCF -> pgsc_calc`
- keep external `reuse.imputed_prs_vcf` support, but only as an explicit stage-boundary import

Concrete next steps:
- done: `analysis.prs.reference_db` now defaults to `pgsc_HGDP+1kGP_v1`
- done: PRS scaffold site generation now starts from the official reference database and layers score-extra harmonized loci on top
- done: PRS scaffold shards are now materialized from whole-contig gVCFs, then subset back to scaffold sites with explicit non-variant calls preserved
- done: `pgsc_calc` scoring now uses the repo-built compatibility surface directly instead of the repo forcing the older dense-target path
- done: local Beagle imputation now acts as a scaffold-imputation layer, including passthrough handling for panel-unavailable or no-overlap shards
- done: direct-scaffold and scaffold+imputed paths have both been exercised on the real sample with `min_overlap` fixed at `0.75`
- next: if PRS work resumes, shift effort from pipeline plumbing to score-level portability diagnostics and target-representation experiments

Decision rules:
- if direct scaffold scoring clears one or more intended scores above `0.75`, keep direct scaffold as a supported PRS path
- if direct scaffold still misses but scaffold+imputation rescues scores, keep imputation as the preferred PRS path
- if neither clears `0.75`, stop widening scope further and reassess the selected scores rather than adding more pipeline complexity by default

Research basis:
- [pgsc_calc prepare input genomes](https://pgsc-calc.readthedocs.io/en/latest/how-to/prepare.html)
- [pgsc_calc samplesheet docs](https://pgsc-calc.readthedocs.io/en/latest/how-to/samplesheet.html)
- [pgsc_calc match-rate docs](https://pgsc-calc.readthedocs.io/en/latest/explanation/match.html)
- [PGS Catalog downloads](https://www.pgscatalog.org/downloads/)
- [pgsc_calc discussion #123](https://github.com/PGScatalog/pgsc_calc/discussions/123)
- [WGS scaffold prototype gist](https://gist.github.com/ckandoth/e8d064685b5012937a75d2b7fca584c5)
- [pgsc_calc reference database docs](https://pgsc-calc.readthedocs.io/en/dev/how-to/database.html)

### 2. Improve throughput without widening scope

Goal:
- reduce wall-clock time while keeping the same canonical pipeline

Path of highest sanity:
- keep WSL-native staging as the default
- make worker selection less conservative where it actually helps
- keep shard reuse intact across reruns
- improve progress visibility before adding new operator complexity

Concrete next steps:
- make worker planning aware of the new scaffold-first PRS stages, not only the older dense-target path
- confirm the stop-adjust-resume path remains reliable after the PRS refactor
- investigate safe live tuning during a run; if true live changes are not reliable, keep the documented graceful stop-and-resume path

Operator-facing target behavior:
- support automatic worker defaults derived from `runtime.max_cpus`, `runtime.max_memory`, current shard plan, and observed free disk headroom
- keep manual override knobs narrow and explicit:
  - `analysis.monogenic.workers`
  - `analysis.prs.workers`
  - `analysis.monogenic.pairhmm_threads`
  - `analysis.prs.pairhmm_threads`
- keep clear operator commands in the README:
  - `bash run.sh manifests/local/<sample>.yaml`
  - `bash run.sh manifests/local/<sample>.yaml -resume`
  - `nextflow log <run-name>`
  - `docker stats`

Explicit non-goals for this phase:
- no new workflow engine
- no Kubernetes/cloud detour
- no new UI

### 3. Harden cache and import behavior

Goal:
- make restarts and imports predictable

Concrete next steps:
- confirm completed shard reuse after interruption for the scaffold path
- tighten cache-completeness checks where needed for scaffold and imputed-scaffold artifacts
- ensure imported artifacts are clearly marked in manifests and reports
- keep `reuse.prs_dense_vcf` temporarily for compatibility, but document it everywhere as the PRS scaffold artifact until runtime naming can be cleaned up safely

### 4. Clarify monogenic vs PRS execution modes

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
- keep the shared preparation path common to both analyses
- make outputs clearly separated by analysis type even when launched together
- keep the unified report for `both`, but allow analysis-specific outputs to stand on their own when only one branch is enabled

### 5. Keep reporting useful and minimal

Goal:
- one report that is readable, cautious, and actually helpful

Concrete next steps:
- keep one primary Markdown artifact
- publish reader-facing stage outputs under `output/results/...`, not only in Linux-side cache
- clearly distinguish:
  - PRS scaffold provenance
  - PRS imputation provenance
  - direct-scaffold vs imputed-scaffold scoring inputs
- preserve technical TSV/VCF artifacts as machine-readable companions
- improve wording only when it increases clarity, not length

### 6. Streamline code without forcing clever abstractions

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
- remove dead helpers that are no longer used
- trim obviously duplicated small stage helpers where a single shared helper is clearer than two copies
- keep deleting vestigial code whenever a refactor makes a helper unnecessary

Ranked worthwhile medium refactors:
- keep PRS scaffold generation separate from PRS scoring logic, so the new scaffold path does not disappear into one giant stage function
- split `workflow/core/reference.py::cmd_prepare_context` into a few clear helper steps
- lightly deduplicate stage orchestration patterns shared by monogenic and PRS only when the shared helper is simpler than the duplicated code

Areas to leave alone unless pain increases:
- `workflow/reporting/report.py`
- the split between `workflow/core`, `workflow/stages`, and `workflow/reporting`
- tracked example manifests
- `runtime.py` as the boring utility bucket

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
- canonical PRS strategy is now local-first hybrid:
  - local scaffold-first path is primary
  - external imputed-VCF import is explicit and supported, but secondary

## Status Log

- Done: repo reduced to one canonical CRAM-first pipeline surface.
- Done: reference probing, validation, and WSL-native staging are implemented.
- Done: monogenic completes end-to-end from the CRAM-first path on the real sample.
- Done: the pipeline now proves the core issue is target-genome representation, not missing genomic data.
- Done: PRS harmonization diagnostics now show that our own score ingestion is no longer the main loss bucket.
- Done: local Beagle-based PRS imputation infrastructure exists and can complete shard-level work against a real panel.
- Done: reader-facing monogenic and PRS artifacts are published into `output/results/...`.
- Done: execution mode is explicit and derived as `monogenic-only`, `prs-only`, or `both`.
- Done: resource-aware worker planning is implemented for monogenic and PRS shard work.
- Done: cache-side shard progress now records reused vs rebuilt shard completions.
- Done: the scaffold-first PRS path is now the implementation backbone, and the older PRS-targeted dense flow is no longer the canonical runtime path.
- Done: the full scaffold-first + imputed PRS path completes end to end on the real sample.
- Done: the remaining PRS blocker is no longer Docker, Beagle, contigs, chrX handling, metadata propagation, or scoring-surface plumbing.
- Current measured result: all requested scores still fail `min_overlap: 0.75` even after scaffold-first imputation, with best overlaps around `PGS002255=67.62%` and `PGS000838=62.95%`.
- Next likely improvement: if PRS work continues, focus on score portability analysis, alternate/commoner cohort-style target representations, or narrower score selection rather than more pipeline format fixes.
