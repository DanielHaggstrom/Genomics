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
- long-running monogenic and PRS stages now emit incremental stage-status updates during shard work
- `run.sh` now forwards extra Nextflow flags, so `bash run.sh manifests/local/<sample>.yaml -resume` is a supported operator path
- `make status RUN=...` now reports the cache/log/stage-status files to watch during a run
- shard progress now distinguishes reused vs rebuilt shard completions in both live status and final stage manifests
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
- keep `min_overlap` fixed at `0.75`; improve matching quality instead of relaxing the acceptance bar

Concrete next steps:
- keep using harmonized PGS Catalog scorefiles (`hm_chr`, `hm_pos`, `effect_allele`, `other_allele`, `hm_inferOtherAllele`) as the canonical score definition surface
  - official basis: the PGS Catalog harmonized files are intended to provide consistent chromosome/position/allele annotations for matching target data
  - source: [PGS Catalog downloads](https://www.pgscatalog.org/downloads/)
- keep using a dense direct-from-CRAM export that includes homozygous-reference sites
  - official basis: pgsc_calc maintainers describe the expected route as genotyping targeted loci and including non-variant sites before scoring
  - source: [pgsc_calc discussion #123](https://github.com/PGScatalog/pgsc_calc/discussions/123)
- add score-by-score matching diagnostics as a first-class artifact before changing more matching logic
  - parse and publish the detailed `pgsc_calc` matching categories that the official report uses in its variant matching section
  - break failures down by categories such as unmatched, strand ambiguous, duplicates, multiallelic, and excluded rows
  - source: [pgsc_calc outputs/report docs](https://pgsc-calc.readthedocs.io/en/latest/explanation/output.html)
- add a pre-match audit step for our own scorefile ingestion, keyed by score ID
  - report how many rows were accepted by:
    - explicit `other_allele`
    - `hm_inferOtherAllele`
  - report how many rows were skipped due to:
    - missing other allele
    - missing harmonized position
    - reference lookup failure
    - reference not present in score alleles
    - unsupported effect allele
  - treat these counts as a gating diagnostic, not just debug output
- tighten target-genome normalization to match the official pgsc_calc expectations as closely as possible
  - keep target build fixed to the score build; do not rely on liftover inside `pgsc_calc`
  - keep the PRS-facing export restricted to chromosomes `1-22`, `X`, `Y`
  - keep the export biallelic and free of `<NON_REF>` placeholders before PLINK conversion
  - verify that reference-call rows survive all the way into the final dense VCF/pfile
  - source: [pgsc_calc samplesheet docs](https://pgsc-calc.readthedocs.io/en/latest/how-to/samplesheet.html)
- prioritize fixes using the current measured failure buckets
  - first: skipped score rows on our side (`missing_other_allele`, `reference_not_in_score_alleles`, `unsupported_effect_allele`)
  - second: whatever the official `pgsc_calc` detailed matching log shows as the dominant excluded/unmatched categories
  - third: small compatibility annoyances that affect traceability rather than matching itself (for example, erased long variant IDs in PLINK)
- after each harmonization change, rerun PRS from cache and compare match rates score-by-score against the previous run
  - success criterion for this phase is not "some PRS output exists"
  - success criterion is "at least one intended score crosses the unchanged `0.75` threshold without weakening acceptance rules"
- if overlap is still insufficient after the direct-harmonization work, evaluate imputation as an explicit scope decision rather than a silent pipeline default

Imputation scope note:
- imputation is not "inventing a second genome"; it is a standard statistical inference step that uses linkage disequilibrium plus a reference panel to estimate genotypes at common variants in a form that many PRS pipelines expect
- official basis: pgsc_calc documentation says the workflow works best with imputed input genomes and treats imputed dosage data as a standard target-genome format
  - source: [pgsc_calc prepare input genomes](https://pgsc-calc.readthedocs.io/en/latest/how-to/prepare.html)
  - source: [pgsc_calc samplesheet docs](https://pgsc-calc.readthedocs.io/en/latest/how-to/samplesheet.html)
- for this repo, imputation should be treated as an optional compatibility layer for PRS, not as a required baseline for all analysis
- do not add imputation until we have measured the best direct-from-WGS harmonized result first
- if added later, keep it narrow:
  - PRS-only
  - clearly marked in manifests and reports
  - reusable as a stage-boundary artifact
  - documented with the exact reference panel, tool, and quality thresholds used

Decision update:
- direct-from-WGS harmonization has now been pushed far enough to show the real boundary:
  - the dense PRS export includes homozygous-reference sites
  - our own score-ingestion audit now skips very few score rows
  - the remaining failures are happening in downstream score matching, not because the target export is still sparse
- because of that, the next sane PRS rescue step is now to add an explicit imputation-compatible path rather than continuing to pretend that direct hard-called WGS alone will necessarily meet the calculator's portability threshold for every intended score
- initial implementation rule:
  - first make imputation a real manifest/stage-boundary concept in the repo
  - then add the local backend against a pinned reference panel
  - do not try to hide an unfinished local backend behind silent fallbacks
- first local-backend caveat:
  - because this repo is single-sample, do not wire `conform-gt` into the canonical local path by default
  - official `conform-gt` docs say the target VCF should have at least 20 individuals so its statistical strand checks are meaningful
  - for the first local backend, require same-build reference-panel resources and use normalized target VCFs directly with Beagle instead

### 4a. Add A Standard PRS-Compatible Imputation Layer

Goal:
- calculate PRS from the genome using a more standard common-variant workflow, while keeping monogenic analysis unchanged

Why this is now justified:
- the `pgsc_calc` docs explicitly say the workflow works best with output from an imputation server such as Michigan or TOPMed
  - source: [pgsc_calc prepare input genomes](https://pgsc-calc.readthedocs.io/en/latest/how-to/prepare.html)
- the `pgsc_calc` samplesheet docs explicitly support dosage-oriented target data
  - source: [pgsc_calc samplesheet docs](https://pgsc-calc.readthedocs.io/en/latest/how-to/samplesheet.html)
- SHAPEIT5 documentation describes common-variant phasing as the first step for large WES/WGS datasets
  - source: [SHAPEIT5 phase_common](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/)
- Beagle documentation explicitly supports phasing plus imputation from a reference panel and states that study and reference markers must have identical `CHROM`, `POS`, `REF`, and `ALT`, with `conform-gt` recommended to align the target to the reference panel
  - source: [Beagle 5.4 page](https://faculty.washington.edu/browning/beagle/b5_4.html)
  - source: [Beagle 5.4 documentation](https://faculty.washington.edu/browning/beagle/beagle_5.4_18Mar22.pdf)
- public imputation servers document the exact target-genome conventions they expect: per-chromosome `vcf.gz`, sorted variants, correct build/chromosome naming, and pre-imputation QC
  - source: [TOPMed data preparation](https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/)
  - source: [Michigan Imputation Server docs](https://imputationserver.readthedocs.io/en/latest/)

High-sanity design choice:
- keep one canonical local PRS path in the repo
- add an optional PRS-only imputation stage between dense PRS genotypes and scoring
- allow external stage-boundary import of an imputed VCF if the user prefers an external imputation server, but do not make remote upload the canonical default

Canonical PRS-imputation path to implement:
1. `CRAM -> dense common-variant PRS scaffold VCF`
   - keep the current dense PRS export logic, but treat it as the scaffold for phasing/imputation rather than the final scoring input
2. `normalize + conform target to reference panel`
   - split/decompose multiallelics where needed
   - left-normalize alleles against the validated FASTA
   - ensure reference alleles are correct
   - conform marker positions / allele order / strand to the chosen imputation reference panel
   - official basis: Beagle requires identical marker definitions between target and reference and recommends `conform-gt`
3. `phase target common variants`
   - preferred first implementation path: Beagle, because it can phase and impute in one local toolchain for a single sample
   - future upgrade path if needed: SHAPEIT5 for a more specialized WGS phasing stage
4. `impute common variants from a phased reference panel`
   - output dosage-aware VCF/BCF where possible
   - keep the imputation backend and reference panel fully recorded in stage manifests
5. `imputation QC`
   - keep only variants that pass explicit quality thresholds from the imputation tool/server
   - record INFO/R2 thresholds in manifests and reports
6. `pgsc_calc scoring from imputed target data`
   - prefer dosage-compatible input if available
   - keep `min_overlap` unchanged at `0.75`
   - compare direct-from-WGS vs imputed match rates score-by-score

Implementation constraints:
- do not change the monogenic path
- do not make the whole pipeline depend on remote upload
- do not create a second competing repo or second orchestration surface
- keep imputation as `analysis.prs.imputation.enabled: true|false`
- keep stage-boundary imports explicit:
  - `reuse.imputed_prs_vcf`
- publish imputation-facing deliverables under `output/results/.../prs/`

Initial implementation order:
1. choose and pin one local imputation backend plus one reference panel family
   - default choice for first implementation: local Beagle-based path
   - rationale: it is documented, local, and directly supports phasing/imputation from a reference panel
2. add a reference-panel manifest type for PRS imputation resources
3. add target normalization + `conform-gt` preparation before imputation
4. add the imputation stage and cache it by:
   - sample fingerprint
   - reference FASTA identity
   - imputation reference panel identity
   - backend version
   - QC thresholds
5. add PRS scoring from the imputed artifact
6. compare overlap against the current direct-from-WGS branch before deciding whether to keep both modes or make imputation the default PRS path

Decision criteria after implementation:
- if imputation lifts one or more intended scores above `0.75`, keep it as the preferred PRS path
- if imputation still does not rescue the intended scores, stop widening scope and reassess the score set itself rather than adding more complexity by default

Direct-WGS compatibility note:
- the pgsc_calc discussion at `PGScatalog/pgsc_calc#123` suggests a standard route for including reference calls at PRS loci:
  - call targeted loci from BAM/CRAM
  - genotype with non-variant sites included
  - emit a PRS-facing VCF that contains both variant and homozygous-reference positions
- for this repo, keep that idea single-sample and cacheable; do not add multisample combine steps unless project scope changes
- that route is now effectively in place in this repo; the current blocker is matching quality, not sparse target export

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
- publish reader-facing stage outputs under `output/results/...`, not only in Linux-side cache
  - monogenic TSVs/VCFs that a user is expected to open should be copied or linked into the run results directory
  - PRS result tables/reports should likewise be published into the run results directory
  - keep heavy cache/work state in Linux-native storage, but do not make the user dig through cache paths for final deliverables
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
- Done: the dense PRS export now includes homozygous-reference sites; the current PRS blocker is score matching at the `0.75` threshold, not a sparse variant-only export.
- Done: monogenic now completes end-to-end from the CRAM-first path on the real sample.
- Done: PRS compatibility export now builds a sex-aware PLINK pfile from the dense PRS VCF.
- Done: PRS union-site building now supports `hm_inferOtherAllele` fallback and keeps valid multiallelic score rows instead of discarding them outright.
- Done: PRS dense-genotype cache lineage now changes when harmonization logic changes, so improved targeting is not masked by stale cached artifacts.
- Done: PRS now emits first-class matching diagnostics, including a per-score pre-match audit from our harmonization step and a published match-threshold failure summary artifact from scoring.
- Done: PRS union-site building now rescues score rows whose alleles only match the reference after strand complement, and records rescued-row counts in the diagnostics/report surface.
- Done: PRS union-site building now supports literal multi-base allele strings against the reference sequence, so simple indel-style score rows are no longer dropped as unsupported SNPs.
- Done: PRS union-site building now includes effect-allele-only score rows as targeted positions instead of skipping them before `pgsc_calc` can apply its own effect-only matching logic.
- Done: the report stage now publishes reader-facing monogenic and PRS artifacts into `output/results/...` instead of leaving them only in Linux-side cache paths.
- Done: PRS manifests and stage-boundary imports now have an explicit imputation surface:
  - `analysis.prs.imputation`
  - `reuse.imputed_prs_vcf`
  - report/manifests now distinguish dense-direct PRS input from imported imputed PRS input
- Done: the first local PRS imputation backend is now wired into the PRS stage:
  - `analysis.prs.imputation.enabled: true` now builds a cached imputed PRS VCF instead of failing fast
  - local Beagle JAR download is pinned and cached under the Linux-side tool cache
  - per-chromosome panel manifests are now a first-class config concept via `manifests/*/imputation_panels/`
  - the report now surfaces local imputation backend/panel details
- In progress: the local PRS imputation backend has not yet been exercised end-to-end against a real GRCh38 reference panel in this repo.
- Done: a real local GRCh38 panel family has now been chosen and staged for the first end-to-end imputed PRS run:
  - IGSR / 1000 Genomes GRCh38 phased panel from `20190312_biallelic_SNV_and_INDEL`
  - Beagle GRCh38 `no_chr_in_chrom_field` genetic maps
  - local ignored manifest path prepared under `manifests/local/imputation_panels/`
- Done: execution mode is now explicit and derived as `monogenic-only`, `prs-only`, or `both`, with tracked example manifests and report output that adapts to the selected mode.
- Done: resource-aware worker planning is now implemented for monogenic and PRS dense-genotype stages, with optional narrow manifest overrides for workers and `pairhmm_threads`.
- Done: cache-side shard progress now records reused vs rebuilt shard completions for monogenic and PRS dense genotyping.
- Done: monogenic annotated-VCF / `CSQ` walking now flows through one shared iterator, with focused helper tests covering multiallelic carried-allele filtering and downstream priority/ClinVar selection.
- Done: the final report now surfaces imported-artifact provenance, source paths, trust level, and source stage manifests more clearly for monogenic and PRS dense inputs.
- Done: stage-level reuse checks are now stricter about complete artifacts, so cached VCF/BCF-style outputs must have their required inferred sidecars before a stage is treated as reusable.
- In progress: PRS direct-from-WGS harmonization; scoring still fails the configured overlap threshold in the last measured run.
- In progress: local PRS imputation backend; the repo now supports imputed-VCF scoring imports, but the local Beagle/reference-panel path is still the next implementation step.
- Next likely improvement: confirm completed shard reuse after interruption now that resume/visibility are easier to use, then decide whether any true live tuning support is worth adding.
