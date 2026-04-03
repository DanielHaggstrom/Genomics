# AGENTS.md

This repository uses a small amount of shared agent guidance on purpose.

## Purpose

Keep this project on the path of highest sanity:
- one canonical pipeline
- one narrow scope
- one clean repo surface
- one reader-facing report

The supported product is:

`single-sample CRAM/CRAI -> monogenic triage + broad PRS scaffold genotypes + optional PRS imputation + PRS scoring -> one report`

## Scope Guardrails

In scope:
- single-sample research-use analysis
- CRAM-first execution
- monogenic SNV/indel triage
- PRS from a broad common-variant scaffold generated from CRAM
- optional PRS-only imputation as a compatibility layer on that scaffold
- explicit stage-boundary imports only
- reference validation, caching, reporting, and performance work directly needed for the above

Out of scope unless the user explicitly changes project scope:
- clinical diagnosis claims
- SV/CNV pipelines
- trio or pedigree workflows
- general-purpose bioinformatics utilities
- extra dashboards, web apps, databases, or orchestration layers
- broad framework abstractions that do not materially simplify this repo

## Use `PLANS.md`

Use [PLANS.md](./PLANS.md) for complex or multi-step work, especially when the task:
- takes longer than a quick fix
- spans multiple files or stages
- changes architecture, pipeline shape, caching, or performance behavior
- needs a handoff-friendly execution record

When using `PLANS.md`:
- keep it concrete and runnable
- update status rather than creating extra planning docs
- prefer one living plan over multiple overlapping notes

## Repo Conventions

- Keep the root lean. Do not add new top-level docs or directories unless clearly necessary.
- Prefer consolidation over many tiny wrappers.
- Keep one canonical execution path. Avoid reintroducing legacy flows.
- Reuse/imports must stay at explicit stage boundaries only.
- Keep tracked examples generic and safe to share.
- Keep personal manifests, reference bundles, presets, raw data, outputs, and secrets out of git.
- Heavy runtime state belongs in Linux-native cache, not the repo.
- Reader-facing outputs belong under `output/results/`.

## Privacy And Safety

Never track:
- `manifests/local/`
- raw genomics files
- personal score presets
- generated reports or runtime state
- machine-specific local paths unless intentionally genericized

If privacy-sensitive data appears in tracked files, fix the current tree and consider history rewrite if needed.

## Change Strategy

Prefer the smallest change that:
- improves correctness
- reduces complexity
- preserves resumability and cache reuse
- keeps the repo understandable to a new reader

When choosing between "more features" and "clearer pipeline", choose the clearer pipeline.

## Validation Expectations

For meaningful changes, run the smallest relevant checks available, typically:
- manifest validation
- Python syntax/tests
- stage-specific smoke checks
- a real pipeline run when the change affects execution behavior

## OpenAI/Codex Docs

For questions about Codex/OpenAI behavior, prefer official OpenAI documentation and official OpenAI domains.
