# Genomics

Single-sample, research-use genome analysis from `CRAM/CRAI` to:
- monogenic triage
- dense PRS genotypes plus PRS scoring
- one final Markdown report

Supported execution modes:
- `both`
- `monogenic-only`
- `prs-only`

## Scope

Supported:
- one sample at a time
- CRAM-first SNV/indel monogenic triage
- PRS from a dense genotype VCF generated from the CRAM
- explicit stage-boundary imports for:
  - called VCF
  - annotated VCF
  - dense PRS VCF

Out of scope:
- clinical diagnosis
- SV/CNV pipelines
- trio/family workflows
- general-purpose bioinformatics utilities unrelated to the above

## Run

Preferred command:

```bash
bash run.sh manifests/local/<sample>.yaml
```

Direct command:

```bash
nextflow run . -params-file manifests/local/<sample>.yaml -profile docker
```

The wrapper resolves the Linux-side cache root from the manifest and keeps Nextflow runtime state there, so the repo only keeps final results.

## Pipeline Shape

Canonical path:

```text
CRAM/CRAI
-> validated manifest
-> prepared context (reference probe + bundle validation + WSL-native staging)
-> monogenic stage
-> dense PRS genotypes
-> PRS scoring
-> one report
```

Reuse is allowed only at explicit stage boundaries:
- `reuse.called_vcf`
- `reuse.annotated_vcf`
- `reuse.prs_dense_vcf`

Those imports do not create alternate execution graphs. They only skip the stage that would normally generate the same artifact.

## Repo Layout

- `main.nf`: the pipeline graph
- `nextflow.config`: Nextflow runtime settings
- `run.sh`: the only user-facing shell entrypoint
- `AGENTS.md`: shared repo instructions for Codex/agents
- `PLANS.md`: shared execution roadmap for larger work
- `manifests/`: tracked examples plus schema; private manifests live under ignored `manifests/local/`
- `workflow/core/`: manifest, runtime, and reference/context logic
- `workflow/stages/`: monogenic and PRS stage implementations
- `workflow/reporting/`: report rendering and outcome enforcement
- `tests/`: unit tests
- `output/results/`: final run outputs only

The heavy runtime/cache state lives under `runtime.linux_cache_root`, not in the repo.

`AGENTS.md` and `PLANS.md` are meant to be tracked. They document shared project guidance and the current narrow roadmap. Personal manifests, private presets, raw data, and generated outputs should remain untracked.

## Manifest

Example tracked manifest:
- `manifests/examples/sample.yaml`
- `manifests/examples/sample.monogenic.yaml`
- `manifests/examples/sample.prs.yaml`

Main fields:
- `sample_id`
- `inputs.cram`
- `inputs.crai`
- `reference.build`
- `reference.bundle_id`
- `analysis.monogenic`
- `analysis.prs`
- `analysis.prs.sample_sex`
- `runtime.linux_cache_root`
- `outputs.results_dir`

PRS score lists can come from either:
- `analysis.prs.score_preset`
- `analysis.prs.pgs_ids`

Execution mode is derived from the enabled analysis branches:
- `analysis.monogenic.enabled: true` and `analysis.prs.enabled: true` -> `both`
- `analysis.monogenic.enabled: true` and `analysis.prs.enabled: false` -> `monogenic-only`
- `analysis.monogenic.enabled: false` and `analysis.prs.enabled: true` -> `prs-only`

The pipeline intentionally rejects manifests where both branches are disabled.

If `analysis.prs.sample_sex` is set to `male` or `female`, the PRS compatibility step keeps sex chromosomes and writes a one-sample `.psam` for downstream tools. If it is omitted, PRS compatibility falls back to autosomes only.

Tracked example presets live in:
- `manifests/examples/score_presets.yaml`

Private sample manifests, private reference bundle manifests, and private score presets belong under ignored `manifests/local/`.

## Reference Bundles

CRAM-based stages require the exact reference used for alignment.

The pipeline now handles that in `PREPARE_CONTEXT` by:
- probing the CRAM header
- extracting contig MD5s and provenance hints
- validating a local FASTA/FAI/dict bundle against the CRAM header
- staging heavy reference files into Linux-native cache when needed

Optional remote MD5 decode probing is available for diagnostics, but it does not replace a validated local bundle.

Tracked example:
- `manifests/examples/reference_bundles/grch38_dragen_alt_masked.yaml`

## Outputs

Final reader-facing result:
- `output/results/<sample>/<run_id>/report/report.md`

The report title and sections adapt to the execution mode:
- `Genome Report` for `both`
- `Monogenic Report` for `monogenic-only`
- `PRS Report` for `prs-only`

Technical manifests:
- `output/results/<sample>/<run_id>/report/artifact_manifest.json`

Heavy cache/work state:
- `<runtime.linux_cache_root>/<project-cache-id>/...`

## Performance Model

The pipeline automatically stages heavy inputs from `/mnt/c/...` or Windows paths into Linux-native cache, then runs heavy work there.

Performance-sensitive stages are structured as scatter/gather:
- monogenic calling is sharded by chromosome plus one `alt` shard
- dense PRS genotyping is sharded by chromosome plus one `alt` shard

That avoids the earlier design where PRS walked chromosomes in one long serial shell loop and monogenic ran as one giant whole-genome caller task.

## Development

```bash
make validate
make test
make clean
make clean-cache
```

`make clean` removes repo-facing results only.  
`make clean-cache` removes Linux-side cache/work state for this project.

## License

MIT. See [LICENSE](LICENSE).
