# Genomics

Single-sample, research-use genome analysis from `CRAM/CRAI` to:
- monogenic triage
- targeted PRS
- one final Markdown report

## Scope

Supported:
- one sample at a time
- SNV/indel monogenic triage
- targeted PRS from CRAM-derived genotype sites
- reuse of trusted existing VCF or annotated VCF artifacts when declared in the manifest

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

Advanced direct command:

```bash
nextflow run . -params-file manifests/local/<sample>.yaml -profile docker
```

The wrapper keeps Nextflow runtime state under `output/` so the repo root stays clean.

## Repo Layout

- `main.nf`: the pipeline graph
- `nextflow.config`: Nextflow runtime settings
- `run.sh`: thin launcher that keeps Nextflow state out of the repo root
- `manifests/`: the only supported input/config surface
- `workflow/`: implementation code for core, monogenic, and PRS stages
- `tests/`: unit tests for the Python parts
- `output/`: generated state and run outputs only

`manifests/` is intentionally separate from `output/`: it is tracked configuration, not generated data.

`manifests/` contains:
- `manifests/examples/sample.yaml`
- `manifests/run_manifest.schema.json`
- `manifests/examples/score_presets.yaml`
- `manifests/reference_bundles/*.example.yaml`
- `manifests/local/` for private sample manifests and private reference bundle manifests

`manifests/run_manifest.schema.json` is the machine-checked contract for a valid run manifest. It exists so the pipeline fails early on bad configuration instead of failing deep into a long run.

`output/` is intentionally split by purpose:
- `output/results/`: reader-facing outputs for a specific run
- `output/cache/`: reusable expensive intermediates keyed by the input fingerprint
- `output/work/`: Nextflow task working directories
- `output/logs/`: launcher and Nextflow logs
- `output/nextflow_state/`: Nextflow runtime state
- `output/nextflow_launcher/`: launch directory used to keep root-level Nextflow files out of the repo

## Reference Bundles

CRAM-based stages require the exact reference used for alignment. The pipeline now separates:
- `reference_probe`: inspect the CRAM header, contig MD5s, and aligner provenance hints
- `validate_reference_bundle`: verify that a local FASTA/FAI/dict bundle matches the CRAM header

Optional remote MD5 decode probing can be enabled with `reference.remote_decode_probe: true`, but that is only a debugging aid. It does not replace a validated local bundle.

Tracked example:
- [`grch38_dragen_alt_masked.example.yaml`](manifests/reference_bundles/grch38_dragen_alt_masked.example.yaml)

## Outputs

Reader-facing result:
- `output/results/<sample>/<run_id>/report/report.md`

Technical artifacts:
- `output/results/<sample>/<run_id>/manifests/`
- `output/results/<sample>/<run_id>/monogenic/`
- `output/results/<sample>/<run_id>/prs/`
- `output/cache/samples/<input_fingerprint>/...`

## Development

```bash
make validate
make test
make clean
```

`make clean` removes generated runtime/output directories under `output/`.

## Note On PRS

This repo does not score a sparse hard-filtered WGS VCF directly. It builds the union of requested PGS loci, genotypes those loci from the CRAM, and then runs `pgsc_calc`.

## License

MIT. See [LICENSE](LICENSE).
