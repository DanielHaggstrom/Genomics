# Genomics

Single-sample research analysis of WGS, centered on calculating polygenic scores
without losing completed work when a long run stops.

`CRAM/CRAI → reference validation → scaffold genotypes → optional imputation → PRS scoring → one report`

PRS is enabled by default. Monogenic SNV/indel triage and PRS-only Beagle imputation
are optional. There are no clinical diagnosis, SV/CNV, or pedigree workflows.
Implementation and synthetic checks are in place; actual caller/scorer integration
and results from personal data remain unvalidated. Do not interpret raw PGS sums
as absolute disease probabilities.

## Configure

Keep original downloads outside the repository. Their directory is protected from
runtime/output writes automatically. Add any larger original-data tree to
`inputs.protected_roots`; the pipeline treats its contents as read-only, including
files whose provenance is uncertain. Missing reference/import indexes are generated
on cache copies. Generated state belongs in `~/.cache/genomics`, separate from originals.

Copy [the example manifest](manifests/examples/sample.yaml) into ignored
`manifests/local/`, and set the CRAM/CRAI paths, real PGS Catalog score IDs
(`analysis.prs.pgs_ids`), resource budget, and exact reference bundle. The demo
preset contains a placeholder ID and is for manifest validation only.

Copy [the reference example](manifests/examples/reference_bundles/grch38_dragen_alt_masked.yaml)
into `manifests/local/reference_bundles/` and supply the original reference FASTA,
FAI, and dictionary paths. Set `reference.bundle_manifest` relative to your run
manifest. A similarly named GRCh38 reference is insufficient: the dictionary must
match the CRAM contig lengths and MD5s. Reference downloading/remote probing is
intentionally absent; an explicit validated local bundle is required for calling.
Remove retired `reference.autodownload` and `reference.remote_decode_probe` keys
from older manifests.

Python dependencies are in [requirements.txt](requirements.txt). Runtime tools are
Docker, samtools, bcftools, and tabix; PRS also needs Java, Nextflow, and curl.
Building the broad scaffold or ancestry calibration needs tar and zstd. Optional
VEP annotation needs curl and tar. Preparation checks tools before staging large
inputs. Linux path/tool migration beyond this launcher remains a separate task;
legacy WSL staging and Docker compatibility helpers are retained for now.

## Run, stop, resume

When ready to authorize a real run:

```bash
bash run.sh manifests/local/sample.yaml
```

Use Ctrl-C or SIGTERM to stop. Run the same command to resume; there is no separate
resume mode (`-resume` is accepted for compatibility). Active child process groups
and directly launched Docker containers are stopped. Completed chromosome/task
checkpoints remain; an interrupted caller task restarts from its last completed
boundary. This is task-level resume, not continuation inside GATK or Beagle.
Power loss/SIGKILL cannot perform cleanup, so verify orphaned tools have stopped
before restarting after an unclean shutdown.

Checkpoints are written atomically after successful commands and fingerprint
nonempty outputs plus required indexes. An interrupted output without a valid
receipt is rebuilt. Legacy outputs lacking receipts may be recomputed once.
Downloads retain `.part` files and request HTTP range continuation on restart.
The server must support range requests. Completed whole-contig gVCFs are shared
across score selections and the optional monogenic branch. Changing data or analysis
settings selects the dependent cache paths; completed work is retained.

One local Python driver runs stages sequentially. A cache lock rejects a second
writer. The pinned upstream `pgsc_calc` v2.3.0 retains its own Nextflow work and
always uses `-resume`; CPU and memory budgets apply to its local executor too.
Scorer containers carry a stable cache-specific label for bounded cleanup after
interruption, using [Nextflow Docker options](https://docs.seqera.io/nextflow/reference/config/docker).
Changing resource settings takes effect on restart. Concurrent analyses against
different caches must also use different result directories/run IDs.

```bash
make status RUN=manifests/local/sample.yaml
```

Status reads saved run/stage metadata without scanning genomic data. Logs stream
to cache files while tools run. After an unclean shutdown, saved status can still
say `running`; it records the last update, not a process liveness monitor.

## Genotypes and results

Direct PRS calling uses a broad common-variant backbone plus score-specific sites,
finalizes gVCFs at those sites, and reconstructs explicit alleles for confident
reference SNP calls. It preserves depth and GQ/RGQ; defaults require DP ≥ 10 and
GQ/RGQ ≥ 20. Missing/uncertain observations are excluded from the scoring copy.
Indel absence is not inferred from an ordinary reference SNP call. Direct calling
and local imputation currently support diploid autosomes; validated imported
sex-chromosome calls require explicit sample sex.

Enable `analysis.prs.imputation.enabled` only after supplying a compatible phased
panel using [the panel example](manifests/examples/imputation_panels/grch38_beagle_panel.yaml).
Imputation augments the scaffold, retaining confident direct observations and
sites absent from the panel. Inspect per-score matching and missing high-weight
variants; the overlap cutoff is a failure gate, not proof of score validity.
Ancestry/reference calibration is separately enabled with `run_ancestry`.

Explicit stage-boundary imports are `reuse.called_vcf`, `reuse.annotated_vcf`,
`reuse.prs_dense_vcf` (scaffold), and `reuse.imputed_prs_vcf`. Avoid importing a
sparse variants-only VCF as if it contained confident reference observations.

The reader-facing report is
`output/results/<sample_id>/<run_id>/report/report.md`, with score/QC tables and
one `artifact_manifest.json` containing complete provenance and canonical cache
paths. Large genotype VCFs stay in the cache. An enabled failed/incomplete stage
makes the run fail, while available failure details remain in the report.

## Maintenance

```bash
make test
make validate
make validate-manifest RUN=manifests/examples/sample.yaml
```

Tests use synthetic data, mocked biological tools, and short Python subprocesses
for interruption checks. They do not run the WGS pipeline. The next acceptance
gate is a bounded synthetic or truth-known-region test of the real external tools,
then a manageable score with independently checked genotypes and arithmetic.

`make clean` removes default published results. `make clean-cache RUN=...` removes
this project's cache and **discards resumable work**; `make clean-all` does both.
Neither is a normal stop/resume step. Private manifests, original data, personal
score presets, results, and runtime state must stay out of Git.
