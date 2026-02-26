# Targeted Genotyping Plan (CRAM -> PRS-Ready VCF)

## Why this plan exists

Our current PRS runs show low overlap for several scores. This is expected when using a sparse, variant-only VCF for PRS-style matching.

Official `pgsc_calc` docs state that WGS data is not currently a fully supported input mode and that creating a compatible genotype representation may require extra work.

## Goal

Create a PRS-focused genotype file from `CRAM` by force-calling only loci needed by selected PGS IDs, then rerun `pgsc_calc` on this denser PRS-targeted VCF.

## Scope and safety

1. Never modify original input files in `C:/Users/<YOUR_USER>/Salud/WGS`.
2. Write all derived files under this repo (`results/prs/targeted_genotypes/`, `logs/`).
3. Keep monogenic analysis on current VEP branch unchanged.

## High-level workflow

1. Build and harmonization checks.
2. Build a union target-site list from chosen PGS score files (GRCh38).
3. Force-call/genotype these target sites from CRAM.
4. Normalize/index/filter to a PRS-ready VCF.
5. Rerun `pgsc_calc` on targeted VCF.
6. Compare overlap and score pass rates versus current run.

## Detailed steps

## 1) Verify reference/build consistency

1. Confirm CRAM and VCF are GRCh38.
2. Confirm PGS scorefiles used are harmonized to GRCh38.
3. Confirm chromosome naming compatibility (`chr1` vs `1`).
4. Confirm `sampleset` uses only safe characters (no `.` or `_`).

Record outputs in `logs/targeted_genotyping_precheck.txt`.

Suggested commands:

```bash
bcftools view -h /mnt/c/Users/<YOUR_USER>/Salud/WGS/DEMO_SAMPLE_001.hard-filtered.vcf.gz \
  | egrep -i "reference|contig|assembly" | head -n 40

bcftools index -n /mnt/c/Users/<YOUR_USER>/Salud/WGS/DEMO_SAMPLE_001.hard-filtered.vcf.gz
```

## 2) Build target-site list from selected PGS IDs

1. Download selected scores in GRCh38.
2. Extract variants into a union table with at least:
   - `CHROM`, `POS`, `REF`, `ALT`
3. Deduplicate to unique target loci.

Output:
- `results/prs/targeted_genotypes/targets/pgs_union_sites.tsv.gz`

## 3) Force-call only target loci from CRAM

Use `bcftools mpileup` + `bcftools call` with constrained calling at supplied alleles/sites.

Template command (to be finalized against local reference paths):

```bash
bcftools mpileup \
  -f <GRCh38.fa> \
  -R results/prs/targeted_genotypes/targets/pgs_union_sites.tsv.gz \
  -a FORMAT/DP,FORMAT/AD \
  /mnt/c/Users/<YOUR_USER>/Salud/WGS/DEMO_SAMPLE_001.cram -Ou \
| bcftools call \
  -m \
  -C alleles \
  -T results/prs/targeted_genotypes/targets/pgs_union_sites.tsv.gz \
  -Oz -o results/prs/targeted_genotypes/DEMO_SAMPLE_001.pgs_targeted.raw.vcf.gz
```

Notes:
1. `-C alleles` performs constrained calling at provided alleles.
2. Site list format must match `bcftools call -C alleles` expectations.

## 4) Normalize and index

1. Left-normalize against GRCh38 reference.
2. Split multiallelics if needed.
3. Index output.

Outputs:
- `results/prs/targeted_genotypes/DEMO_SAMPLE_001.pgs_targeted.vcf.gz`
- `results/prs/targeted_genotypes/DEMO_SAMPLE_001.pgs_targeted.vcf.gz.tbi`

## 5) Run `pgsc_calc` on targeted VCF

1. Point `config/prs_samplesheet.csv` to targeted VCF prefix.
2. Keep Docker profile.
3. Start with same batch strategy and same threshold used in successful run (`PGSC_MIN_OVERLAP=0.50`), then evaluate if threshold can be raised.

## 6) Acceptance criteria

1. More scores pass `min_overlap` than with hard-filtered VCF.
2. Median match rate increases materially.
3. No increase in obvious allele/build mismatch warnings.
4. Pipeline completes reproducibly across reruns.

## Risks and mitigations

1. Wrong reference FASTA/build:
   - Mitigation: strict GRCh38 validation before calling.
2. Bad targets file schema:
   - Mitigation: test on small subset first.
3. Overinterpreting low-overlap scores:
   - Mitigation: keep match-rate reporting and quality flags in final summary.

## Official references used

1. PGS Catalog Calculator docs: input preparation (WGS note, imputation preference)
   - https://pgsc-calc.readthedocs.io/en/latest/how-to/prepare.html
2. PGS Catalog Calculator docs: match-rate errors and `--min_overlap` guidance
   - https://pgsc-calc.readthedocs.io/en/latest/explanation/match.html
3. PGS Catalog Calculator docs: samplesheet constraints (reserved characters)
   - https://pgsc-calc.readthedocs.io/en/latest/how-to/samplesheet.html
   - https://pgsc-calc.readthedocs.io/en/latest/reference/input.html
4. PGS Catalog Calculator docs: outputs layout (`match`, `score`, `pipeline_info`)
   - https://pgsc-calc.readthedocs.io/en/latest/explanation/output.html
5. bcftools manual: targets/regions and constrained calling with `-C alleles`
   - https://samtools.github.io/bcftools/bcftools.html

## Next execution steps (project)

1. Freeze current PRS outputs as baseline.
2. Implement a scripted target-site extractor from selected PGS IDs.
3. Run a pilot force-call on 1 batch worth of sites.
4. Compare overlap metrics and decide full rollout.
