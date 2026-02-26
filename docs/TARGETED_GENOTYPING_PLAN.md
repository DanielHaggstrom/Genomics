# Targeted Genotyping Plan (CRAM -> PRS-Ready VCF)

## Why this plan exists

PRS overlap can be low when using sparse, variant-only VCF files.
`pgsc_calc` may require extra input preparation for WGS workflows.

## Goal

Create a PRS-focused genotype file from `CRAM` by force-calling loci needed by selected PGS IDs, then rerun `pgsc_calc`.

## Scope and safety

1. Never modify original input files in your source WGS directory.
2. Write all derived files under this repo (`results/prs/targeted_genotypes/`, `logs/`).
3. Keep monogenic analysis unchanged.

## High-level workflow

1. Build and harmonization checks.
2. Build a union target-site list from selected PGS score files (GRCh38).
3. Force-call/genotype these target sites from CRAM.
4. Normalize/index/filter to a PRS-ready VCF.
5. Rerun `pgsc_calc` on targeted VCF.
6. Compare overlap and pass rates versus current baseline.

## Detailed steps

## 1) Verify reference/build consistency

1. Confirm CRAM and VCF are GRCh38.
2. Confirm PGS scorefiles are harmonized to GRCh38.
3. Confirm chromosome naming compatibility (`chr1` vs `1`).
4. Confirm `sampleset` uses safe characters.

Suggested commands:

```bash
bcftools view -h /path/to/sample.hard-filtered.vcf.gz \
  | egrep -i "reference|contig|assembly" | head -n 40

bcftools index -n /path/to/sample.hard-filtered.vcf.gz
```

## 2) Build target-site list from selected PGS IDs

1. Download selected scores in GRCh38.
2. Extract variants into a union table with at least:
   - `CHROM`, `POS`, `REF`, `ALT`
3. Deduplicate to unique loci.

Output:
- `results/prs/targeted_genotypes/targets/pgs_union_sites.tsv.gz`

## 3) Force-call only target loci from CRAM

Template command:

```bash
bcftools mpileup \
  -f <GRCh38.fa> \
  -R results/prs/targeted_genotypes/targets/pgs_union_sites.tsv.gz \
  -a FORMAT/DP,FORMAT/AD \
  /path/to/sample.cram -Ou \
| bcftools call \
  -m \
  -C alleles \
  -T results/prs/targeted_genotypes/targets/pgs_union_sites.tsv.gz \
  -Oz -o results/prs/targeted_genotypes/sample.pgs_targeted.raw.vcf.gz
```

## 4) Normalize and index

1. Left-normalize against GRCh38 reference.
2. Split multiallelics if needed.
3. Index output.

Outputs:
- `results/prs/targeted_genotypes/sample.pgs_targeted.vcf.gz`
- `results/prs/targeted_genotypes/sample.pgs_targeted.vcf.gz.tbi`

## 5) Run `pgsc_calc` on targeted VCF

1. Point `config/prs_samplesheet.csv` to targeted VCF prefix.
2. Keep Docker profile where available.
3. Start with the same baseline `PGSC_MIN_OVERLAP`, then tune.

## 6) Acceptance criteria

1. More scores pass `min_overlap` than baseline.
2. Median match rate increases materially.
3. No increase in clear allele/build mismatch warnings.
4. Pipeline reruns reproducibly.

## References

1. https://pgsc-calc.readthedocs.io/en/latest/how-to/prepare.html
2. https://pgsc-calc.readthedocs.io/en/latest/explanation/match.html
3. https://pgsc-calc.readthedocs.io/en/latest/how-to/samplesheet.html
4. https://pgsc-calc.readthedocs.io/en/latest/reference/input.html
5. https://pgsc-calc.readthedocs.io/en/latest/explanation/output.html
6. https://samtools.github.io/bcftools/bcftools.html