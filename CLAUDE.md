# CLAUDE.md — Sorghum Multi-omics Thesis

## Project Overview

Multi-omics analysis of four *Sorghum bicolor* samples (BTx623 reference, `GCF_000003195.3`) with contrasting TAA (trans-aconitic acid) production phenotypes. Other phenotypes, i.e., Sorgoleone production, are also included to the story. Sequenced with ONT R10.4.1 duplex 400 bps SUP (Dorado v1.3.0).

Samples:

| Library | Sample | TAA conc. in juice | Juice Production | Callus Formation |
|---------|--------|---------------|---------------|-----------------|
| r0074 | SBC4 | High | ++ | Mid |
| r0066 | SBC10 | Low | +++ | Good |
| r0075, r0078, r0078-2 | SBC11 | High | - | Mid |
| r0076 | SBC23 | High | ++ | Good |

SBC11 is special: its three libraries must be merged with `samtools merge` before use (see DATA.md).

## Running Commands

**All bioinformatics tools run inside Docker.** Never run pipeline tools directly on the host.

```shell
# Build the image (once)
./docker/build.sh

# Run a Snakemake workflow target
./docker/run.sh snakemake <target> --cores <N>

# Dry-run
./docker/run.sh snakemake <target> -n

# Run an ad-hoc command
./docker/run.sh <command>
```

## Pipeline Targets (in order)

| Step | Target | Snakefile rule file |
|------|--------|---------------------|
| 1 — Read preprocessing | `reads_all` | `workflow/rules/reads_preprocessing.smk` |
| 2 — Variant calling | `variants_all` | `workflow/rules/variant_analysis.smk` |
| 3 — VCF filtering & phasing | `vcf_all` | `workflow/rules/vcf_processing.smk` |
| 4 — SnpEff annotation | `annotate_vcf` | `workflow/rules/snpeff_annotation.smk` |
| 5 — SIFT4G annotation | `annotate_sift` | `workflow/rules/sift_annotation.smk` |
| 6 — Methylation calling | `methylation_all` | `workflow/rules/methylation.smk` |
| 7 — DMR analysis (TAA) | `annotate_dmr` | `workflow/rules/dmr_analysis.smk` |
| 8 — Multi-omics gene ranking | `ranked_genes` | `workflow/rules/ranked_genes.smk` |

Key outputs:
- `results/vcf_processing/{sample}.phased.vcf.gz` — phased VCFs
- `results/snpeff/{sample}.private.snpeff.vcf.gz` — SnpEff-annotated private variants
- `results/sift4g/{sample}.private.sift4g.vcf.gz` — SIFT4G-annotated private variants
- `resources/bedmethyl/{sample}.bed` — raw 5mC pileup
- `resources/bedmethyl/{sample}.filtered.bed` — positions with ≥ 10 valid reads
- `results/DMR/{pair}.5mC.DMR.tsv` — annotated DMRs per TAA contrast pair
- `results/ranked_genes_lists/SBC10.multiomics_ranked.tsv` — final multi-omics gene ranking

## Repository Layout

```
workflow/           # Snakemake workflows
  Snakefile         # main entry point; includes all rule files
  rules/            # .smk rule files (one per pipeline step)
  scripts/          # all helper scripts (used by rules and analysis)
analysis/           # downstream analysis notebooks and figures
  figures/          # output figures
resources/          # input data and intermediates
  align_bam/        # per-library BAMs (r0066, r0074, …)
  align_bam_sample/ # per-sample BAMs (SBC4, SBC10, SBC11, SBC23)
  bedmethyl/        # modkit 5mC pileup outputs
  depth/            # per-sample depth files and plots
  vcf/              # symlinks to Clair3 VCF outputs
results/            # final pipeline outputs
  variant_calling/  # Clair3 output directories
  vcf_processing/   # filtered, phased VCFs
  snpeff/           # SnpEff-annotated VCFs
  sift4g/           # SIFT4G-annotated VCFs
  private_variants/ # bcftools isec outputs
  DSS/              # per-pair DSS input files
  DMR/              # DSS DMR calls and annotations
  ranked_genes_lists/ # multi-omics gene ranking outputs
docker/             # Dockerfile, environment.yml, helper shell scripts
dag/                # Snakemake DAG PDFs
```

## Toolchain (inside Docker image `thesis-tools:latest`)

minimap2 2.30, samtools 1.21, bcftools 1.21, htslib 1.21, bedtools 2.31.1, whatshap 2.8, modkit 0.2.6, snpEff (bioconda), SIFT4G, DSS (R/Bioconductor), clair3 (models at `/opt/models/`), snakemake 9.16.3, Python 3.11 (pandas, matplotlib, seaborn, cyvcf2, matplotlib-venn).

## Key Gotchas

- **SBC11 BAM must be manually merged** from three libraries before the variant calling and methylation steps. See DATA.md Step 1.
- **SnpEff chromosome synonyms**: VCF contig IDs (e.g. `NC_012870.2`) differ from SnpEff names (`1`, `2`, …). The synonym file at `workflow/scripts/creating_synonyms_chr.py` bridges this.
- **Methylation requires MM/ML tags**: BAMs must be basecalled with Dorado `--modified-bases`. Verify with `modkit summary`; if 0 modified bases reported, re-basecall.
- **SnpEff database**: must be downloaded once with `./docker/run.sh snpEff download Sorghum_bicolor` before running `annotate_vcf`.

## Analysis Scripts

All scripts live in `workflow/scripts/`. Run via `./docker/run.sh python3 workflow/scripts/<script>.py`.

| Script | Purpose |
|--------|---------|
| `variant_landscape.py` | Generates variant landscape figures from SnpEff stats CSVs |
| `methylation_landscape.py` | Methylation landscape plots |
| `vcf_benchmark.py` | VCF benchmark figure |
| `annot_vcf_to_tsv.py` | Converts annotated VCF → TSV for notebooks |
| `private_variants.sh` | `bcftools isec` to find sample-private variants |
| `merge_vcf.sh` | Merges per-sample phased VCFs into multi-sample VCF |
| `rank_dmr_genes.py` | Ranks genes by DMR proximity for multi-omics integration |
| `build_ranked_genes.py` | Builds final multi-omics ranked gene list |

## Reference Docs

- `DATA.md` — full data generation pipeline with exact commands
