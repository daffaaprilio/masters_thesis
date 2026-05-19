# CLAUDE.md — Sorghum Multi-omics Thesis

## Project Overview

Multi-omics analysis of four *Sorghum bicolor* samples (BTx623 reference, `GCF_000003195.3`) with contrasting TAA (trans-aconitic acid) production phenotypes. Other phenotypes, i.e., Sorgoleone production, are also included to the story. Sequenced with ONT R10.4.1 duplex 400 bps SUP (Dorado v1.3.0).

Samples:

| Sample | Libraries | TAA Production | TAA Secretion | Callus Formation |
|--------|-----------|---------------|---------------|--------|
| SBC4 | r0074 | ++ | High | Mid |
| SBC10 | r0066 | +++ | Low | Good |
| SBC11 | r0075, r0078, r0078-2 | − | High | Mid |
| SBC23 | r0076 | ++ | High | Good |

SBC11 is special: its three libraries must be merged with `samtools merge` before use (see DATA.md).

## Running Commands

**All bioinformatics tools run inside Docker.** Never run pipeline tools directly on the host.

```shell
# Build the image (once)
./docker/build.sh

# Run a Snakemake workflow target
./docker/snakemake.sh <target> --cores <N>

# Dry-run
./docker/snakemake.sh <target> -n

# Run an ad-hoc command
./docker/run.sh <command>
```

Override defaults via env vars:
```shell
THESIS_IMAGE=thesis-tools:latest CORES=24 SNPEFF_DIR=/path/to/snpEff ./docker/snakemake.sh
```

## Pipeline Targets (in order)

| Step | Target | Snakefile rule file |
|------|--------|---------------------|
| 1 — Read preprocessing | `reads_all` | `workflow/rules/reads_preprocessing.smk` |
| 2 — Variant calling | `variants_all` | `workflow/rules/variant_analysis.smk` |
| 3 — VCF postprocessing | `vcf_all` | `workflow/rules/vcf_processing.smk` |
| 4 — Methylation calling | `methylation_all` | `workflow/rules/methylation.smk` |

Key outputs:
- `results/vcf_processing/{sample}.phased.vcf.gz` — phased, annotated VCFs
- `resources/bedmethyl/{sample}.bed` — raw 5mC pileup
- `resources/bedmethyl/{sample}.filtered.bed` — positions with ≥ 10 valid reads

## Repository Layout

```
workflow/           # Snakemake workflows
  rules/            # .smk rule files (one per pipeline step)
  scripts/          # helper scripts called by rules
analysis/           # downstream analysis
  00_data_quality/  # read depth QC, VCF benchmarks
  01_variant_landscape/ # SnpEff summary figures (fig01–fig10)
  02_methylation_landscape/ # methylation figures
  04_sorgoleone/    # sorgoleone biosynthetic pathway variant analysis
  scripts/          # standalone Python/shell scripts for analysis
  data/tsv/         # VCF converted to TSV for notebooks
resources/          # input data and intermediates
  align_bam/        # per-library BAMs (r0066, r0074, …)
  align_bam_sample/ # per-sample BAMs (SBC4, SBC10, SBC11, SBC23)
  bedmethyl/        # modkit 5mC pileup outputs
  depth/            # per-sample depth files and plots
  vcf/              # symlinks to Clair3 VCF outputs
results/            # final pipeline outputs
  variant_calling/  # Clair3 output directories
  vcf_processing/   # filtered, phased, annotated VCFs
docker/             # Dockerfile, environment.yml, helper shell scripts
dag/                # Snakemake DAG PDFs
```

## Toolchain (inside Docker image `thesis-tools:latest`)

minimap2 2.30, samtools 1.21, bcftools 1.21, htslib 1.21, bedtools 2.31.1, whatshap 2.8, modkit 0.2.6, snpEff (bioconda), clair3 (models at `/opt/models/`), snakemake 9.16.3, Python 3.11 (pandas, matplotlib, seaborn, cyvcf2, matplotlib-venn).

## Key Gotchas

- **SBC11 BAM must be manually merged** from three libraries before the variant calling and methylation steps. See DATA.md Step 1.
- **SnpEff chromosome synonyms**: VCF contig IDs (e.g. `NC_012870.2`) differ from SnpEff names (`1`, `2`, …). The synonym file at `workflow/scripts/creating_synonyms_chr.py` bridges this.
- **Methylation requires MM/ML tags**: BAMs must be basecalled with Dorado `--modified-bases`. Verify with `modkit summary`; if 0 modified bases reported, re-basecall.
- **SnpEff database**: must be downloaded once with `./docker/run.sh snpEff download Sorghum_bicolor` before running `vcf_all`.

## Analysis Scripts

Run all analysis scripts via `./docker/run.sh python3 analysis/scripts/<script>.py`.

| Script | Purpose |
|--------|---------|
| `variant_landscape.py` | Generates fig01–fig10 from SnpEff stats CSVs |
| `methylation_landscape.py` | Methylation landscape plots |
| `vcf_benchmark.py` | VCF benchmark figure |
| `annot_vcf_to_tsv.py` | Converts annotated VCF → TSV for notebooks |
| `private_variants.sh` | `bcftools isec` to find sample-private variants |
| `merge_vcf.sh` | Merges per-sample phased VCFs into multi-sample VCF |

Jupyter notebooks in `analysis/04_sorgoleone/` explore sorgoleone pathway variants.

## Reference Docs

- `DATA.md` — full data generation pipeline with exact commands
- `analysis/ANALYSIS.md` — downstream analysis overview
