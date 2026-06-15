# Multi-omics Integration for Sorghum Metabolic Variation Analysis

Genomics (variant calling), transcriptomics (gene co-expression), and epigenomics (methylation) of four *Sorghum bicolor* samples with contrasting TAA production phenotypes.

## Samples

| Library | Sample | TAA conc. in juice | Juice Production | Callus Formation |
|---------|--------|---------------|---------------|-----------------|
| r0074 | SBC4 | High | ++ | Mid |
| r0066 | SBC10 | Low | +++ | Good |
| r0075, r0078, r0078-2 | SBC11 | High | - | Mid |
| r0076 | SBC23 | High | ++ | Good |

Reference genome: `GCF_000003195.3` (BTx623, NCBIv3)
Sequencing: ONT R10.4.1 duplex, 400 bps SUP (Dorado v1.3.0)

## Repository Structure

```
thesis/
├── workflow/           # Snakemake workflows and scripts
│   ├── Snakefile
│   ├── rules/
│   │   ├── reads_preprocessing.smk   # alignment, indexing, depth
│   │   ├── variant_analysis.smk      # Clair3 variant calling
│   │   ├── vcf_processing.smk        # filtering, phasing
│   │   ├── snpeff_annotation.smk     # SnpEff annotation, private variants
│   │   ├── sift_annotation.smk       # SIFT4G functional annotation
│   │   ├── methylation.smk           # modkit 5mC pileup
│   │   ├── dmr_analysis.smk          # DSS DMR calling (TAA contrast)
│   │   └── ranked_genes.smk          # multi-omics gene ranking
│   └── scripts/        # helper scripts used by workflows and analysis
├── analysis/           # downstream analysis notebooks and figures
├── resources/          # input data and intermediate files
├── results/            # final output files
├── docker/             # Dockerfile and helper scripts
└── dag/                # Snakemake DAG visualizations
```

## Quick Start

All bioinformatics tools run inside a Docker image. See [DATA.md](DATA.md) for the full data generation pipeline.

```shell
# Build the Docker image
./docker/build.sh

# Run a Snakemake workflow target
./docker/run.sh snakemake <target> --cores 16

# Dry-run
./docker/run.sh snakemake <target> -n

# Run an ad-hoc command
./docker/run.sh samtools view -H resources/align_bam_sample/SBC10.bam
```

## Pipeline Targets

| Step | Target | Rule file |
|------|--------|-----------|
| 1 — Read preprocessing | `reads_all` | `reads_preprocessing.smk` |
| 2 — Variant calling | `variants_all` | `variant_analysis.smk` |
| 3 — VCF filtering & phasing | `vcf_all` | `vcf_processing.smk` |
| 4 — SnpEff annotation | `annotate_vcf` | `snpeff_annotation.smk` |
| 5 — SIFT4G annotation | `annotate_sift` | `sift_annotation.smk` |
| 6 — Methylation calling | `methylation_all` | `methylation.smk` |
| 7 — DMR analysis (TAA) | `annotate_dmr` | `dmr_analysis.smk` |
| 8 — Multi-omics gene ranking | `ranked_genes` | `ranked_genes.smk` |
