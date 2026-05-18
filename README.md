# Multi-omics Integration for Sorghum Metabolic Variation Analysis

Genomics (variant calling), transcriptomics (gene co-expression), and epigenomics (methylation) of four *Sorghum bicolor* samples with contrasting TAA production phenotypes.

## Samples

| Library | Sample | TAA Production | TAA Secretion | Callus Formation |
|---------|--------|---------------|---------------|-----------------|
| r0074 | SBC4 | ++ | High | Mid |
| r0066 | SBC10 | +++ | Low | Good |
| r0075, r0078, r0078-2 | SBC11 | - | High | Mid |
| r0076 | SBC23 | ++ | High | Good |

Reference genome: `GCF_000003195.3` (BTx623, NCBIv3)
Sequencing: ONT R10.4.1 duplex, 400 bps SUP (Dorado v1.3.0)

## Repository Structure

```
thesis/
├── workflow/           # Snakemake workflows for data generation
│   ├── rules/
│   │   ├── reads_preprocessing.smk   # alignment, indexing, depth
│   │   ├── variant_analysis.smk      # Clair3 variant calling
│   │   ├── vcf_processing.smk        # filtering, phasing, annotation
│   │   └── methylation.smk           # modkit 5mC pileup
│   └── scripts/        # helper scripts used by workflows
├── analysis/           # downstream analysis notebooks and scripts
├── resources/          # input data and intermediate files
├── results/            # final output files
├── docker/             # Dockerfile and helper scripts
└── dag/                # Snakemake DAG visualizations
```

## Quick Start

All bioinformatics tools run inside a Docker image. See [DATA.md](DATA.md) for the full data generation pipeline and [ANALYSIS.md](analysis/ANALYSIS.md) for downstream analysis.

```shell
# Build the Docker image
./docker/build.sh

# Run a Snakemake workflow
./docker/snakemake.sh reads_all --cores 16

# Run an ad-hoc command
./docker/run.sh samtools view -H resources/align_bam_sample/SBC10.bam
```
