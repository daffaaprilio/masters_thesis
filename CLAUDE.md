# CLAUDE.md — Sorghum Multi-omics Thesis

## Project Overview

Multi-omics analysis of four *Sorghum bicolor* samples (BTx623 reference, `GCF_000003195.3`) with contrasting TAA (trans-aconitic acid) production phenotypes. Another phenotype, stem juiciness, that is regulated by D-gene transcription factor, is also included to the story, mainly as the positive control. Samples are sequenced with ONT R10.4.1 duplex 400 bps SUP (Dorado v1.3.0). For gene co-expression analysis, use data from ATTED-II database.

Samples:

| Library | Sample | TAA conc. in juice | Juice Production | Callus Formation |
|---------|--------|---------------|---------------|-----------------|
| r0074 | SBC4 | High | ++ | Mid |
| r0066 | SBC10 | Low | +++ | Good |
| r0075, r0078, r0078-2 | SBC11 | High | - | Mid |
| r0076 | SBC23 | High | ++ | Good |

SBC11 is special: its three libraries must be merged with `samtools merge` before use (see README.md Step 1).

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
| 9 — Structural variant calling (Sniffles2, vanilla) | `sv_all` | `workflow/rules/sv_analysis.smk` |

Key outputs:
- `results/vcf_processing/{sample}.phased.vcf.gz` — phased VCFs
- `results/variant_groups/{group}.vcf.gz` — variants per sample-sharing group. `{group}` is an underscore-joined subset of samples in canonical order (`SBC4`, `SBC4_SBC11_SBC23`, `SBC4_SBC10_SBC11_SBC23`, …); all 15 non-empty subsets. Single-sample groups (e.g. `SBC10`) are the old sample-private variants.
- `results/snpeff/{group}.annotated.vcf.gz` — SnpEff-annotated group variants
- `results/sift4g/{group}.sift4g.vcf.gz` — SIFT4G-annotated group variants
- `resources/bedmethyl/{sample}.bed` — raw 5mC pileup
- `resources/bedmethyl/{sample}.filtered.bed` — positions with ≥ 10 valid reads
- `results/DMR/{pair}.5mC.DMR.tsv` — annotated DMRs per TAA contrast pair
- `results/ranked_genes_lists/SBC10.multiomics_ranked.tsv` — final multi-omics gene ranking
- `results/sv_calling/{sample}.sniffles.vcf.gz` — per-sample Sniffles2 SV calls (+ `.snf`)
- `results/sv_calling/combined.sniffles.vcf.gz` — combined multi-sample SV VCF (Sniffles2 `.snf` merge); vanilla, unfiltered and unannotated

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
  variant_groups/   # bcftools isec outputs, one VCF per sample-sharing group
  snpeff/           # SnpEff-annotated VCFs
  sift4g/           # SIFT4G-annotated VCFs
  DSS/              # per-pair DSS input files
  DMR/              # DSS DMR calls and annotations
  ranked_genes_lists/ # multi-omics gene ranking outputs
  sv_calling/       # vanilla Sniffles2 SV VCFs (per-sample + combined) + .snf files
docker/             # Dockerfile, environment.yml, helper shell scripts
dag/                # Snakemake DAG PDFs
```

## Toolchain (inside Docker image `thesis-tools:latest`)

minimap2 2.30, samtools 1.21, bcftools 1.21, htslib 1.21, bedtools 2.31.1, whatshap 2.8, modkit 0.2.6, snpEff (bioconda), SIFT4G, DSS (R/Bioconductor), clair3 (models at `/opt/models/`), sniffles2 2.4 (structural variants), snakemake 9.16.3, Python 3.11 (pandas, matplotlib, seaborn, cyvcf2, matplotlib-venn).

## Key Gotchas

- **SBC11 BAM must be manually merged** from three libraries before the variant calling and methylation steps. See README.md Step 1.
- **SnpEff chromosome synonyms**: VCF contig IDs (e.g. `NC_012870.2`) differ from SnpEff names (`1`, `2`, …). The synonym file at `workflow/scripts/creating_synonyms_chr.py` bridges this.
- **Methylation requires MM/ML tags**: BAMs must be basecalled with Dorado `--modified-bases`. Verify with `modkit summary`; if 0 modified bases reported, re-basecall.
- **SnpEff database**: must be downloaded once with `./docker/run.sh snpEff download Sorghum_bicolor` before running `annotate_vcf`.
- **SVs use a separate merge path**: structural variants are merged across samples with Sniffles2's `.snf` population mode, **not** `bcftools isec` (exact POS/REF/ALT matching fragments SV breakpoints). The Sniffles2 output is kept **vanilla** (raw per-sample + combined VCFs); filtering/annotation/downstream processing is decided after reviewing the raw calls.

## Analysis Scripts

All scripts live in `workflow/scripts/`. Run via `./docker/run.sh python3 workflow/scripts/<script>.py`.

| Script | Purpose |
|--------|---------|
| `variant_landscape.py` | Generates variant landscape figures from SnpEff stats CSVs |
| `variant_upset.py` | UpSet plot of variant set intersections, reading the `results/variant_groups/{group}.vcf.gz` counts (does not run `bcftools isec`) |
| `methylation_landscape.py` | Methylation landscape plots |
| `vcf_benchmark.py` | VCF benchmark figure |
| `annot_vcf_to_tsv.py` | Converts annotated VCF → TSV for notebooks |
| `merge_vcf.sh` | Merges per-sample phased VCFs into multi-sample VCF |
| `rank_dmr_genes.py` | Ranks genes by DMR proximity for multi-omics integration |
| `build_ranked_genes.py` | Builds final multi-omics ranked gene list |

## Reference Docs

- `README.md` — full data generation pipeline with exact commands

# Thesis Narrative

Sorghum (*Sorghum bicolor*) is a resilient, multipurpose crop with agronomic value beyond its grain. Among its secondary metabolites is trans-aconitic acid (TAA), a tricarboxylic acid that accumulates in the stem juice and has attracted industrial interest as a bio-based plasticizer precursor, a sustainable alternative to chemically hazardous synthetic compounds such as DEHP. Compared to other TAA sources, including bioengineered microbial strains and sugarcane, sorghum stands out as the most practical feedstock: its grain and stem serve independent purposes, allowing TAA extraction from the stem without competing with food production. Within the plant, TAA accumulates as a stable carbon pool derived from the TCA cycle intermediate cis-aconitate, with its biosynthesis coupled to photosynthetic activity. TAA also serves as a chemical defence compound, functioning as an antifungal phytoalexin precursor in wheat and as an antifeedant against insect herbivores in grasses, and exerts allelopathic suppression of neighbouring plants when released into soil. Despite these documented roles, the gene encoding the primary biosynthetic enzyme, aconitate isomerase (EC 5.3.3.7), remains unidentified in any plant species, and the downstream transport and methylation steps are equally uncharacterized.

Closing this gap calls for a strategy that does not depend on sequence similarity to the characterized microbial and fungal orthologs. Forward genomics offers such a route: by exploiting the natural variation in TAA accumulation and secretion across *Sorghum bicolor* accessions, candidate genes can be prioritized directly from genotype–phenotype associations rather than from homology. To this end, we construct a multi-omics platform that integrates genomic, epigenomic, and transcriptomic evidence into a ranked list of candidate genes underlying metabolic variation in sorghum. TAA biosynthesis and secretion serve as the primary application of the platform. In parallel, the D-gene, a well-characterized transcription factor controlling stem juiciness, serves as a positive control: recovering this known gene validates that the platform behaves as intended before it is trusted on the unresolved TAA pathway.

The platform integrates three omics layers, all anchored on four sorghum accessions of contrasting phenotypes sequenced with Oxford Nanopore (ONT) long-read sequencing. At the genomic layer, variants are called for each accession against the reference genome, yielding accession-specific variant sites. At the epigenomic layer, DNA methylation is profiled and used to identify differentially methylated regions (DMRs) between each pair of accessions, producing six pairwise comparisons. At the transcriptomic layer, a gene co-expression network is constructed following the ATTED-II framework from publicly available expression data. Evidence from these three layers is then combined to rank candidate genes, with the goal of identifying the gene encoding aconitate isomerase — the first committed step of TAA biosynthesis — in *Sorghum bicolor*.