# CLAUDE.md — Sorghum Multi-omics Thesis

## Project Overview

Multi-omics analysis of four *Sorghum bicolor* samples (BTx623 reference, `GCF_000003195.3`) with contrasting TAA (trans-aconitic acid) production phenotypes. Another phenotype, stem juiciness, that is regulated by D-gene transcription factor, is also included to the story, mainly as the positive control. Samples are sequenced with ONT R10.4.1 duplex 400 bps SUP (Dorado v1.3.0). For gene co-expression analysis, use data from ATTED-II database.

Samples:

| Cultivar name (Region) | Library | Sample | TAA conc. in juice | Juice Production | Callus Formation |
|--|--------|--------|---------------|---------------|-----------------|
| IS32787 (Somalia) | r0074 | SBC4 | High | ++ | Mid |
| IS20956 (Indonesia) | r0066 | SBC10 | Low | +++ | Good |
| S. VULGARE 72-726-7 (Uganda) | r0075, r0078, r0078-2 | SBC11 | High | - | Mid |
| 240 WAD UMM BENEIN (Sudan) | r0076 | SBC23 | High | ++ | Good |

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
| 8 — SV calling + SnpEff SV→gene table | `sv_all` | `workflow/rules/sv_analysis.smk` |

> Multi-omics gene ranking is **not** part of the Snakemake pipeline — it is performed manually in a notebook environment from the layer outputs above (SnpEff/SIFT variant annotations, DMRs, co-expression).

Key outputs:
- `results/vcf_processing/{sample}.phased.vcf.gz` — phased VCFs
- `results/variant_groups/{group}.vcf.gz` — variants per sample-sharing group. `{group}` is an underscore-joined subset of samples in canonical order (`SBC4`, `SBC4_SBC11_SBC23`, `SBC4_SBC10_SBC11_SBC23`, …); all 15 non-empty subsets. Single-sample groups (e.g. `SBC10`) are the old sample-private variants.
- `results/snpeff/{group}.annotated.vcf.gz` — SnpEff-annotated group variants
- `results/sift4g/{group}.sift4g.vcf.gz` — SIFT4G-annotated group variants
- `resources/bedmethyl/{sample}.bed` — raw 5mC pileup
- `resources/bedmethyl/{sample}.filtered.bed` — positions with ≥ 10 valid reads
- `results/DMR/{pair}.5mC.DMR.tsv` — annotated DMRs per TAA contrast pair
- `results/sv_calling/{sample}.sniffles.vcf.gz` — per-sample Sniffles2 SV calls (+ `.snf`)
- `results/sv_calling/combined.sniffles.vcf.gz` — combined multi-sample SV VCF (Sniffles2 `.snf` merge); vanilla, unfiltered and unannotated
- `results/sv_groups/{group}.vcf.gz` — combined SV VCF split into the 15 sample-sharing groups (the SV counterpart of `variant_groups/`); GT-based exact membership, same `{group}` labels
- `results/snpeff_sv/{group}.annotated.vcf.gz` — SnpEff-annotated SV groups (SV counterpart of `results/snpeff/`)
- `results/sv_genes/{group}.sv_genes.tsv` — SV→gene table from the SnpEff `ANN` field; long format, one row per SV–gene with effect/impact, no scoring

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
  sv_calling/       # vanilla Sniffles2 SV VCFs (per-sample + combined) + .snf files
  sv_groups/        # combined SV VCF split into the 15 sample-sharing groups
  snpeff_sv/        # SnpEff-annotated SV groups
  sv_genes/         # SV→gene tables (one row per SV–gene, from SnpEff ANN)
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
- **SV sample-sharing groups** (`results/sv_groups/`) are produced by splitting the combined VCF on **GT presence** (`GT[i]="alt"` / `!="alt"`), not `bcftools isec` and not `SUPP_VEC` (which is read-support based and looser). Same 15 `{group}` labels as `variant_groups/`, so SV and SNP groups line up 1:1.
- **SV→gene uses SnpEff `ANN`, not GFF overlap**: `annotate_sv` runs the same `annotate_vcf.sh` + DB as the SNP groups, then `sv_gene_mapping.py` parses the `ANN` field (one entry per gene, worst impact wins). SnpEff classifies the per-gene consequence (`transcript_ablation`, `feature_ablation`, `exon_loss_variant`, `feature_fusion` for BND, `duplication`/`inversion` for DUP/INV) and also tags upstream/downstream genes within ~5 kb as MODIFIER. It does **not** cap multi-gene SVs — a multi-megabase DUP/INV still emits one ANN per spanned gene, so **mega-SV artifacts (|SVLEN| in the tens of Mb, ≈ whole chromosome arms) dominate `sv_genes/` row counts (~90%+)**; filter by `svlen`/`impact`/`effect` downstream. SnpEff also emits a gene-less chromosome-level summary entry per DUP/INV (dropped by the mapper).

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
| `rank_dmr_genes.py` | Ranks genes by DMR proximity (standalone helper for the manual notebook-based gene ranking; not wired into Snakemake) |
| `genomics_scoring.py`, `genomics_scoring_stable.py` | Per-gene genomic scoring from annotated VCFs (standalone helpers for the manual notebook-based gene ranking; not wired into Snakemake) |
| `sv_gene_mapping.py` | Maps each SV in a **SnpEff-annotated** SV group VCF (`results/snpeff_sv/{group}.annotated.vcf.gz`) to the gene(s) it affects, by parsing the `ANN` field (one ANN per gene, worst impact wins); long-format one-row-per-SV–gene TSV with effect/impact, **no scoring** → `results/sv_genes/{group}.sv_genes.tsv` |

## Reference Docs

- `README.md` — full data generation pipeline with exact commands

# Thesis Narrative

> Writing guide
> Use a simple language, while maintaining academic tone

Sorghum (*Sorghum bicolor*) is a resilient, multipurpose crop whose natural genetic diversity offers a route to understanding the molecular basis of agronomically important metabolic and developmental traits. This thesis constructs a multi-omics platform that combines genomic variant analysis, epigenomic (DNA methylation) analysis, and transcriptomic co-expression analysis to nominate candidate causative genes for phenotypic variation in sorghum, anchored on four accessions of contrasting phenotypes sequenced with Oxford Nanopore (ONT) long-read sequencing. Rather than collapsing evidence from the three omics layers into a single combined score, each layer independently produces its own list of candidate causative genes for a given phenotype; a phenotype-relevant gene module --- curated from the literature or inferred from co-expression with known regulators --- is then used to test whether each candidate list is enriched for genes of established functional relevance.

The platform is applied to two case studies. The first targets the genetic basis of stem juiciness, a trait governed by developmentally programmed cell death (dPCD) of the stem pith parenchyma and controlled by the D-gene, a transcription factor with an established causal role, making it a tractable positive control. A module of dPCD executor genes is first curated from the literature. Variant analysis lists genes associated with the phenotype across the juicy accessions relative to the dry accession, and DMR analysis lists genes with promoter methylation changes consistent with transcriptional repression of the dPCD module; both candidate lists are tested for enrichment of the module, with variant-derived candidates showing modest but significant enrichment and methylation-derived candidates showing none. A gene co-expression network built from the ATTED-II database, using the dPCD module genes as query genes, failed to identify a cluster unifying them, and did not recover the D-gene itself --- a mixed result that qualifies, rather than cleanly validates, the platform's positive control.

The second case study targets the genetic basis of trans-aconitic acid (TAA) biosynthesis, a plant secondary metabolite of interest as a bio-based plasticizer precursor whose central biosynthetic enzyme, aconitate isomerase (EC 5.3.3.7), remains unidentified in any plant species. Variant and DMR analyses, applied analogously to contrast the low-TAA accession against the three high-TAA accessions, each produce a candidate causative gene list that may include the aconitate isomerase gene. Because no literature-curated module exists for TAA biosynthesis, the co-expression layer instead ranks every isomerase-family gene annotated in the sorghum genome by its co-expression with genes regulated by phytochrome, reflecting the hypothesis that aconitate isomerase activity is light-regulated; the highest-ranking isomerase genes by this criterion constitute the prioritised aconitate isomerase candidates. Together, the two case studies demonstrate the platform's capacity to connect multi-omics evidence to candidate genes in sorghum, while the stem-juiciness case study also exposes its current limits for cleanly validating a known genetic basis.