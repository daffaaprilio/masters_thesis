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

## Running Commands

All tools run inside the `thesis-tools` Docker image. Build it once before running any workflow.

```shell
./docker/build.sh
```

Use `rsync` to transfer files between devices (MacStudio ↔ lab server):

```shell
rsync -avn daffa@matsu:/home/daffa/Work/2026/thesis/resources/bedmethyl ./resources/
```

Use `docker/snakemake.sh` for Snakemake targets and `docker/run.sh` for one-off commands:

```shell
# Run a workflow target
./docker/snakemake.sh <target> --cores <N>

# Dry-run
./docker/snakemake.sh <target> -n

# Ad-hoc command
./docker/run.sh <command>
```

Override image, SnpEff path, or core count via env vars:

```shell
THESIS_IMAGE=thesis-tools:latest CORES=24 SNPEFF_DIR=/path/to/snpEff ./docker/snakemake.sh
```

Print the DAG:

```shell
./docker/run.sh bash -c "cd workflow && snakemake --dag | dot -Tpdf > ../dag/dag.pdf"
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
| 9 — Structural variant calling (parallel track) | `sv_all` | `sv_analysis.smk` |

---

## Step 1 — Read Preprocessing

**Snakefile:** `workflow/rules/reads_preprocessing.smk` | **Target:** `reads_all`

Aligns each library to the reference genome (minimap2), indexes the BAMs (samtools), and computes per-position read depth.

```shell
./docker/snakemake.sh reads_all --cores 24
```

### Per-sample BAMs

SBC4, SBC10, and SBC23 each map to a single library; symlink to avoid duplication:

```shell
declare -A samples=([SBC10]="r0066" [SBC4]="r0074" [SBC23]="r0076")
for sample in ${(k)samples}; do
    ln -s "$(pwd)/resources/align_bam/${samples[$sample]}.bam" "resources/align_bam_sample/${sample}.bam"
    ln -s "$(pwd)/resources/align_bam/${samples[$sample]}.bam.bai" "resources/align_bam_sample/${sample}.bam.bai"
done
```

SBC11 has three libraries that must be merged and re-sorted:

```shell
samtools merge -r -@ 6 - \
    resources/align_bam/r0075.bam \
    resources/align_bam/r0078.bam \
    resources/align_bam/r0078-2.bam \
  | samtools sort -@ 6 -o resources/align_bam_sample/SBC11.bam
samtools index resources/align_bam_sample/SBC11.bam
```

### Depth files for SBC11

The workflow handles single-library samples automatically. For SBC11, generate the depth file and plots manually after merging:

```shell
samtools depth -a resources/align_bam_sample/SBC11.bam -o resources/depth/SBC11.depth

# Plot — saves SBC11_depth.{png,svg,pdf} and SBC11_depth.pkl
python workflow/scripts/visualize_depth.py \
    --input resources/depth/SBC11.depth \
    --output resources/depth/SBC11_depth \
    --library "r0078, r0075, r0078-2" \
    --save-pickle resources/depth/SBC11_depth.pkl

# Re-plot from pickle (skips depth recalculation)
python workflow/scripts/visualize_depth.py \
    --load-pickle resources/depth/SBC11_depth.pkl \
    --output resources/depth/SBC11_depth \
    --library "r0078, r0075, r0078-2"
```

---

## Step 2 — Variant Calling

**Snakefile:** `workflow/rules/variant_analysis.smk` | **Target:** `variants_all`

Calls variants with Clair3 using the ONT R10.4.1 SUP v500 model. Requires per-sample BAMs from Step 1.

```shell
./docker/snakemake.sh variants_all --cores 24
```

Run on a subset of samples:

```shell
./docker/snakemake.sh \
    results/variant_calling/SBC10/merge_output.vcf.gz \
    results/variant_calling/SBC11/merge_output.vcf.gz \
    --cores 16
```

### Link Clair3 output to the VCF directory

```shell
for sample in SBC4 SBC10 SBC11 SBC23; do
    ln -s results/variant_calling/${sample}/merge_output.vcf.gz results/vcf/${sample}.vcf.gz
    ./docker/run.sh bcftools index results/vcf/${sample}.vcf.gz
done
```

---

## Step 3 — VCF Postprocessing

**Snakefile:** `workflow/rules/vcf_processing.smk` | **Target:** `vcf_all`

Filters variants (QUAL ≥ 20, DP 10–100, 10 nuclear chromosomes + MT + chloroplast), phases haplotypes with WhatsHap, and annotates variant consequences with SnpEff.

### One-time SnpEff setup

```shell
./docker/run.sh snpEff download Sorghum_bicolor
```

The database is stored in `${SNPEFF_DIR}/data/Sorghum_bicolor/`. A chromosome synonym file is required because VCF contig IDs (e.g. `NC_012870.2`) differ from SnpEff's numeric names (`1`, `2`, …). Pre-generated files are at `workflow/scripts/vcf_chr_list.txt` and `workflow/scripts/snpeff_db_chr_list.txt`. To regenerate them:

```shell
# VCF contig names (from any sample VCF header)
bcftools view -h results/vcf/SBC4.vcf.gz | grep "^##contig" > workflow/scripts/vcf_chr_list.txt

# SnpEff chromosome list (trim file after the last contig entry to avoid bloat)
./docker/run.sh snpEff dump Sorghum_bicolor > workflow/scripts/snpeff_db_chr_list.txt

# Build the synonym key
python3 workflow/scripts/creating_synonyms_chr.py
```

### Run VCF postprocessing

```shell
./docker/snakemake.sh vcf_all --cores 24
```

Rules applied in order:
1. `bcftools reheader` — embed SAMPLE tag in each VCF
2. Filter to PASS variants on the 10 nuclear chromosomes + organelles
3. `whatshap phase` — add haplotype information
4. `snpEff ann` — annotate variant consequences

---

## Step 4 — Methylation Calling

**Snakefile:** `workflow/rules/methylation.smk` | **Target:** `methylation_all`

Produces per-sample 5mC bedMethyl files using modkit. Requires aligned per-sample BAMs from Step 1 that carry `MM`/`ML` tags (present when Dorado is run with `--modified-bases`).

### Verify MM/ML tags before running

```shell
samtools view resources/align_bam_sample/SBC10.bam | head -1000 | \
    awk '{found=0; for(i=12;i<=NF;i++) if($i~/^MM:/) found=1; print found}' | \
    sort | uniq -c
# Expected: 1000 1   (all reads carry MM tag)

./docker/run.sh modkit summary resources/align_bam_sample/SBC10.bam
```

If `modkit summary` reports 0 modified bases, the Dorado BAMs were basecalled without `--modified-bases` and must be re-basecalled.

### Run methylation calling

```shell
./docker/snakemake.sh methylation_all --cores 24
```

Outputs:
- `resources/bedmethyl/{sample}.bed` — raw pileup (CpG, strands combined)
- `resources/bedmethyl/{sample}.filtered.bed` — positions with ≥ 10 valid reads

---

## Step 9 — Structural Variant Calling (parallel track)

**Snakefile:** `workflow/rules/sv_analysis.smk` | **Target:** `sv_all`

Calls structural variants (large deletions, insertions, duplications, inversions,
breakends) that Clair3 cannot see, using **Sniffles2** on the per-sample BAMs from
Step 1. This is a standalone track: it does **not** feed the SNV+methylation
`genomic_score`. SVs are merged across samples with Sniffles2's own `.snf`
population mode (not `bcftools isec`, which fragments SV breakpoints), annotated
with the existing SnpEff database (reusing `annotate_vcf.sh`), and assembled into a
parallel candidate table that cross-references the ranked gene lists.

SVs intentionally **skip SIFT4G** (substitution-scoring only) and **WhatsHap
phasing**.

### Run SV calling

```shell
./docker/run.sh snakemake sv_all --cores 8
```

Pipeline: `sniffles_call` (per-sample VCF + `.snf`) → `sniffles_combine`
(force-genotyped multi-sample VCF) → `filter_sv` (PASS, |SVLEN| ≥ 50 bp, 12
retained contigs) → `annotate_sv` (SnpEff) → `build_sv_table`.

Outputs:
- `results/sv_calling/{sample}.sniffles.vcf.gz`, `{sample}.snf` — per-sample SV calls
- `results/sv_calling/combined.annotated.vcf.gz` — SnpEff-annotated multi-sample SV VCF
- `results/sv_candidates/sv_candidate_table.tsv` — gene-associated SVs with per-sample
  genotypes, a TAA-segregation pattern (`high-specific` / `low-specific` / `mixed`),
  and cross-reference columns (`in_ranked_list`, `best_rank`, `best_variant_score`,
  `best_methylation_score`) into `results/ranked_genes_lists/`

The candidate table is meant to be read alongside the multi-omics ranking: an SV
whose gene also ranks highly in the SNV+methylation list is a strong candidate.
