# Multi-omics Integration for Sorghum Metabolic Variation Analysis

Genomics (variant calling), transcriptomics (gene co-expression), and epigenomics (methylation) of four *Sorghum bicolor* samples with contrasting TAA production phenotypes.

## Samples

| Cultivar name (Region) | Library | Sample | TAA conc. in juice | Juice Production | Callus Formation |
|--|--------|--------|---------------|---------------|-----------------|
| IS32787 (Somalia) | r0074 | SBC4 | High | ++ | Mid |
| IS20956 (Indonesia) | r0066 | SBC10 | Low | +++ | Good |
| S. VULGARE 72-726-7 (Uganda) | r0075, r0078, r0078-2 | SBC11 | High | - | Mid |
| 240 WAD UMM BENEIN (Sudan) | r0076 | SBC23 | High | ++ | Good |

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
│   │   └── dmr_analysis.smk          # DSS DMR calling (TAA contrast)
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
| 8 — SV calling + SnpEff SV→gene table | `sv_all` | `sv_analysis.smk` |

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

## Step 3 — VCF Filtering & Phasing

**Snakefile:** `workflow/rules/vcf_processing.smk` | **Target:** `vcf_all`

Reheaders each Clair3 VCF with its sample name, filters variants (QUAL ≥ 20, DP 10–100), restricts to the 10 nuclear chromosomes + MT + chloroplast, and phases haplotypes read-backed with WhatsHap.

```shell
./docker/snakemake.sh vcf_all --cores 24
```

Rules applied in order:
1. `reheader_vcf` — replace the generic `SAMPLE` column name with the real sample name
2. `filter_vcf` — keep PASS variants within the QUAL/DP bounds
3. `filter_chromosomes` — retain the 10 nuclear chromosomes + MT + chloroplast
4. `phase_vcf` — WhatsHap read-backed haplotype phasing (per sample, against the reference BAM)

Output: `results/vcf_processing/{sample}.phased.vcf.gz` (+ `.csi`).

---

## Step 4 — SnpEff Annotation

**Snakefile:** `workflow/rules/snpeff_annotation.smk` | **Target:** `annotate_vcf`

Splits the four phased VCFs into the 15 **sample-sharing groups** (every non-empty subset of samples) with `bcftools isec`, then annotates each group's variant consequences with SnpEff.

A `{group}` label is the underscore-joined subset of samples in canonical order (`SBC4`, `SBC4_SBC23`, `SBC4_SBC10_SBC11_SBC23`, …). Single-sample groups (e.g. `SBC10`) are that sample's private variants.

### One-time SnpEff database build

SnpEff uses a **custom** database built from the same NCBI GTF as SIFT4G (Step 5), so both tools share NCBI chromosome IDs (`NC_012870.2`, …) — no chromosome-synonym/rename step is needed.

The `snpeff_prep` rule builds the database automatically the first time `annotate_vcf` runs. To build it explicitly:

```shell
./docker/run.sh bash workflow/scripts/snpeff_prep.sh
```

The database lands in `resources/snpeff/data/Sorghum_bicolor_NCBIv3/`.

### Run annotation

```shell
./docker/snakemake.sh annotate_vcf --cores 24
```

Pipeline: `snpeff_prep` (one-time DB) → `intersect_group` (15 sample-sharing group VCFs) → `annotate_vcf` (SnpEff `ANN` field per group).

Outputs:
- `results/variant_groups/{group}.vcf.gz` — variants per sample-sharing group
- `results/snpeff/{group}.annotated.vcf.gz` — SnpEff-annotated group variants (+ per-group `.stats.html`/`.stats.csv`)

---

## Step 5 — SIFT4G Annotation

**Snakefile:** `workflow/rules/sift_annotation.smk` | **Target:** `annotate_sift`

Adds amino-acid–level functional predictions (SIFT score; < 0.05 = DAMAGING) on top of the SnpEff-annotated group VCFs from Step 4.

### One-time SIFT4G database build

The prediction database is built once from the genome FASTA, the NCBI GTF, and the UniRef90 protein set (`prepare_sift_db` → `build_sift_db`; several hours). It is rebuilt only if you remove the sentinel `resources/sift4g/sorghum_db/.setup_done`. The finished DB lands at `resources/sift4g/sorghum_db/NCBIv3/`.

### Run annotation

```shell
./docker/snakemake.sh annotate_sift --cores 24
```

For each group, `SIFT4G_Annotator.jar` consumes the SnpEff-annotated VCF and emits a SIFT-annotated VCF plus a `_SIFTannotations.xls` table; the VCF is then sorted, bgzipped, and indexed.

Outputs:
- `results/sift4g/{group}.sift4g.vcf.gz` — SnpEff + SIFT4G–annotated group variants
- `results/sift4g/{group}.annotated_SIFTannotations.xls` — per-variant SIFT table

---

## Step 6 — Methylation Calling

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

## Step 7 — DMR Analysis (TAA)

**Snakefile:** `workflow/rules/dmr_analysis.smk` | **Target:** `annotate_dmr`

Calls differentially methylated regions (5mC) between samples with DSS, then annotates them to genes. Requires the filtered bedMethyl files from Step 6.

```shell
# DSS saturates all cores per pair; run with --cores 32 (or match your host) so
# Snakemake serialises pairs and avoids oversubscription.
./docker/snakemake.sh annotate_dmr --cores 32
```

Pipeline: `prepare_dss` (bedMethyl → DSS input, collapsing Watson/Crick CpG pairs) → `run_dss_pair` (DML/DMR per pairwise comparison) → `summarise_dmr` (combine all pairs) → `annotate_dmr` (assign DMRs to genes via strand-aware 2 kbp promoter regions).

Six pairwise comparisons are run: `SBC10_vs_SBC4`, `SBC10_vs_SBC11`, `SBC10_vs_SBC23`, `SBC11_vs_SBC4`, `SBC11_vs_SBC23`, `SBC4_vs_SBC23`.

Outputs:
- `results/DMR/{pair}.5mC.DMR.tsv`, `{pair}.5mC.DML.tsv` — per-pair DMR/DML calls
- `results/DMR/DMR_all_combined.tsv`, `DMR_summary.tsv` — combined tables
- `results/DMR/DMR_annotated.tsv` — DMRs annotated to genes

---

## Step 8 — Structural Variants (Sniffles2 → SnpEff SV→gene table)

**Snakefile:** `workflow/rules/sv_analysis.smk` | **Target:** `sv_all`

Calls structural variants (large deletions, insertions, duplications, inversions,
breakends) that Clair3 cannot see, using **Sniffles2** on the per-sample BAMs from
Step 1, then splits, annotates, and maps them to genes.

Cross-sample merging uses Sniffles2's own `.snf` population mode (not
`bcftools isec`, whose exact POS/REF/ALT matching fragments SV breakpoints).

### Run the SV pipeline

```shell
# SnpEff is a memory-hungry JVM; on a low-RAM host run sequentially (--cores 1).
./docker/run.sh snakemake sv_all --cores 1
```

Pipeline: `sniffles_call` (per-sample VCF + `.snf`) → `sniffles_combine`
(force-genotyped multi-sample VCF) → `sv_group` (split into sample-sharing groups)
→ `annotate_sv` (SnpEff) → `sv_gene_table` (SV→gene mapping).

`sv_group` splits the combined VCF into the 15 sample-sharing groups — the SV
counterpart of `variant_groups/`, with the same `{group}` labels. Membership is
**GT-based** (`GT[i]="alt"` / `!="alt"`), the SV analog of the SNP `intersect_group`
exact-membership split — not `bcftools isec` (SV breakpoints wobble) and not
`SUPP_VEC` (read-support based, looser).

`annotate_sv` reuses the same `annotate_vcf.sh` + SnpEff DB as the SNP groups, so SV
ANN gene IDs are `LOC*` and line up with the SNP track. `sv_gene_table` then runs
`sv_gene_mapping.py`, which parses the SnpEff `ANN` field (one entry per gene, worst
impact wins) and emits a long-format table — **one row per SV–gene**, with the SV's
per-gene effect and impact, **no scoring**. SnpEff does not cap multi-gene SVs, so a
multi-megabase DUP/INV (≈ whole chromosome arm) still produces one row per spanned
gene; these mega-SV artifacts dominate the row counts and should be filtered by
`svlen`/`impact`/`effect` downstream.

Outputs:
- `results/sv_calling/{sample}.sniffles.vcf.gz`, `{sample}.snf` — per-sample SV calls
- `results/sv_calling/combined.sniffles.vcf.gz` — combined multi-sample SV VCF (raw)
- `results/sv_groups/{group}.vcf.gz` — combined SV VCF split into the 15 sample-sharing groups
- `results/snpeff_sv/{group}.annotated.vcf.gz` — SnpEff-annotated SV groups
- `results/sv_genes/{group}.sv_genes.tsv` — SV→gene table (one row per SV–gene, effect/impact, no scoring)

---

## Multi-omics Gene Ranking (manual)

The final integration — ranking candidate genes across the genomic (Steps 4–5), epigenomic (Step 7), and SV (Step 8) layers plus the ATTED-II co-expression network — is **not** part of the Snakemake pipeline. It is performed manually in a notebook environment from the layer outputs above. The standalone helper scripts `workflow/scripts/rank_dmr_genes.py`, `genomics_scoring.py`, and `genomics_scoring_stable.py` are available as references.
