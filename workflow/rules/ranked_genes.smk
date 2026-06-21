# Ranked gene list rule:
# Parse SIFT-annotated VCF (SnpEff ANN + SIFTINFO) + DMR_annotated.tsv
# → per-sample TSV with genomic_score and epigenomic_score per gene.
#
# DMR_DIR is defined in rules/dmr_analysis.smk (included before this file).

RANKED_GENES_DIR = f"{WDIR}/results/ranked_genes_lists"
GENE_INFO        = f"{WDIR}/resources/NCBI_FTP/gene_info_4558"


rule ranked_genes:
    """Merge SnpEff+SIFT variant annotations with DMR epigenomics → ranked gene TSV.

    Genomic score  = impact_score + 0.5 × sift_disruption
      impact_score:    HIGH=1.0  MODERATE=0.67  LOW=0.33  MODIFIER=0.0
      sift_disruption: 1 − min(SIFT_score) for scored CDS variants (0 if absent)

    Epigenomic score = max(|diff.Methy|) across all DMR rows involving the sample.

    Output is sorted by genomic_score DESC, epigenomic_score DESC.
    """
    input:
        # size-1 variant group == sample-private variants
        vcf       = f"{SIFT_DIR}/{{sample}}.sift4g.vcf.gz",
        tbi       = f"{SIFT_DIR}/{{sample}}.sift4g.vcf.gz.tbi",
        dmr       = f"{DMR_DIR}/DMR_annotated.tsv",
        gene_info = GENE_INFO,
    output:
        f"{RANKED_GENES_DIR}/{{sample}}.multiomics_ranked.tsv",
    log:
        f"{WDIR}/workflow/logs/ranked_genes/ranked_genes.{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        mkdir -p {RANKED_GENES_DIR} $(dirname {log})
        python3 {WDIR}/workflow/scripts/build_ranked_genes.py \
            --vcf       {input.vcf} \
            --dmr       {input.dmr} \
            --gene-info {input.gene_info} \
            --sample    {wildcards.sample} \
            --out       {output} \
            > {log} 2>&1
        """
