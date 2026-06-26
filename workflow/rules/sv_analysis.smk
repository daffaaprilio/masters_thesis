# Structural variant analysis rules (Sniffles2):
# per-sample SV calling → population .snf merge → filter → SnpEff annotate →
# parallel SV candidate table cross-referenced against the ranked gene lists.
#
# This track is deliberately SEPARATE from the SNV+methylation ranking: it does
# NOT feed genomic_score. SVs are merged with Sniffles2's own .snf population
# mode (NOT bcftools isec, whose exact POS/REF/ALT matching fragments the wobbly
# breakpoints of structural variants).
#
# Shared constants reused from earlier-included rule files (see Snakefile include
# order): KEEP_CHROMS (vcf_processing.smk), SNPEFF_DIR (snpeff_annotation.smk),
# GENE_INFO + RANKED_GENES_DIR (ranked_genes.smk).

ALIGN_BAM_DIR = f"{WDIR}/resources/align_bam_sample"
SV_DIR        = f"{WDIR}/results/sv_calling"
SV_CAND_DIR   = f"{WDIR}/results/sv_candidates"
SV_LOG_DIR    = f"{WDIR}/workflow/logs/sv_analysis"

# Minimum |SVLEN| to retain in the filtered set (bp).
MIN_SVLEN = 50


rule sv_all:
    input:
        f"{SV_CAND_DIR}/sv_candidate_table.tsv",


rule sniffles_call:
    """Per-sample structural variant calling: emit both a VCF and a .snf for
    population merging."""
    input:
        bam = ancient(f"{ALIGN_BAM_DIR}/{{sample}}.bam"),
        bai = ancient(f"{ALIGN_BAM_DIR}/{{sample}}.bam.bai"),
        ref = REF,
    output:
        vcf = f"{SV_DIR}/{{sample}}.sniffles.vcf.gz",
        snf = f"{SV_DIR}/{{sample}}.snf",
    log:
        f"{SV_LOG_DIR}/sniffles_call/{{sample}}.{TIMESTAMP}.log",
    threads: 8
    shell:
        """
        mkdir -p {SV_DIR} $(dirname {log})
        sniffles \
            --input {input.bam} \
            --reference {input.ref} \
            --vcf {output.vcf} \
            --snf {output.snf} \
            --sample-id {wildcards.sample} \
            --threads {threads} \
            > {log} 2>&1
        """


rule sniffles_combine:
    """Force-genotyped multi-sample SV VCF from per-sample .snf files."""
    input:
        snfs = expand(f"{SV_DIR}/{{sample}}.snf", sample=SAMPLES),
    output:
        vcf = f"{SV_DIR}/combined.sniffles.vcf.gz",
        tbi = f"{SV_DIR}/combined.sniffles.vcf.gz.tbi",
    log:
        f"{SV_LOG_DIR}/sniffles_combine/combine.{TIMESTAMP}.log",
    threads: 4
    shell:
        """
        (
            sniffles --input {input.snfs} \
                --vcf {SV_DIR}/combined.raw.vcf.gz \
                --threads {threads}
            # Re-sort + bgzip defensively so htslib indexing always succeeds.
            bcftools sort -Oz -o {output.vcf} {SV_DIR}/combined.raw.vcf.gz
            bcftools index -t {output.vcf}
            rm -f {SV_DIR}/combined.raw.vcf.gz {SV_DIR}/combined.raw.vcf.gz.tbi
        ) > {log} 2>&1
        """


rule filter_sv:
    """Keep PASS SVs with |SVLEN| >= MIN_SVLEN on the 12 retained contigs.
    Records with no SVLEN (e.g. BND breakends) are kept."""
    input:
        vcf = f"{SV_DIR}/combined.sniffles.vcf.gz",
        tbi = f"{SV_DIR}/combined.sniffles.vcf.gz.tbi",
    output:
        vcf = f"{SV_DIR}/combined.filtered.vcf.gz",
        tbi = f"{SV_DIR}/combined.filtered.vcf.gz.tbi",
    log:
        f"{SV_LOG_DIR}/filter_sv/filter.{TIMESTAMP}.log",
    params:
        regions = KEEP_CHROMS,
        expr    = f"SVLEN>-{MIN_SVLEN} && SVLEN<{MIN_SVLEN}",
    shell:
        """
        (
            bcftools view -f PASS -r {params.regions} {input.vcf} \
                | bcftools filter -e '{params.expr}' -O z -o {output.vcf}
            bcftools index -t {output.vcf}
        ) > {log} 2>&1
        """


rule annotate_sv:
    """Annotate the filtered multi-sample SV VCF with SnpEff (reuses
    annotate_vcf.sh, which inherits the contig-synonym bridge from the built DB).
    SIFT4G is intentionally not applied (substitution-scoring only)."""
    input:
        vcf     = f"{SV_DIR}/combined.filtered.vcf.gz",
        tbi     = f"{SV_DIR}/combined.filtered.vcf.gz.tbi",
        db_flag = f"{SNPEFF_DIR}/data/Sorghum_bicolor_NCBIv3/sequence.NC_012870.2.bin",
    output:
        vcf = f"{SV_DIR}/combined.annotated.vcf.gz",
        csi = f"{SV_DIR}/combined.annotated.vcf.gz.csi",
    log:
        f"{SV_LOG_DIR}/annotate_sv/annotate.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/workflow/scripts/annotate_vcf.sh \
            combined \
            {input.vcf} \
            {SV_DIR} \
            > {log} 2>&1
        """


rule build_sv_table:
    """Build the parallel SV candidate table: gene-associated SVs with per-sample
    genotypes, a TAA-segregation pattern, and one-directional cross-reference
    columns into the multi-omics ranked gene lists."""
    input:
        vcf       = f"{SV_DIR}/combined.annotated.vcf.gz",
        csi       = f"{SV_DIR}/combined.annotated.vcf.gz.csi",
        gene_info = GENE_INFO,
        # Concrete dependency guarantees the ranking step has run; the script
        # additionally globs RANKED_GENES_DIR for any other available samples.
        ranked    = f"{RANKED_GENES_DIR}/SBC10.multiomics_ranked.tsv",
    output:
        table = f"{SV_CAND_DIR}/sv_candidate_table.tsv",
    log:
        f"{SV_LOG_DIR}/build_sv_table/build.{TIMESTAMP}.log",
    shell:
        """
        mkdir -p {SV_CAND_DIR} $(dirname {log})
        python3 {WDIR}/workflow/scripts/build_sv_table.py \
            --vcf        {input.vcf} \
            --gene-info  {input.gene_info} \
            --ranked-dir {RANKED_GENES_DIR} \
            --out        {output.table} \
            > {log} 2>&1
        """
