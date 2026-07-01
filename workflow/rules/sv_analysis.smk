#   SV track: structural variant calling (Sniffles2) — vanilla output only.
#
#   Per-sample SV VCFs (+ .snf) and a combined multi-sample VCF from the .snf
#   population merge. Annotation and downstream integration live in vcf_annotation.smk.
#   SV_DIR is defined as a shared path in the Snakefile.

ALIGN_BAM_DIR = f"{WDIR}/resources/align_bam_sample"
SV_LOG_DIR    = f"{WDIR}/workflow/logs/sv_analysis"

# Mega-SVs: |SVLEN| >= 100 kb are whole-chromosome-arm DUP/INV/DEL artifacts
# (long-read mapping noise across repeats) that SnpEff explodes into thousands of
# per-gene ANN rows. 100 kb is the distribution elbow: SVs below it average ~1.8
# spanned genes, those >= 1 Mb average 260-1220. Such calls are ~0.65% of SVs but
# ~74% of the SV->gene row bloat. BND breakpoints carry no SVLEN and are kept.
MEGA_SVLEN = 100_000


rule sv_all:
    input:
        expand(f"{SV_DIR}/{{sample}}.sniffles.vcf.gz", sample=SAMPLES),
        expand(f"{SV_DIR}/filtered/{{sample}}.sniffles.filtered.vcf.gz", sample=SAMPLES),
        f"{SV_DIR}/combined.sniffles.vcf.gz",
        f"{SV_DIR}/combined.sniffles.filtered.vcf.gz",


rule sniffles_call:
    """Per-sample structural variant calling: VCF + .snf for population merging."""
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


rule filter_sv:
    """Drop mega-SVs (|SVLEN| >= MEGA_SVLEN) from each per-sample Sniffles2 VCF.
    SV-track counterpart of the SNP-track filter_vcf (params.expr + bcftools filter
    -e). BND records have no SVLEN and are always kept.

    NOTE: this is a standalone per-sample product. sniffles_combine population-merges
    the .snf files (not these VCFs), so the combined/annotated output is de-mega'd
    separately by filter_combined_sv below — that is the path feeding the notebook.
    """
    input:
        vcf = f"{SV_DIR}/{{sample}}.sniffles.vcf.gz",
    output:
        vcf = f"{SV_DIR}/filtered/{{sample}}.sniffles.filtered.vcf.gz",
        csi = f"{SV_DIR}/filtered/{{sample}}.sniffles.filtered.vcf.gz.csi",
    log:
        f"{SV_LOG_DIR}/filter_sv/{{sample}}.{TIMESTAMP}.log",
    params:
        expr = f"SVLEN>={MEGA_SVLEN} || SVLEN<=-{MEGA_SVLEN}",
    shell:
        """
        (
            mkdir -p $(dirname {output.vcf}) $(dirname {log})
            bcftools filter -e '{params.expr}' -O z -o {output.vcf} {input.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


rule filter_combined_sv:
    """Drop mega-SVs (|SVLEN| >= MEGA_SVLEN) from the combined multi-sample SV VCF.
    This is the version that feeds SnpEff annotation and, downstream, the notebook's
    results/combined/all.annotated.vcf.gz — so this is the filter that actually removes
    the per-gene mega-SV row explosion from the analysis. BND breakpoints (no SVLEN)
    are kept. Same MEGA_SVLEN cutoff and expression as the per-sample filter_sv."""
    input:
        vcf = f"{SV_DIR}/combined.sniffles.vcf.gz",
        tbi = f"{SV_DIR}/combined.sniffles.vcf.gz.tbi",
    output:
        vcf = f"{SV_DIR}/combined.sniffles.filtered.vcf.gz",
        tbi = f"{SV_DIR}/combined.sniffles.filtered.vcf.gz.tbi",
    log:
        f"{SV_LOG_DIR}/filter_combined_sv/combined.{TIMESTAMP}.log",
    params:
        expr = f"SVLEN>={MEGA_SVLEN} || SVLEN<=-{MEGA_SVLEN}",
    shell:
        """
        (
            mkdir -p $(dirname {output.vcf}) $(dirname {log})
            bcftools filter -e '{params.expr}' -O z -o {output.vcf} {input.vcf}
            bcftools index -t {output.vcf}
        ) > {log} 2>&1
        """


rule sniffles_combine:
    """Force-genotyped multi-sample SV VCF from the per-sample .snf files."""
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
