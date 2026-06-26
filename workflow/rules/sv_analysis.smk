# Structural variant calling (Sniffles2) — vanilla output only.
#
# Produces standard Sniffles2 results for review: per-sample SV VCFs (+ .snf) and
# a combined multi-sample VCF from the .snf population merge. No filtering,
# annotation, or downstream integration — how to process these is decided after
# reviewing the raw calls.

ALIGN_BAM_DIR = f"{WDIR}/resources/align_bam_sample"
SV_DIR        = f"{WDIR}/results/sv_calling"
SV_LOG_DIR    = f"{WDIR}/workflow/logs/sv_analysis"


rule sv_all:
    input:
        expand(f"{SV_DIR}/{{sample}}.sniffles.vcf.gz", sample=SAMPLES),
        f"{SV_DIR}/combined.sniffles.vcf.gz",


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
