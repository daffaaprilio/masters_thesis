# Structural variant calling (Sniffles2) — vanilla output only.
#
# Produces standard Sniffles2 results for review: per-sample SV VCFs (+ .snf) and
# a combined multi-sample VCF from the .snf population merge. No filtering,
# annotation, or downstream integration — how to process these is decided after
# reviewing the raw calls.

ALIGN_BAM_DIR = f"{WDIR}/resources/align_bam_sample"
SV_DIR        = f"{WDIR}/results/sv_calling"
SV_GROUP_DIR  = f"{WDIR}/results/sv_groups"
SV_LOG_DIR    = f"{WDIR}/workflow/logs/sv_analysis"


def _sv_group_expr(pattern):
    """bcftools -i expression selecting SVs with this exact present/absent pattern.

    `pattern` is positional over SAMPLES (SBC4,SBC10,SBC11,SBC23) — the same order
    as the combined VCF's sample columns and as the SNP VARGROUPS patterns, so an SV
    group lines up 1:1 with its SNP variant_groups counterpart. Presence is GT-based
    (carries an ALT allele), which is stricter and cleaner than SUPP_VEC (read-support
    based). This is the SV analog of the SNP `intersect_group` exact-membership split,
    GT-aware instead of `bcftools isec` (SV breakpoints wobble, so isec over-fragments).
    """
    terms = [
        f'GT[{i}]="alt"' if bit == "1" else f'GT[{i}]!="alt"'
        for i, bit in enumerate(pattern)
    ]
    return " && ".join(terms)


rule sv_all:
    input:
        expand(f"{SV_DIR}/{{sample}}.sniffles.vcf.gz", sample=SAMPLES),
        f"{SV_DIR}/combined.sniffles.vcf.gz",
        expand(f"{SV_GROUP_DIR}/{{group}}.vcf.gz", group=VARGROUP_LABELS),


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


rule sv_group:
    """Split the combined Sniffles2 VCF into the 15 sample-sharing groups — the SV
    counterpart of results/variant_groups/. One VCF per group label, selected by the
    exact GT presence pattern across the four samples (see _sv_group_expr)."""
    input:
        vcf = f"{SV_DIR}/combined.sniffles.vcf.gz",
        tbi = f"{SV_DIR}/combined.sniffles.vcf.gz.tbi",
    output:
        vcf = f"{SV_GROUP_DIR}/{{group}}.vcf.gz",
        tbi = f"{SV_GROUP_DIR}/{{group}}.vcf.gz.tbi",
    log:
        f"{SV_LOG_DIR}/sv_group/{{group}}.{TIMESTAMP}.log",
    params:
        expr = lambda wc: _sv_group_expr(VARGROUPS[wc.group]["pattern"]),
    shell:
        """
        (
            mkdir -p {SV_GROUP_DIR}
            bcftools view -i '{params.expr}' -O z -o {output.vcf} {input.vcf}
            bcftools index -t {output.vcf}
        ) > {log} 2>&1
        """
