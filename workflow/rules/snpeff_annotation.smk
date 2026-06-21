# SnpEff annotation rules:
# per-group variant extraction (every sample-sharing subset) → SnpEff annotation.

SNPEFF_DIR = f"{WDIR}/resources/snpeff"

rule snpeff_prep:
    """Build custom SnpEff database for Sorghum bicolor from the NCBI GTF (one-time setup)."""
    output:
        f"{SNPEFF_DIR}/data/Sorghum_bicolor_NCBIv3/sequence.NC_012870.2.bin",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/snpeff_prep.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/workflow/scripts/snpeff_prep.sh > {log} 2>&1
        """


rule intersect_group:
    """Extract variants whose presence pattern across the four samples == this group
    (one UpSet cell) using bcftools isec. Size-1 groups == sample-private variants.

    -n ~PATTERN selects sites matching the exact present/absent pattern (positional
    over SAMPLES); -w writes the records from the first member's file.
    """
    input:
        vcfs = expand(f"{PHASED_VCF_DIR}/{{sample}}.phased.vcf.gz", sample=SAMPLES),
        csis = expand(f"{PHASED_VCF_DIR}/{{sample}}.phased.vcf.gz.csi", sample=SAMPLES),
    output:
        vcf = f"{VARGROUP_DIR}/{{group}}.vcf.gz",
        tbi = f"{VARGROUP_DIR}/{{group}}.vcf.gz.tbi",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/intersect.{{group}}.{TIMESTAMP}.log",
    params:
        pattern = lambda wc: VARGROUPS[wc.group]["pattern"],
        w       = lambda wc: VARGROUPS[wc.group]["w"],
    shell:
        """
        (
            mkdir -p {VARGROUP_DIR}
            bcftools isec -n ~{params.pattern} -w{params.w} -O z \
                -o {output.vcf} {input.vcfs}
            bcftools index -t {output.vcf}
        ) > {log} 2>&1
        """


rule annotate_vcf:
    """Annotate a variant group with SnpEff and index the output."""
    input:
        vcf     = f"{VARGROUP_DIR}/{{group}}.vcf.gz",
        tbi     = f"{VARGROUP_DIR}/{{group}}.vcf.gz.tbi",
        db_flag = f"{SNPEFF_DIR}/data/Sorghum_bicolor_NCBIv3/sequence.NC_012870.2.bin",
    output:
        vcf = f"{SNPEFF_ANNOT_DIR}/{{group}}.annotated.vcf.gz",
        csi = f"{SNPEFF_ANNOT_DIR}/{{group}}.annotated.vcf.gz.csi",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/annotate_vcf.{{group}}.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/workflow/scripts/annotate_vcf.sh \
            {wildcards.group} \
            {input.vcf} \
            {SNPEFF_ANNOT_DIR} \
            > {log} 2>&1
        """


# rule vcf_to_tsv:
#     """Convert an annotated variant group VCF to TSV for notebook exploration."""
#     input:
#         vcf = f"{SNPEFF_ANNOT_DIR}/{{group}}.annotated.vcf.gz",
#         csi = f"{SNPEFF_ANNOT_DIR}/{{group}}.annotated.vcf.gz.csi",
#     output:
#         f"{TSV_DIR}/{{group}}.annotated.tsv",
#     log:
#         f"{WDIR}/workflow/logs/vcf_annotation/vcf_to_tsv.{{group}}.{TIMESTAMP}.log",
#     shell:
#         """
#         python3 {WDIR}/workflow/scripts/annot_single_vcf_to_tsv.py \
#             -v {input.vcf} \
#             -o {TSV_DIR} \
#             > {log} 2>&1
#         """
