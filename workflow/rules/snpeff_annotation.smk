# SnpEff annotation rules:
# private variant extraction per sample → SnpEff annotation → VCF-to-TSV conversion.

rule snpeff_prep:
    """Download prebuilt SnpEff database for Sorghum bicolor (one-time setup)."""
    output:
        f"{SNPEFF_DIR}/data/Sorghum_bicolor/sequence.1.bin",
    log:
        f"{LOG_DIR}/snpeff_prep.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/workflow/scripts/snpeff_prep.sh > {log} 2>&1
        """


rule private_variants:
    """Find variants private to each sample using bcftools isec (-n =1)."""
    input:
        vcfs = expand(f"{PHASED_VCF_DIR}/{{sample}}.phased.vcf.gz", sample=SAMPLES),
        csis = expand(f"{PHASED_VCF_DIR}/{{sample}}.phased.vcf.gz.csi", sample=SAMPLES),
    output:
        vcfs = expand(f"{PRIVATE_DIR}/{{sample}}.private.vcf.gz", sample=SAMPLES),
        tbis = expand(f"{PRIVATE_DIR}/{{sample}}.private.vcf.gz.tbi", sample=SAMPLES),
    log:
        f"{LOG_DIR}/private_variants.{TIMESTAMP}.log",
    params:
        out_dir = PRIVATE_DIR,
        samples = " ".join(SAMPLES),
    shell:
        """
        (
            mkdir -p {params.out_dir}
            bcftools isec -n =1 -p {params.out_dir} {input.vcfs}

            # Rename 0000.vcf … to sample-named files (order matches SAMPLES)
            samples=({params.samples})
            for i in "${{!samples[@]}}"; do
                idx=$(printf "%04d" "$i")
                src="{params.out_dir}/${{idx}}.vcf"
                dst="{params.out_dir}/${{samples[$i]}}.private.vcf"
                if [[ -f "$src" ]]; then
                    mv "$src" "$dst"
                    bgzip -f "$dst"
                    bcftools index -t "${{dst}}.gz"
                fi
            done
        ) > {log} 2>&1
        """


rule annotate_vcf:
    """Annotate private variants with SnpEff and index the output."""
    input:
        vcf     = f"{PRIVATE_DIR}/{{sample}}.private.vcf.gz",
        tbi     = f"{PRIVATE_DIR}/{{sample}}.private.vcf.gz.tbi",
        db_flag = f"{SNPEFF_DIR}/data/Sorghum_bicolor/sequence.1.bin",
    output:
        vcf = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz",
        csi = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz.csi",
    log:
        f"{LOG_DIR}/annotate_vcf.{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/workflow/scripts/annotate_vcf.sh \
            {wildcards.sample}.private \
            {input.vcf} \
            {PRIVATE_DIR} \
            > {log} 2>&1
        """


rule vcf_to_tsv:
    """Convert annotated private VCF to TSV for notebook exploration."""
    input:
        vcf = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz",
        csi = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz.csi",
    output:
        f"{TSV_DIR}/{{sample}}.private.annotated.tsv",
    log:
        f"{LOG_DIR}/vcf_to_tsv.{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        python3 {WDIR}/workflow/scripts/annot_single_vcf_to_tsv.py \
            -v {input.vcf} \
            -o {TSV_DIR} \
            > {log} 2>&1
        """
