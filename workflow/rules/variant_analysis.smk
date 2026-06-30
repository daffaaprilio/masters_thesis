#   SNP track: variant calling + per-sample VCF processing
#   - Calling:    Clair3 (ONT R10.4.1 duplex SUP model)
#   - Processing: reheader -> QUAL/DP filter -> chromosome filter -> WhatsHap phasing

from datetime import datetime
from pathlib import Path

WDIR      = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
SAMPLES   = ["SBC4", "SBC10", "SBC11", "SBC23"]
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR   = f"{WDIR}/workflow/logs/variant_analysis"

REF     = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
# MODEL   = "/opt/models/r1041_e82_400bps_sup_v520_with_mv"
MODEL   = "/opt/models/r1041_e82_400bps_sup_v500"

# ── VCF processing config (per-sample filter + phasing) ────────────────────────
VCF_DIR  = f"{WDIR}/results/vcf"
OUT_DIR  = f"{WDIR}/results/vcf_processing"
BAM_DIR  = f"{WDIR}/resources/align_bam_sample"

QUAL_MIN = 20
DP_MIN   = 10
DP_MAX   = 100

# Chromosomes to retain: 10 main nuclear + MT + chloroplast
KEEP_CHROMS = ",".join(
    [f"NC_01287{i}.2" for i in range(10)] + ["NC_008360.1", "NC_008602.1"]
)


rule variants_all:
    input:
        expand(f"{VCF_DIR}/{{sample}}.vcf.gz", sample=SAMPLES)


rule vcf_all:
    input:
        expand(f"{OUT_DIR}/{{sample}}.phased.vcf.gz", sample=SAMPLES),
        expand(f"{OUT_DIR}/{{sample}}.phased.vcf.gz.csi", sample=SAMPLES),


rule clair3_cpu:
    input:
        bam = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam",
        bai = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam.bai",
        ref = REF,
    output:
        vcf = f"{WDIR}/results/variant_calling/{{sample}}/merge_output.vcf.gz",
    log:
        f"{LOG_DIR}/clair3_cpu/{{sample}}.{TIMESTAMP}.log",
    threads: 8
    shell:
        """
        run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.ref} \
            --threads={threads} \
            --platform=ont \
            --model_path={MODEL} \
            --output={WDIR}/results/variant_calling/{wildcards.sample}/ \
            --include_all_ctgs \
            > {log} 2>&1
        """


rule publish_vcf:
    """Symlink Clair3 merge output to the canonical results/vcf/ location."""
    input:
        vcf = ancient(f"{WDIR}/results/variant_calling/{{sample}}/merge_output.vcf.gz"),
    output:
        vcf = f"{VCF_DIR}/{{sample}}.vcf.gz",
        csi = f"{VCF_DIR}/{{sample}}.vcf.gz.csi",
    log:
        f"{LOG_DIR}/publish_vcf/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            out_dir=$(dirname {output.vcf})
            ln -sf "$(realpath --relative-to="$out_dir" {input.vcf})" {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


rule reheader_vcf:
    """Rename the generic SAMPLE column header to the actual sample name."""
    input:
        vcf = f"{VCF_DIR}/{{sample}}.vcf.gz",
        csi = f"{VCF_DIR}/{{sample}}.vcf.gz.csi",
    output:
        vcf = f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz",
        csi = f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz.csi",
    log:
        f"{LOG_DIR}/reheader_vcf/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            echo "{wildcards.sample}" | bcftools reheader -s /dev/stdin {input.vcf} -o {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


rule filter_vcf:
    """Keep PASS variants within configurable QUAL and DP bounds."""
    input:
        vcf = f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz",
        csi = f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz.csi",
    output:
        vcf = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz",
        csi = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz.csi",
    log:
        f"{LOG_DIR}/filter_vcf/{{sample}}.{TIMESTAMP}.log",
    params:
        expr = f"QUAL<{QUAL_MIN} || FORMAT/DP<{DP_MIN} || FORMAT/DP>{DP_MAX}",
    shell:
        """
        (
            bcftools view -f PASS {input.vcf} \
                | bcftools filter -e '{params.expr}' -O z -o {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


rule filter_chromosomes:
    """Retain only the 10 main nuclear chromosomes, MT, and chloroplast."""
    input:
        vcf = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz",
        csi = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz.csi",
    output:
        vcf = f"{OUT_DIR}/chr_filtered/{{sample}}.chr_filtered.vcf.gz",
        csi = f"{OUT_DIR}/chr_filtered/{{sample}}.chr_filtered.vcf.gz.csi",
    log:
        f"{LOG_DIR}/filter_chromosomes/{{sample}}.{TIMESTAMP}.log",
    params:
        regions = KEEP_CHROMS,
    shell:
        """
        (
            bcftools view -r {params.regions} -O z -o {output.vcf} {input.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


rule phase_vcf:
    """Read-backed haplotype phasing with WhatsHap per sample."""
    input:
        vcf = f"{OUT_DIR}/chr_filtered/{{sample}}.chr_filtered.vcf.gz",
        csi = f"{OUT_DIR}/chr_filtered/{{sample}}.chr_filtered.vcf.gz.csi",
        bam = f"{BAM_DIR}/{{sample}}.bam",
        bai = f"{BAM_DIR}/{{sample}}.bam.bai",
        ref = REF,
    output:
        vcf = f"{OUT_DIR}/{{sample}}.phased.vcf.gz",
        csi = f"{OUT_DIR}/{{sample}}.phased.vcf.gz.csi",
    log:
        f"{LOG_DIR}/phase_vcf/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            whatshap phase \
                --reference {input.ref} \
                --ignore-read-groups \
                --output /dev/stdout \
                {input.vcf} {input.bam} \
                | bgzip -c > {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """
