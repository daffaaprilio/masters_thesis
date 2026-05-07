from datetime import datetime
from pathlib import Path

WDIR      = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
SAMPLES   = ["SBC4", "SBC10", "SBC11", "SBC23"]
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR   = f"{WDIR}/workflow/logs/variant_analysis"

REF     = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
# MODEL   = "/opt/models/r1041_e82_400bps_sup_v520_with_mv"
MODEL   = "/opt/models/r1041_e82_400bps_sup_v500"


rule variants_all:
    input:
        expand(
            f"{WDIR}/results/vcf/{{sample}}.vcf.gz",
            sample=SAMPLES
        )


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
        podman run --rm \
          -v {WDIR}/resources/:{WDIR}/resources/:z \
          -v {WDIR}/results/:{WDIR}/results/:z \
          docker.io/hkubal/clair3:latest \
          /opt/bin/run_clair3.sh \
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
        vcf = f"{WDIR}/results/variant_calling/{{sample}}/merge_output.vcf.gz",
    output:
        vcf = f"{WDIR}/results/vcf/{{sample}}.vcf.gz",
        csi = f"{WDIR}/results/vcf/{{sample}}.vcf.gz.csi",
    log:
        f"{LOG_DIR}/publish_vcf/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            ln -sf {input.vcf} {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """
