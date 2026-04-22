WDIR    = "/home/daffa/Work/2026/thesis"
SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]

REF     = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
# MODEL   = "/opt/models/r1041_e82_400bps_sup_v520_with_mv"
MODEL   = "/opt/models/r1041_e82_400bps_sup_v500"


rule all:
    input:
        expand(
            f"{WDIR}/results/variant_calling/{{sample}}/merge_output.vcf.gz",
            sample=SAMPLES
        )


rule clair3_cpu:
    input:
        bam = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam",
        bai = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam.bai",
        ref = REF,
    output:
        vcf = f"{WDIR}/results/variant_calling/{{sample}}/merge_output.vcf.gz",
    threads: 8
    shell:
        """
        podman run --rm \
          -v {WDIR}/resources/:{WDIR}/resources/ \
          -v {WDIR}/results/:{WDIR}/results/ \
          docker.io/hkubal/clair3:latest \
          /opt/bin/run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.ref} \
            --threads={threads} \
            --platform=ont \
            --model_path={MODEL} \
            --output={WDIR}/results/variant_calling/{wildcards.sample}/ \
            --include_all_ctgs
        """
