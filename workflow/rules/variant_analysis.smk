WDIR        = "/home/daffa/Work/2026/thesis"

rule all:
    input:
        expand(f"{WDIR}/results/gvcf/{{sample}}.g.vcf.gz", sample=["SBC4", "SBC10", "SBC11", "SBC23"])

rule call_variants_per_sample:
    input:
        bam = f"{WDIR}/resources/align_bam/{{sample}}.aligned.bam",
        ref = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
    output:
        f"{WDIR}/results/gvcf/{{sample}}.g.vcf.gz"
    shell:
        '''
        podman run --rm \
            -v {WDIR}:{WDIR}:Z -w {WDIR} \
            broadinstitute/gatk:latest gatk --java-options "-Xmx8g" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            -ERC GVCF
        '''