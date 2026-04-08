WDIR        = "/home/daffa/Work/2026/thesis"

rule all:
    input:
        expand(f"{WDIR}/resources/align_bam/{{sample}}.aligned.bam", sample=["SBC4", "SBC10", "SBC11", "SBC23"]),
        expand(f"{WDIR}/resources/align_sam/{{sample}}.aligned.sam", sample=["SBC4", "SBC10", "SBC11", "SBC23"])


rule merge_bam_files:
    # merge bam files for SBC11 only
    input:
        f"{WDIR}/resources/trim_bam/r0075.bam",
        f"{WDIR}/resources/trim_bam/r0078.bam",
        f"{WDIR}/resources/trim_bam/r0078-2.bam"
    output:
        f"{WDIR}/resources/bam/SBC11.bam"
    shell:
        '''
        samtools merge -u {output} {input}
        '''

rule align_bam_files:
    input:
        bam = f"{WDIR}/resources/bam/{{sample}}.bam",
        ref = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
    output:
        f"{WDIR}/resources/align_bam/{{sample}}.aligned.bam"
    shell:
        '''
        samtools fastq -T '*' {input.bam} | minimap2 -ax map-ont -y --MD -t 8 {input.ref} - | samtools sort -@ 8 -o {output}
        samtools index {output}
        '''

rule create_sam_from_bam:
    input:
        bam = f"{WDIR}/resources/align_bam/{{sample}}.aligned.bam"
    output:
        f"{WDIR}/resources/align_sam/{{sample}}.aligned.sam"
    shell:
        '''
        samtools view -o {output} {input}
        '''