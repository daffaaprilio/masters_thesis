WDIR        = "/home/daffa/Work/2026/thesis"

rule all:
    input:
        # aligned bam files
        expand(
            f"{WDIR}/resources/align_bam/{{library}}.bam", 
            library=glob_wildcards(f"{WDIR}/resources/fastq/{{library}}.fq").library
            ),
        # bam files with additional metadata (retrieved from trim_bam, the original bam_file)
        # expand(
        #     f"{WDIR}/resources/metadata_bam/{{library}}.bam", 
        #     library=glob_wildcards(f"{WDIR}/resources/fastq/{{library}}.fq").library
        #     ),
        # final bam ready for haplotype calling
        # expand(
        #     f"{WDIR}/results/bam/{{sample}}.bam", sample=["SBC4", "SBC10", "SBC11", "SBC23"]
        #     )


rule align_reads:
    input:
        fastq = f"{WDIR}/resources/fastq/{{library}}.fq",
        ref = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna",
        fai = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.fai",
    output:
        align_bam = f"{WDIR}/resources/align_bam/{{library}}.bam",
    threads: 2
    shell:
        '''
        minimap2 -ax map-ont -t {threads} {input.ref} {input.fastq}
        '''