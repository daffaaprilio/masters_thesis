WDIR        = "/home/daffa/Work/2026/thesis"

repetition = config.get('repetition', 1)

rule all:
    input:
        # coverage depth files
        expand(
            f"{WDIR}/resources/depth/{{library}}.depth",
            library=glob_wildcards(f"{WDIR}/resources/fastq/{{library}}.fq").library
            ),
        # coverage depth plots & pickle files
        expand(
            f"{WDIR}/resources/depth/{{library}}_depth.pdf",
            library=glob_wildcards(f"{WDIR}/resources/fastq/{{library}}.fq").library
            ),
        expand(
            f"{WDIR}/resources/depth/{{library}}_depth.pkl",
            library=glob_wildcards(f"{WDIR}/resources/fastq/{{library}}.fq").library
            ),
        # aligned bam files
        expand(
            f"{WDIR}/resources/align_bam/{{library}}.bam",
            library=glob_wildcards(f"{WDIR}/resources/fastq/{{library}}.fq").library
            ),

rule just_plot:
    input:
        # coverage depth plots with repetition
        expand(
            f"{WDIR}/resources/depth/{{library}}_depth_{repetition}.pdf",
            library=glob_wildcards(f"{WDIR}/resources/fastq/{{library}}.fq").library
            ),


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
        minimap2 -ax map-ont -t {threads} {input.ref} {input.fastq} \
            | samtools sort -@ {threads} -o {output.align_bam}
        '''


rule index_bam:
    input:
        bam = f"{WDIR}/resources/align_bam/{{library}}.bam",
    output:
        bai = f"{WDIR}/resources/align_bam/{{library}}.bam.bai",
    shell:
        "samtools index {input.bam}"


rule coverage_depth:
    input:
        bam = f"{WDIR}/resources/align_bam/{{library}}.bam",
        bai = f"{WDIR}/resources/align_bam/{{library}}.bam.bai",
    output:
        depth = f"{WDIR}/resources/depth/{{library}}.depth",
    shell:
        "samtools depth -a {input.bam} -o {output.depth}"


rule visualize_depthfile:
    input:
        depth = f"{WDIR}/resources/depth/{{library}}.depth",
    output:
        plot = f"{WDIR}/resources/depth/{{library}}_depth.pdf",
        pickle = f"{WDIR}/resources/depth/{{library}}_depth.pkl",
    shell:
        """
        python {WDIR}/workflow/scripts/visualize_depth.py \
            --input {input.depth} \
            --output {output.plot} \
            --library {wildcards.library} \
            --save-pickle {output.pickle}
        """

rule bypass_pickle:
    input:
        pickle = f"{WDIR}/resources/depth/{{library}}_depth.pkl",
    output:
        plot = f"{WDIR}/resources/depth/{{library}}_depth_{repetition}.pdf",
    shell:
        """
        python {WDIR}/workflow/scripts/visualize_depth.py \
        --load-pickle {input.pickle} \
        --output {output.plot} \
        --library {wildcards.library}
        """

