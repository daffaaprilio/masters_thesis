from datetime import datetime

WDIR        = "/home/daffa/Work/2026/thesis"
TIMESTAMP   = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR     = f"{WDIR}/workflow/logs/reads_preprocessing"

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
            f"{WDIR}/resources/depth/{{library}}_depth{{ext}}",
            library=glob_wildcards(f"{WDIR}/resources/fastq/{{library}}.fq").library,
            ext=[".png", ".svg", ".pdf"]
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
            f"{WDIR}/resources/depth/{{library}}_depth_{repetition}{{ext}}",
            library=glob_wildcards(f"{WDIR}/resources/fastq/{{library}}.fq").library,
            ext=[".png", ".svg", ".pdf"]
            ),


rule align_reads:
    input:
        fastq = f"{WDIR}/resources/fastq/{{library}}.fq",
        ref = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna",
        fai = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.fai",
    output:
        align_bam = f"{WDIR}/resources/align_bam/{{library}}.bam",
    log:
        f"{LOG_DIR}/align_reads/{{library}}.{TIMESTAMP}.log",
    threads: 2
    shell:
        '''
        (
            minimap2 -ax map-ont -t {threads} {input.ref} {input.fastq} \
                | samtools sort -@ {threads} -o {output.align_bam}
        ) > {log} 2>&1
        '''


rule index_bam:
    input:
        bam = f"{WDIR}/resources/align_bam/{{library}}.bam",
    output:
        bai = f"{WDIR}/resources/align_bam/{{library}}.bam.bai",
    log:
        f"{LOG_DIR}/index_bam/{{library}}.{TIMESTAMP}.log",
    shell:
        "samtools index {input.bam} > {log} 2>&1"


rule coverage_depth:
    input:
        bam = f"{WDIR}/resources/align_bam/{{library}}.bam",
        bai = f"{WDIR}/resources/align_bam/{{library}}.bam.bai",
    output:
        depth = f"{WDIR}/resources/depth/{{library}}.depth",
    log:
        f"{LOG_DIR}/coverage_depth/{{library}}.{TIMESTAMP}.log",
    shell:
        "samtools depth -a {input.bam} -o {output.depth} > {log} 2>&1"


rule visualize_depthfile:
    input:
        depth = f"{WDIR}/resources/depth/{{library}}.depth",
    output:
        plot_png = f"{WDIR}/resources/depth/{{library}}_depth.png",
        plot_svg = f"{WDIR}/resources/depth/{{library}}_depth.svg",
        plot_pdf = f"{WDIR}/resources/depth/{{library}}_depth.pdf",
        pickle   = f"{WDIR}/resources/depth/{{library}}_depth.pkl",
    log:
        f"{LOG_DIR}/visualize_depthfile/{{library}}.{TIMESTAMP}.log",
    shell:
        """
        python {WDIR}/workflow/scripts/visualize_depth.py \
            --input {input.depth} \
            --output {WDIR}/resources/depth/{wildcards.library}_depth \
            --library {wildcards.library} \
            --save-pickle {output.pickle} \
            > {log} 2>&1
        """

rule bypass_pickle:
    input:
        pickle = f"{WDIR}/resources/depth/{{library}}_depth.pkl",
    output:
        plot_png = f"{WDIR}/resources/depth/{{library}}_depth_{repetition}.png",
        plot_svg = f"{WDIR}/resources/depth/{{library}}_depth_{repetition}.svg",
        plot_pdf = f"{WDIR}/resources/depth/{{library}}_depth_{repetition}.pdf",
    log:
        f"{LOG_DIR}/bypass_pickle/{{library}}.{TIMESTAMP}.log",
    shell:
        """
        python {WDIR}/workflow/scripts/visualize_depth.py \
            --load-pickle {input.pickle} \
            --output {WDIR}/resources/depth/{wildcards.library}_depth_{repetition} \
            --library {wildcards.library} \
            > {log} 2>&1
        """
