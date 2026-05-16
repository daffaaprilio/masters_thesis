from datetime import datetime
from pathlib import Path

WDIR        = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
TIMESTAMP   = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR     = f"{WDIR}/workflow/logs/reads_preprocessing"

repetition = config.get('repetition', 1)

SAMPLE_LIBRARIES = {
    "SBC4":  ["r0074"],
    "SBC10": ["r0066"],
    "SBC11": ["r0075", "r0078", "r0078-2"],
    "SBC23": ["r0076"],
}
SAMPLES             = list(SAMPLE_LIBRARIES.keys())
SINGLE_LIB_SAMPLES = {k: v[0] for k, v in SAMPLE_LIBRARIES.items() if len(v) == 1}
MULTI_LIB_SAMPLES  = {k: v    for k, v in SAMPLE_LIBRARIES.items() if len(v) > 1}

LIBRARIES = glob_wildcards(f"{WDIR}/resources/trim_bam/{{library}}.bam").library

rule reads_all:
    input:
        # coverage depth files
        expand(
            f"{WDIR}/resources/depth/{{library}}.depth",
            library=LIBRARIES,
            ),
        # coverage depth plots & pickle files
        expand(
            f"{WDIR}/resources/depth/{{library}}_depth{{ext}}",
            library=LIBRARIES,
            ext=[".png", ".svg", ".pdf"]
            ),
        expand(
            f"{WDIR}/resources/depth/{{library}}_depth.pkl",
            library=LIBRARIES,
            ),
        # aligned bam files (per library)
        expand(
            f"{WDIR}/resources/align_bam/{{library}}.bam",
            library=LIBRARIES,
            ),
        # merged sample-level bam files
        expand(
            f"{WDIR}/resources/align_bam_sample/{{sample}}.bam.bai",
            sample=SAMPLES
            ),
        # QC: verify MM/ML tags survived alignment
        expand(
            f"{WDIR}/resources/qc/modkit_summary/{{sample}}.txt",
            sample=SAMPLES,
            ),

rule just_plot:
    input:
        # coverage depth plots with repetition
        expand(
            f"{WDIR}/resources/depth/{{library}}_depth_{repetition}{{ext}}",
            library=LIBRARIES,
            ext=[".png", ".svg", ".pdf"]
            ),


rule align_reads:
    input:
        bam = f"{WDIR}/resources/trim_bam/{{library}}.bam",
        ref = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna",
    output:
        align_bam = f"{WDIR}/resources/align_bam/{{library}}.bam",
    log:
        f"{LOG_DIR}/align_reads/{{library}}.{TIMESTAMP}.log",
    threads: 6
    shell:
        '''
        (
            samtools fastq -T '*' {input.bam} \
                | minimap2 -ax map-ont -t {threads} -y --secondary=no \
                    {input.ref} - \
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


rule symlink_sample_bam:
    """Symlink single-library sample BAMs to avoid duplication."""
    wildcard_constraints:
        sample = "|".join(SINGLE_LIB_SAMPLES.keys()),
    input:
        bam = lambda wc: f"{WDIR}/resources/align_bam/{SINGLE_LIB_SAMPLES[wc.sample]}.bam",
        bai = lambda wc: f"{WDIR}/resources/align_bam/{SINGLE_LIB_SAMPLES[wc.sample]}.bam.bai",
    output:
        bam = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam",
        bai = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam.bai",
    log:
        f"{LOG_DIR}/symlink_sample_bam/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            ln -sf {input.bam} {output.bam}
            ln -sf {input.bai} {output.bai}
        ) > {log} 2>&1
        """


rule merge_sample_bam:
    """Merge multiple library BAMs into one sorted sample BAM."""
    wildcard_constraints:
        sample = "|".join(MULTI_LIB_SAMPLES.keys()),
    input:
        bams = lambda wc: expand(
            f"{WDIR}/resources/align_bam/{{library}}.bam",
            library=MULTI_LIB_SAMPLES[wc.sample]
        ),
        bais = lambda wc: expand(
            f"{WDIR}/resources/align_bam/{{library}}.bam.bai",
            library=MULTI_LIB_SAMPLES[wc.sample]
        ),
    output:
        bam = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam",
        bai = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam.bai",
    log:
        f"{LOG_DIR}/merge_sample_bam/{{sample}}.{TIMESTAMP}.log",
    threads: 6
    shell:
        """
        (
            samtools merge -f -r -@ {threads} - {input.bams} \
                | samtools sort -@ {threads} -o {output.bam}
            samtools index -@ {threads} {output.bam}
        ) > {log} 2>&1
        """


rule validate_mod_tags:
    """QC: verify MM/ML modification tags survived alignment (Step 5)."""
    input:
        bam = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam",
        bai = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam.bai",
    output:
        summary = f"{WDIR}/resources/qc/modkit_summary/{{sample}}.txt",
    log:
        f"{LOG_DIR}/validate_mod_tags/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            echo "=== MM/ML tag presence (first 1000 reads) ===" > {output.summary}
            samtools view {input.bam} | head -1000 \
                | awk '{{found=0; for(i=12;i<=NF;i++) if($i~/^MM:/) found=1; print found}}' \
                | sort | uniq -c >> {output.summary}

            echo "" >> {output.summary}
            echo "=== modkit summary ===" >> {output.summary}
            modkit summary {input.bam} >> {output.summary}
        ) > {log} 2>&1
        """
