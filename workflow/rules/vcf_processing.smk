from datetime import datetime
from pathlib import Path

WDIR      = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
SAMPLES   = ["SBC10", "SBC11", "SBC4", "SBC23"]
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR   = f"{WDIR}/workflow/logs/vcf_processing"

VCF_DIR     = f"{WDIR}/results/vcf"
OUT_DIR     = f"{WDIR}/results/vcf_processing"
ANN_DIR     = f"{WDIR}/results/vcf_for_annotation"
BAM_DIR     = f"{WDIR}/resources/align_bam_sample"
REF         = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
GFF_RAW     = f"{WDIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.gz"


# Filtering thresholds
QUAL_MIN = 20
DP_MIN   = 10
DP_MAX   = 100

# Chromosomes to retain: 10 main nuclear + MT + chloroplast
KEEP_CHROMS = ",".join(
    [f"NC_01287{i}.2" for i in range(10)] + ["NC_008360.1", "NC_008602.1"]
)


rule vcf_all:
    input:
        expand(f"{OUT_DIR}/{{sample}}.phased.vcf.gz", sample=SAMPLES),
        expand(f"{OUT_DIR}/{{sample}}.phased.vcf.gz.csi", sample=SAMPLES),        
        expand(f"{ANN_DIR}/{{sample}}.renamed.vcf.gz", sample=SAMPLES),
        expand(f"{ANN_DIR}/{{sample}}.renamed.vcf.gz.csi", sample=SAMPLES),

rule reheader_vcf:
    """Rename the generic SAMPLE column header to the actual sample name."""
    input:
        vcf = f"{VCF_DIR}/{{sample}}.vcf.gz",
        csi = f"{VCF_DIR}/{{sample}}.vcf.gz.csi",
    output:
        vcf = temp(f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz"),
        csi = temp(f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz.csi"),
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
        vcf = temp(f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz"),
        csi = temp(f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz.csi"),
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
        vcf = temp(f"{OUT_DIR}/chr_filtered/{{sample}}.chr_filtered.vcf.gz"),
        csi = temp(f"{OUT_DIR}/chr_filtered/{{sample}}.chr_filtered.vcf.gz.csi"),
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


rule rename_chromosomes:
    """Rename VCF contig names to match SnpEff database chromosome naming."""
    input:
        vcf        = f"{OUT_DIR}/{{sample}}.phased.vcf.gz",
        csi        = f"{OUT_DIR}/{{sample}}.phased.vcf.gz.csi",
        rename_map = f"{WDIR}/workflow/scripts/synonyms.txt",
    output:
        vcf = f"{ANN_DIR}/{{sample}}.renamed.vcf.gz",
        csi = f"{ANN_DIR}/{{sample}}.renamed.vcf.gz.csi",
    log:
        f"{LOG_DIR}/rename_chromosomes/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            bcftools annotate \
                --rename-chrs {input.rename_map} \
                -O z -o {output.vcf} {input.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """
