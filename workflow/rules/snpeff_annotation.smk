from datetime import datetime
from pathlib import Path

WDIR      = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
SAMPLES   = ["SBC10", "SBC11", "SBC4", "SBC23"]
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR   = f"{WDIR}/workflow/logs/vcf_processing"

VCF_DIR     = f"{WDIR}/results/vcf"
OUT_DIR     = f"{WDIR}/results/vcf_processing"
BAM_DIR     = f"{WDIR}/resources/align_bam_sample"
REF         = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
GFF_RAW     = f"{WDIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.gz"
SNPEFF_DB   = "Sorghum_bicolor"
SNPEFF_DIR  = "/Users/daffa/local/lib/snpEff"
SNPEFF_DATA = f"{SNPEFF_DIR}/data/{SNPEFF_DB}"


# Filtering thresholds — override at runtime with --config key=value
QUAL_MIN = 20
DP_MIN   = 10
DP_MAX   = 100

# Discover intersected VCFs — must precede any rule that uses these lists.
# [0-9]+ constraint excludes _EGI-suffixed files.
STRAT, TYPE, INTERSECTION = glob_wildcards(
    f"{WDIR}/discussions/{{strat}}/{{type}}/{{intersection,[0-9]+}}.vcf"
)

rule annotate_all:
    input:
        expand(f"{OUT_DIR}/annotated/{{sample}}.stats.html", sample=SAMPLES),
        expand(f"{OUT_DIR}/annotated/{{sample}}.stats.csv", sample=SAMPLES),
        expand(
            f"{WDIR}/discussions/{{strat}}/{{type}}_annotated/{{intersection}}.annotated.vcf.gz",
            zip, strat=STRAT, type=TYPE, intersection=INTERSECTION,
        ),

rule annotate_vcf:
    """Annotate variants with gene features and predicted consequences via snpEff."""
    input:
        vcf  = f"{OUT_DIR}/renamed/{{sample}}.renamed.vcf.gz",
        csi  = f"{OUT_DIR}/renamed/{{sample}}.renamed.vcf.gz.csi",
        done = f"{SNPEFF_DATA}/snpEffectPredictor.bin",
    output:
        vcf   = f"{OUT_DIR}/annotated/{{sample}}.annotated.vcf.gz",
        csi   = f"{OUT_DIR}/annotated/{{sample}}.annotated.vcf.gz.csi",
        stats = f"{OUT_DIR}/annotated/{{sample}}.stats.html",
        csv   = f"{OUT_DIR}/annotated/{{sample}}.stats.csv",
    log:
        f"{LOG_DIR}/annotate_vcf/{{sample}}.{TIMESTAMP}.log",
    params:
        db      = SNPEFF_DB,
        datadir = f"{SNPEFF_DIR}/data",
        config  = f"{SNPEFF_DIR}/snpEff.config",
    shell:
        """
        (
            snpEff ann \
                -config {params.config} \
                -dataDir {params.datadir} \
                -v \
                -nodownload \
                -stats {output.stats} \
                -csvStats {output.csv} \
                {params.db} \
                {input.vcf} \
                | bgzip -c > {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """

rule annotate_intersected:
    """Annotate intersected variants with gene features and predicted consequences via snpEff."""
    input:
        vcf  = f"{WDIR}/discussions/{{strat}}/{{type}}/{{intersection}}.vcf",
        done = f"{SNPEFF_DATA}/snpEffectPredictor.bin",
    output:
        vcf   = f"{WDIR}/discussions/{{strat}}/{{type}}_annotated/{{intersection}}.annotated.vcf.gz",
        csi   = f"{WDIR}/discussions/{{strat}}/{{type}}_annotated/{{intersection}}.annotated.vcf.gz.csi",
        stats = f"{WDIR}/discussions/{{strat}}/{{type}}_annotated/{{intersection}}.stats.html",
        csv   = f"{WDIR}/discussions/{{strat}}/{{type}}_annotated/{{intersection}}.stats.csv",
    log:
        f"{LOG_DIR}/annotate_intersected/{{strat}}/{{type}}_annotated/{{intersection}}.{TIMESTAMP}.log",
    params:
        db      = SNPEFF_DB,
        datadir = f"{SNPEFF_DIR}/data",
        config  = f"{SNPEFF_DIR}/snpEff.config",
    shell:
        """
        (
            snpEff ann \
                -config {params.config} \
                -dataDir {params.datadir} \
                -v \
                -nodownload \
                -stats {output.stats} \
                -csvStats {output.csv} \
                {params.db} \
                {input.vcf} \
                | bgzip -c > {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """

# ── Cleanup ───────────────────────────────────────────────────────────────────

rule clean:
    shell:
        """
        rm -rf {OUT_DIR}/annotated {SNPEFF_DIR}
        """
