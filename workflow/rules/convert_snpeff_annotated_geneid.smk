from datetime import datetime
from pathlib import Path

WDIR        = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
TIMESTAMP   = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR     = f"{WDIR}/workflow/logs/convert_annotated_geneid"

TAXID     = config.get("taxid", 4558)
STRAT_CFG = config.get("strat", "all_samples")
SAMPLES   = ["SBC4", "SBC10", "SBC11", "SBC23"]
OUT_DIR   = f"{WDIR}/results/vcf_processing"

rule convert_all:
    input:
        f"{WDIR}/resources/NCBI_FTP/gene_info_{TAXID}",
        expand(
            f"{OUT_DIR}/annotated/{{sample}}.annotated_EGI.vcf",
            sample=SAMPLES,
        ),
        expand(
            f"{WDIR}/discussions/{STRAT_CFG}/{{type}}_annotated/.egi_done",
            type=["private", "shared"],
        ),

rule filter_gene_info:
    """Filter NCBI gene_info for a specific taxon ID using dockerized PySpark."""
    input:
        gene_info = f"{WDIR}/resources/NCBI_FTP/gene_info",
        script    = f"{WDIR}/workflow/scripts/extract_sbi_ncbi_gene_info.py",
    output:
        filtered  = f"{WDIR}/resources/NCBI_FTP/gene_info_{TAXID}",
    log:
        f"{LOG_DIR}/filter_gene_info_{TIMESTAMP}.log",
    params:
        taxid         = TAXID,
        data_dir      = f"{WDIR}/resources/NCBI_FTP",
        driver_memory = "12g",
    shell:
        """
        docker run --rm \
            -v {input.script}:/opt/spark/work-dir/extract_sbi_ncbi_gene_info.py \
            -v {params.data_dir}:/data \
            apache/spark-py \
            /opt/spark/bin/spark-submit \
                /opt/spark/work-dir/extract_sbi_ncbi_gene_info.py \
                --input /data/gene_info \
                --output /data/gene_info_{params.taxid} \
                --tax-id {params.taxid} \
                --driver-memory {params.driver_memory} \
        2>&1 | tee {log}
        """

rule convert_sample_vcf_annot_geneid:
    """Convert SnpEff LocusTags to NCBI GeneIDs for a per-sample annotated VCF."""
    input:
        vcf                = f"{OUT_DIR}/annotated/{{sample}}.annotated.vcf.gz",
        script             = f"{WDIR}/workflow/scripts/snpeff_output_convert2EGI.py",
        filtered_gene_info = f"{WDIR}/resources/NCBI_FTP/gene_info_{TAXID}",
    output:
        output_vcf         = f"{OUT_DIR}/annotated/{{sample}}.annotated_EGI.vcf",
    log:
        f"{LOG_DIR}/convert_annot_vcf_{{sample}}_{TIMESTAMP}.log",
    shell:
        """
        python3 {input.script} \
            --input {input.vcf} \
            --output {output.output_vcf} \
            --key {input.filtered_gene_info} \
        2>&1 | tee {log}
        """

rule convert_vcf_annot_geneid:
    """Convert SnpEff LocusTags to NCBI GeneIDs for an annotated intersected VCF."""
    input:
        vcf                = f"{WDIR}/discussions/{{strat}}/{{type}}_annotated/{{intersection}}.annotated.vcf.gz",
        script             = f"{WDIR}/workflow/scripts/snpeff_output_convert2EGI.py",
        filtered_gene_info = f"{WDIR}/resources/NCBI_FTP/gene_info_{TAXID}",
    output:
        output_vcf         = f"{WDIR}/discussions/{{strat}}/{{type}}_annotated/{{intersection}}.annotated_EGI.vcf",
    log:
        f"{LOG_DIR}/convert_annot_vcf_{{strat}}_{{type}}_{{intersection}}_{TIMESTAMP}.log",
    shell:
        """
        python3 {input.script} \
            --input {input.vcf} \
            --output {output.output_vcf} \
            --key {input.filtered_gene_info} \
        2>&1 | tee {log}
        """

def get_annotated_intersected_vcfs(wc):
    """Wait for the relevant checkpoint, then return the actual EGI-converted VCF paths."""
    if wc.type == "private":
        checkpoints.private_variants.get()
    else:
        checkpoints.shared_variants.get()
    intersections, = glob_wildcards(
        f"{WDIR}/discussions/{wc.strat}/{wc.type}/{{intersection,[0-9]+}}.vcf"
    )
    return expand(
        f"{WDIR}/discussions/{wc.strat}/{wc.type}_annotated/{{intersection}}.annotated_EGI.vcf",
        intersection=intersections,
    )

rule collect_egi_converted:
    """Sentinel: all intersected VCFs for a given strat/type have been EGI-converted."""
    input: get_annotated_intersected_vcfs
    output: touch(f"{WDIR}/discussions/{{strat}}/{{type}}_annotated/.egi_done")