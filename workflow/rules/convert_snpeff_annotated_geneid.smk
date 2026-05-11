from datetime import datetime
from pathlib import Path

WDIR        = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
TIMESTAMP   = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR     = f"{WDIR}/workflow/logs/convert_annotated_geneid"

TAXID = config.get("taxid", 4558)

# Glob on plain .vcf (always present) to derive wildcard combinations.
# annotate_intersected must run first; its output is the actual input below.
STRAT, TYPE, INTERSECTION = glob_wildcards(
    f"{WDIR}/discussions/{{strat}}/{{type}}/{{intersection,[0-9]+}}.vcf"
)

rule convert_all:
    input:
        f"{WDIR}/resources/NCBI_FTP/gene_info_{TAXID}",
        expand(
            f"{WDIR}/discussions/{{strat}}/{{type}}_annotated/{{intersection}}.annotated_EGI.vcf",
            zip, strat=STRAT, type=TYPE, intersection=INTERSECTION,
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