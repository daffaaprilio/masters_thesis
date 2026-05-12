from datetime import datetime
from pathlib import Path

WDIR      = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR   = f"{WDIR}/workflow/logs/filter_candidates"


rule filter_private_candidates:
    """Filter private annotated EGI VCF for high-confidence candidate variants."""
    wildcard_constraints:
        intersection = r"[0-9]+",
    input:
        vcf    = f"{WDIR}/discussions/{{strat}}/private_annotated/{{intersection}}.annotated_EGI.vcf",
        script = f"{WDIR}/workflow/scripts/filter_candidate_variants.py",
    output:
        tsv    = f"{WDIR}/discussions/{{strat}}/private_annotated/{{intersection}}.candidate_variants.tsv",
    log:
        f"{LOG_DIR}/{{strat}}/{{intersection}}_{TIMESTAMP}.log",
    shell:
        """
        python3 {input.script} \
            --input {input.vcf} \
            --output {output.tsv} \
        2>&1 | tee {log}
        """


def get_private_egi_vcfs_for_filter(wc):
    """Wait for the private_variants checkpoint, return candidate TSV paths."""
    checkpoints.private_variants.get()
    intersections, = glob_wildcards(
        f"{WDIR}/discussions/{wc.strat}/private/{{intersection,[0-9]+}}.vcf"
    )
    return expand(
        f"{WDIR}/discussions/{wc.strat}/private_annotated/{{intersection}}.candidate_variants.tsv",
        intersection=intersections,
    )


rule collect_filtered_private:
    """Sentinel: all private annotated VCFs have been filtered to candidate TSVs."""
    input:  get_private_egi_vcfs_for_filter
    output: touch(f"{WDIR}/discussions/{{strat}}/private_annotated/.filter_done")
