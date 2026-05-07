from datetime import datetime
from pathlib import Path

WDIR      = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR   = f"{WDIR}/workflow/logs/intersect_variant"

STRAT     = config.get("strat", "all_samples")

# Add new named strategies here
STRAT_SAMPLES = {
    "all_samples": ["SBC4", "SBC10", "SBC11", "SBC23"],
}

# Direct override: --config isec_samples=SBC10,SBC4
ISEC_SAMPLES = [
    s.strip()
    for s in config.get("isec_samples", ",".join(STRAT_SAMPLES[STRAT])).split(",")
]

OUT_DIR = f"{WDIR}/discussions/{STRAT}"


ANN_DIR = f"{WDIR}/results/vcf_processing/annotated"


rule isec_all:
    input:
        f"{OUT_DIR}/private/README.txt",
        f"{OUT_DIR}/shared/README.txt",


rule private_variants:
    """Variants private to each individual sample (present in exactly one sample)."""
    input:
        vcfs = expand(f"{ANN_DIR}/{{sample}}.annotated.vcf.gz",     sample=ISEC_SAMPLES),
        csis = expand(f"{ANN_DIR}/{{sample}}.annotated.vcf.gz.csi", sample=ISEC_SAMPLES),
    output:
        readme = f"{OUT_DIR}/private/README.txt",
    log:
        f"{LOG_DIR}/{STRAT}.private.{TIMESTAMP}.log",
    params:
        outdir = f"{OUT_DIR}/private",
    shell:
        """
        bcftools isec -n =1 -p {params.outdir} {input.vcfs} > {log} 2>&1
        """


rule shared_variants:
    """Variants shared across all samples in the strategy."""
    input:
        vcfs = expand(f"{ANN_DIR}/{{sample}}.annotated.vcf.gz",     sample=ISEC_SAMPLES),
        csis = expand(f"{ANN_DIR}/{{sample}}.annotated.vcf.gz.csi", sample=ISEC_SAMPLES),
    output:
        readme = f"{OUT_DIR}/shared/README.txt",
    log:
        f"{LOG_DIR}/{STRAT}.shared.{TIMESTAMP}.log",
    params:
        outdir   = f"{OUT_DIR}/shared",
        n_shared = len(ISEC_SAMPLES),
    shell:
        """
        bcftools isec -n ={params.n_shared} -p {params.outdir} {input.vcfs} > {log} 2>&1
        """