from datetime import datetime
from pathlib import Path

WDIR      = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR   = f"{WDIR}/workflow/logs/methylation"

REF = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"

SAMPLES = ["SBC4", "SBC10", "SBC11", "SBC23"]

rule methylation_all:
    input:
        expand(
            f"{WDIR}/results/methylation/{{sample}}.bedmethyl.gz",
            sample=SAMPLES
        ),
        expand(
            f"{WDIR}/results/methylation/{{sample}}.bedmethyl.gz.tbi",
            sample=SAMPLES
        ),

rule extract_methylation:
    input:
        bam = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam",
        bai = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam.bai",
        ref = REF,
    output:
        bedmethyl = f"{WDIR}/results/methylation/{{sample}}.bedmethyl.gz",
        tbi       = f"{WDIR}/results/methylation/{{sample}}.bedmethyl.gz.tbi",
    log:
        f"{LOG_DIR}/{{sample}}.{TIMESTAMP}.log",
    threads: 1
    shell:
        """
        set -euo pipefail
        modkit pileup \
            --ref {input.ref} \
            --threads {threads} \
            {input.bam} - \
            2> {log} \
            | bgzip -@ {threads} > {output.bedmethyl}
        tabix -p bed {output.bedmethyl} >> {log} 2>&1
        """
