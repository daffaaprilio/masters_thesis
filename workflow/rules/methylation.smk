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
            f"{WDIR}/resources/bedmethyl/{{sample}}.bed",
            sample=SAMPLES,
        ),
        expand(
            f"{WDIR}/resources/bedmethyl/{{sample}}.filtered.bed",
            sample=SAMPLES,
        ),


rule modkit_pileup:
    """Step 6 — aggregate per-read mod probabilities into per-position bedMethyl."""
    input:
        bam = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam",
        bai = f"{WDIR}/resources/align_bam_sample/{{sample}}.bam.bai",
        ref = REF,
    output:
        bed = f"{WDIR}/resources/bedmethyl/{{sample}}.bed",
    log:
        snakemake = f"{LOG_DIR}/modkit_pileup/{{sample}}.{TIMESTAMP}.log",
        modkit    = f"{WDIR}/logs/modkit/{{sample}}.log",
    threads: 24
    shell:
        """
        modkit pileup \
            {input.bam} \
            {output.bed} \
            --ref {input.ref} \
            --threads {threads} \
            --log-filepath {log.modkit} \
            > {log.snakemake} 2>&1
        """


rule filter_bedmethyl:
    """Step 7 — keep positions with N_valid >= 10 and write summary statistics."""
    input:
        bed = f"{WDIR}/resources/bedmethyl/{{sample}}.bed",
    output:
        filtered = f"{WDIR}/resources/bedmethyl/{{sample}}.filtered.bed",
        summary  = f"{WDIR}/resources/bedmethyl/{{sample}}.summary.txt",
    log:
        f"{LOG_DIR}/filter_bedmethyl/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            awk '$10 >= 10' {input.bed} > {output.filtered}

            echo "=== {wildcards.sample} bedMethyl summary (N_valid >= 10) ===" > {output.summary}

            echo "Positions passing filter:" >> {output.summary}
            wc -l < {output.filtered} >> {output.summary}

            echo "Mean pct_mod and total positions:" >> {output.summary}
            awk 'BEGIN{{n=0; s=0}} {{n++; s+=$11}} END{{print "Mean pct_mod:", s/n, "  N positions:", n}}' \
                {output.filtered} >> {output.summary}

            echo "Highly methylated positions (pct_mod > 80):" >> {output.summary}
            awk '$11 > 80' {output.filtered} | wc -l >> {output.summary}
        ) > {log} 2>&1
        """
