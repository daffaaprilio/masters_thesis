# DMR analysis rules:
# bedMethyl → DSS format conversion → DSS DMR calling → DMR annotation.

GFF           = f"{WDIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff"
BEDMETHYL_DIR = f"{WDIR}/resources/bedmethyl"
DSS_DIR       = f"{WDIR}/results/DSS"
DMR_DIR       = f"{WDIR}/results/DMR"

rule prepare_dss:
    """Convert filtered bedMethyl files to DSS input format (collapses Watson/Crick CpG pairs)."""
    input:
        expand(f"{BEDMETHYL_DIR}/{{sample}}.filtered.bed", sample=SAMPLES),
    output:
        dss     = expand(f"{DSS_DIR}/{{sample}}.5mC.dss.txt", sample=SAMPLES),
        summary = f"{DSS_DIR}/conversion_summary.tsv",
    log:
        f"{WDIR}/workflow/logs/prepare_dss.{TIMESTAMP}.log",
    shell:
        """
        python3 {WDIR}/workflow/scripts/prepare_dss.py > {log} 2>&1
        """


rule run_dss_pair:
    """Run DSS DML/DMR calling for one pairwise comparison (5mC only).

    threads: 32 — each DSS run saturates all cores internally; this prevents
    Snakemake from launching more than one pair in parallel on a 32-core host.
    Run with --cores 32 (or match to your machine) to enforce serialisation.
    """
    input:
        expand(f"{DSS_DIR}/{{sample}}.5mC.dss.txt", sample=SAMPLES),
    output:
        dml = f"{DMR_DIR}/{{pair}}.5mC.DML.tsv",
        dmr = f"{DMR_DIR}/{{pair}}.5mC.DMR.tsv",
    log:
        f"{WDIR}/workflow/logs/run_dss_{{pair}}.{TIMESTAMP}.log",
    threads: 32
    shell:
        """
        cd {WDIR} && Rscript workflow/scripts/run_dss_dmr.R {wildcards.pair} > {log} 2>&1
        """


rule summarise_dmr:
    """Combine per-pair DMR tables into DMR_summary.tsv and DMR_all_combined.tsv."""
    input:
        expand(f"{DMR_DIR}/{{pair}}.5mC.DMR.tsv", pair=PAIRS),
    output:
        summary  = f"{DMR_DIR}/DMR_summary.tsv",
        combined = f"{DMR_DIR}/DMR_all_combined.tsv",
    log:
        f"{WDIR}/workflow/logs/summarise_dmr.{TIMESTAMP}.log",
    shell:
        """
        cd {WDIR} && Rscript workflow/scripts/summarise_dss_dmr.R > {log} 2>&1
        """


rule annotate_dmr:
    """Annotate DMRs with strand-aware 2 kbp promoter regions across all genes."""
    input:
        dmr = f"{DMR_DIR}/DMR_all_combined.tsv",
        gff = GFF,
    output:
        f"{DMR_DIR}/DMR_annotated.tsv",
    log:
        f"{WDIR}/workflow/logs/annotate_dmr.{TIMESTAMP}.log",
    shell:
        """
        python3 {WDIR}/workflow/scripts/annotate_DMR.py \
            --dmr    {input.dmr} \
            --gff    {input.gff} \
            --outdir {DMR_DIR} \
            --all-genes \
            > {log} 2>&1
        """
