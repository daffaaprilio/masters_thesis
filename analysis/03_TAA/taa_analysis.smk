# Snakefile for TAA candidate gene identification:
# covers both genomics (private variants → SnpEff annotation) and
# epigenomics (bedMethyl → DSS DMR calling → promoter annotation).

from datetime import datetime
from pathlib import Path

# workflow.basedir = analysis/03_TAA/ → parent.parent = repo root
WDIR      = config.get("WDIR", str(Path(workflow.basedir).parent.parent))
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR   = f"{WDIR}/analysis/logs"

SAMPLES      = ["SBC4", "SBC10", "SBC11", "SBC23"]
PAIRS        = ["SBC10_vs_SBC4", "SBC10_vs_SBC11", "SBC10_vs_SBC23"]

# ── Genomics paths ─────────────────────────────────────────────────────────────
VCF_DIR     = f"{WDIR}/results/vcf_processing"
PRIVATE_DIR = f"{WDIR}/analysis/data/vcf/private_variants"
TSV_DIR     = f"{WDIR}/analysis/data/tsv"
SNPEFF_DIR  = f"{WDIR}/resources/snpeff"

# ── Epigenomics paths ──────────────────────────────────────────────────────────
BEDMETHYL_DIR = f"{WDIR}/resources/bedmethyl"
DSS_DIR       = f"{WDIR}/analysis/data/taa_DSS"
DMR_DIR       = f"{WDIR}/analysis/data/taa_DMR"
GFF           = f"{WDIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff"

# ── SIFT4G paths ───────────────────────────────────────────────────────────────
# All SIFT4G tools live inside Docker. Override via --config KEY=value.
SIFT4G_SCRIPTS   = config.get("SIFT4G_SCRIPTS",   "/opt/SIFT4G_Create_Genomic_DB")
SIFT4G_JAR       = config.get("SIFT4G_JAR",       "/opt/SIFT4G_Annotator.jar")
SIFT4G_BUILD_DIR = config.get("SIFT4G_BUILD_DIR", f"{WDIR}/resources/sift4g/sorghum_db")
SIFT4G_DB        = config.get("SIFT4G_DB",        f"{SIFT4G_BUILD_DIR}/NCBIv3")
UNIPROT_FASTA    = config.get("UNIPROT_FASTA",    f"{WDIR}/resources/sift4g/uniprot_sprot.fasta")
SIFT_DIR         = f"{WDIR}/analysis/data/vcf/sift4g"
REF              = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"

SORGHUM_CHROMS = [f"NC_01287{i}.2" for i in range(10)] + ["NC_008360.1", "NC_008602.1"]


rule taa_all:
    """Run all TAA steps: genomics (private variants + SIFT4G) + epigenomics (DMR annotation)."""
    input:
        expand(f"{TSV_DIR}/{{sample}}.private.annotated.tsv", sample=SAMPLES),
        expand(f"{SIFT_DIR}/{{sample}}/{{sample}}.private.annotated_SIFTannotations.txt", sample=SAMPLES),
        expand(f"{SIFT_DIR}/{{sample}}/{{sample}}.private.sift4g.vcf.gz", sample=SAMPLES),
        f"{DMR_DIR}/DMR_annotated.tsv",


rule snpeff_prep:
    """Download prebuilt SnpEff database for Sorghum bicolor (one-time setup)."""
    output:
        f"{SNPEFF_DIR}/data/Sorghum_bicolor/sequence.1.bin",
    log:
        f"{LOG_DIR}/snpeff_prep.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/analysis/scripts/snpeff_prep.sh > {log} 2>&1
        """


rule private_variants:
    """Find variants private to each sample using bcftools isec (-n =1)."""
    input:
        vcfs = expand(f"{VCF_DIR}/{{sample}}.phased.vcf.gz", sample=SAMPLES),
        csis = expand(f"{VCF_DIR}/{{sample}}.phased.vcf.gz.csi", sample=SAMPLES),
    output:
        vcfs = expand(f"{PRIVATE_DIR}/{{sample}}.private.vcf.gz", sample=SAMPLES),
        tbis = expand(f"{PRIVATE_DIR}/{{sample}}.private.vcf.gz.tbi", sample=SAMPLES),
    log:
        f"{LOG_DIR}/private_variants.{TIMESTAMP}.log",
    params:
        out_dir = PRIVATE_DIR,
        samples = " ".join(SAMPLES),
    shell:
        """
        (
            mkdir -p {params.out_dir}
            bcftools isec -n =1 -p {params.out_dir} {input.vcfs}

            # Rename 0000.vcf … to sample-named files (order matches SAMPLES)
            samples=({params.samples})
            for i in "${{!samples[@]}}"; do
                idx=$(printf "%04d" "$i")
                src="{params.out_dir}/${{idx}}.vcf"
                dst="{params.out_dir}/${{samples[$i]}}.private.vcf"
                if [[ -f "$src" ]]; then
                    mv "$src" "$dst"
                    bgzip -f "$dst"
                    bcftools index -t "${{dst}}.gz"
                fi
            done
        ) > {log} 2>&1
        """


rule annotate_vcf:
    """Annotate private variants with SnpEff and index the output."""
    input:
        vcf     = f"{PRIVATE_DIR}/{{sample}}.private.vcf.gz",
        tbi     = f"{PRIVATE_DIR}/{{sample}}.private.vcf.gz.tbi",
        db_flag = f"{SNPEFF_DIR}/data/Sorghum_bicolor/sequence.1.bin",
    output:
        vcf = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz",
        csi = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz.csi",
    log:
        f"{LOG_DIR}/annotate_vcf.{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/analysis/scripts/annotate_vcf.sh \
            {wildcards.sample}.private \
            {input.vcf} \
            {PRIVATE_DIR} \
            > {log} 2>&1
        """


rule vcf_to_tsv:
    """Convert annotated private VCF to TSV for notebook exploration."""
    input:
        vcf = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz",
        csi = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz.csi",
    output:
        f"{TSV_DIR}/{{sample}}.private.annotated.tsv",
    log:
        f"{LOG_DIR}/vcf_to_tsv.{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        python3 {WDIR}/analysis/scripts/annot_single_vcf_to_tsv.py \
            -v {input.vcf} \
            -o {TSV_DIR} \
            > {log} 2>&1
        """


# ── SIFT4G rules ──────────────────────────────────────────────────────────────

rule sift4g_download_uniprot:
    """Download UniProt/Swiss-Prot protein FASTA (one-time, ~820 MB uncompressed)."""
    output:
        UNIPROT_FASTA,
    log:
        f"{LOG_DIR}/sift4g_download_uniprot.{TIMESTAMP}.log",
    shell:
        """
        (
            mkdir -p $(dirname {output})
            wget -q \
                https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
                -O {output}.gz
            gunzip {output}.gz
        ) > {log} 2>&1
        """


rule sift4g_db_setup:
    """Set up SIFT4G input directory structure and write the sorghum config file.

    Splits the genome FASTA into one .fa per main chromosome, compresses the
    NCBI GFF3 into gene-annotation-src/, and writes the SIFT4G config file.
    """
    input:
        fasta   = REF,
        fai     = f"{REF}.fai",
        gff     = GFF,
        uniprot = UNIPROT_FASTA,
    output:
        config   = f"{SIFT4G_BUILD_DIR}/sorghum_bicolor.config",
        sentinel = touch(f"{SIFT4G_BUILD_DIR}/.setup_done"),
    log:
        f"{LOG_DIR}/sift4g_db_setup.{TIMESTAMP}.log",
    params:
        build_dir = SIFT4G_BUILD_DIR,
        chroms    = " ".join(SORGHUM_CHROMS),
    shell:
        """
        (
            mkdir -p {params.build_dir}/chr-src {params.build_dir}/gene-annotation-src

            for chr in {params.chroms}; do
                samtools faidx {input.fasta} "$chr" \
                    > {params.build_dir}/chr-src/"$chr".fa
            done

            # Convert NCBI GFF3 → GTF (gff_gene_format_to_ucsc.pl only reads .gtf/.gtf.gz)
            python3 - {input.gff} \
                {params.build_dir}/gene-annotation-src/sorghum_bicolor.gtf.gz << 'PYEOF'
import sys, re, gzip
feat_map = {{'mRNA': 'transcript', 'exon': 'exon', 'CDS': 'CDS',
             'start_codon': 'start_codon', 'stop_codon': 'stop_codon'}}
def parse_attrs(s):
    return dict(m.groups() for m in re.finditer(r'([^=;\s]+)=([^;]*)', s))
def strip_pfx(s):
    return re.sub(r'^(gene-|rna-|cds-|exon-)', '', s)
gff_in, gtf_out = sys.argv[1], sys.argv[2]
with open(gff_in) as fin, gzip.open(gtf_out, 'wt') as fout:
    for line in fin:
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\\n').split('\\t')
        if len(f) != 9 or f[2] not in feat_map:
            continue
        a = parse_attrs(f[8])
        if f[2] == 'mRNA':
            tid = strip_pfx(a.get('ID', ''))
            gid = strip_pfx(a.get('gene', a.get('Parent', tid)))
        else:
            tid = strip_pfx(a.get('Parent', '').split(',')[0])
            gid = strip_pfx(a.get('gene', tid))
        f[2] = feat_map[f[2]]
        f[8] = 'gene_id "' + gid + '"; transcript_id "' + tid + '";'
        fout.write('\\t'.join(f) + '\\n')
PYEOF

            cat > {output.config} << 'CONFIG'
SIFT4G_PATH=/usr/local/bin/sift4g
PROTEIN_DB={input.uniprot}
PARENT_DIR={params.build_dir}
ORG=Sorghum_bicolor
ORG_VERSION=NCBIv3
GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard
MITO_GENETIC_CODE_TABLE=1
MITO_GENETIC_CODE_TABLENAME=Standard
PLASTID_GENETIC_CODE_TABLE=11
PLASTID_GENETIC_CODE_TABLENAME=Bacterial,_Archaea_and_Plant_Plastid
GENE_DOWNLOAD_DEST=gene-annotation-src
CHR_DOWNLOAD_DEST=chr-src
LOGFILE=Log.txt
ZLOGFILE=Log2.txt
FASTA_DIR=fasta
SUBST_DIR=subst
ALIGN_DIR=SIFT_alignments
SIFT_SCORE_DIR=SIFT_predictions
SINGLE_REC_BY_CHR_DIR=singleRecords
SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores
DBSNP_DIR=dbSNP
FASTA_LOG=fasta.log
INVALID_LOG=invalid.log
PEPTIDE_LOG=peptide.log
ENS_PATTERN=ENS
SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord
CONFIG
        ) > {log} 2>&1
        """


rule sift4g_db_build:
    """Build SIFT4G prediction database for Sorghum bicolor NCBIv3 (one-time, several hours).

    Must cd into SIFT4G_SCRIPTS because make-SIFT-db-all.pl calls all sub-scripts
    by bare filename. The final database lives in resources/sift4g/sorghum_db/NCBIv3/.
    """
    input:
        config   = f"{SIFT4G_BUILD_DIR}/sorghum_bicolor.config",
        sentinel = f"{SIFT4G_BUILD_DIR}/.setup_done",
    output:
        touch(f"{SIFT4G_BUILD_DIR}/.db_built"),
    log:
        f"{LOG_DIR}/sift4g_db_build.{TIMESTAMP}.log",
    params:
        scripts = SIFT4G_SCRIPTS,
    threads: 8
    shell:
        """
        (
            cd {params.scripts}
            perl make-SIFT-db-all.pl -config {input.config}
        ) > {log} 2>&1
        """


rule sift4g_annotate:
    """Chain SIFT4G onto SnpEff-annotated private variants (per sample).

    Decompresses the SnpEff VCF, runs SIFT4G_Annotator.jar, then bgzips and
    indexes the annotated VCF. Output carries both ANN= (SnpEff) and
    SIFT_SCORE/SIFT_PRED in INFO. SIFT score < 0.05 = DAMAGING.

    Override paths at runtime:
        --config SIFT4G_JAR=/path/SIFT4G_Annotator.jar SIFT4G_DB=/path/db
    """
    input:
        vcf     = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz",
        csi     = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz.csi",
        db_flag = f"{SIFT4G_BUILD_DIR}/.db_built",
    output:
        txt = f"{SIFT_DIR}/{{sample}}/{{sample}}.private.annotated_SIFTannotations.txt",
        vcf = f"{SIFT_DIR}/{{sample}}/{{sample}}.private.sift4g.vcf.gz",
        tbi = f"{SIFT_DIR}/{{sample}}/{{sample}}.private.sift4g.vcf.gz.tbi",
    log:
        f"{LOG_DIR}/sift4g_annotate.{{sample}}.{TIMESTAMP}.log",
    params:
        jar = SIFT4G_JAR,
        db  = SIFT4G_DB,
    shell:
        """
        (
            out_dir={SIFT_DIR}/{wildcards.sample}
            stem={wildcards.sample}.private.annotated
            mkdir -p "$out_dir"

            # SIFT4G Annotator requires an uncompressed VCF
            bcftools view {input.vcf} -O v -o "$out_dir/$stem.vcf"

            # Annotate — SIFT_SCORE/SIFT_PRED added to INFO alongside SnpEff ANN=
            java -jar {params.jar} -c \
                -i "$out_dir/$stem.vcf" \
                -d {params.db} \
                -r "$out_dir"

            # SIFT4G names its output VCF alongside the annotations txt.
            # Filename varies by version; find it, compress to a known name.
            sift_vcf=$(ls "$out_dir"/*SIFT*.vcf 2>/dev/null | head -1)
            if [[ -z "$sift_vcf" ]]; then
                # Some versions rewrite the input VCF in the output dir
                sift_vcf="$out_dir/$stem.vcf"
            fi
            bcftools sort "$sift_vcf" | bgzip -c > {output.vcf}
            bcftools index -t {output.vcf}

            # Clean up uncompressed intermediates
            find "$out_dir" -name "*.vcf" -delete
        ) > {log} 2>&1
        """


# ── Epigenomics rules ──────────────────────────────────────────────────────────

rule prepare_dss:
    """Convert filtered bedMethyl files to DSS input format (collapses Watson/Crick CpG pairs)."""
    input:
        expand(f"{BEDMETHYL_DIR}/{{sample}}.filtered.bed", sample=SAMPLES),
    output:
        dss     = expand(f"{DSS_DIR}/{{sample}}.5mC.dss.txt", sample=SAMPLES),
        summary = f"{DSS_DIR}/conversion_summary.tsv",
    log:
        f"{LOG_DIR}/prepare_dss.{TIMESTAMP}.log",
    shell:
        """
        python3 {WDIR}/analysis/scripts/prepare_taa_dss.py > {log} 2>&1
        """


rule run_dss_pair:
    """Run DSS DML/DMR calling for one SBC10 vs other pair (5mC only).

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
        f"{LOG_DIR}/run_dss_{{pair}}.{TIMESTAMP}.log",
    threads: 32
    shell:
        """
        cd {WDIR} && Rscript analysis/scripts/run_dss_dmr_taa.R {wildcards.pair} > {log} 2>&1
        """


rule summarise_dmr:
    """Combine per-pair DMR tables into DMR_summary.tsv and DMR_all_combined.tsv."""
    input:
        expand(f"{DMR_DIR}/{{pair}}.5mC.DMR.tsv", pair=PAIRS),
    output:
        summary  = f"{DMR_DIR}/DMR_summary.tsv",
        combined = f"{DMR_DIR}/DMR_all_combined.tsv",
    log:
        f"{LOG_DIR}/summarise_dmr.{TIMESTAMP}.log",
    shell:
        """
        cd {WDIR} && Rscript analysis/scripts/summarise_dss_dmr_taa.R > {log} 2>&1
        """


rule annotate_dmr:
    """Annotate DMRs with strand-aware 2 kbp promoter regions across all genes."""
    input:
        dmr = f"{DMR_DIR}/DMR_all_combined.tsv",
        gff = GFF,
    output:
        f"{DMR_DIR}/DMR_annotated.tsv",
    log:
        f"{LOG_DIR}/annotate_dmr.{TIMESTAMP}.log",
    shell:
        """
        python3 {WDIR}/analysis/scripts/annotate_DMR.py \
            --dmr    {input.dmr} \
            --gff    {input.gff} \
            --outdir {DMR_DIR} \
            --all-genes \
            > {log} 2>&1
        """
