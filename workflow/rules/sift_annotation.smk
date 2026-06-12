# SIFT4G annotation rules:
# build the Sorghum bicolor prediction database (one-time) → annotate per-sample variants.

rule prepare_sift_db:
    """Prepare SIFT4G input directory structure and config for Sorghum bicolor NCBIv3.

    - Splits the genome FASTA into one .fa per chromosome under chr-src/
    - Gzips the GTF annotation into gene-annotation-src/
    - Writes the config file consumed by build_sift_db
    """
    input:
        fasta = REF,
        fai   = f"{REF}.fai",
        gtf   = GTF,
        prot  = UNIREF90,
    output:
        config   = f"{SIFT4G_BUILD_DIR}/sorghum_bicolor.config",
        sentinel = f"{SIFT4G_BUILD_DIR}/.setup_done",
    log:
        f"{LOG_DIR}/prepare_sift_db.{TIMESTAMP}.log",
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

            cp {input.gtf} {params.build_dir}/gene-annotation-src/sorghum_bicolor.gtf
            gzip -f {params.build_dir}/gene-annotation-src/sorghum_bicolor.gtf

            cat > {output.config} << 'CONFIG'
SIFT4G_PATH=/usr/local/bin/sift4g
PROTEIN_DB={input.prot}
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
            touch {output.sentinel}
        ) > {log} 2>&1
        """


rule build_sift_db:
    """Build the SIFT4G prediction database (one-time, several hours).

    Must cd into the scripts dir because make-SIFT-db-all.pl calls sub-scripts
    by bare filename. The finished DB lands at resources/sift4g/sorghum_db/NCBIv3/.
    """
    input:
        config   = f"{SIFT4G_BUILD_DIR}/sorghum_bicolor.config",
        sentinel = f"{SIFT4G_BUILD_DIR}/.setup_done",
    output:
        f"{SIFT4G_BUILD_DIR}/.db_built",
    log:
        f"{LOG_DIR}/build_sift_db.{TIMESTAMP}.log",
    threads: 8
    shell:
        """
        (
            cd /opt/SIFT4G_Create_Genomic_DB
            perl make-SIFT-db-all.pl -config {input.config}
            touch {output}
        ) > {log} 2>&1
        """


rule annotate_sift:
    """Annotate SnpEff-private variants with SIFT4G functional scores (per sample).

    SIFT4G_Annotator.jar outputs {stem}_SIFTannotations.xls alongside a
    SIFT-annotated VCF. The VCF is sorted, bgzipped, and indexed.
    SIFT score < 0.05 = DAMAGING.
    """
    input:
        vcf     = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz",
        csi     = f"{PRIVATE_DIR}/{{sample}}.private.annotated.vcf.gz.csi",
        db_flag = f"{SIFT4G_BUILD_DIR}/.db_built",
    output:
        xls = f"{SIFT_DIR}/{{sample}}/{{sample}}.private.annotated_SIFTannotations.xls",
        vcf = f"{SIFT_DIR}/{{sample}}/{{sample}}.private.sift4g.vcf.gz",
        tbi = f"{SIFT_DIR}/{{sample}}/{{sample}}.private.sift4g.vcf.gz.tbi",
    log:
        f"{LOG_DIR}/annotate_sift.{{sample}}.{TIMESTAMP}.log",
    params:
        jar = "/opt/SIFT4G_Annotator.jar",
        db  = SIFT4G_DB,
    shell:
        """
        (
            out_dir={SIFT_DIR}/{wildcards.sample}
            stem={wildcards.sample}.private.annotated
            mkdir -p "$out_dir"

            bcftools view {input.vcf} -O v -o "$out_dir/$stem.vcf"

            java -jar {params.jar} -c \
                -i "$out_dir/$stem.vcf" \
                -d {params.db} \
                -r "$out_dir"

            sift_vcf=$(ls "$out_dir"/*_SIFTpredictions*.vcf 2>/dev/null | head -1)
            if [[ -z "$sift_vcf" ]]; then
                sift_vcf=$(ls "$out_dir"/*SIFT*.vcf 2>/dev/null | head -1)
            fi
            bcftools sort "$sift_vcf" | bgzip -c > {output.vcf}
            bcftools index -t {output.vcf}

            find "$out_dir" -name "*.vcf" -delete
        ) > {log} 2>&1
        """
