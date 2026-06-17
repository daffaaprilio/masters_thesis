# SIFT4G annotation rules:
# build the Sorghum bicolor prediction database (one-time) → annotate per-sample variants.

# if you want to redo the process, remove the sentinel file first (f"{SIFT4G_BUILD_DIR}/.setup_done")

GTF              = f"{WDIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gtf"
UNIREF90         = f"{WDIR}/resources/UNIPROT_FTP/uniref90.fasta"
SIFT4G_BUILD_DIR = f"{WDIR}/resources/sift4g/sorghum_db"
SIFT4G_DB        = f"{SIFT4G_BUILD_DIR}/NCBIv3"
SORGHUM_CHROMS   = [f"NC_01287{i}.2" for i in range(10)] + ["NC_008360.1", "NC_008602.1"]

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
        f"{WDIR}/workflow/logs/vcf_annotation/prepare_sift_db.{TIMESTAMP}.log",
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
        f"{WDIR}/workflow/logs/vcf_annotation/build_sift_db.{TIMESTAMP}.log",
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

    Chromosome renaming: SnpEff outputs simple names (1–10); the SIFT4G DB was
    built with NCBI IDs (NC_012870.2 etc.).  We reverse-map before annotating
    so the VCF chromosomes match the DB .regions files.
    """
    input:
        vcf      = f"{SNPEFF_ANNOT_DIR}/{{sample}}.private.annotated.vcf.gz",
        csi      = f"{SNPEFF_ANNOT_DIR}/{{sample}}.private.annotated.vcf.gz.csi",
        db_flag  = f"{SIFT4G_BUILD_DIR}/.db_built",
        synonyms = f"{WDIR}/workflow/scripts/synonyms.txt",
    output:
        xls = f"{SIFT_DIR}/{{sample}}.private.annotated_SIFTannotations.xls",
        vcf = f"{SIFT_DIR}/{{sample}}.private.sift4g.vcf.gz",
        tbi = f"{SIFT_DIR}/{{sample}}.private.sift4g.vcf.gz.tbi",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/annotate_sift.{{sample}}.{TIMESTAMP}.log",
    params:
        jar = "/opt/SIFT4G_Annotator.jar",
        db  = SIFT4G_DB,
    shell:
        """
        (
            out_dir={SIFT_DIR}
            tmp_dir={SIFT_DIR}/{wildcards.sample}_tmp
            stem={wildcards.sample}.private.annotated
            mkdir -p "$out_dir" "$tmp_dir"

            # Reverse the synonym map: SnpEff names (col 2) → NCBI IDs (col 1)
            awk '{{print $2"\t"$1}}' {input.synonyms} > "$tmp_dir/chrom_rename_rev.txt"

            # Rename chromosomes back to NCBI IDs so they match the SIFT DB
            bcftools annotate --rename-chrs "$tmp_dir/chrom_rename_rev.txt" \
                {input.vcf} -O v -o "$tmp_dir/$stem.vcf"

            java -jar {params.jar} -c \
                -i "$tmp_dir/$stem.vcf" \
                -d {params.db} \
                -r "$out_dir"

            sift_vcf=$(ls "$out_dir"/${{stem}}_SIFTpredictions*.vcf 2>/dev/null | head -1)
            if [[ -z "$sift_vcf" ]]; then
                sift_vcf=$(ls "$out_dir"/${{stem}}*SIFT*.vcf 2>/dev/null | head -1)
            fi
            bcftools sort "$sift_vcf" | bgzip -c > {output.vcf}
            bcftools index -t {output.vcf}

            rm -f "$sift_vcf"
            rm -rf "$tmp_dir"
        ) > {log} 2>&1
        """
