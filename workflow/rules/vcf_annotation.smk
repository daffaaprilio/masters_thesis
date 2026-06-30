#   vcf_annotation
#   (SNP calling) Annotate merged VCF with Snpeff, SIFT
#   (SV calling) Annotate combined SV VCF with SnpEff only,
#   Then, combine these annotated VCF into one VCF

SNPEFF_DIR = f"{WDIR}/resources/snpeff"

rule snpeff_prep:
    """Build custom SnpEff database for Sorghum bicolor from the NCBI GTF (one-time setup)."""
    output:
        f"{SNPEFF_DIR}/data/Sorghum_bicolor_NCBIv3/sequence.NC_012870.2.bin",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/snpeff_prep.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/workflow/scripts/snpeff_prep.sh > {log} 2>&1
        """


rule merge_samples:
    """Merge the four phased per-sample VCFs into one multi-sample VCF in canonical
    SAMPLES order, splitting multiallelics + left-aligning indels (bcftools norm) so
    each ALT is its own record. This merged VCF is the pipeline's single source of
    truth; any group/contrast slicing is a downstream `bcftools view -i` step, not
    part of the pipeline."""
    input:
        vcfs = expand(f"{PHASED_VCF_DIR}/{{sample}}.phased.vcf.gz", sample=SAMPLES),
        csis = expand(f"{PHASED_VCF_DIR}/{{sample}}.phased.vcf.gz.csi", sample=SAMPLES),
        ref  = REF,
    output:
        vcf = MERGED_PHASED_VCF,
        tbi = f"{MERGED_PHASED_VCF}.tbi",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/merge_samples.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/workflow/scripts/merge_vcf.sh > {log} 2>&1
        """


rule annotate_merged:
    """Annotate the merged multi-sample VCF with SnpEff — the pipeline's only SnpEff
    run. SnpEff is per-allele/sample-agnostic, so the multi-sample annotated VCF
    carries each variant's consequence once, with all four genotype columns intact."""
    input:
        vcf     = MERGED_PHASED_VCF,
        tbi     = f"{MERGED_PHASED_VCF}.tbi",
        db_flag = f"{SNPEFF_DIR}/data/Sorghum_bicolor_NCBIv3/sequence.NC_012870.2.bin",
    output:
        vcf = f"{SNPEFF_ANNOT_DIR}/merged.annotated.vcf.gz",
        csi = f"{SNPEFF_ANNOT_DIR}/merged.annotated.vcf.gz.csi",
        csv = f"{SNPEFF_ANNOT_DIR}/merged.stats.csv",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/annotate_merged.{TIMESTAMP}.log",
    shell:
        """
        bash {WDIR}/workflow/scripts/annotate_vcf.sh \
            merged \
            {input.vcf} \
            {SNPEFF_ANNOT_DIR} \
            > {log} 2>&1
        """

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


rule sift_merged:
    """Annotate the merged SnpEff VCF with SIFT4G functional scores — the pipeline's
    only SIFT4G run. SIFT4G is per-allele, so the multi-sample VCF is scored once.

    SIFT4G_Annotator.jar outputs {stem}_SIFTannotations.xls alongside a
    SIFT-annotated VCF. The VCF is sorted, bgzipped, and indexed.
    SIFT score < 0.05 = DAMAGING.

    SnpEff uses a custom DB built from the same NCBI GTF as SIFT4G, so both tools
    share NCBI chromosome IDs (NC_012870.2 etc.) — no rename step needed.
    """
    input:
        vcf      = f"{SNPEFF_ANNOT_DIR}/merged.annotated.vcf.gz",
        csi      = f"{SNPEFF_ANNOT_DIR}/merged.annotated.vcf.gz.csi",
        db_flag  = f"{SIFT4G_BUILD_DIR}/.db_built",
    output:
        xls = f"{SIFT_DIR}/merged.annotated_SIFTannotations.xls",
        vcf = f"{SIFT_DIR}/merged.sift4g.vcf.gz",
        tbi = f"{SIFT_DIR}/merged.sift4g.vcf.gz.tbi",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/sift_merged.{TIMESTAMP}.log",
    params:
        jar = "/opt/SIFT4G_Annotator.jar",
        db  = SIFT4G_DB,
    shell:
        """
        (
            out_dir={SIFT_DIR}
            tmp_dir={SIFT_DIR}/merged_tmp
            stem=merged.annotated
            mkdir -p "$out_dir" "$tmp_dir"

            # Decompress for SIFT4G_Annotator (requires uncompressed VCF)
            bcftools view {input.vcf} -O v -o "$tmp_dir/$stem.vcf"

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


# ── SV annotation + combine ────────────────────────────────────────────────────

rule annotate_sv:
    """SnpEff-annotate the combined multi-sample SV VCF — SnpEff only (SIFT4G is not
    meaningful for structural variants). Reuses annotate_vcf.sh + the same SnpEff DB as
    the SNP track, so SV ANN gene IDs (LOC*) line up with the SNP annotations. SnpEff
    classifies each SV's per-gene consequence (transcript_ablation, exon_loss,
    feature_fusion for BND, duplication/inversion for DUP/INV, …)."""
    input:
        vcf     = f"{SV_DIR}/combined.sniffles.vcf.gz",
        tbi     = f"{SV_DIR}/combined.sniffles.vcf.gz.tbi",
        db_flag = f"{SNPEFF_DIR}/data/Sorghum_bicolor_NCBIv3/sequence.NC_012870.2.bin",
    output:
        vcf = f"{SV_SNPEFF_DIR}/combined.annotated.vcf.gz",
        csi = f"{SV_SNPEFF_DIR}/combined.annotated.vcf.gz.csi",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/annotate_sv.{TIMESTAMP}.log",
    shell:
        """
        mkdir -p {SV_SNPEFF_DIR} $(dirname {log})
        bash {WDIR}/workflow/scripts/annotate_vcf.sh \
            combined \
            {input.vcf} \
            {SV_SNPEFF_DIR} \
            > {log} 2>&1
        """


rule concat_vcf:
    """Combine the two final annotated VCFs into one multi-sample VCF carrying both
    variant classes: the SNP track (SnpEff+SIFT, merged.sift4g) and the SV track
    (SnpEff, combined). Both have the same four samples in canonical order, so
    `bcftools concat -a` stacks their records; the result is sorted, bgzipped, indexed."""
    input:
        snp = f"{SIFT_DIR}/merged.sift4g.vcf.gz",
        sv  = f"{SV_SNPEFF_DIR}/combined.annotated.vcf.gz",
    output:
        vcf = f"{COMBINED_DIR}/all.annotated.vcf.gz",
        tbi = f"{COMBINED_DIR}/all.annotated.vcf.gz.tbi",
    log:
        f"{WDIR}/workflow/logs/vcf_annotation/concat_vcf.{TIMESTAMP}.log",
    shell:
        """
        (
            mkdir -p {COMBINED_DIR}
            bcftools concat -a {input.snp} {input.sv} \
                | bcftools sort -Oz -o {output.vcf}
            bcftools index -t {output.vcf}
        ) > {log} 2>&1
        """


rule annotate_all:
    """Genomic-layer deliverable: SNP (SnpEff+SIFT) and SV (SnpEff) annotations
    concatenated into one multi-sample VCF (results/combined/all.annotated.vcf.gz)."""
    input:
        f"{COMBINED_DIR}/all.annotated.vcf.gz",
        
