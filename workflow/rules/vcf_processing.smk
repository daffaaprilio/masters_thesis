from datetime import datetime

WDIR       = "/home/daffa/Work/2026/thesis"
SAMPLES    = ["SBC10", "SBC11", "SBC4", "SBC23"]
PAIR_FOCUS = ["SBC10", "SBC11"]   # the biologically significant pair (high vs low TAA)
TIMESTAMP  = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR    = f"{WDIR}/workflow/logs/vcf_processing"

VCF_DIR = f"{WDIR}/results/vcf"
OUT_DIR = f"{WDIR}/results/vcf_processing"
BAM_DIR = f"{WDIR}/resources/align_bam_sample"
REF     = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
GFF_RAW = f"{WDIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.gz"
GFF_BGZ = f"{WDIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.bgz"

# 10 nuclear chromosomes (sorted by size); excludes unplaced scaffolds and organelles
CHROMS = [
    "NC_012870.2", "NC_012871.2", "NC_012872.2", "NC_012873.2", "NC_012874.2",
    "NC_012875.2", "NC_012876.2", "NC_012877.2", "NC_012878.2", "NC_012879.2",
]

# Filtering thresholds # — override at runtime with --config key=value
QUAL_MIN = 20   # config.get("qual_min", 20)
DP_MIN   = 10   # config.get("dp_min",   10)
DP_MAX   = 100  # config.get("dp_max",  100)

# All isec output labels consumed by annotate_vcf
DIFF_LABELS = (
    [f"{s}_private" for s in SAMPLES]   # private to each of the 4 samples
    + ["all_shared"]                    # present in all 4
    + ["pair_shared"]                   # shared between SBC10 and SBC11 only
)


rule all:
    input:
        expand(f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz",  sample=SAMPLES),
        expand(f"{OUT_DIR}/phased/{{sample}}.phased.vcf.gz",      sample=SAMPLES),
        f"{OUT_DIR}/merged/all_samples.merged.vcf.gz",
        expand(f"{OUT_DIR}/isec/{{sample}}_private.vcf.gz",       sample=SAMPLES),
        f"{OUT_DIR}/isec/all_shared.vcf.gz",
        f"{OUT_DIR}/isec/pair_shared.vcf.gz",
        expand(f"{OUT_DIR}/annotated/{{label}}.annotated.vcf.gz", label=DIFF_LABELS),


# ── GFF preparation ───────────────────────────────────────────────────────────

rule prepare_gff:
    """Re-compress NCBI GFF3 with bgzip and create tabix index for bcftools csq."""
    input:
        GFF_RAW,
    output:
        bgz = GFF_BGZ,
        tbi = GFF_BGZ + ".tbi",
    log:
        f"{LOG_DIR}/prepare_gff.{TIMESTAMP}.log",
    shell:
        """
        (
            TMPGFF=$(mktemp)
            zcat {input} > "$TMPGFF"
            (grep "^#" "$TMPGFF"; grep -v "^#" "$TMPGFF" | grep -v "trans-splicing" | sort -k1,1 -k4,4n) | bgzip -c > {output.bgz}
            rm -f "$TMPGFF"
            tabix -p gff {output.bgz}
        ) > {log} 2>&1
        """

# add SAMPLE information to the VCF file

rule reheader_vcf:
    """Rename the generic SAMPLE column header to the actual sample name."""
    input:
        vcf = f"{VCF_DIR}/{{sample}}.vcf.gz",
        csi = f"{VCF_DIR}/{{sample}}.vcf.gz.csi",
    output:
        vcf = f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz",
        csi = f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz.csi",
    log:
        f"{LOG_DIR}/reheader_vcf/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            echo "{wildcards.sample}" | bcftools reheader -s /dev/stdin {input.vcf} -o {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


# ── Per-sample filtering ──────────────────────────────────────────────────────

rule filter_vcf:
    """Keep PASS variants within configurable QUAL and DP bounds."""
    input:
        vcf = f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz",
        csi = f"{OUT_DIR}/reheadered/{{sample}}.reheadered.vcf.gz.csi",
    output:
        vcf = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz",
        csi = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz.csi",
    log:
        f"{LOG_DIR}/filter_vcf/{{sample}}.{TIMESTAMP}.log",
    params:
        expr    = f"QUAL<{QUAL_MIN} || FORMAT/DP<{DP_MIN} || FORMAT/DP>{DP_MAX}",
        regions = ",".join(CHROMS),
    shell:
        """
        (
            bcftools view -f PASS -r {params.regions} {input.vcf} \
                | bcftools filter -e '{params.expr}' -O z -o {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


# ── Haplotype phasing ─────────────────────────────────────────────────────────

rule phase_vcf:
    """Read-backed haplotype phasing with WhatsHap per sample."""
    input:
        vcf = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz",
        csi = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz.csi",
        bam = f"{BAM_DIR}/{{sample}}.bam",
        bai = f"{BAM_DIR}/{{sample}}.bam.bai",
        ref = REF,
    output:
        vcf = f"{OUT_DIR}/phased/{{sample}}.phased.vcf.gz",
        csi = f"{OUT_DIR}/phased/{{sample}}.phased.vcf.gz.csi",
    log:
        f"{LOG_DIR}/phase_vcf/{{sample}}.{TIMESTAMP}.log",
    shell:
        """
        (
            whatshap phase \
                --reference {input.ref} \
                --ignore-read-groups \
                --output /dev/stdout \
                {input.vcf} {input.bam} \
                | bgzip -c > {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


# ── Multi-sample merge ────────────────────────────────────────────────────────

rule merge_vcf:
    """Merge per-sample filtered VCFs into one multi-sample VCF."""
    input:
        vcfs = expand(f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz", sample=SAMPLES),
        csis = expand(f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz.csi", sample=SAMPLES),
    output:
        vcf = f"{OUT_DIR}/merged/all_samples.merged.vcf.gz",
        csi = f"{OUT_DIR}/merged/all_samples.merged.vcf.gz.csi",
    log:
        f"{LOG_DIR}/merge_vcf.{TIMESTAMP}.log",
    shell:
        """
        (
            bcftools merge --missing-to-ref {input.vcfs} -O z -o {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


# ── Differential variants (isec) ─────────────────────────────────────────────
#
# isec_private    — all 4 samples, -n =1
#   0000→SBC10_private, 0001→SBC11_private, 0002→SBC4_private, 0003→SBC23_private
#
# isec_all_shared — all 4 samples, -n =4
#   0000 used as representative (same loci across all outputs)
#
# isec_pair_shared — SBC10 + SBC11 only (2-sample default)
#   0002 = shared loci using SBC10 calls

rule isec_private:
    """Variants present in exactly one sample — truly private to each of the 4 samples."""
    input:
        vcfs = expand(f"{OUT_DIR}/phased/{{sample}}.phased.vcf.gz", sample=SAMPLES),
        csis = expand(f"{OUT_DIR}/phased/{{sample}}.phased.vcf.gz.csi", sample=SAMPLES),
    output:
        sbc10     = f"{OUT_DIR}/isec/SBC10_private.vcf.gz",
        sbc10_csi = f"{OUT_DIR}/isec/SBC10_private.vcf.gz.csi",
        sbc11     = f"{OUT_DIR}/isec/SBC11_private.vcf.gz",
        sbc11_csi = f"{OUT_DIR}/isec/SBC11_private.vcf.gz.csi",
        sbc4      = f"{OUT_DIR}/isec/SBC4_private.vcf.gz",
        sbc4_csi  = f"{OUT_DIR}/isec/SBC4_private.vcf.gz.csi",
        sbc23     = f"{OUT_DIR}/isec/SBC23_private.vcf.gz",
        sbc23_csi = f"{OUT_DIR}/isec/SBC23_private.vcf.gz.csi",
    log:
        f"{LOG_DIR}/isec_private.{TIMESTAMP}.log",
    params:
        raw_dir = f"{OUT_DIR}/isec/raw_private",
    shell:
        """
        (
            mkdir -p {params.raw_dir}
            bcftools isec {input.vcfs} -p {params.raw_dir} -O z -n =1
            cp {params.raw_dir}/0000.vcf.gz {output.sbc10}
            cp {params.raw_dir}/0001.vcf.gz {output.sbc11}
            cp {params.raw_dir}/0002.vcf.gz {output.sbc4}
            cp {params.raw_dir}/0003.vcf.gz {output.sbc23}
            bcftools index {output.sbc10}
            bcftools index {output.sbc11}
            bcftools index {output.sbc4}
            bcftools index {output.sbc23}
        ) > {log} 2>&1
        """


rule isec_all_shared:
    """Variants present in all 4 samples."""
    input:
        vcfs = expand(f"{OUT_DIR}/phased/{{sample}}.phased.vcf.gz", sample=SAMPLES),
        csis = expand(f"{OUT_DIR}/phased/{{sample}}.phased.vcf.gz.csi", sample=SAMPLES),
    output:
        vcf = f"{OUT_DIR}/isec/all_shared.vcf.gz",
        csi = f"{OUT_DIR}/isec/all_shared.vcf.gz.csi",
    log:
        f"{LOG_DIR}/isec_all_shared.{TIMESTAMP}.log",
    params:
        raw_dir    = f"{OUT_DIR}/isec/raw_all_shared",
        nfiles_arg = f"={len(SAMPLES)}",   # resolves to "=4"
    shell:
        """
        (
            mkdir -p {params.raw_dir}
            bcftools isec {input.vcfs} -p {params.raw_dir} -O z -n {params.nfiles_arg}
            cp {params.raw_dir}/0000.vcf.gz {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


rule isec_pair_shared:
    """Variants shared between SBC10 and SBC11, regardless of SBC4/SBC23 status."""
    input:
        vcfs = expand(f"{OUT_DIR}/phased/{{sample}}.phased.vcf.gz", sample=PAIR_FOCUS),
        csis = expand(f"{OUT_DIR}/phased/{{sample}}.phased.vcf.gz.csi", sample=PAIR_FOCUS),
    output:
        vcf = f"{OUT_DIR}/isec/pair_shared.vcf.gz",
        csi = f"{OUT_DIR}/isec/pair_shared.vcf.gz.csi",
    log:
        f"{LOG_DIR}/isec_pair_shared.{TIMESTAMP}.log",
    params:
        raw_dir = f"{OUT_DIR}/isec/raw_pair_shared",
    shell:
        """
        (
            mkdir -p {params.raw_dir}
            bcftools isec {input.vcfs} -p {params.raw_dir} -O z
            cp {params.raw_dir}/0002.vcf.gz {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


# ── Gene feature annotation ───────────────────────────────────────────────────

rule annotate_vcf:
    """
    Annotate variants with gene features and predicted consequences via bcftools csq.
    Adds BCSQ INFO tag: gene|transcript|biotype|strand|consequence|amino-acid-change.
    """
    input:
        vcf = f"{OUT_DIR}/isec/{{label}}.vcf.gz",
        csi = f"{OUT_DIR}/isec/{{label}}.vcf.gz.csi",
        gff = GFF_BGZ,
        tbi = GFF_BGZ + ".tbi",
        ref = REF,
    output:
        vcf = f"{OUT_DIR}/annotated/{{label}}.annotated.vcf.gz",
        csi = f"{OUT_DIR}/annotated/{{label}}.annotated.vcf.gz.csi",
    log:
        f"{LOG_DIR}/annotate_vcf/{{label}}.{TIMESTAMP}.log",
    shell:
        """
        (
            bcftools csq \
                --phase a \
                --ncsq 30 \
                -f {input.ref} \
                -g {input.gff} \
                {input.vcf} -O z -o {output.vcf}
            bcftools index {output.vcf}
        ) > {log} 2>&1
        """


# ── Cleanup ───────────────────────────────────────────────────────────────────

rule clean:
    shell:
        """
        rm -rf {OUT_DIR}/reheadered {OUT_DIR}/filtered {OUT_DIR}/phased {OUT_DIR}/merged {OUT_DIR}/isec {OUT_DIR}/annotated
        rm -f {GFF_BGZ} {GFF_BGZ}.tbi
        """
