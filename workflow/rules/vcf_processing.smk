# VCF post-processing: filter → merge → differential → annotate
#
# Run standalone:
#   snakemake --snakefile workflow/rules/vcf_processing.smk [--config qual_min=20 dp_min=10 dp_max=100]
#
# Annotation uses bcftools csq (consequence calling) rather than bcftools annotate
# because the only annotation resource is a GFF3 file.  csq reads GFF3 natively
# and writes a BCSQ INFO tag with gene name + predicted consequence.

WDIR    = "/home/daffa/Work/2026/thesis"
SAMPLES = ["SBC10", "SBC11"]   # SBC4, SBC23 oncoming — extend this list when ready

VCF_DIR = f"{WDIR}/results/vcf"
OUT_DIR = f"{WDIR}/results/vcf_processing"
REF     = f"{WDIR}/resources/ref/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna"
GFF_RAW = f"{WDIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.gz"
GFF_BGZ = f"{WDIR}/resources/annot/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.bgz"

# Filtering thresholds — override at runtime with --config key=value
QUAL_MIN = 20   # config.get("qual_min", 20)
DP_MIN   = 10   # config.get("dp_min",   10)
DP_MAX   = 100  # config.get("dp_max",  100)

# Differential sets produced by bcftools isec (2-sample case)
DIFF_LABELS = ["SBC10_private", "SBC11_private", "shared"]


rule all:
    input:
        expand(f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz",      sample=SAMPLES),
        f"{OUT_DIR}/merged/all_samples.merged.vcf.gz",
        expand(f"{OUT_DIR}/isec/{{label}}.vcf.gz",                    label=DIFF_LABELS),
        expand(f"{OUT_DIR}/annotated/{{label}}.annotated.vcf.gz",     label=DIFF_LABELS),


# ── GFF preparation ───────────────────────────────────────────────────────────

rule prepare_gff:
    """Re-compress NCBI GFF3 with bgzip and create tabix index for bcftools csq."""
    input:
        GFF_RAW,
    output:
        bgz = GFF_BGZ,
        tbi = GFF_BGZ + ".tbi",
    shell:
        """
        zcat {input} | bgzip -c > {output.bgz}
        tabix -p gff {output.bgz}
        """


# ── Per-sample filtering ──────────────────────────────────────────────────────

rule filter_vcf:
    """Keep PASS variants within configurable QUAL and DP bounds."""
    input:
        vcf = f"{VCF_DIR}/{{sample}}.vcf.gz",
        csi = f"{VCF_DIR}/{{sample}}.vcf.gz.csi",
    output:
        vcf = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz",
        csi = f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz.csi",
    params:
        expr = f"QUAL<{QUAL_MIN} || FORMAT/DP<{DP_MIN} || FORMAT/DP>{DP_MAX}",
    shell:
        """
        bcftools view -f PASS {input.vcf} \
            | bcftools filter -e '{params.expr}' -O z -o {output.vcf}
        bcftools index {output.vcf}
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
    shell:
        """
        bcftools merge --missing-to-ref {input.vcfs} -O z -o {output.vcf}
        bcftools index {output.vcf}
        """


# ── Differential variants (isec) ─────────────────────────────────────────────
#
# bcftools isec on 2 files writes:
#   0000.vcf.gz  → private to SBC10 (high-TAA candidate variants)
#   0001.vcf.gz  → private to SBC11 (low-TAA candidate variants)
#   0002.vcf.gz  → positions shared, SBC10 calls  ─┐ same loci,
#   0003.vcf.gz  → positions shared, SBC11 calls  ─┘ SBC10 used as representative

rule isec_variants:
    input:
        vcfs = expand(f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz", sample=SAMPLES),
        csis = expand(f"{OUT_DIR}/filtered/{{sample}}.filtered.vcf.gz.csi", sample=SAMPLES),
    output:
        sbc10_private = f"{OUT_DIR}/isec/SBC10_private.vcf.gz",
        sbc10_priv_csi= f"{OUT_DIR}/isec/SBC10_private.vcf.gz.csi",
        sbc11_private = f"{OUT_DIR}/isec/SBC11_private.vcf.gz",
        sbc11_priv_csi= f"{OUT_DIR}/isec/SBC11_private.vcf.gz.csi",
        shared        = f"{OUT_DIR}/isec/shared.vcf.gz",
        shared_csi    = f"{OUT_DIR}/isec/shared.vcf.gz.csi",
    params:
        raw_dir = f"{OUT_DIR}/isec/raw",
    shell:
        """
        mkdir -p {params.raw_dir}
        bcftools isec {input.vcfs} -p {params.raw_dir} -O z
        cp {params.raw_dir}/0000.vcf.gz {output.sbc10_private}
        cp {params.raw_dir}/0001.vcf.gz {output.sbc11_private}
        cp {params.raw_dir}/0002.vcf.gz {output.shared}
        bcftools index {output.sbc10_private}
        bcftools index {output.sbc11_private}
        bcftools index {output.shared}
        """


# ── Gene feature annotation ───────────────────────────────────────────────────

rule annotate_vcf:
    """
    Annotate variants with gene features and predicted consequences via bcftools csq.
    Adds BCSQ INFO tag: gene|transcript|biotype|strand|consequence|amino-acid-change.
    --local-csq handles unphased/partially-phased ONT calls.
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
    shell:
        """
        bcftools csq \
            -f {input.ref} \
            -g {input.gff} \
            --local-csq \
            {input.vcf} -O z -o {output.vcf}
        bcftools index {output.vcf}
        """
