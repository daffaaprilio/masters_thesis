# ─────────────────────────────────────────────
# Stage 1: extract Clair3 model files
# ─────────────────────────────────────────────
FROM --platform=linux/amd64 docker.io/hkubal/clair3:latest AS clair3-source

# ─────────────────────────────────────────────
# Stage 2: build sift4g C++ binary
# ─────────────────────────────────────────────
FROM --platform=linux/amd64 ubuntu:22.04 AS sift4g-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    gcc-10 \
    g++-10 \
    make \
    ca-certificates \
 && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 \
                         --slave   /usr/bin/g++ g++ /usr/bin/g++-10 \
 && rm -rf /var/lib/apt/lists/*

RUN git clone --recursive https://github.com/rvaser/sift4g.git /sift4g
WORKDIR /sift4g
RUN make

# ─────────────────────────────────────────────
# Stage 3: main thesis-tools image
# ─────────────────────────────────────────────
FROM --platform=linux/amd64 docker.io/condaforge/miniforge3:latest

LABEL org.opencontainers.image.title="thesis-tools"
LABEL org.opencontainers.image.description="ONT sorghum bioinformatics pipeline (minimap2, samtools, bcftools, modkit, whatshap, snpeff, bedtools, clair3, snakemake, sift4g)"

# ── System deps for SIFT4G_Create_Genomic_DB (Perl/BioPerl) + Java ──────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    perl \
    python3 \
    wget \
    openjdk-11-jdk \
 && rm -rf /var/lib/apt/lists/*

# ── conda/bioconda tools ─────────────────────────────────────────────────────
COPY docker/environment.yml /tmp/environment.yml
RUN mamba env update -n base -f /tmp/environment.yml && \
    conda clean -afy && \
    cpanm --quiet --notest Switch

# ── Clair3 model files ───────────────────────────────────────────────────────
COPY --from=clair3-source /opt/models /opt/models

# ── sift4g binary (compiled in Stage 2) ─────────────────────────────────────
COPY --from=sift4g-builder /sift4g/bin/sift4g /usr/local/bin/sift4g
RUN chmod +x /usr/local/bin/sift4g

# ── SIFT4G_Create_Genomic_DB Perl scripts ───────────────────────────────────
RUN git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git \
        /opt/SIFT4G_Create_Genomic_DB \
 && chmod 775 /opt/SIFT4G_Create_Genomic_DB/common-utils.pl

# ── SIFT4G_Annotator.jar ─────────────────────────────────────────────────────
RUN wget -q https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar \
        -O /opt/SIFT4G_Annotator.jar \
 && chmod 644 /opt/SIFT4G_Annotator.jar

# ── Environment variable for convenience ────────────────────────────────────
ENV SIFT4G_DB_SCRIPTS=/opt/SIFT4G_Create_Genomic_DB
ENV SIFT4G_ANNOTATOR=/opt/SIFT4G_Annotator.jar

# ── Verify all tools ─────────────────────────────────────────────────────────
RUN minimap2 --version && \
    samtools --version | head -1 && \
    bcftools --version | head -1 && \
    modkit --version && \
    whatshap --version && \
    bedtools --version | head -1 && \
    run_clair3.sh --help > /dev/null 2>&1 && \
    snakemake --version && \
    Rscript -e "library(DSS); library(data.table); cat('R packages OK\n')" && \
    sift4g --version && \
    java -jar /opt/SIFT4G_Annotator.jar --help > /dev/null 2>&1 || true && \
    perl /opt/SIFT4G_Create_Genomic_DB/make-SIFT-db-all.pl --help > /dev/null 2>&1 || true

WORKDIR /data