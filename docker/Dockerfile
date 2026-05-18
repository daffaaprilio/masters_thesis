# Stage 1: extract Clair3 model files from the official image
FROM --platform=linux/amd64 docker.io/hkubal/clair3:latest AS clair3-source

# Stage 2: main thesis-tools image
FROM --platform=linux/amd64 docker.io/condaforge/miniforge3:latest

LABEL org.opencontainers.image.title="thesis-tools"
LABEL org.opencontainers.image.description="ONT sorghum bioinformatics pipeline (minimap2, samtools, bcftools, modkit, whatshap, snpeff, bedtools, clair3, snakemake)"

# Install all tools via conda/bioconda
COPY docker/environment.yml /tmp/environment.yml
RUN mamba env update -n base -f /tmp/environment.yml && \
    conda clean -afy

# Copy Clair3 model files from the official image (keeps models at the same /opt/models path)
COPY --from=clair3-source /opt/models /opt/models

# Verify key tools are on PATH
RUN minimap2 --version && \
    samtools --version | head -1 && \
    bcftools --version | head -1 && \
    modkit --version && \
    whatshap --version && \
    bedtools --version | head -1 && \
    run_clair3.sh --help > /dev/null 2>&1 && \
    snakemake --version

WORKDIR /data
