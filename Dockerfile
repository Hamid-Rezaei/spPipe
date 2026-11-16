# spPipe runtime image: the same pinned conda environment used by -profile conda.
# Build (from the repo root):  docker build -t sppipe:1.0 .
# Nextflow then uses it with:  nextflow run main.nf -profile docker
FROM mambaorg/micromamba:1.5.10

# procps provides `ps`, which Nextflow needs to collect task metrics.
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes

# Make the env usable without activation (Nextflow runs `bash -c ...`)
# Caches must be writable when Nextflow runs the container as the host user
ENV PATH=/opt/conda/bin:$PATH \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/mplconfig \
    NUMBA_CACHE_DIR=/tmp/numba_cache

LABEL org.opencontainers.image.title="spPipe" \
      org.opencontainers.image.description="Spatial transcriptomics (Visium DLPFC) pipeline environment: scanpy, squidpy, leidenalg"
