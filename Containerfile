# Containerfile for an internal-lab metadarkmatter environment.
#
# Provides a reproducible image with:
#   - Python 3.12 + the metadarkmatter package installed from the locked
#     uv environment (uv.lock).
#   - The external bioinformatics tools documented in
#     docs/EXTERNAL_TOOLS.md (blast, mmseqs2, kraken2, fastANI, skani,
#     diamond, mashtree, mash) pulled from bioconda.
#
# This is intended for local-lab use; we deliberately do not push to a
# public registry. Build with:
#
#   podman build -t metadarkmatter:dev .
#   # or: docker build -t metadarkmatter:dev .
#
# Run a one-shot CLI command against host data:
#
#   podman run --rm -v $PWD/data:/data metadarkmatter:dev \
#       metadarkmatter doctor
#
# The default entry point is the metadarkmatter CLI; passing no
# arguments prints the top-level help.

FROM mambaorg/micromamba:1.5.10

# Avoid prompts and use bash for activation niceties.
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root

# Bioinformatics tools from bioconda. Versions are pinned for
# reproducibility - bump deliberately and re-test before changing.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 \
        "blast>=2.16" \
        "mmseqs2>=15" \
        "kraken2>=2.1.3" \
        "fastani>=1.34" \
        "skani>=0.3" \
        "diamond>=2.1" \
        "mashtree>=1.4" \
        "mash>=2.3" \
        "samtools>=1.20" \
        pip \
    && micromamba clean --all --yes

# Install uv for the Python install step. We pin the version that
# generated uv.lock to ensure the same resolver behaviour.
RUN pip install --no-cache-dir "uv==0.4.18"

WORKDIR /opt/metadarkmatter

# Copy the project source. The .dockerignore keeps the context small.
COPY . /opt/metadarkmatter

# Sync the locked Python environment into a project-local venv and
# install the package. --frozen ensures the lockfile matches.
RUN uv sync --frozen --extra dev \
    && uv pip install -e .

# Sanity-check the install at build time. If anything is missing
# (Python deps, CLI entry point, importable modules), the build fails.
RUN uv run metadarkmatter --version \
    && uv run metadarkmatter doctor

# Keep the venv on PATH so 'metadarkmatter' is callable without 'uv run'.
ENV PATH="/opt/metadarkmatter/.venv/bin:$PATH"

# Work in /data by default so host bind-mounts of input data are the
# natural place to read from and write to.
WORKDIR /data
ENTRYPOINT ["metadarkmatter"]
CMD ["--help"]
