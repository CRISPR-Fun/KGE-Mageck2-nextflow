# ============================================================================
# Dockerfile for MAGeCK2-KGE: UMI-aware CRISPR screen analysis with
# interactive HTML reports
# ============================================================================
# mageck2 is pip-installed from a pinned git commit (the environment.yml pip
# section) — its setup.py compiles the bundled C++ helpers (RRA, mageckGSEA)
# at install time, so the base image needs a C/C++ compiler + make (provided
# by the env spec: c-compiler/cxx-compiler/make).
#
# Build:
#   docker build -t mobasherbarsi/mageck2-kge:0.4.0 .
# Push:
#   docker push mobasherbarsi/mageck2-kge:0.4.0
# ============================================================================

FROM continuumio/miniconda3:25.3.1-1 AS base

LABEL org.opencontainers.image.title="MAGeCK2-KGE"
LABEL org.opencontainers.image.description="MAGeCK2 (UMI-aware) with KGE interactive Plotly HTML reports for CRISPR screen analysis"
LABEL org.opencontainers.image.version="0.4.0"
LABEL org.opencontainers.image.authors="Mubasher Mohammed"
LABEL org.opencontainers.image.source="https://github.com/MubasherMohammed/KGE-Mageck2-nextflow"
LABEL org.opencontainers.image.base="mageck2 git+https://github.com/davidliwei/mageck2.git@fb2f9083"

ENV CONDA_ENV=mageckenv2
ENV CONDA_DIR=/opt/conda
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# ---------------------------------------------------------------------------
# Create conda environment with MAGeCK2 + dependencies
# (the pip step installs mageck2 from the pinned git commit and compiles the
# RRA + mageckGSEA C++ helpers)
# ---------------------------------------------------------------------------
COPY environment.yml /tmp/environment.yml

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda env create -f /tmp/environment.yml -n ${CONDA_ENV} && \
    conda clean -afy

# ---------------------------------------------------------------------------
# Overlay KGE modifications (enhanced HTML report, UMI count figure,
# robust PDF report templates)
# ---------------------------------------------------------------------------
COPY mageck2/ /tmp/kge/mageck2/

RUN MAGECK2_DIR=$(${CONDA_DIR}/envs/${CONDA_ENV}/bin/python -c "import mageck2; print(mageck2.__file__)" | xargs dirname) && \
    for ext in py Rnw; do \
        cp /tmp/kge/mageck2/*.${ext} "${MAGECK2_DIR}/" 2>/dev/null || true; \
    done && \
    rm -rf "${MAGECK2_DIR}/__pycache__" /tmp/kge && \
    echo "KGE modifications installed to ${MAGECK2_DIR}"

# ---------------------------------------------------------------------------
# Make conda env the default PATH
# ---------------------------------------------------------------------------
ENV CONDA_DEFAULT_ENV=${CONDA_ENV}
ENV PATH=/opt/conda/envs/${CONDA_ENV}/bin:${PATH}

# Verify installation
RUN mageck2 --version && \
    python -c "from mageck2.htmlReport import mageck_report_main; print('KGE HTML report: OK')" && \
    python -c "import plotly; print(f'Plotly {plotly.__version__}: OK')" && \
    test -x "$(command -v RRA)" && echo "RRA: OK (compiled C++ helper)" && \
    test -x "$(command -v mageckGSEA)" && echo "mageckGSEA: OK (compiled C++ helper)" && \
    mageck2 count --help | grep -q -- '--umi' && echo "UMI support: OK"

# Default working directory for pipeline runs
WORKDIR /data

# NOTE: No ENTRYPOINT set intentionally.
# Nextflow runs process scripts via `bash -c` inside containers, so an ENTRYPOINT
# like `mageck2` would prepend itself to every command (e.g. `mageck2 bash -c "fastqc ..."`)
# and break non-mageck2 processes. The PATH already points to the conda env with mageck2.

# Default command — shows help when container is run interactively
CMD ["mageck2", "--help"]