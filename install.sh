#!/bin/bash
#
# KGE-Mageck2 Installer — MAGeCK2 with interactive HTML reports (UMI-aware)
#
# Creates the mageckenv2 conda environment (mageck2 is pip-installed from a
# pinned git commit — its setup.py compiles the RRA + mageckGSEA C++
# helpers at install time), then overlays the KGE-modified source
# (interactive Plotly HTML report + UMI count figure).
#
# Usage:
#   bash install.sh            # Create new 'mageckenv2' conda environment
#   bash install.sh myenv      # Create environment with custom name
#
set -euo pipefail

ENV_NAME="${1:-mageckenv2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================"
echo "  KGE-Mageck2 Installer — MAGeCK2 + HTML"
echo "============================================"
echo ""

# ---- Step 1: Check / Install conda ----
if ! command -v conda &> /dev/null; then
    echo "[0/4] conda not found — installing Miniconda..."

    # Detect OS and architecture
    OS="$(uname -s)"
    ARCH="$(uname -m)"
    case "${OS}" in
        Linux)  PLATFORM="Linux" ;;
        Darwin) PLATFORM="MacOSX" ;;
        *)      echo "ERROR: Unsupported OS: ${OS}"; exit 1 ;;
    esac
    case "${ARCH}" in
        x86_64)  MARCH="x86_64" ;;
        aarch64|arm64) MARCH="arm64" ;;
        *)       echo "ERROR: Unsupported architecture: ${ARCH}"; exit 1 ;;
    esac

    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-${PLATFORM}-${MARCH}.sh"
    MINICONDA_INSTALLER="/tmp/miniconda_installer.sh"

    echo "  Downloading: ${MINICONDA_URL}"
    curl -fsSL -o "${MINICONDA_INSTALLER}" "${MINICONDA_URL}"
    chmod +x "${MINICONDA_INSTALLER}"

    MINICONDA_DIR="${HOME}/miniconda3"
    echo "  Installing to: ${MINICONDA_DIR}"
    bash "${MINICONDA_INSTALLER}" -b -p "${MINICONDA_DIR}"
    rm -f "${MINICONDA_INSTALLER}"

    # Initialize conda for both bash and zsh (macOS defaults to zsh)
    eval "$("${MINICONDA_DIR}/bin/conda" shell.bash hook)"
    conda init bash 2>/dev/null || true
    conda init zsh  2>/dev/null || true

    echo "  Miniconda installed successfully."
    echo "  Shell integration configured for bash and zsh."
    echo "  (Restart your terminal for conda to be available in new sessions.)"
    echo ""
fi

# The env spec (environment.yml) is cross-platform — it resolves on
# linux-64, linux-aarch64, osx-64 and osx-arm64 (mageck2 is a pip git
# install, the R stack is noarch), so no CONDA_SUBDIR emulation is needed.
CONDA_PLATFORM="$(conda info --json 2>/dev/null | "$(conda info --base)/bin/python" -c "import sys,json; print(json.load(sys.stdin)['platform'])" 2>/dev/null || echo unknown)"
echo "  Detected conda platform: ${CONDA_PLATFORM}"

# ---- Step 2: Create conda environment ----
echo "[1/4] Creating conda environment '${ENV_NAME}' (mageck2 + R/LaTeX toolchain + C/C++ compiler)..."
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "  Environment '${ENV_NAME}' already exists."
    read -p "  Remove and recreate? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n "${ENV_NAME}" -y
    else
        echo "  Keeping existing environment. Will (re-)apply the KGE overlay."
    fi
fi

if ! conda env list | grep -q "^${ENV_NAME} "; then
    # The pip step of the env spec installs mageck2 from the pinned git
    # commit (davidliwei/mageck2 @ fb2f9083) and compiles the bundled
    # RRA + mageckGSEA C++ helpers — this can take a few minutes.
    conda env create -f "${SCRIPT_DIR}/environment.yml" -n "${ENV_NAME}"
fi

# ---- Step 3: Find site-packages path ----
echo "[2/4] Locating site-packages..."
SITE_PACKAGES=$(conda run -n "${ENV_NAME}" python -c "import site; print(site.getsitepackages()[0])" 2>/dev/null || true)

# Fallback: activation (some shells/conda versions break `conda run`)
if [ -z "${SITE_PACKAGES}" ] || [ ! -d "${SITE_PACKAGES}" ]; then
    eval "$(conda shell.bash hook)"
    set +u   # some conda-forge activation hooks reference bare vars (conda/conda#9966)
    conda activate "${ENV_NAME}" 2>/dev/null || true
    set -u
    SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])" 2>/dev/null || true)
fi

MAGECK2_DIR="${SITE_PACKAGES:-/nonexistent}/mageck2"

if [ ! -d "${MAGECK2_DIR}" ]; then
    echo "ERROR: mageck2 package not found at ${MAGECK2_DIR}"
    echo "  The pip git install (step 1) may have failed — check network access"
    echo "  to github.com and the env's C/C++ toolchain (make, g++/clang++)."
    exit 1
fi
echo "  Found: ${MAGECK2_DIR}"

# ---- Step 4: Overlay the KGE-modified source ----
echo "[3/4] Installing KGE modifications (HTML report + UMI figure)..."
for ext in py gmt txt Rnw Rmd RTemplate; do
    for f in "${SCRIPT_DIR}"/mageck2/*.${ext}; do
        [ -e "${f}" ] || continue
        # Delete-then-copy (NOT a plain cp): env files are hard links into
        # the conda package cache; copying over them in place would rewrite
        # the shared cache inode and corrupt it. Done explicitly because
        # macOS BSD cp has no --remove-destination flag.
        rm -f "${MAGECK2_DIR}/$(basename "${f}")"
        cp "${f}" "${MAGECK2_DIR}/"
    done
done

# Clear bytecode cache so Python picks up new source
rm -rf "${MAGECK2_DIR}/__pycache__"

echo "  Copied modified files to ${MAGECK2_DIR}"

# ---- Step 5: Verify ----
echo "[4/4] Verifying installation..."
conda run -n "${ENV_NAME}" python -c "
from mageck2.version import __version__
from mageck2.htmlReport import mageck_report_main
import plotly
import decoupler
print(f'  MAGeCK2 version:   {__version__}')
print(f'  Plotly version:    {plotly.__version__}')
print(f'  Decoupler version: {decoupler.__version__}')
print(f'  HTML report:       OK (mageck2.htmlReport.mageck_report_main)')
" 2>/dev/null

# The compiled C++ helpers live on the env's PATH (installed by mageck2's setup.py)
conda run -n "${ENV_NAME}" bash -c '
for bin in RRA mageckGSEA; do
    if command -v ${bin} >/dev/null 2>&1; then
        echo "  ${bin}: present (compiled C++ helper)"
    else
        echo "  ${bin}: MISSING — mageck2 count/test will fail; re-run the installer"
    fi
done
' 2>/dev/null || true

echo ""
echo "============================================"
echo "  Installation complete!"
echo "============================================"
echo ""
echo "  Activate:  conda activate ${ENV_NAME}"
echo ""
echo "  Usage (UMI-aware counting, R1 = guide / R2 = UMI):"
echo "    mageck2 count -l library.csv -n demo \\"
echo "        --fastq  S1_R1 S2_R1 ... \\"
echo "        --fastq-2 S1_R2 S2_R2 ... \\"
echo "        --sample-label ctrl,treat1,treat2 \\"
echo "        --umi secondpair --umi-start-2 0 --umi-end-2 10 \\"
echo "        --pdf-report --html-report"
echo ""
echo "    mageck2 test  -k demo.count.txt -n demo -t treat1,treat2 -c ctrl \\"
echo "        --norm-method none --pdf-report --html-report"
echo "    mageck2 mle   -k demo.count.txt -d design.txt -n demo --html-report"
echo ""
echo "  Note: mageck2 has no standalone 'report' subcommand — the pipeline's"
echo "  KGE_REPORT step calls mageck2.htmlReport.mageck_report_main directly."
echo ""