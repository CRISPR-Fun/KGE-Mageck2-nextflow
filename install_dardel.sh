#!/bin/bash
#
# KGE-Mageck2 Installer for Dardel (PDC/KTH HPC)
#
# Loads required modules, creates a conda environment under Klemming
# (mageck2 is pip-installed from a pinned git commit — its setup.py
# compiles the RRA + mageckGSEA C++ helpers), then overlays the KGE-
# modified source (interactive Plotly HTML report + UMI count figure).
#
# Usage:
#   bash install_dardel.sh              # Create 'mageckenv2' environment
#   bash install_dardel.sh myenv        # Create environment with custom name
#
set -euo pipefail

ENV_NAME="${1:-mageckenv2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================"
echo "  KGE-Mageck2 Installer — Dardel (PDC/KTH)"
echo "============================================"
echo ""

# ---- Step 1: Load Dardel modules ----
echo "[1/5] Loading Dardel modules..."
if ! command -v ml &> /dev/null && ! command -v module &> /dev/null; then
    echo "ERROR: module system not found. Are you on Dardel?"
    echo "  This installer is for PDC Dardel. Use install.sh for local machines."
    exit 1
fi

# Load required Dardel modules
# Adjust versions if needed for your Dardel allocation
ml PDC/23.12
ml bioinfo-tools
ml java/11
ml miniconda3/25.3.1-1-cpeGNU-24.11
ml nextflow/24.04.2
ml singularity/4.1.1-cpeGNU-23.12

echo "  Loaded PDC, bioinfo-tools, Java 11, miniconda3, Nextflow 24.04, Singularity 4.1.1"

# Verify conda is available after module load
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not available after loading modules."
    echo "  Check available versions with: ml spider miniconda3"
    exit 1
fi

# ---- Step 2: Create conda environment ----
echo "[2/5] Creating conda environment '${ENV_NAME}' (mageck2 + R/LaTeX toolchain + C/C++ compiler)..."
if conda env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
    echo "  Environment '${ENV_NAME}' already exists."
    read -p "  Remove and recreate? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n "${ENV_NAME}" -y
    else
        echo "  Keeping existing environment. Will (re-)apply the KGE overlay."
    fi
fi

if ! conda env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
    # The pip step of the env spec installs mageck2 from the pinned git
    # commit (davidliwei/mageck2 @ fb2f9083) and compiles the bundled
    # RRA + mageckGSEA C++ helpers — this can take a few minutes.
    conda env create -f "${SCRIPT_DIR}/environment.yml" -n "${ENV_NAME}"
fi

# ---- Step 3: Activate and find site-packages ----
echo "[3/5] Locating site-packages..."

# On Dardel, 'conda run' may not work reliably — use source activate
# Save current state so we can provide activation instructions
SITE_PACKAGES=$(conda run -n "${ENV_NAME}" python -c "import site; print(site.getsitepackages()[0])" 2>/dev/null || true)

# Fallback: try source activate if conda run failed
if [ -z "${SITE_PACKAGES}" ] || [ ! -d "${SITE_PACKAGES}" ]; then
    echo "  conda run failed, trying source activate..."
    source activate "${ENV_NAME}" 2>/dev/null || conda activate "${ENV_NAME}" 2>/dev/null || true
    SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])" 2>/dev/null)
fi

MAGECK2_DIR="${SITE_PACKAGES:-/nonexistent}/mageck2"

if [ ! -d "${MAGECK2_DIR}" ]; then
    echo "ERROR: mageck2 package not found at ${MAGECK2_DIR}"
    echo "  The pip git install (step 2) may have failed."
    echo "  Try activating manually: source activate ${ENV_NAME}"
    exit 1
fi
echo "  Found: ${MAGECK2_DIR}"

# ---- Step 4: Overlay the KGE-modified source ----
echo "[4/5] Installing KGE modifications (HTML report + UMI figure)..."
for ext in py gmt txt Rnw Rmd RTemplate; do
    for f in "${SCRIPT_DIR}"/mageck2/*.${ext}; do
        [ -e "${f}" ] || continue
        # Delete-then-copy (NOT a plain cp): env files are hard links into
        # the conda package cache; copying over them in place would rewrite
        # the shared cache inode and corrupt it.
        rm -f "${MAGECK2_DIR}/$(basename "${f}")"
        cp "${f}" "${MAGECK2_DIR}/"
    done
done

# Clear bytecode cache so Python picks up new source
rm -rf "${MAGECK2_DIR}/__pycache__"

echo "  Copied modified files to ${MAGECK2_DIR}"

# ---- Step 5: Verify ----
echo "[5/5] Verifying installation..."
conda run -n "${ENV_NAME}" python -c "
from mageck2.version import __version__
from mageck2.htmlReport import mageck_report_main
import plotly
import decoupler
print(f'  MAGeCK2 version:   {__version__}')
print(f'  Plotly version:    {plotly.__version__}')
print(f'  Decoupler version: {decoupler.__version__}')
print(f'  HTML report:       OK (mageck2.htmlReport.mageck_report_main)')
" 2>/dev/null || {
    # Fallback verification with source activate
    source activate "${ENV_NAME}" 2>/dev/null || true
    python -c "
from mageck2.version import __version__
from mageck2.htmlReport import mageck_report_main
import plotly
import decoupler
print(f'  MAGeCK2 version:   {__version__}')
print(f'  Plotly version:    {plotly.__version__}')
print(f'  Decoupler version: {decoupler.__version__}')
print(f'  HTML report:       OK (mageck2.htmlReport.mageck_report_main)')
"
}

echo ""
echo "============================================"
echo "  Installation complete!"
echo "============================================"
echo ""
echo "  To use on Dardel, load modules then activate:"
echo "  (The pipeline requires Nextflow >= 24.04 — Nextflow/22.10.1 is NOT compatible)"
echo ""
echo "    ml PDC/23.12"
echo "    ml bioinfo-tools"
echo "    ml java/11"
echo "    ml miniconda3/25.3.1-1-cpeGNU-24.11"
echo "    ml nextflow/24.04.2"
echo "    ml singularity/4.1.1-cpeGNU-23.12"
echo "    source activate ${ENV_NAME}"
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