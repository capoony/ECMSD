#!/bin/bash

# install.sh - Installs ECMSD dependencies into the current conda environment

set -euo pipefail

WD="$(dirname "$(dirname "$(realpath "$0")")")"


# Check conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install conda first."
    exit 1
fi

# Check we're inside a conda environment
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "Error: No conda environment is active. Please activate one first:"
    echo "  conda activate <your-env>"
    exit 1
fi

echo "Installing ECMSD dependencies into: $CONDA_PREFIX"


conda install -y \
    -c bioconda \
    -c conda-forge \
    --no-channel-priority \
    bioconda::minimap2 \
    bioconda::bbmap \
    "conda-forge::r-base>=4.0" \
    "conda-forge::r-tidyverse>=2.0.0" \
    "conda-forge::r-dplyr>=1.1.0" \
    conda-forge::pigz \
    conda-forge::wget \
    "conda-forge::python>=3.6" \
    conda-forge::numpy

# Remove existing installation if present
if [[ -d "$CONDA_PREFIX/lib/ecmsd" ]]; then
    echo "Existing ECMSD installation found, replacing..."
    rm -rf "$CONDA_PREFIX/lib/ecmsd"
fi
if [[ -f "$CONDA_PREFIX/bin/ECMSD" ]]; then
    rm -f "$CONDA_PREFIX/bin/ECMSD"
fi

# Install scripts into the conda environment
mkdir -p "$CONDA_PREFIX/lib/ecmsd/shell"
mkdir -p "$CONDA_PREFIX/lib/ecmsd/scripts"

cp "$WD/shell/"*.sh "$CONDA_PREFIX/lib/ecmsd/shell/"
cp "$WD/scripts/"*.py "$CONDA_PREFIX/lib/ecmsd/scripts/"
cp "$WD/scripts/"*.R "$CONDA_PREFIX/lib/ecmsd/scripts/"

# Create the ECMSD entry point in the environment's bin
cp "$WD/shell/ECMSD.sh" "$CONDA_PREFIX/bin/ECMSD"
chmod +x "$CONDA_PREFIX/bin/ECMSD"

echo ""
echo "Done! You can now run:"
echo "  ECMSD --help"