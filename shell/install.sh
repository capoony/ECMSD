#!/bin/bash

# install.sh - Installs ECMSD and its dependencies
#
# By default the dependencies and the ECMSD command are installed into the
# currently active conda environment. With --prefix a new conda environment is
# created at the given path and ECMSD is installed there instead.

set -euo pipefail

WD="$(dirname "$(dirname "$(realpath "$0")")")"

usage() {
    cat <<EOF
Usage: bash install.sh [-p|--prefix PATH]

Installs ECMSD and its dependencies (minimap2, bbmap, R, Python, pigz, wget).

Options:
  -p, --prefix PATH   Create a new local conda environment at PATH (conda -p)
                      and install ECMSD into it. PATH may be relative or
                      absolute. If PATH is already a conda environment, ECMSD
                      is (re)installed into it. A non-empty directory that is
                      not a conda environment is refused.
  -h, --help          Show this help message and exit.

Without --prefix, ECMSD is installed into the active conda environment
(conda activate <your-env> first).

Examples:
  bash install.sh                              # into the active environment
  bash install.sh --prefix /path/to/ecmsd_env  # into a new local environment
EOF
}

PREFIX_ARG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
    -p | --prefix)
        if [[ $# -lt 2 || -z "${2:-}" || "$2" == -* ]]; then
            echo "Error: $1 requires a PATH argument."
            exit 1
        fi
        PREFIX_ARG="$2"
        shift 2
        ;;
    --prefix=*)
        PREFIX_ARG="${1#--prefix=}"
        if [[ -z "$PREFIX_ARG" ]]; then
            echo "Error: --prefix requires a PATH argument."
            exit 1
        fi
        shift
        ;;
    -h | --help)
        usage
        exit 0
        ;;
    *)
        echo "Error: unknown option '$1'"
        usage
        exit 1
        ;;
    esac
done

# Check conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install conda first."
    exit 1
fi

CHANNELS=(-c bioconda -c conda-forge --no-channel-priority)
PACKAGES=(
    bioconda::minimap2
    bioconda::bbmap
    "conda-forge::r-base>=4.0"
    "conda-forge::r-tidyverse>=2.0.0"
    "conda-forge::r-dplyr>=1.1.0"
    conda-forge::r-gridextra
    conda-forge::r-data.table
    conda-forge::pigz
    conda-forge::wget
    "conda-forge::python>=3.6=*cpython*"
)

if [[ -n "$PREFIX_ARG" ]]; then
    # Resolve to an absolute path; the directory does not have to exist yet
    TARGET="$(realpath -m "$PREFIX_ARG")"

    # conda does not support spaces in environment paths
    if [[ "$TARGET" == *[[:space:]]* ]]; then
        echo "Error: conda environment paths must not contain spaces: '$TARGET'"
        exit 1
    fi

    if [[ -d "$TARGET/conda-meta" ]]; then
        echo "Existing conda environment found at: $TARGET"
        echo "Installing ECMSD dependencies into it..."
        conda install -y -p "$TARGET" "${CHANNELS[@]}" "${PACKAGES[@]}"
    elif [[ -e "$TARGET" && ( ! -d "$TARGET" || -n "$(ls -A "$TARGET")" ) ]]; then
        echo "Error: '$TARGET' exists and is not a conda environment. Choose another path."
        exit 1
    else
        echo "Creating new conda environment at: $TARGET"
        conda create -y -p "$TARGET" "${CHANNELS[@]}" "${PACKAGES[@]}"
    fi
    ACTIVATE_HINT="conda activate \"$TARGET\""
else
    # Check we're inside a conda environment
    if [[ -z "${CONDA_PREFIX:-}" ]]; then
        echo "Error: No conda environment is active. Please activate one first:"
        echo "  conda activate <your-env>"
        echo "or create a new one with:"
        echo "  bash install.sh --prefix /path/to/ecmsd_env"
        exit 1
    fi
    TARGET="$CONDA_PREFIX"

    echo "Installing ECMSD dependencies into: $TARGET"
    conda install -y -p "$TARGET" "${CHANNELS[@]}" "${PACKAGES[@]}"
    ACTIVATE_HINT=""
fi

# Remove existing installation if present
if [[ -d "$TARGET/lib/ecmsd" ]]; then
    echo "Existing ECMSD installation found, replacing..."
    rm -rf "$TARGET/lib/ecmsd"
fi
if [[ -f "$TARGET/bin/ECMSD" ]]; then
    rm -f "$TARGET/bin/ECMSD"
fi

# Install scripts into the conda environment
mkdir -p "$TARGET/bin"
mkdir -p "$TARGET/lib/ecmsd/shell"
mkdir -p "$TARGET/lib/ecmsd/scripts"

cp "$WD/shell/"*.sh "$TARGET/lib/ecmsd/shell/"
cp "$WD/scripts/"*.py "$TARGET/lib/ecmsd/scripts/"
cp "$WD/scripts/"*.R "$TARGET/lib/ecmsd/scripts/"

# Create the ECMSD entry point in the environment's bin
cp "$WD/shell/ECMSD.sh" "$TARGET/bin/ECMSD"
chmod +x "$TARGET/bin/ECMSD"

echo ""
echo "Done! ECMSD was installed into: $TARGET"
echo "You can now run:"
if [[ -n "$ACTIVATE_HINT" ]]; then
    echo "  $ACTIVATE_HINT"
fi
echo "  ECMSD --help"
