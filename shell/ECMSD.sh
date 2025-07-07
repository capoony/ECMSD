#!/bin/bash

###############################################################################
# ECMSD.sh - Efficient Comprehensive Mitochondrial Sequence Detection Pipeline #
###############################################################################

set -euo pipefail

###############################################################################
#                               Usage Function                                #
###############################################################################
usage() {
    cat <<EOF
Usage:
    ECMSD.sh -f|--fwd FWD_FASTQ -o|--out OUTPUT_FOLDER [options]

REQUIRED ARGUMENTS:
  -f | --fwd         Path to the forward FASTQ file
  -o | --out         Path to the output folder

OPTIONAL ARGUMENTS:
  -r | --rev         Path to the reverse FASTQ file (default: no)
  -m | --merged      Path to merged FASTQ file (default: no)
  -b | --Binsize     Bin size for analysis (default: 1000)
  -u | --RMUS-threshold  RMUS threshold for analysis (default: 0.15)
  -q | --mapping_quality Mapping quality threshold (default: 20)
  -t | --threads     Number of threads to use (default: 10)
  -c | --force       Force overwrite of existing output files
  -x | --taxonomic-hierarchy Taxonomic hierarchy (default: species)
  -v | --version     Show version and exit
  -h | --help        Show this help message and exit

Example:
  ECMSD.sh -f reads_R1.fastq -r reads_R2.fastq -o results/
EOF
    exit 1
}

###############################################################################
#                              Default Values                                 #
###############################################################################
fwd=""
rev=""
merged=""
binsize=1000
rmus_threshold=0.15
quality=20
threads=10
force="no"
version="1.0"
taxonomic_hierarchy="species"
output=""

###############################################################################
#                             Argument Parsing                                #
###############################################################################
while [[ $# -gt 0 ]]; do
    case $1 in
    -f | --fwd)
        fwd="$2"
        shift 2
        ;;
    -r | --rev)
        rev="$2"
        shift 2
        ;;
    -m | --merged)
        merged="$2"
        shift 2
        ;;
    -b | --Binsize)
        binsize="$2"
        shift 2
        ;;
    -u | --RMUS-threshold)
        rmus_threshold="$2"
        shift 2
        ;;
    -q | --mapping_quality)
        quality="$2"
        shift 2
        ;;
    -t | --threads)
        threads="$2"
        shift 2
        ;;
    -c | --force)
        force="yes"
        shift
        ;;
    -x | --taxonomic-hierarchy)
        taxonomic_hierarchy="$2"
        shift 2
        ;;
    -o | --out)
        output="$2"
        shift 2
        ;;
    -v | --version)
        echo "ECMSD version ${version}"
        exit 0
        ;;
    -h | --help)
        usage
        ;;
    *)
        echo "Unknown option: $1"
        usage
        ;;
    esac
done

###############################################################################
#                        Check Required Arguments                             #
###############################################################################
if [[ -z "${fwd}" || -z "${output}" ]]; then
    echo "Error: --fwd and --out are required."
    usage
fi

###############################################################################
#                            Check Input Files                                #
###############################################################################
[[ ! -f "${fwd}" ]] && {
    echo "Error: Forward FASTQ file '${fwd}' does not exist."
    exit 1
}
[[ -n "${rev}" && ! -f "${rev}" ]] && {
    echo "Error: Reverse FASTQ file '${rev}' does not exist."
    exit 1
}
[[ -n "${merged}" && ! -f "${merged}" ]] && {
    echo "Error: Merged FASTQ file '${merged}' does not exist."
    exit 1
}

###############################################################################
#                        Prepare Output Directory                             #
###############################################################################
Output="${output%/}"
[[ -z "${Output}" ]] && {
    echo "Error: Output directory is not specified."
    exit 1
}

if [[ "${force}" == "yes" && -d "${Output}" ]]; then
    echo "Removing existing output directory: ${Output}"
    rm -rf "${Output}"
elif [[ -d "${Output}" ]]; then
    echo "Output directory '${Output}' already exists. Use -c or --force to overwrite."
    exit 1
fi

###############################################################################
#                        Find Base Directory                                  #
###############################################################################
WD="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

###############################################################################
#                        Clean Logs Directory                                 #
###############################################################################
if [[ -d "${Output}/logs" ]]; then
    rm -rf "${Output}/logs"
fi

###############################################################################
#                        Create Output Directories                            #
###############################################################################
echo "Creating output directories..."
mkdir -p "${Output}/mapping"
mkdir "${Output}/logs"

###############################################################################
#                        Dependency Setup                                     #
###############################################################################
echo "Checking dependencies and reference data..."
[[ ! -d "${WD}/scripts/conda_env" ]] && bash "${WD}/shell/requirements.sh" "${WD}" "${Output}"
[[ ! -d "${WD}/data/refseq_mito" ]] && bash "${WD}/shell/MakeRef.sh" "${WD}" "${threads}"

###############################################################################
#                        Activate Conda Environment                           #
###############################################################################
echo "Activating conda environment..."
eval "$(conda shell.bash hook)" || {
    echo "Conda shell hook could not be initialized"
    exit 1
}
conda activate "${WD}/scripts/conda_env" || {
    echo "Conda environment could not be activated"
    exit 1
}

###############################################################################
#                        Mapping Function                                     #
###############################################################################
run_mapping() {
    # Function to run minimap2 mapping and log the progress in a log file
    local reads="$1"
    echo "Running minimap2 mapping for: ${reads}"
    minimap2 \
        -x sr \
        --secondary=no \
        -t "${threads}" "${REF}" ${reads} \
        2>"${Output}/logs/logfile.log" |
        awk -v Q="${quality}" '$12 >= Q {print}'
}

###############################################################################
#                        Mapping Step                                         #
###############################################################################
PAF="${Output}/mapping/Mito.paf"
REF="${WD}/data/refseq_mito/mitochondrion_refseq_taxid_masked.fna.gz"

if [[ -z "${rev}" || "${rev}" == "no" ]]; then
    echo "Running single-end mapping..."
    run_mapping "${fwd}" | gzip >"${PAF}.gz"
else
    if [[ -z "${fwd}" || -z "${rev}" ]]; then
        echo "Error: Both forward and reverse reads must be provided for paired-end mapping."
        exit 1
    fi
    echo "Running paired-end mapping..."
    run_mapping "${fwd} ${rev}" | gzip >"${PAF}.gz"
fi

###############################################################################
#                        Mapping on Merged Reads                              #
###############################################################################
if [[ -n "${merged}" && "${merged}" != "no" ]]; then
    if [[ ! -f "${merged}" ]]; then
        echo "Error: Merged reads file '${merged}' does not exist."
        exit 1
    fi
    echo "Running mapping on merged reads..."
    run_mapping "${merged}" | gzip >>"${PAF}.gz"
fi

###############################################################################
#                        Parse PAF and Calculate RMUS                         #
###############################################################################
echo "Parsing PAF and calculating RMUS/taxonomic proportions..."
python "${WD}/scripts/LinkTaxonomy.py" \
    --Nodes "${WD}/data/refseq_mito/NCBI_taxdump/nodes.dmp" \
    --Names "${WD}/data/refseq_mito/NCBI_taxdump/names.dmp" \
    --PAF "${PAF}.gz" \
    --Bins "${binsize}" \
    --RMUS "${rmus_threshold}" \
    --output "${Output}/mapping/Mito_summary"

###############################################################################
#                        Plotting Results                                     #
###############################################################################
echo "Plotting results..."
Rscript "${WD}/scripts/process_files.R" ${Output} ${taxonomic_hierarchy}

###############################################################################
#                        Pipeline Completed                                   #
###############################################################################
echo "ECMSD pipeline completed successfully."
