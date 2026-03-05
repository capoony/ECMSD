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
  -f | --fwd or -m | --merged  FASTQ_FILE  Path to the forward or merged FASTQ file
  -o | --out OUTPUT_FOLDER  Path to the output folder

OPTIONAL ARGUMENTS:
  -f | --fwd FWD_FASTQ                  Path to the forward FASTQ file (default: None)
  -r | --rev REV_FASTQ                  Path to the reverse FASTQ file (default: None)
  -m | --merged MERGED_FASTQ            Path to merged FASTQ file (default: None)
  -b | --binsize BINSIZE                Bin size for analysis (default: BINSIZE = 1000)
  -u | --RMUS-threshold THRESHHOLD      RMUS threshold for analysis (default: THRESHHOLD = 0.15)
  -q | --mapping_quality QUALITY        Mapping quality threshold (default: QUALITY = 20)
  -p | --prefix PREFIX                  Prefix for output files (default: None)
  -s | --skip_environment               Skip conda environment setup (default: false)
  -t | --threads THREADS                Number of threads to use (default: THREADS = 10)
  -c | --force                          Force overwrite of existing output files (default: false)
  -x | --taxonomic-hierarchy HIERARCHY  Taxonomic hierarchy (default: HIERARCHY = species)
  -v | --version                        Show version and exit
  -h | --help                           Show this help message and exit

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
force=false
version="1.0"
taxonomic_hierarchy="species"
skip_environment=false
output=""
prefix=""

echo "Starting ECMSD pipeline..."

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
    -b | --binsize)
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
        force=true
        shift
        ;;
    -x | --taxonomic-hierarchy)
        taxonomic_hierarchy="$2"
        shift 2
        ;;
    -p | --prefix)
        prefix="$2"
        shift 2
        ;;
    -s | --skip_environment)
        skip_environment=true
        shift
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

# if rev is provided, also fwd is required
if [[ -n "${rev}" && -z "${fwd}" ]]; then
    echo "Error: Forward FASTQ file must be provided if reverse FASTQ file is specified."
    exit 1
fi

# if merged is provided, fwd and rev should be skipped
if [[ -n "${merged}" && ( -n "${fwd}" || -n "${rev}" ) ]]; then
    echo "Error: Merged FASTQ file cannot be provided with forward or reverse FASTQ files."
    exit 1
fi

###############################################################################
#                        Prepare Output Directory                             #
###############################################################################
Output="${output%/}"
[[ -z "${Output}" ]] && {
    echo "Error: Output directory is not specified."
    exit 1
}

if [[ "${force}" == true && -d "${Output}" ]]; then
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
if [[ "${skip_environment}" == false ]]; then
    echo "Setting up conda environment..."
    [[ ! -d "${WD}/scripts/conda_env" ]] && bash "${WD}/shell/requirements.sh" "${WD}" "${Output}"
else
    echo "Skipping conda environment setup as per user request."
fi

[[ ! -d "${WD}/data/refseq_mito" ]] && bash "${WD}/shell/MakeRef.sh" "${WD}" "${threads}"

###############################################################################
#                        Activate Conda Environment                           #
###############################################################################
if [[ ! -d "${WD}/scripts/conda_env" ]] && [[ "${skip_environment}" == false ]]; then
    echo "Error: Conda environment not found."
    exit 1
fi

if [[ "${skip_environment}" == false ]]; then
    echo "Conda environment found at ${WD}/scripts/conda_env"

    echo "Activating conda environment..."
    eval "$(conda shell.bash hook)" || {
        echo "Conda shell hook could not be initialized"
        exit 1
    }
    conda activate "${WD}/scripts/conda_env" || {
        echo "Conda environment could not be activated"
        exit 1
    }
else
    echo "Skipping conda environment activation as per user request."
fi

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
        2> "${Output}/logs/minimap2.log" |
        awk -v Q="${quality}" '$12 >= Q {print}'

    #prtint log file content
    echo "Minimap2 log:"
    cat "${Output}/logs/minimap2.log"
}

###############################################################################
#                        Mapping Step                                         #
###############################################################################
echo "Starting mapping step..."

PAF="${Output}/mapping/Mito.paf"

if [[ -n "${prefix}" ]]; then
    PAF="${Output}/mapping/${prefix}_Mito.paf"
fi

REF="${WD}/data/refseq_mito/mitochondrion_refseq_taxid_masked.fna.gz"

# check if fwd is provided
if [[ -n "${fwd}" ]]; then
    if [[ -z "${rev}" || "${rev}" == "no" ]]; then
        echo "Running single-end mapping..."
        run_mapping "${fwd}" | gzip >"${PAF}.gz"
    else
        echo "Running paired-end mapping..."
        run_mapping "${fwd} ${rev}" | gzip >"${PAF}.gz"
    fi
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

#check if PAF file is created
if [[ ! -f "${PAF}.gz" ]]; then
    echo "Error: PAF file '${PAF}.gz' was not created."
    exit 1
fi

echo "Parsing PAF and calculating RMUS/taxonomic proportions..."

output_base="${Output}/mapping/Mito_summary"
# if prefix is provided, modify output file name
if [[ -n "${prefix}" ]]; then
    output_base="${Output}/mapping/${prefix}_Mito_summary"
fi

python "${WD}/scripts/LinkTaxonomy.py" \
    --Nodes "${WD}/data/refseq_mito/NCBI_taxdump/nodes.dmp" \
    --Names "${WD}/data/refseq_mito/NCBI_taxdump/names.dmp" \
    --PAF "${PAF}.gz" \
    --Bins "${binsize}" \
    --RMUS "${rmus_threshold}" \
    --output "${output_base}"

###############################################################################
#                        Plotting Results                                     #
###############################################################################

summary_file="${output_base}.txt"

if [[ ! -f "${summary_file}" ]]; then
    echo "Error:Mito summary file '${summary_file}' was not created."
    exit 1
fi

#check if summary file has more than one line (header + at least one data line)
if [[ $(wc -l <"${summary_file}") -le 1 ]]; then
    echo "Error: Mito summary file '${summary_file}' is empty or has no data."
    exit 1
fi

echo "Plotting results..."
Rscript "${WD}/scripts/process_files.R" ${summary_file} ${Output} ${taxonomic_hierarchy} ${prefix}

###############################################################################
#                        Pipeline Completed                                   #
###############################################################################
echo "ECMSD pipeline completed successfully."
