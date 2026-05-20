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
  ECMSD -m MERGED_FASTQ -o OUTPUT_FOLDER -d DB_FOLDER [options]

REQUIRED ARGUMENTS FOR ANALYSIS (choose one):
  -f | --fwd or -m | --merged  FASTQ_FILE   Path to the forward FASTQ file or merged FASTQ file
  -o | --out OUTPUT_FOLDER                  Path to the output folder
  -d | --db-folder DB_FOLDER                Folder for the database

REQUIRED ARGUMENTS FOR BUILDING DATABASE:
  -z | --create-db                          Create a new database
  -d | --db-folder DB_FOLDER                Folder for the database

ALL ARGUMENTS:
  -h | --help                               Show this help message and exit
  -v | --version                            Show version and exit
  -f | --fwd FWD_FASTQ                      Path to the forward FASTQ file (default: None)
  -r | --rev REV_FASTQ                      Path to the reverse FASTQ file (default: None)
  -m | --merged MERGED_FASTQ                Path to merged FASTQ file (default: None)
  -b | --binsize BINSIZE                    Bin size for analysis (default: BINSIZE = 1000)
  -u | --RMUS-threshold THRESHHOLD          RMUS threshold for analysis (default: THRESHHOLD = 0.15)
  -q | --mapping_quality QUALITY            Mapping quality threshold (default: QUALITY = 20)
  -p | --prefix PREFIX                      Prefix for output files (default: None)
  -t | --threads THREADS                    Number of threads to use (default: THREADS = 10)
  -x | --taxonomic-hierarchy HIERARCHY      Taxonomic hierarchy (default: HIERARCHY = species)
  -c | --force                              Force overwrite of existing output files (default: false)
  -z | --create-db                          Create a new database
  -d | --db-folder DB_FOLDER                Folder for the database

Example:
  ECMSD -f reads_R1.fastq -r reads_R2.fastq -o results/
  ECMSD --create-db --db-folder /path/to/db_folder
  ECMSD -f reads_R1.fastq -o results/ --db-folder /path/to/db_folder
EOF
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
version="1.1.0"
taxonomic_hierarchy="species"
# skip_environment=false
output=""
prefix=""
db_folder=""
create_db=false

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
        exit 0
        ;;
    -z | --create-db)
        create_db=true
        shift
        ;;
    -d | --db-folder)
        db_folder="$2"
        shift 2
        ;;
    *)
        echo "Unknown option: $1"
        usage
        exit 1
        ;;
    esac
done

echo "Starting ECMSD pipeline..."

###############################################################################
#                        Resolve Script Directory                            #
###############################################################################
# Detect if running from a conda install or directly from the repo
# Once installed via conda, scripts live in $CONDA_PREFIX/lib/ecmsd/scripts/
if [[ -d "$(dirname "$(realpath "$0")")/../lib/ecmsd/scripts" ]]; then
    # Running from conda install
    SCRIPT_DIR="$(dirname "$(realpath "$0")")/../lib/ecmsd/scripts"
    SHELL_DIR="$(dirname "$(realpath "$0")")/../lib/ecmsd/shell"
else
    # Running directly from repo
    SCRIPT_DIR="$(dirname "$(realpath "$0")")/../scripts"
    SHELL_DIR="$(dirname "$(realpath "$0")")/../shell"
fi

# Convert db_folder to absolute path without requiring the directory to exist yet
if [[ -n "${db_folder}" && "${db_folder}" != /* ]]; then
    db_folder="$(pwd)/${db_folder}"
fi

###############################################################################
#                        Handle Database Creation Early                      #
###############################################################################
if [[ "${create_db}" == true ]]; then
    if [[ -z "${db_folder}" ]]; then
        echo "Error: --db-folder must be specified with --create-db."
        exit 1
    fi

    echo "Creating database in folder: ${db_folder}"
    mkdir -p "${db_folder}"

    bash "${SHELL_DIR}/MakeRef.sh" "${db_folder}" "${SCRIPT_DIR}" "${threads}"

    echo "Database created successfully in ${db_folder}."
    exit 0 
fi

###############################################################################
#                        Check Database Folder                               #
###############################################################################
if [[ -n "${db_folder}" ]]; then
    # Convert db_folder to an absolute path (works even if it doesn't yet exist)

    echo "Using database folder: ${db_folder}"

    if [[ ! -f "${db_folder}/NCBI_taxdump/nodes.dmp" || \
          ! -f "${db_folder}/NCBI_taxdump/names.dmp" || \
          ! -f "${db_folder}/mitochondrion_refseq_taxid_masked.fna.gz" ]]; then
        echo "Database files missing in ${db_folder}. Generating missing files..."
        bash "${SHELL_DIR}/MakeRef.sh" "${db_folder}" "${SCRIPT_DIR}" "${threads}"
    fi

    REF="${db_folder}/mitochondrion_refseq_taxid_masked.fna.gz"
    NODES="${db_folder}/NCBI_taxdump/nodes.dmp"
    NAMES="${db_folder}/NCBI_taxdump/names.dmp"
else
    echo "No database folder provided. Please provide a database folder using --db-folder or create a new database using --create-db."
    usage
    exit 1
    # DEFAULT_DB="${WD}/data/refseq_mito"
    # REF="${DEFAULT_DB}/mitochondrion_refseq_taxid_masked.fna.gz"
    # NODES="${DEFAULT_DB}/NCBI_taxdump/nodes.dmp"
    # NAMES="${DEFAULT_DB}/NCBI_taxdump/names.dmp"

    # if [[ ! -f "${REF}" || ! -f "${NODES}" || ! -f "${NAMES}" ]]; then
    #     echo "Default database files are missing. Generating database..."
    #     mkdir -p "${DEFAULT_DB}"
    #     bash "${SHELL_DIR}/MakeRef.sh" "${DEFAULT_DB}" "${threads}"
    # fi
fi

###############################################################################
#                        Check Required Arguments                             #
###############################################################################
if [[ ( -z "${fwd}" && -z "${merged}" ) || -z "${output}" ]]; then
    echo "Error: A read input (--fwd or --merged) and --out are required."
    usage
    exit 1
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
#                        Mapping Function                                     #
###############################################################################
run_mapping() {
    # Function to run minimap2 mapping and log the progress in a log file
    # Accepts one or more read file paths as arguments
    echo "Running minimap2 mapping for: $*"
    minimap2 \
        -x sr \
        --secondary=no \
        -t "${threads}" "${REF}" "$@" \
        2> "${Output}/logs/minimap2.log" |
        awk -v Q="${quality}" '$12 >= Q {print}'

    # print log file content
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

# check if fwd is provided
if [[ -n "${fwd}" ]]; then
    if [[ -z "${rev}" || "${rev}" == "no" ]]; then
        echo "Running single-end mapping..."
        run_mapping "${fwd}" | gzip >"${PAF}.gz"
    else
        echo "Running paired-end mapping..."
        run_mapping "${fwd}" "${rev}" | gzip >"${PAF}.gz"
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

###############################################################################
#                        Update Script References                             #
###############################################################################
python "${SCRIPT_DIR}/LinkTaxonomy.py" \
    --Nodes "${NODES}" \
    --Names "${NAMES}" \
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
Rscript "${SCRIPT_DIR}/process_files.R" "${summary_file}" "${Output}" "${taxonomic_hierarchy}" "${prefix}"

###############################################################################
#                        Pipeline Completed                                   #
###############################################################################
echo "ECMSD pipeline completed successfully."
    