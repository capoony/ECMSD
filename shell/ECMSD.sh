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
  -b | --force-rebuild                      Optional: discard an existing database and rebuild it from scratch

ALL ARGUMENTS:
  -h | --help                               Show this help message and exit
  -v | --version                            Show version and exit
  -f | --fwd FWD_FASTQ                      Path to the forward FASTQ file (default: None)
  -r | --rev REV_FASTQ                      Path to the reverse FASTQ file (default: None)
  -m | --merged MERGED_FASTQ                Path to merged FASTQ file (default: None)
  -u | --cov-threshold THRESHOLD            Minimum % of reference covered by reads to retain it (default: 5)
  -n | --top-n N                            Number of top references to generate alignment plots for (default: 25)
  -q | --mapping_quality QUALITY            Mapping quality threshold (default: QUALITY = 20)
  -p | --prefix PREFIX                      Prefix for output files (default: None)
  -t | --threads THREADS                    Number of threads to use (default: THREADS = 10)
  -x | --taxonomic-hierarchy HIERARCHY      Taxonomic hierarchy: subspecies, species, genus, family, order, class, phylum, kingdom or domain (default: HIERARCHY = species)
  -c | --force                              Force overwrite of existing output files (default: false)
  -k | --skip-mapping                       Skip mapping and coverage steps; reuse existing PAF and coverage files (implies --force)
  -z | --create-db                          Create a new database
  -d | --db-folder DB_FOLDER                Folder for the database
  -b | --force-rebuild                      Discard an existing database and rebuild it from scratch (implies --create-db)

Example:
  ECMSD -f reads_R1.fastq -r reads_R2.fastq -o results/ -d /path/to/db
  ECMSD --create-db --db-folder /path/to/db_folder
  ECMSD --create-db --force-rebuild --db-folder /path/to/db_folder
  ECMSD -f reads_R1.fastq -o results/ --db-folder /path/to/db_folder
EOF
}

###############################################################################
#                              Default Values                                 #
###############################################################################
fwd=""
rev=""
merged=""
cov_threshold=5
top_n=25
quality=20
threads=10
force="no"
skip_mapping="no"
version="1.3.2"
taxonomic_hierarchy="species"
output=""
prefix=""
db_folder=""
create_db=false
force_rebuild="no"

###############################################################################
#                             Argument Parsing                                #
###############################################################################
need_value() {
    [[ $# -ge 2 ]] || { echo "Error: option '$1' requires a value." >&2; exit 1; }
}

while [[ $# -gt 0 ]]; do
    case $1 in
    -f | --fwd)
        need_value "$@"
        fwd="$2"
        shift 2
        ;;
    -r | --rev)
        need_value "$@"
        rev="$2"
        shift 2
        ;;
    -m | --merged)
        need_value "$@"
        merged="$2"
        shift 2
        ;;
    -u | --cov-threshold)
        need_value "$@"
        cov_threshold="$2"
        shift 2
        ;;
    -n | --top-n)
        need_value "$@"
        top_n="$2"
        shift 2
        ;;
    -q | --mapping_quality)
        need_value "$@"
        quality="$2"
        shift 2
        ;;
    -t | --threads)
        need_value "$@"
        threads="$2"
        shift 2
        ;;
    -c | --force)
        force="yes"
        shift
        ;;
    -k | --skip-mapping)
        skip_mapping="yes"
        force="yes"
        shift
        ;;
    -x | --taxonomic-hierarchy | --taxonomic_hierarchy)
        need_value "$@"
        taxonomic_hierarchy="$2"
        shift 2
        ;;
    -p | --prefix)
        need_value "$@"
        prefix="$2"
        shift 2
        ;;
    -o | --out)
        need_value "$@"
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
    -b | --force-rebuild)
        force_rebuild="yes"
        shift
        ;;
    -d | --db-folder)
        need_value "$@"
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

for pair in "-u:${cov_threshold}" "-q:${quality}" "-n:${top_n}" "-t:${threads}"; do
    [[ "${pair#*:}" =~ ^[0-9]+(\.[0-9]+)?$ ]] || { echo "Error: ${pair%%:*} requires a number, got '${pair#*:}'." >&2; exit 1; }
done
for pair in "-q:${quality}" "-n:${top_n}" "-t:${threads}"; do
    [[ "${pair#*:}" =~ ^[0-9]+$ ]] || { echo "Error: ${pair%%:*} requires an integer, got '${pair#*:}'." >&2; exit 1; }
done
taxonomic_hierarchy="${taxonomic_hierarchy,,}"
case "${taxonomic_hierarchy}" in
subspecies | species | genus | family | order | class | phylum | kingdom | domain) ;;
*)
    echo "Error: -x must be one of subspecies, species, genus, family, order, class, phylum, kingdom, domain (got '${taxonomic_hierarchy}')." >&2
    exit 1
    ;;
esac

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
# --force-rebuild only ever applies to a database build, so treat it as implying
# --create-db rather than rejecting the combination.
if [[ "${force_rebuild}" == "yes" && "${create_db}" != true ]]; then
    echo "Note: --force-rebuild implies --create-db."
    create_db=true
fi

if [[ "${create_db}" == true ]]; then
    if [[ -z "${db_folder}" ]]; then
        echo "Error: --db-folder must be specified with --create-db."
        exit 1
    fi

    # Database creation is a standalone operation that exits when finished
    if [[ -n "${fwd}" || -n "${rev}" || -n "${merged}" || -n "${output}" ]]; then
        echo "Warning: read and output arguments are ignored during database creation."
        echo "         Run ECMSD again without --create-db/--force-rebuild to analyse samples."
    fi

    echo "Creating database in folder: ${db_folder}"
    mkdir -p "${db_folder}"

    bash "${SHELL_DIR}/MakeRef.sh" "${db_folder}" "${SCRIPT_DIR}" "${threads}" "${force_rebuild}"

    echo "Database ready in ${db_folder}."
    exit 0
fi

###############################################################################
#                        Check Database Folder                               #
###############################################################################
if [[ -n "${db_folder}" ]]; then
    echo "Using database folder: ${db_folder}"

    if [[ ! -f "${db_folder}/NCBI_taxdump/nodes.dmp" || \
          ! -f "${db_folder}/NCBI_taxdump/names.dmp" || \
          ! -f "${db_folder}/mitochondrion_refseq_taxid_masked.fna.gz" ]]; then
        echo "Database not found. Run \`ECMSD --create-db\` first" >&2
        exit 1
    fi

    REF="${db_folder}/mitochondrion_refseq_taxid_masked.fna.gz"
    NODES="${db_folder}/NCBI_taxdump/nodes.dmp"
    NAMES="${db_folder}/NCBI_taxdump/names.dmp"
else
    echo "No database folder provided. Please provide a database folder using --db-folder or create a new database using --create-db."
    usage
    exit 1
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
[[ -n "${fwd}" && ! -f "${fwd}" ]] && {
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

if [[ "${skip_mapping}" == "yes" && -d "${Output}" ]]; then
    echo "Skipping mapping — reusing existing PAF in: ${Output}/mapping"
    rm -rf "${Output}/logs"
elif [[ "${force}" == "yes" && -d "${Output}" ]]; then
    echo "Removing existing ECMSD output in: ${Output}"
    rm -rf "${Output}/mapping" "${Output}/logs"
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
mkdir -p "${Output}/logs"

###############################################################################
#                        Mapping Function                                     #
###############################################################################
run_mapping() {
    # Accepts one or more read file paths as arguments
    echo "Running minimap2 mapping for: $*" >&2
    minimap2 \
        -x sr \
        --secondary=no \
        -t "${threads}" "${REF}" "$@" \
        2> "${Output}/logs/minimap2.log" |
        awk -v Q="${quality}" '$12 >= Q {print}'

    echo "Minimap2 log:" >&2
    cat "${Output}/logs/minimap2.log" >&2
}

###############################################################################
#                        Mapping Step                                         #
###############################################################################
echo "Starting mapping step..."

PAF="${Output}/mapping/Mito.paf"
if [[ -n "${prefix}" ]]; then
    PAF="${Output}/mapping/${prefix}_Mito.paf"
fi

if [[ "${skip_mapping}" == "yes" ]]; then
    echo "Skipping mapping step — using existing ${PAF}.gz"
    if [[ ! -f "${PAF}.gz" ]]; then
        echo "Error: Expected PAF file '${PAF}.gz' not found. Cannot skip mapping."
        exit 1
    fi
else
    if [[ -n "${merged}" ]]; then
        echo "Running mapping on merged reads..."
        run_mapping "${merged}" | gzip >"${PAF}.gz"
    elif [[ -z "${rev}" ]]; then
        echo "Running single-end mapping..."
        run_mapping "${fwd}" | gzip >"${PAF}.gz"
    else
        echo "Running paired-end mapping..."
        run_mapping "${fwd}" "${rev}" | gzip >"${PAF}.gz"
    fi
fi

if [[ ! -f "${PAF}.gz" ]]; then
    echo "Error: PAF file '${PAF}.gz' was not created."
    exit 1
fi

###############################################################################
#                        Coverage Calculation                                 #
###############################################################################
COVERAGE_OUT="${Output}/mapping/Mito_coverage.txt"
if [[ -n "${prefix}" ]]; then
    COVERAGE_OUT="${Output}/mapping/${prefix}_Mito_coverage.txt"
fi

if [[ "${skip_mapping}" == "yes" && -f "${COVERAGE_OUT}" ]]; then
    echo "Skipping coverage calculation — reusing existing ${COVERAGE_OUT}"
else
    echo "Calculating per-reference coverage..."
    python3 "${SCRIPT_DIR}/coverage.py" "${PAF}.gz" "${COVERAGE_OUT}" "${quality}"
fi

###############################################################################
#                        Parse PAF and Assign Taxonomy                        #
###############################################################################
output_base="${Output}/mapping/Mito_summary"
if [[ -n "${prefix}" ]]; then
    output_base="${Output}/mapping/${prefix}_Mito_summary"
fi

echo "Parsing PAF and assigning taxonomy..."
python3 "${SCRIPT_DIR}/LinkTaxonomy.py" \
    --Nodes "${NODES}" \
    --Names "${NAMES}" \
    --PAF "${PAF}.gz" \
    --Coverage "${COVERAGE_OUT}" \
    --CovThreshold "${cov_threshold}" \
    --TaxonomicHierarchy "${taxonomic_hierarchy}" \
    --MapQuality "${quality}" \
    --output "${output_base}"

###############################################################################
#                        Plotting Results                                     #
###############################################################################
summary_file="${output_base}.txt"

if [[ ! -f "${summary_file}" ]]; then
    echo "Error: Mito summary file '${summary_file}' was not created."
    exit 1
fi

if [[ $(wc -l <"${summary_file}") -le 1 ]]; then
    echo "Warning: No references passed the coverage threshold (--cov-threshold ${cov_threshold}). Plots will be empty."
fi

echo "Plotting results..."
Rscript "${SCRIPT_DIR}/process_files.R" "${Output}" "${taxonomic_hierarchy}" "${prefix}"

###############################################################################
#                        PAF Alignment Coverage Plots                         #
###############################################################################
echo "Generating per-reference PAF alignment plots (top ${top_n} references)..."
Rscript "${SCRIPT_DIR}/plot_paf_alignments.R" \
    "${Output}" \
    "${PAF}.gz" \
    "${output_base}.ref_summary.txt" \
    "${top_n}" \
    "${taxonomic_hierarchy}"

###############################################################################
#                        Pipeline Completed                                   #
###############################################################################
echo "ECMSD pipeline completed successfully."
