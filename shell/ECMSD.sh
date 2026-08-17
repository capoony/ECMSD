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
  -u | --cov-threshold THRESHOLD            Minimum % of reference covered by reads to retain it (default: 25)
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
cov_threshold=25
top_n=25
quality=20
threads=10
force="no"
skip_mapping="no"
version="1.3.0"
taxonomic_hierarchy="species"
output=""
prefix=""
db_folder=""
create_db=false
force_rebuild="no"

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
    -u | --cov-threshold)
        cov_threshold="$2"
        shift 2
        ;;
    -n | --top-n)
        top_n="$2"
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
    -k | --skip-mapping)
        skip_mapping="yes"
        force="yes"
        shift
        ;;
    -x | --taxonomic-hierarchy | --taxonomic_hierarchy)
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
    -b | --force-rebuild)
        force_rebuild="yes"
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
        echo "Database files missing in ${db_folder}. Generating missing files..."
        # Never force a rebuild here: an analysis run must not discard a database
        # that parallel jobs may be reading from.
        bash "${SHELL_DIR}/MakeRef.sh" "${db_folder}" "${SCRIPT_DIR}" "${threads}" "no"
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
mkdir -p "${Output}/logs"

###############################################################################
#                        Mapping Function                                     #
###############################################################################
run_mapping() {
    # Accepts one or more read file paths as arguments
    echo "Running minimap2 mapping for: $*"
    minimap2 \
        -x sr \
        --secondary=no \
        -t "${threads}" "${REF}" "$@" \
        2> "${Output}/logs/minimap2.log" |
        awk -v Q="${quality}" '$12 >= Q {print}'

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

if [[ "${skip_mapping}" == "yes" ]]; then
    echo "Skipping mapping step — using existing ${PAF}.gz"
    if [[ ! -f "${PAF}.gz" ]]; then
        echo "Error: Expected PAF file '${PAF}.gz' not found. Cannot skip mapping."
        exit 1
    fi
else
    if [[ -n "${fwd}" ]]; then
        if [[ -z "${rev}" || "${rev}" == "no" ]]; then
            echo "Running single-end mapping..."
            run_mapping "${fwd}" | gzip >"${PAF}.gz"
        else
            echo "Running paired-end mapping..."
            run_mapping "${fwd}" "${rev}" | gzip >"${PAF}.gz"
        fi
    fi

    if [[ -n "${merged}" && "${merged}" != "no" ]]; then
        if [[ ! -f "${merged}" ]]; then
            echo "Error: Merged reads file '${merged}' does not exist."
            exit 1
        fi
        echo "Running mapping on merged reads..."
        run_mapping "${merged}" | gzip >>"${PAF}.gz"
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
    python3 - <<'PYEOF' "${PAF}.gz" "${COVERAGE_OUT}" "${quality}"
import sys
import gzip
import math
from collections import defaultdict

paf_file = sys.argv[1]
out_file = sys.argv[2]
min_mapq = int(sys.argv[3])

diff_arrays = defaultdict(lambda: defaultdict(int))
ref_lengths = {}
seen_read_ref = set()

opener = gzip.open if paf_file.endswith(".gz") else open

with opener(paf_file, "rt") as fh:
    for line in fh:
        parts = line.strip().split("\t")
        if len(parts) < 12:
            continue
        read_name = parts[0]
        ref_name  = parts[5]
        ref_len   = int(parts[6])
        ref_start = int(parts[7])
        ref_end   = int(parts[8])
        mapq      = int(parts[11])
        ref_lengths[ref_name] = ref_len

        if mapq < min_mapq:
            continue

        key = (read_name, ref_name)
        if key in seen_read_ref:
            continue
        seen_read_ref.add(key)

        diff_arrays[ref_name][ref_start] += 1
        diff_arrays[ref_name][ref_end]   -= 1

with open(out_file, "w") as out:
    out.write("reference\tref_length\tmean_coverage\tstd_coverage\tpct_covered\n")
    for ref_name, ref_len in sorted(ref_lengths.items()):
        diff  = diff_arrays[ref_name]
        depth = 0
        covered = 0
        sum_d   = 0
        sum_sq_d = 0
        for pos in range(ref_len):
            depth    += diff.get(pos, 0)
            sum_d    += depth
            sum_sq_d += depth * depth
            if depth > 0:
                covered += 1
        mean_cov    = sum_d / ref_len if ref_len > 0 else 0
        variance    = (sum_sq_d / ref_len - mean_cov ** 2) if ref_len > 0 else 0
        std_cov     = math.sqrt(max(variance, 0.0))
        pct_covered = (covered / ref_len * 100) if ref_len > 0 else 0
        out.write(f"{ref_name}\t{ref_len}\t{mean_cov:.2f}\t{std_cov:.2f}\t{pct_covered:.2f}\n")

print(f"Coverage statistics written to: {out_file}")
PYEOF
fi

###############################################################################
#                        Parse PAF and Assign Taxonomy                        #
###############################################################################
output_base="${Output}/mapping/Mito_summary"
if [[ -n "${prefix}" ]]; then
    output_base="${Output}/mapping/${prefix}_Mito_summary"
fi

echo "Parsing PAF and assigning taxonomy..."
python "${SCRIPT_DIR}/LinkTaxonomy.py" \
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
    echo "Warning: No references passed the coverage threshold (--cov-threshold ${cov_threshold}). Skipping plots."
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
