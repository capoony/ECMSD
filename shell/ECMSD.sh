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
  -u | --cov-threshold   Minimum percentage of reference covered by reads (0-100) for a reference to be retained (default: 50)
  -n | --top-n           Number of top references to generate alignment plots for (default: 25)
  -q | --mapping_quality Mapping quality threshold (default: 20)
  -t | --threads     Number of threads to use (default: 10)
  -c | --force       Force overwrite of existing output files
  -k | --skip-mapping  Skip mapping step and reuse existing PAF file (implies --force)
  -x | --taxonomic-hierarchy Taxonomic hierarchy for classification; one of: species, genus, family, order, phylum, kingdom (default: species)
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
cov_threshold=50
top_n=100
quality=10
threads=10
force="no"
skip_mapping="no"
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
        force="yes"   # must be yes so the existing-dir check doesn't abort
        shift
        ;;
    -x | --taxonomic-hierarchy | --taxonomic_hierarchy)
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

if [[ "${skip_mapping}" == "yes" && -d "${Output}" ]]; then
    echo "Skipping mapping — reusing existing PAF in: ${Output}/mapping"
    rm -rf "${Output}/logs"
elif [[ "${force}" == "yes" && -d "${Output}" ]]; then
    echo "Clearing mapping and log directories in: ${Output}"
    rm -rf "${Output}/mapping" "${Output}/logs"
elif [[ -d "${Output}" ]]; then
    echo "Output directory '${Output}' already exists. Use -c or --force to overwrite."
    exit 1
fi

###############################################################################
#                        Find Base Directory                                  #
###############################################################################
WD="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

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
if [[ ! -d "${WD}/scripts/conda_env" ]]; then 
    bash "${WD}/shell/requirements.sh" "${WD}" "${Output}"
fi
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
    # Function to run minimap2 mapping (primary alignments only) and log the progress
    local reads="$1"
    echo "Running minimap2 mapping for: ${reads}"
    minimap2 \
        -x sr \
        --secondary=no \
        -t "${threads}" "${REF}" ${reads} \
        2>>"${Output}/logs/logfile.log" |
        awk -v Q="${quality}" '$12 >= Q {print}'
}

###############################################################################
#                        Mapping Step                                         #
###############################################################################
PAF="${Output}/mapping/Mito.paf"
REF="${WD}/data/refseq_mito/mitochondrion_refseq_taxid_masked.fna.gz"

if [[ "${skip_mapping}" == "yes" ]]; then
    echo "Skipping mapping step — using existing ${PAF}.gz"
    if [[ ! -f "${PAF}.gz" ]]; then
        echo "Error: Expected PAF file '${PAF}.gz' not found. Cannot skip mapping."
        exit 1
    fi
else
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

    ###########################################################################
    #                    Mapping on Merged Reads                              #
    ###########################################################################
    if [[ -n "${merged}" && "${merged}" != "no" ]]; then
        if [[ ! -f "${merged}" ]]; then
            echo "Error: Merged reads file '${merged}' does not exist."
            exit 1
        fi
        echo "Running mapping on merged reads..."
        run_mapping "${merged}" | gzip >>"${PAF}.gz"
    fi
fi



###############################################################################
#                        Coverage Calculation                                 #
###############################################################################
COVERAGE_OUT="${Output}/mapping/Mito_coverage.txt"

if [[ "${skip_mapping}" == "yes" && -f "${COVERAGE_OUT}" ]]; then
    echo "Skipping coverage calculation — reusing existing ${COVERAGE_OUT}"
else
    echo "Calculating coverage of mapped reads on reference..."
    python3 - <<'PYEOF' "${PAF}.gz" "${COVERAGE_OUT}" "${quality}"
import sys
import gzip
import math
from collections import defaultdict

paf_file    = sys.argv[1]
out_file    = sys.argv[2]
min_mapq    = int(sys.argv[3])

# Use difference arrays for O(n + ref_len) coverage calculation instead of O(n * read_len)
# diff_arrays[ref_name][pos] stores the increment/decrement for the running depth
diff_arrays = defaultdict(lambda: defaultdict(int))
ref_lengths = {}

# Track (read_name, ref_name) pairs already counted to avoid double-counting
# a read that has both a primary and a secondary alignment to the same reference
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

        # Enforce MQ threshold on all (primary-only) alignments
        if mapq < min_mapq:
            continue

        # Skip if this read has already been counted on this reference
        key = (read_name, ref_name)
        if key in seen_read_ref:
            continue
        seen_read_ref.add(key)

        # Difference array: +1 at start, -1 at end
        diff_arrays[ref_name][ref_start] += 1
        diff_arrays[ref_name][ref_end]   -= 1

with open(out_file, "w") as out:
    out.write("reference\tref_length\tmean_coverage\tstd_coverage\tpct_covered\n")
    for ref_name, ref_len in sorted(ref_lengths.items()):
        diff = diff_arrays[ref_name]
        # Reconstruct depth via running sum over difference array
        depth      = 0
        covered    = 0
        sum_d      = 0
        sum_sq_d   = 0
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
fi  # end skip_mapping coverage guard

###############################################################################
#                        Parse PAF and Calculate Coverage                         #
###############################################################################
echo "Parsing PAF and calculating RMUS/taxonomic proportions..."
python "${WD}/scripts/LinkTaxonomy.py" \
    --Nodes "${WD}/data/refseq_mito/NCBI_taxdump/nodes.dmp" \
    --Names "${WD}/data/refseq_mito/NCBI_taxdump/names.dmp" \
    --PAF "${PAF}.gz" \
    --Coverage "${COVERAGE_OUT}" \
    --CovThreshold "${cov_threshold}" \
    --TaxonomicHierarchy "${taxonomic_hierarchy}" \
    --MapQuality "${quality}" \
    --output "${Output}/mapping/Mito_summary"
###############################################################################
#                        Plotting Results                                     #
###############################################################################
echo "Plotting results..."
Rscript "${WD}/scripts/process_files.R" ${Output} ${taxonomic_hierarchy} 2>&1 | grep -v 'Fontconfig error'

###############################################################################
#                        PAF Alignment Coverage Plots                         #
###############################################################################
echo "Generating per-reference PAF alignment plots (top ${top_n} references)..."
Rscript "${WD}/scripts/plot_paf_alignments.R" \
    "${Output}" \
    "${PAF}.gz" \
    "${Output}/mapping/Mito_summary.ref_summary.txt" \
    "${top_n}" 2>&1 | grep -v 'Fontconfig error'

###############################################################################
#                        Pipeline Completed                                   #
###############################################################################
echo "ECMSD pipeline completed successfully."
