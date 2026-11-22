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
  -b | --force-rebuild                      Optional: rebuild the database from scratch even if it is current

MANAGING THE DATABASE:
     | --update-db                          Rebuild against the current NCBI release if it is newer; extra references are kept
     | --extend-db --fasta FASTA            Add mitogenomes with headers >taxid|id
     | --extend-db --accessions FILE        Add mitogenomes by NCBI accession (one per line), fetched with their taxids
ALL ARGUMENTS:
  -h | --help                               Show this help message and exit
  -v | --version                            Show version and exit
  -f | --fwd FWD_FASTQ                      Path to the forward FASTQ file (default: None)
  -r | --rev REV_FASTQ                      Path to the reverse FASTQ file (default: None)
  -m | --merged MERGED_FASTQ                Path to merged FASTQ file (default: None)
  -u | --cov-threshold THRESHOLD            Minimum % of reference covered by reads to retain it (default: 25)
  -n | --top-n N                            Number of top references to generate alignment plots for (default: 25)
  -q | --mapping_quality QUALITY            Drop reads whose best alignments all have a lower MAPQ (default: 0 = no MAPQ filtering)
     | --seed-k K                           minimap2 k-mer size (default: 17)
     | --seed-w W                           minimap2 minimizer window size (default: 7)
     | --min-chain-score M                   minimap2 minimal chaining score -m; higher values align fewer weak reads and map faster (default: 20)
     | --max-alns N                          Maximum number of secondary alignments per read, minimap2 -N (default: 20, or 20 + the largest number of references of one taxid in an extended database)
     | --min-disc-reads N                    A reference is only reported when at least N reads discriminate it from all other references (default: 20)
     | --min-disc-breadth PCT                ... and these reads cover at least PCT % of the reference (default: 5)
     | --min-disc-evenness E                 ... and are spread over it with an evenness of at least E, 0-1 (default: 0.5)
     | --min-disc-window-ratio R             ... and hit at least R times the 1 kb windows that as many random reads would hit; 0 disables (default: 0.5)
     | --min-decline D                       Label a reported reference poor_fit when the edit distance decline of its discriminatory reads is below D (default: 0.9)
     | --min-breadth-ratio R                 ... or when their breadth is below R times the breadth expected from their depth; 0 disables (default: 0.5)
     | --max-overlap J                       Label reported references that share more than J (Jaccard) of their compatible reads as one ambiguous group (default: 0.5)
     | --min-aln-len N                       Minimum alignment length in bp of a read that counts as discriminatory evidence; lower it for very short reads (default: 20)
     | --no-dedup                            Count PCR duplicates (same references, start, end and strand) as separate reads in the evidence table
     | --group-by-taxid MODE                 auto, yes or no: count all references of one taxid (FASTA headers <taxid>|<accession>, e.g. several haplotypes of a species) as one unit for evidence, credibility and assignment; auto = yes when the database has several references of a taxid (default: auto)
     | --keep-all-alignments                 Also keep every alignment (all_alignments PAF and a table with tp, AS, NM, MAPQ); large and slower
  -p | --prefix PREFIX                      Prefix for output files (default: None)
  -t | --threads THREADS                    Number of threads to use (default: THREADS = 10)
  -x | --taxonomic-hierarchy HIERARCHY      Taxonomic hierarchy: subspecies, species, genus, family, order, class, phylum, kingdom or domain (default: HIERARCHY = species)
  -c | --force                              Force overwrite of existing output files (default: false)
  -k | --skip-mapping                       Skip mapping and coverage steps; reuse existing PAF and coverage files (implies --force)
  -z | --create-db                          Create a new database
  -d | --db-folder DB_FOLDER                Folder for the database
  -b | --force-rebuild                      Rebuild the database from scratch even if it is current (implies --create-db)
     | --update-db                          Update the database to the current NCBI release (built in DB_FOLDER/update.new, then swapped in)
     | --extend-db                          Add references from --fasta FASTA (headers >taxid|id) or --accessions FILE

All reads are mapped once with base-level alignment (-c) and secondary alignments, and all
best-scoring alignments of each read are kept by alignment score (AS), not by MAPQ. A
reference is only reported when reads that discriminate it from the other references
support it (see --min-disc-*), not when it merely shares reads with a reference that is
present.

Example:
  ECMSD -f reads_R1.fastq -r reads_R2.fastq -o results/ -d /path/to/db
  ECMSD --create-db --db-folder /path/to/db_folder
  ECMSD --create-db --force-rebuild --db-folder /path/to/db_folder
  ECMSD --update-db --db-folder /path/to/db_folder
  ECMSD --extend-db --accessions accessions.txt --db-folder /path/to/db_folder
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
quality=0
seed_k=17
seed_w=7
max_alns=20
max_alns_set="no"
chain_score=20
keep_all="no"
min_disc_reads=20
min_disc_breadth=5
min_disc_evenness=0.5
min_disc_window_ratio=0.5
min_decline=0.9
min_breadth_ratio=0.5
max_overlap=0.5
min_aln_len=20
dedup="yes"
group_by_taxid="auto"
threads=10
force="no"
skip_mapping="no"
version="1.3.3"
taxonomic_hierarchy="species"
output=""
prefix=""
db_folder=""
create_db=false
force_rebuild="no"
update_db=false
extend_db=false
ext_kind=""
ext_input=""

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
    --seed-k)
        need_value "$@"
        seed_k="$2"
        shift 2
        ;;
    --seed-w)
        need_value "$@"
        seed_w="$2"
        shift 2
        ;;
    --max-alns)
        need_value "$@"
        max_alns="$2"
        max_alns_set="yes"
        shift 2
        ;;
    --min-chain-score)
        need_value "$@"
        chain_score="$2"
        shift 2
        ;;
    --min-disc-reads)
        need_value "$@"
        min_disc_reads="$2"
        shift 2
        ;;
    --min-disc-breadth)
        need_value "$@"
        min_disc_breadth="$2"
        shift 2
        ;;
    --min-disc-evenness)
        need_value "$@"
        min_disc_evenness="$2"
        shift 2
        ;;
    --min-disc-window-ratio)
        need_value "$@"
        min_disc_window_ratio="$2"
        shift 2
        ;;
    --min-decline)
        need_value "$@"
        min_decline="$2"
        shift 2
        ;;
    --min-breadth-ratio)
        need_value "$@"
        min_breadth_ratio="$2"
        shift 2
        ;;
    --max-overlap)
        need_value "$@"
        max_overlap="$2"
        shift 2
        ;;
    --min-aln-len)
        need_value "$@"
        min_aln_len="$2"
        shift 2
        ;;
    --no-dedup)
        dedup="no"
        shift
        ;;
    --group-by-taxid)
        need_value "$@"
        group_by_taxid="$2"
        shift 2
        ;;
    --keep-all-alignments)
        keep_all="yes"
        shift
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
    --update-db)
        update_db=true
        shift
        ;;
    --extend-db)
        extend_db=true
        shift
        ;;
    --fasta | --accessions)
        need_value "$@"
        ext_kind="${1#--}"
        ext_input="$2"
        shift 2
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

for pair in "-u:${cov_threshold}" "-q:${quality}" "-n:${top_n}" "-t:${threads}" "--min-disc-reads:${min_disc_reads}" "--min-disc-breadth:${min_disc_breadth}" "--min-disc-evenness:${min_disc_evenness}" "--min-disc-window-ratio:${min_disc_window_ratio}" "--min-decline:${min_decline}" "--min-breadth-ratio:${min_breadth_ratio}" "--max-overlap:${max_overlap}"; do
    [[ "${pair#*:}" =~ ^[0-9]+(\.[0-9]+)?$ ]] || { echo "Error: ${pair%%:*} requires a number, got '${pair#*:}'." >&2; exit 1; }
done
for pair in "-q:${quality}" "-n:${top_n}" "-t:${threads}" "--seed-k:${seed_k}" "--seed-w:${seed_w}" "--max-alns:${max_alns}" "--min-chain-score:${chain_score}" "--min-disc-reads:${min_disc_reads}" "--min-aln-len:${min_aln_len}"; do
    [[ "${pair#*:}" =~ ^[0-9]+$ ]] || { echo "Error: ${pair%%:*} requires an integer, got '${pair#*:}'." >&2; exit 1; }
done
case "${group_by_taxid}" in
auto | yes | no) ;;
*)
    echo "Error: --group-by-taxid must be auto, yes or no (got '${group_by_taxid}')." >&2
    exit 1
    ;;
esac
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

    mode="build"
    [[ "${force_rebuild}" == "yes" ]] && mode="rebuild"
    ECMSD_VERSION="${version}" bash "${SHELL_DIR}/MakeRef.sh" "${db_folder}" "${SCRIPT_DIR}" "${threads}" "${mode}"

    echo "Database ready in ${db_folder}."
    exit 0
fi

# Updating and extending change an existing database and exit when finished
if [[ "${update_db}" == true || "${extend_db}" == true ]]; then
    if [[ "${update_db}" == true && "${extend_db}" == true ]]; then
        echo "Error: run --update-db and --extend-db separately." >&2
        exit 1
    fi
    if [[ -z "${db_folder}" || ! -d "${db_folder}" ]]; then
        echo "Error: --db-folder must name an existing database folder." >&2
        exit 1
    fi
    if [[ "${update_db}" == true ]]; then
        ECMSD_VERSION="${version}" bash "${SHELL_DIR}/MakeRef.sh" "${db_folder}" "${SCRIPT_DIR}" "${threads}" update
    else
        if [[ -z "${ext_input}" || ! -f "${ext_input}" ]]; then
            echo "Error: --extend-db needs an existing --fasta or --accessions file." >&2
            exit 1
        fi
        ECMSD_VERSION="${version}" bash "${SHELL_DIR}/MakeRef.sh" "${db_folder}" "${SCRIPT_DIR}" "${threads}" extend \
            "${ext_kind}" "$(realpath "${ext_input}")"
    fi
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
    # Reference FASTA files: RefSeq, plus the references added with --extend-db
    REFS=("${REF}")
    [[ -f "${db_folder}/extra_refs.fna.gz" ]] && REFS+=("${db_folder}/extra_refs.fna.gz")
    DB_INFO="${db_folder}/db_info.txt"
    if [[ -f "${DB_INFO}" ]]; then
        echo "Database version:"
        sed 's/^/  /' "${DB_INFO}"
    else
        echo "Database version: unknown (built before ECMSD recorded it; see --update-db)"
    fi
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
#                        Mapping Functions                                    #
###############################################################################
# Chaining thresholds. The sr preset requires a chaining score of 25 (-m) and a peak
# score of 40 (-s), which discards most reads below ~30 bp. -m 20 (default,
# --min-chain-score) keeps short historical DNA fragments while reads that can only
# form a very weak chain are not aligned; -n 1 leaves -m as the only chaining limit
# and -s 20 does not cut short but perfect alignments.
MM2_SHORT=(-n 1 -m "${chain_score}" -s 20)

# Single mapping pass of the original read files against the full database.
# Base-level alignment (-c) and secondary alignments; nothing is removed by MAPQ.
# With extra references, both FASTA files are streamed to minimap2 as one target on
# stdin (joined gzip files are one gzip stream); the reads are always files.
run_mapping() {
    echo "Mapping reads: $*" >&2
    local mm2=(minimap2 -x sr -k "${seed_k}" -w "${seed_w}" "${MM2_SHORT[@]}"
        --secondary=yes -N "${max_alns}" -c -t "${threads}")
    if (( ${#REFS[@]} > 1 )); then
        cat "${REFS[@]}" | "${mm2[@]}" - "$@" 2> "${Output}/logs/minimap2.log"
    else
        "${mm2[@]}" "${REF}" "$@" 2> "${Output}/logs/minimap2.log"
    fi

    echo "Minimap2 log:" >&2
    cat "${Output}/logs/minimap2.log" >&2
}

###############################################################################
#                        References per Taxid                                 #
###############################################################################
# Several references of one taxid only exist in databases with <taxid>|<accession>
# headers (e.g. haplotypes added with --extend-db). They turn on --group-by-taxid auto,
# and --max-alns must leave room for all of them plus 20 competitors, or the
# discriminatory reads of the taxid hit the minimap2 -N limit.
max_per_taxid=$(gzip -dc "${REFS[@]}" | grep '^>' | cut -d'|' -f1 | sort | uniq -c | awk 'm < $1 { m = $1 } END { print m + 0 }')
if [[ "${group_by_taxid}" == "auto" ]]; then
    group_by_taxid="no"
    (( max_per_taxid > 1 )) && group_by_taxid="yes"
fi
echo "Grouping references by taxid: ${group_by_taxid}"
if [[ "${max_alns_set}" == "no" ]] && (( max_per_taxid > 1 )); then
    max_alns=$(( 20 + max_per_taxid ))
    echo "Up to ${max_per_taxid} references per taxid: --max-alns set to ${max_alns}"
fi

###############################################################################
#                        Mapping Step                                         #
###############################################################################
echo "Starting mapping step..."

PAF="${Output}/mapping/Mito.paf"
if [[ -n "${prefix}" ]]; then
    PAF="${Output}/mapping/${prefix}_Mito.paf"
fi
# Unfiltered mapping output and alignment table (only with --keep-all-alignments)
ALL_ALNS="${PAF%.paf}_all_alignments.paf.gz"
ALN_TABLE="${PAF%.paf}_alignments.tsv.gz"

if [[ "${skip_mapping}" == "yes" ]]; then
    echo "Skipping mapping step — using existing ${PAF}.gz"
    if [[ ! -f "${PAF}.gz" ]]; then
        echo "Error: Expected PAF file '${PAF}.gz' not found. Cannot skip mapping."
        exit 1
    fi
else
    # Original read files (merged and single-end: one file; paired-end: two)
    if [[ -n "${merged}" ]]; then
        echo "Mapping merged reads..."
        reads_in=("${merged}")
    elif [[ -z "${rev}" ]]; then
        echo "Mapping single-end reads..."
        reads_in=("${fwd}")
    else
        echo "Mapping paired-end reads..."
        reads_in=("${fwd}" "${rev}")
    fi

    # Keep all best-scoring (by AS) alignments of each read. By default the mapping
    # output (tens of millions of lines for deep samples) is streamed straight into
    # the best-alignment selection and never written to disk.
    if [[ "${keep_all}" == "yes" ]]; then
        run_mapping "${reads_in[@]}" | gzip >"${ALL_ALNS}"
        python3 "${SCRIPT_DIR}/select_best_alignments.py" --min-mapq "${quality}" \
            "${ALL_ALNS}" "${PAF}.gz" --table "${ALN_TABLE}"
    else
        run_mapping "${reads_in[@]}" |
            python3 "${SCRIPT_DIR}/select_best_alignments.py" --min-mapq "${quality}" - "${PAF}.gz"
    fi
    if [[ -z "$(gzip -cd "${PAF}.gz" | head -c 1)" ]]; then
        echo "Warning: no read aligned to the mitochondrial database."
    fi
fi

if [[ ! -f "${PAF}.gz" ]]; then
    echo "Error: PAF file '${PAF}.gz' was not created."
    exit 1
fi

###############################################################################
#                        Reference Evidence                                   #
###############################################################################
# A reference is only reported when reads discriminate it from the other references,
# not when it merely shares reads with a reference that is present (Step A). The
# coverage of the credible references is then calculated from all compatible
# alignments, including reads that tie with other references (Step B).
EVIDENCE_OUT="${Output}/mapping/Mito_reference_evidence.txt"
if [[ -n "${prefix}" ]]; then
    EVIDENCE_OUT="${Output}/mapping/${prefix}_Mito_reference_evidence.txt"
fi
echo "Assessing discriminatory evidence per reference..."
python3 "${SCRIPT_DIR}/assess_references.py" \
    --PAF "${PAF}.gz" \
    --Names "${NAMES}" \
    --output "${EVIDENCE_OUT}" \
    --max-alns "${max_alns}" \
    --min-disc-reads "${min_disc_reads}" \
    --min-disc-breadth "${min_disc_breadth}" \
    --min-disc-evenness "${min_disc_evenness}" \
    --min-window-ratio "${min_disc_window_ratio}" \
    --min-decline "${min_decline}" \
    --min-breadth-ratio "${min_breadth_ratio}" \
    --max-overlap "${max_overlap}" \
    --min-aln-len "${min_aln_len}" \
    $([[ "${dedup}" == "no" ]] && echo --no-dedup) \
    $([[ "${group_by_taxid}" == "yes" ]] && echo --group-by-taxid)

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
    echo "Calculating per-reference coverage (all compatible alignments)..."
    # -q was already applied to whole reads by select_best_alignments.py (best MAPQ of the
    # read's compatible alignments); applying it again per alignment would drop tied
    # alignments of reads that were kept, so no further MAPQ filter here.
    python3 "${SCRIPT_DIR}/coverage.py" "${PAF}.gz" "${COVERAGE_OUT}" 0
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
    --Evidence "${EVIDENCE_OUT}" \
    --CovThreshold "${cov_threshold}" \
    --TaxonomicHierarchy "${taxonomic_hierarchy}" \
    --MapQuality "${quality}" \
    $([[ "${group_by_taxid}" == "yes" ]] && echo --GroupByTaxid) \
    --output "${output_base}"

###############################################################################
#                        Plotting Results                                     #
###############################################################################
summary_file="${output_base}.txt"

if [[ ! -f "${summary_file}" ]]; then
    echo "Error: Mito summary file '${summary_file}' was not created."
    exit 1
fi

# One row per call with a status (confident, ambiguous, poor_fit, low_evidence).
python3 "${SCRIPT_DIR}/call_species.py" \
    --Evidence "${EVIDENCE_OUT}" \
    --Summary "${output_base}.txt" \
    --output "${output_base%_summary}_calls.txt" \
    $([[ -f "${DB_INFO}" ]] && echo --DbInfo "${DB_INFO}")

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
    "${taxonomic_hierarchy}" \
    "${output_base%_summary}_calls.txt"

###############################################################################
#                        Pipeline Completed                                   #
###############################################################################
echo "ECMSD pipeline completed successfully."
