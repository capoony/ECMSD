#!/bin/bash
set -euo pipefail

#set working directory
# This script builds a minimap2 reference database for mitochondrial genomes using RefSeq data.
WD=$1
SCRIPT_DIR=$2
threads=$3
# Mode: build (default), rebuild (discard and build from scratch), extend (add the
# references of $6, a FASTA or accession list as named by $5) or update (rebuild against
# the current NCBI release in a staging folder and swap the files in at the end).
# "yes"/"no" are the old force-rebuild values.
mode=${4:-build}
case "${mode}" in
yes) mode=rebuild ;;
no) mode=build ;;
esac

# The only three files the analysis pipeline needs from the database folder
# (see ECMSD.sh: REF, NODES, NAMES; extra_refs.fna.gz and db_info.txt are optional). Everything else produced below is a build
# intermediate and is removed by cleanup_intermediates() once the build succeeds.
REF_FILE="mitochondrion_refseq_taxid_masked.fna.gz"
NODES_FILE="NCBI_taxdump/nodes.dmp"
NAMES_FILE="NCBI_taxdump/names.dmp"
# Optional extra references added with --extend-db (headers <taxid>|<accession>)
EXTRA_FILE="extra_refs.fna.gz"
EXT_LOG="extensions.log"
INFO_FILE="db_info.txt"
# The full NCBI taxonomy, kept gzipped so extensions can be re-trimmed offline
FULL_DIR="NCBI_taxdump/full"
FULL_NODES="${FULL_DIR}/nodes.dmp.gz"
FULL_NAMES="${FULL_DIR}/names.dmp.gz"
FULL_MERGED="${FULL_DIR}/merged.dmp.gz"
RELEASE_URL="https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"
# Staging folder of --update-db inside the database folder (same file system, so the
# finished files are moved into place, never copied)
STAGE_DIR="update.new"

# Build intermediates, removed once the database is complete. The compressed
# unmasked FASTA is listed for databases built by earlier script versions.
INTERMEDIATE_FILES=(
    mitochondrion_refseq.fa.gz
    nucl_gb.accession2taxid.gz
    nucl_refseq.accession2taxid.tsv
    mitochondrion_refseq_taxid.fna
    mitochondrion_refseq_taxid.fna.gz
    mitochondrion_refseq_taxid_masked.fna
    refseq_release.txt
    refseq_downloaded.txt
)

echo "*********************"
echo "Setting up RefSeq mitochondrial genome database... "
echo "Have a cup of coffee, this may take a while... "
echo '''
   ( (
    ) )
  ........
  |      |]
  \      /
   `----´ '''

echo ""

# Returns success only if all files required by an analysis run are present
database_complete() {
    [ -f "${REF_FILE}" ] && [ -f "${NODES_FILE}" ] && [ -f "${NAMES_FILE}" ]
}

# Download a gzip file to a temp name and promote it only once it is intact.
# A killed or damaged download therefore never leaves a file that the
# "already downloaded" guards would accept on the next run.
fetch_gz() {
    local url=$1 dest=$2
    wget -O "${dest}.tmp" "${url}"
    echo "Checking ${dest} for damage (gzip -t), this can take several minutes for large files... "
    if ! gzip -t "${dest}.tmp"; then
        rm -f "${dest}.tmp"
        echo "Error: ${dest} was damaged in transit. Re-run to download it again."
        exit 1
    fi
    mv "${dest}.tmp" "${dest}"
}

# Function to download mitochondrial genomes
# Uses HTTPS, not the recursive FTP fetch used before: that wrote every file into
# one -O target and was seen to deliver a corrupt archive with extra bytes.
download_genomes() {
    if [ -f "mitochondrion_refseq.fa.gz" ]; then
        echo "Mitochondrial genomes already downloaded. Skipping download step."
        return
    fi

    echo "Downloading mitochondrial genomes from NCBI RefSeq... "
    # Release number of the files below, recorded in db_info.txt
    wget -qO refseq_release.txt.tmp "${RELEASE_URL}" && mv refseq_release.txt.tmp refseq_release.txt ||
        { rm -f refseq_release.txt.tmp; echo "Warning: could not read the RefSeq release number."; }
    local base="https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion"
    # The release may be split into several numbered parts; fetch all of them
    local parts
    parts=$(wget -qO- "${base}/" | grep -o 'mitochondrion\.[0-9.]*genomic\.fna\.gz' | sort -u) || true
    if [ -z "${parts}" ]; then
        echo "Error: no genomic FASTA files listed at ${base}/"
        exit 1
    fi

    for part in ${parts}; do
        fetch_gz "${base}/${part}" "${part}"
    done
    # Concatenated gzip members are still a valid gzip stream
    cat ${parts} >mitochondrion_refseq.fa.gz.tmp
    mv mitochondrion_refseq.fa.gz.tmp mitochondrion_refseq.fa.gz
    rm -f ${parts}
    # Download date of the mitogenomes, recorded in db_info.txt (the build date can be later)
    date +%F >refseq_downloaded.txt

    if [ ! -f "mitochondrion_refseq.fa.gz" ]; then
        echo "Error: Mitochondrial genome download failed."
        exit 1
    fi

    echo "Mitochondrial genomes downloaded successfully."
}

# Function to download accession to taxid mapping file
download_accession_to_taxid() {
    if [ -f "nucl_gb.accession2taxid.gz" ]; then
        echo "Accession to taxid mapping file already downloaded. Skipping download step."
        return
    fi

    echo "Downloading accession to taxid mapping file... "
    fetch_gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz \
        nucl_gb.accession2taxid.gz

    if [ ! -f "nucl_gb.accession2taxid.gz" ]; then
        echo "Error: Accession to taxid mapping file download failed."
        exit 1
    fi

    echo "Accession to taxid mapping file downloaded successfully."
}

# Function to extract taxids
# Only NC_ accessions are kept. NC_ marks a curated, complete genomic molecule —
# for this dataset, a finished and reviewed mitogenome. NW_ marks an unplaced
# WGS scaffold that automated screening flagged as mitochondrial: uncurated and
# frequently partial. Since renameFASTA_taxid.py keeps only sequences present in
# this mapping, restricting it here is what drops NW_ records from the database.
extract_taxids() {
    if [ -f "nucl_refseq.accession2taxid.tsv" ]; then
        echo "Taxid extraction already completed. Skipping step."
        return
    fi

    echo "Extracting taxids for mitochondrial genomes (NC_ accessions only)... "
    gunzip -c nucl_gb.accession2taxid.gz | awk '$1~/^NC_/' | cut -f 2,3 >nucl_refseq.accession2taxid.tsv.tmp
    mv nucl_refseq.accession2taxid.tsv.tmp nucl_refseq.accession2taxid.tsv

    if [ ! -f "nucl_refseq.accession2taxid.tsv" ]; then
        echo "Error: Taxid extraction failed."
        exit 1
    fi

    echo "Taxid extraction successful."
}

# Function to rename FASTA headers
rename_fasta_headers() {
    if [ -f "mitochondrion_refseq_taxid.fna" ]; then
        echo "FASTA headers already renamed. Skipping step."
        return
    fi

    echo "Renaming FASTA headers to include taxid... "
    python3 "${SCRIPT_DIR}/renameFASTA_taxid.py" \
        --Taxid "${WD}/nucl_refseq.accession2taxid.tsv" \
        --input "${WD}/mitochondrion_refseq.fa.gz" \
        --output "${WD}/mitochondrion_refseq_taxid.fna.tmp"
    mv mitochondrion_refseq_taxid.fna.tmp mitochondrion_refseq_taxid.fna

    if [ ! -f "mitochondrion_refseq_taxid.fna" ]; then
        echo "Error: FASTA header renaming failed."
        exit 1
    fi

    echo "FASTA header renaming successful."
}

# Function to mask low-complexity regions
# Guard checks both the raw and the compressed output, since compress_fasta()
# replaces the uncompressed file and this step is expensive to repeat.
mask_low_complexity_regions() {
    if [ -f "mitochondrion_refseq_taxid_masked.fna" ] || [ -f "${REF_FILE}" ]; then
        echo "Low-complexity regions already masked. Skipping step."
        return
    fi

    echo "Masking low-complexity regions in the mitochondrial genomes... "
    bbmask.sh \
        in=mitochondrion_refseq_taxid.fna \
        out=mitochondrion_refseq_taxid_masked.tmp.fna \
        overwrite=t
    mv mitochondrion_refseq_taxid_masked.tmp.fna mitochondrion_refseq_taxid_masked.fna

    if [ ! -f "mitochondrion_refseq_taxid_masked.fna" ]; then
        echo "Error: Low-complexity region masking failed."
        exit 1
    fi

    echo "Low-complexity region masking successful."
}

# Function to compress masked FASTA file
# Only the masked FASTA is kept; the unmasked one is a build intermediate and is
# deleted by cleanup_intermediates() rather than compressed.
compress_fasta() {
    if [ -f "${REF_FILE}" ]; then
        echo "Masked FASTA file already compressed. Skipping step."
        return
    fi

    echo "Zipping the masked FASTA file... "
    pigz -p "${threads}" -c mitochondrion_refseq_taxid_masked.fna >"${REF_FILE}.tmp"
    mv "${REF_FILE}.tmp" "${REF_FILE}"
    rm -f mitochondrion_refseq_taxid_masked.fna

    if [ ! -f "${REF_FILE}" ]; then
        echo "Error: FASTA file compression failed."
        exit 1
    fi

    echo "FASTA file compression successful."
}

# Function to download the NCBI taxonomy
# Keeps the full nodes.dmp, names.dmp and merged.dmp gzipped in NCBI_taxdump/full/, so
# trim_taxdump() and --extend-db work offline against the release the database was
# built with. Extracted into a temp dir and moved into place, so an interrupted
# download never leaves partial files that the guard would accept.
download_ncbi_taxonomy() {
    if [ -f "${FULL_NODES}" ] && [ -f "${FULL_NAMES}" ] && [ -f "${FULL_MERGED}" ]; then
        echo "NCBI taxonomy files already downloaded. Skipping step."
        return
    fi

    echo "Downloading NCBI taxonomy files... "
    mkdir -p NCBI_taxdump
    fetch_gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz NCBI_taxdump/taxdump.tar.gz

    rm -rf NCBI_taxdump.tmp
    mkdir NCBI_taxdump.tmp
    tar -xzf NCBI_taxdump/taxdump.tar.gz -C NCBI_taxdump.tmp nodes.dmp names.dmp merged.dmp
    pigz -p "${threads}" NCBI_taxdump.tmp/*.dmp
    date +%F >NCBI_taxdump.tmp/downloaded
    rm -rf "${FULL_DIR}"
    mv NCBI_taxdump.tmp "${FULL_DIR}"
    rm -f NCBI_taxdump/taxdump.tar.gz

    echo "NCBI taxonomy files downloaded successfully."
}

# Function to remove every build intermediate that no analysis run reads
# Runs only once the database is verified complete, so an interrupted build keeps
# its intermediates and can be resumed without re-downloading several GB.
cleanup_intermediates() {
    if ! database_complete; then
        echo "Skipping cleanup: database is incomplete, keeping intermediates for resume."
        return
    fi

    echo "Removing build intermediates not used by the analysis pipeline... "

    rm -f "${INTERMEDIATE_FILES[@]}"

    # The NCBI taxdump ships ~9 files but LinkTaxonomy.py only reads nodes.dmp
    # and names.dmp. Whitelist those two so unused files are dropped here
    # (the gzipped full taxonomy in NCBI_taxdump/full/ is a folder and is kept).
    if [ -d "NCBI_taxdump" ]; then
        find NCBI_taxdump -maxdepth 1 -type f \
            ! -name "nodes.dmp" ! -name "names.dmp" -delete
    fi

    echo "Cleanup successful."
}

# Function to trim nodes.dmp and names.dmp to the lineages of taxids in the reference
# FASTA, the extra references and any further FASTA given as arguments (headers are
# taxids, or <taxid>|<accession>). Keeps every ancestor up to root 1, so
# LinkTaxonomy.py lookups still resolve. Reads the full taxonomy in NCBI_taxdump/full/,
# so re-trimming is idempotent; files are replaced via .tmp + mv.
trim_taxdump() {
    echo "Trimming taxonomy files to taxids present in the database... "
    local fastas=("${REF_FILE}") missing
    [ -f "${EXTRA_FILE}" ] && fastas+=("${EXTRA_FILE}")
    missing=$(python3 "${SCRIPT_DIR}/trim_taxdump.py" \
        --nodes "${FULL_NODES}" --names "${FULL_NAMES}" \
        --out-nodes "${NODES_FILE}" --out-names "${NAMES_FILE}" \
        "${fastas[@]}" "$@") || { echo "Error: taxonomy trimming failed."; exit 1; }
    if [ -n "${missing}" ]; then
        echo "Warning: $(wc -l <<<"${missing}") database taxids are not in the taxonomy; their reads get no lineage: $(tr '\n' ' ' <<<"${missing}")"
    fi
    echo "Taxonomy trimming successful."
}

# Value of a key in db_info.txt (empty when missing)
info_get() {
    [ -f "${INFO_FILE}" ] && awk -F'\t' -v k="$1" '$1 == k { print $2 }' "${INFO_FILE}" || true
}

# Function to write db_info.txt: the version of every part of the database.
# Arguments: build date, date of the last extension.
write_db_info() {
    local n_extra=0
    [ -f "${EXTRA_FILE}" ] && n_extra=$(gzip -dc "${EXTRA_FILE}" | grep -c '^>' || true)
    {
        printf 'ecmsd_version\t%s\n' "${ECMSD_VERSION:-unknown}"
        printf 'built\t%s\n' "${1:-unknown}"
        printf 'refseq_release\t%s\n' "$(cat refseq_release.txt 2>/dev/null || info_get refseq_release | grep . || echo unknown)"
        printf 'refseq_downloaded\t%s\n' "$(cat refseq_downloaded.txt 2>/dev/null || info_get refseq_downloaded | grep . || echo unknown)"
        printf 'taxdump_downloaded\t%s\n' "$(cat "${FULL_DIR}/downloaded" 2>/dev/null || echo unknown)"
        printf 'n_references\t%s\n' "$(gzip -dc "${REF_FILE}" | grep -c '^>')"
        printf 'n_extra_references\t%s\n' "${n_extra}"
        printf 'extended\t%s\n' "${2:-never}"
    } >"${INFO_FILE}.tmp"
    mv "${INFO_FILE}.tmp" "${INFO_FILE}"
}

# Function to discard an existing database so it is rebuilt from scratch
# Only files this script is known to create are removed by name, never the
# database folder itself, since --db-folder may point at a shared directory.
reset_database() {
    echo "*** --force-rebuild: discarding the existing database in ${WD} ***"

    for f in "${REF_FILE}" "${INTERMEDIATE_FILES[@]}"; do
        if [ -f "${f}" ]; then
            echo "  removing ${f}"
            rm -f "${f}"
        fi
    done

    if [ -d "NCBI_taxdump" ]; then
        echo "  removing NCBI_taxdump/"
        rm -rf "NCBI_taxdump"
    fi

    # Extra references and their log are kept; trim_taxdump() includes them again.
    echo "Existing database discarded. Rebuilding from scratch... "
    echo ""
}

# Function to download and prepare RefSeq and the taxonomy in the current folder
prepare_refseq() {
    download_genomes
    download_accession_to_taxid
    extract_taxids
    rename_fasta_headers
    mask_low_complexity_regions
    compress_fasta
    download_ncbi_taxonomy
}

# Function to add extra references (--extend-db). Arguments: input kind (fasta or
# accessions) and input file. Each file is replaced via .tmp + mv, and the taxonomy
# before the references, so a running analysis sees either the old or the new set.
extend_database() {
    local kind=$1 input=$2
    if ! database_complete; then
        echo "Error: no complete database in ${WD}. Run ECMSD --create-db first."
        exit 1
    fi
    # Databases built before the full taxonomy was kept get it once here
    download_ncbi_taxonomy

    echo "Checking the new references... "
    python3 "${SCRIPT_DIR}/extend_refs.py" "--${kind}" "${input}" \
        --nodes "${FULL_NODES}" --merged "${FULL_MERGED}" --refs "${REF_FILE}" \
        $([ -f "${EXTRA_FILE}" ] && echo --existing "${EXTRA_FILE}") \
        --output extra_new.fna --log extra_new.tsv
    if [ ! -s extra_new.fna ]; then
        rm -f extra_new.fna extra_new.tsv
        echo "No new references to add; the database is unchanged."
        return
    fi

    echo "Masking low-complexity regions in the new references... "
    bbmask.sh in=extra_new.fna out=extra_new_masked.fna overwrite=t
    # Concatenated gzip members are still a valid gzip stream
    { [ -f "${EXTRA_FILE}" ] && cat "${EXTRA_FILE}"; pigz -p "${threads}" -c extra_new_masked.fna; } >"${EXTRA_FILE}.tmp"

    trim_taxdump "${EXTRA_FILE}.tmp"
    mv "${EXTRA_FILE}.tmp" "${EXTRA_FILE}"

    local today
    today=$(date +%F)
    { [ -f "${EXT_LOG}" ] && cat "${EXT_LOG}" || printf 'date\tsource\tid\ttaxid\tlength\tnote\n'
      awk -v d="${today}" -v s="$(basename "${input}")" '{ print d "\t" s "\t" $0 }' extra_new.tsv; } >"${EXT_LOG}.tmp"
    mv "${EXT_LOG}.tmp" "${EXT_LOG}"
    write_db_info "$(info_get built)" "${today}"
    echo "Added $(wc -l <extra_new.tsv) references ($(info_get n_extra_references) extra references in total); see ${EXT_LOG}."
    rm -f extra_new.fna extra_new_masked.fna extra_new.tsv
}

# Function to update the database (--update-db). Builds the new release in a staging
# folder while the old database stays in use, then moves the finished files in. The
# new taxonomy also covers the old references, so a job running during the swap that
# mapped against the old FASTA still finds its lineages.
update_database() {
    if ! database_complete; then
        echo "Error: no complete database in ${WD}. Run ECMSD --create-db first."
        exit 1
    fi
    local current latest built extended
    current=$(info_get refseq_release)
    latest=$(wget -qO- "${RELEASE_URL}") || latest=""
    if [ -n "${latest}" ] && [ "${latest}" = "${current}" ]; then
        echo "Database is current (RefSeq release ${current})."
        return
    fi
    echo "Updating the database: RefSeq release ${current:-unknown} -> ${latest:-unknown}"
    extended=$(info_get extended)

    mkdir -p "${STAGE_DIR}"
    cd "${STAGE_DIR}"
    prepare_refseq
    if [ -f "../${EXTRA_FILE}" ]; then
        echo "Checking the extra references against the new taxonomy... "
        gzip -dc "../${EXTRA_FILE}" >extra_old.fna
        python3 "${SCRIPT_DIR}/extend_refs.py" --fasta extra_old.fna \
            --nodes "${FULL_NODES}" --merged "${FULL_MERGED}" \
            --output extra_checked.fna --log extra_checked.tsv
        # Already masked when they were added
        pigz -p "${threads}" -c extra_checked.fna >"${EXTRA_FILE}"
        rm -f extra_old.fna extra_checked.fna extra_checked.tsv
    fi
    mkdir -p NCBI_taxdump
    trim_taxdump "../${REF_FILE}"
    if ! database_complete; then
        echo "Error: the update finished but required files are missing in ${WD}/${STAGE_DIR}."
        exit 1
    fi
    built=$(date +%F)
    write_db_info "${built}" "${extended:-never}"
    cleanup_intermediates

    echo "Moving the new database files into place... "
    cd "${WD}"
    rm -rf "${FULL_DIR}"
    mv "${STAGE_DIR}/${FULL_DIR}" "${FULL_DIR}"
    mv -f "${STAGE_DIR}/${NODES_FILE}" "${NODES_FILE}"
    mv -f "${STAGE_DIR}/${NAMES_FILE}" "${NAMES_FILE}"
    if [ -f "${STAGE_DIR}/${EXTRA_FILE}" ]; then
        mv -f "${STAGE_DIR}/${EXTRA_FILE}" "${EXTRA_FILE}"
    fi
    mv -f "${STAGE_DIR}/${REF_FILE}" "${REF_FILE}"
    mv -f "${STAGE_DIR}/${INFO_FILE}" "${INFO_FILE}"
    rm -rf "${STAGE_DIR}"
    echo "Database updated to RefSeq release $(info_get refseq_release)."
}

cd "${WD}"

# One database job at a time per folder. Waits, so parallel workflow jobs that all call
# --create-db at their start simply find the finished database.
# shortcut: no lock where flock is missing (macOS); add a mkdir-based lock if needed there
if command -v flock >/dev/null; then
    exec 9>"${WD}/.ecmsd_db.lock"
    if ! flock -n 9; then
        echo "Another ECMSD database job is running in ${WD}; waiting for it to finish... "
        flock 9
    fi
fi

###############################################################################
#                              Main execution flow                            #
###############################################################################

case "${mode}" in
extend)
    extend_database "${5:?input kind missing}" "${6:?input file missing}"
    exit 0
    ;;
update)
    update_database
    exit 0
    ;;
rebuild)
    # Discard whatever is there and fall through to a full rebuild
    reset_database
    ;;
build)
    # A complete database needs no rebuilding. Cleanup still runs so folders built
    # by earlier versions of this script shed their leftover intermediates.
    if database_complete; then
        echo "Database already complete in ${WD}. Nothing to build."
        cleanup_intermediates
        echo "To rebuild it from scratch, re-run with --create-db --force-rebuild."
        exit 0
    fi
    ;;
*)
    echo "Error: unknown database mode '${mode}'."
    exit 1
    ;;
esac

prepare_refseq
mkdir -p NCBI_taxdump
trim_taxdump

if ! database_complete; then
    echo "Error: database build finished but required files are missing in ${WD}."
    exit 1
fi

write_db_info "$(date +%F)" "$(info_get extended)"
cleanup_intermediates

echo ""
echo "Database ready in ${WD}:"
du -sh "${REF_FILE}" "${NODES_FILE}" "${NAMES_FILE}"
du -sh "${WD}"
