#!/bin/bash
set -euo pipefail

#set working directory
# This script builds a minimap2 reference database for mitochondrial genomes using RefSeq data.
WD=$1
SCRIPT_DIR=$2
threads=$3

# activate the conda environment
#eval "$(conda shell.bash hook)"
#conda activate ${WD}/scripts/conda_env

# The only three files the analysis pipeline reads from the database folder
# (see ECMSD.sh: REF, NODES, NAMES). Everything else produced below is a build
# intermediate and is removed by cleanup_intermediates() once the build succeeds.
REF_FILE="mitochondrion_refseq_taxid_masked.fna.gz"
NODES_FILE="NCBI_taxdump/nodes.dmp"
NAMES_FILE="NCBI_taxdump/names.dmp"

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

# Function to download mitochondrial genomes
download_genomes() {
    if [ -f "mitochondrion_refseq.fa.gz" ]; then
        echo "Mitochondrial genomes already downloaded. Skipping download step."
        return
    fi

    echo "Downloading mitochondrial genomes from NCBI RefSeq... "
    wget -r -np -nd -A "*.genomic.fna.gz" \
        ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/ \
        -O mitochondrion_refseq.fa.gz

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
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

    if [ ! -f "nucl_gb.accession2taxid.gz" ]; then
        echo "Error: Accession to taxid mapping file download failed."
        exit 1
    fi

    echo "Accession to taxid mapping file downloaded successfully."
}

# Function to extract taxids
extract_taxids() {
    if [ -f "nucl_refseq.accession2taxid.tsv" ]; then
        echo "Taxid extraction already completed. Skipping step."
        return
    fi

    echo "Extracting taxids for mitochondrial genomes... "
    gunzip -c nucl_gb.accession2taxid.gz | awk '$1~/^N[CW]/' | cut -f 2,3 >nucl_refseq.accession2taxid.tsv

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
        --output "${WD}/mitochondrion_refseq_taxid.fna"

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
        out=mitochondrion_refseq_taxid_masked.fna

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
    pigz -p "${threads}" mitochondrion_refseq_taxid_masked.fna

    if [ ! -f "${REF_FILE}" ]; then
        echo "Error: FASTA file compression failed."
        exit 1
    fi

    echo "FASTA file compression successful."
}

# Function to download and extract NCBI taxonomy files
download_ncbi_taxonomy() {
    if [ -f "${NODES_FILE}" ] && [ -f "${NAMES_FILE}" ]; then
        echo "NCBI taxonomy files already downloaded and extracted. Skipping step."
        return
    fi

    echo "Downloading NCBI taxonomy files... "
    mkdir -p NCBI_taxdump
    wget -P NCBI_taxdump ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

    if [ ! -f "NCBI_taxdump/taxdump.tar.gz" ]; then
        echo "Error: NCBI taxonomy file download failed."
        exit 1
    fi

    (cd NCBI_taxdump && tar -xzf taxdump.tar.gz && rm taxdump.tar.gz)

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

    # Downloads and intermediates of the reference FASTA. The compressed
    # unmasked FASTA is listed for databases built by earlier script versions.
    rm -f \
        mitochondrion_refseq.fa.gz \
        nucl_gb.accession2taxid.gz \
        nucl_refseq.accession2taxid.tsv \
        mitochondrion_refseq_taxid.fna \
        mitochondrion_refseq_taxid.fna.gz \
        mitochondrion_refseq_taxid_masked.fna

    # The NCBI taxdump ships ~9 files but LinkTaxonomy.py only reads nodes.dmp
    # and names.dmp. Whitelist those two so unused files are dropped here.
    if [ -d "NCBI_taxdump" ]; then
        find NCBI_taxdump -maxdepth 1 -type f \
            ! -name "nodes.dmp" ! -name "names.dmp" -delete
    fi

    echo "Cleanup successful."
}

cd "${WD}"

###############################################################################
#                              Main execution flow                            #
###############################################################################

# A complete database needs no rebuilding. Cleanup still runs so folders built
# by earlier versions of this script shed their leftover intermediates.
if database_complete; then
    echo "Database already complete in ${WD}. Nothing to build."
    cleanup_intermediates
    echo "To rebuild from scratch, delete ${WD} and re-run with --create-db."
    exit 0
fi

download_genomes
download_accession_to_taxid
extract_taxids
rename_fasta_headers
mask_low_complexity_regions
compress_fasta
download_ncbi_taxonomy

if ! database_complete; then
    echo "Error: database build finished but required files are missing in ${WD}."
    exit 1
fi

cleanup_intermediates

echo ""
echo "Database ready in ${WD}:"
du -sh "${REF_FILE}" "${NODES_FILE}" "${NAMES_FILE}"
du -sh "${WD}"
