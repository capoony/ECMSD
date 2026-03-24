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
mask_low_complexity_regions() {
    if [ -f "mitochondrion_refseq_taxid_masked.fna" ]; then
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
compress_fasta() {
    if [ -f "mitochondrion_refseq_taxid_masked.fna.gz" ]; then
        echo "Masked FASTA file already compressed. Skipping step."
        return
    fi

    echo "Zipping the masked FASTA file... "
    pigz -p "${threads}" mitochondrion_refseq_taxid_masked.fna
    pigz -p "${threads}" mitochondrion_refseq_taxid.fna

    if [ ! -f "mitochondrion_refseq_taxid_masked.fna.gz" ]; then
        echo "Error: FASTA file compression failed."
        exit 1
    fi

    echo "FASTA file compression successful."
}

# Function to download and extract NCBI taxonomy files
download_ncbi_taxonomy() {
    if [ -f "NCBI_taxdump/nodes.dmp" ] && [ -f "NCBI_taxdump/names.dmp" ]; then
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

cd "${WD}"

# Main execution flow
download_genomes
download_accession_to_taxid
extract_taxids
rename_fasta_headers
mask_low_complexity_regions
compress_fasta
download_ncbi_taxonomy
