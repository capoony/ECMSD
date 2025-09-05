#set working directory
# This script builds a Kraken2 database for mitochondrial genomes using RefSeq data.
WD=$1
threads=$2

# activate the conda environment
eval "$(conda shell.bash hook)"
conda activate ${WD}/scripts/conda_env

## use Kraken2 to build a database for mitochondrial genomes
## and download all available RefSeq mitochondrial genomes
mkdir -p ${WD}/data/refseq_mito/NCBI_taxdump && cd ${WD}/data/refseq_mito

# download mitochondrial genomes from NCBI RefSeq and store with name "mitochondrion_refseq.fa.gz"
wget -r -np -nd -A "*.genomic.fna.gz" \
    ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/ \
    -O mitochondrion_refseq.fa.gz

# download the accession to taxid mapping file for RefSeq nucleotide sequences
# and filter for mitochondrial genomes (accessions starting with NC or NW)
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip -c nucl_gb.accession2taxid.gz | awk '$1~/^N[C,W]/' | cut -f 2,3 >nucl_refseq.accession2taxid.tsv

## rename the FASTA headers to include the taxid
# This script will read the accession2taxid mapping file and rename the FASTA headers
python3 ${WD}/scripts/renameFASTA_taxid.py \
    --Taxid ${WD}/data/refseq_mito/nucl_refseq.accession2taxid.tsv \
    --input ${WD}/data/refseq_mito/mitochondrion_refseq.fa.gz \
    --output ${WD}/data/refseq_mito/mitochondrion_refseq_taxid.fna

## mask low-complexity regions in the mitochondrial genomes
# which java
# java -version

bbmask.sh \
    in=${WD}/data/refseq_mito/mitochondrion_refseq_taxid.fna \
    out=${WD}/data/refseq_mito/mitochondrion_refseq_taxid_masked.fna

## gzip the masked FASTA file
pigz -p ${threads} \
    ${WD}/data/refseq_mito/mitochondrion_refseq_taxid_masked.fna \
    ${WD}/data/refseq_mito/mitochondrion_refseq_taxid.fna

## download nodes.dmp and names.dmp from NCBI taxonomy
cd NCBI_taxdump
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
rm taxdump.tar.gz
