#set working directory
# This script builds a Kraken2 database for mitochondrial genomes using RefSeq data.
WD=$1

# activate the conda environment
conda activate ${WD}/scripts/conda_env

## use Kraken2 to build a database for mitochondrial genomes
## and download all available RefSeq mitochondrial genomes
mkdir ${WD}/data/refseq_mito
cd ${WD}/data/refseq_mito

# download mitochondrial genomes from NCBI RefSeq
wget -r -np -nd -A "*.genomic.fna.gz" \
    ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/

# download the accession to taxid mapping file for RefSeq nucleotide sequences
# and filter for mitochondrial genomes (accessions starting with NC or NW)
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip -c nucl_gb.accession2taxid.gz | awk '$1~/^N[C,W]/' | cut -f 2,3 >nucl_refseq.accession2taxid.tsv

## rename the FASTA headers to include the taxid
# This script will read the accession2taxid mapping file and rename the FASTA headers
python3 ${WD}/scripts/renameFASTA_taxid.py \
    --Taxid ${WD}/data/refseq_mito/nucl_refseq.accession2taxid.tsv \
    --input ${WD}/data/refseq_mito/mitochondrion.1.1.genomic.fna.gz \
    --output ${WD}/data/refseq_mito/mitochondrion.1.1.genomic_taxid.fna

## mask low-complexity regions in the mitochondrial genomes
bbmask.sh \
    in=${WD}/data/refseq_mito/mitochondrion.1.1.genomic_taxid.fna \
    out=${WD}/data/refseq_mito/mitochondrion.1.1.genomic_taxid_masked.fna

## build kraken2 database
mkdir -p ${WD}/data/mt_kraken_db/library
mkdir -p ${WD}/data/mt_kraken_db/taxonomy

# Download the latest taxonomy database and add the mitochondrial genomes to the Kraken2 database
kraken2-build \
    --download-taxonomy \
    --db ${WD}/data/mt_kraken_db \
    --use-ftp

# Add the mitochondrial genomes to the Kraken2 database
kraken2-build \
    --add-to-library ${WD}/data/refseq_mito/mitochondrion.1.1.genomic_taxid_masked.fna \
    --db ${WD}/data/mt_kraken_db \
    --no-masking

# Build the Kraken2 database with the mitochondrial genomes
kraken2-build \
    --build \
    --db ${WD}/data/mt_kraken_db

### now rund Kraken2 on the trimmed reads

## 3) run Kraken2 on the trimmed reads
conda activate ${WD}/scripts/conda_env

mkdir -p ${WD}/results/kraken2

while IFS=$"," read -r Library Name Age City Country Wolb Type SRA; do
    ## skip the header line
    if [[ "${Library}" != "Library" ]]; then
        echo "Processing library ${Library}"

        kraken2 \
            --db ${WD}/data/mt_kraken_db \
            --threads 150 \
            --report ${WD}/results/kraken2/${Name}_1_report.txt \
            --paired \
            ${WD}/results/trimmed_reads/${Name}_1_trimmed.fastq.gz \
            ${WD}/results/trimmed_reads/${Name}_2_trimmed.fastq.gz >/dev/null

        python ${WD}/scripts/ParseKraken2.py \
            --input ${WD}/results/kraken2/${Name}_1_report.txt \
            --output ${WD}/results/kraken2/${Name}_1_summary.tsv \
            --Level S

        python ${WD}/scripts/ParseKraken2.py \
            --input ${WD}/results/kraken2/${Name}_1_report.txt \
            --output ${WD}/results/kraken2/${Name}_1_summary_genus.tsv \
            --Level S

    fi

done <${WD}/data/datasets.csv

for i in 25 50 75 100 125; do
    while IFS=$"," read -r Library Name Age City Country Wolb Type SRA; do
        ## skip the header line
        if [[ "${Library}" != "Library" && ${Type} == "recent" ]]; then
            echo "Processing library ${Name}"

            kraken2 \
                --db ${WD}/data/mt_kraken_db \
                --threads 150 \
                --report ${WD}/results/kraken2/${Name}_1_${i}_report.txt \
                --paired \
                ${WD}/results/trimmed_reads/${Name}_1_trimmed_${i}.fastq.gz \
                ${WD}/results/trimmed_reads/${Name}_2_trimmed_${i}.fastq.gz >/dev/null

            python ${WD}/scripts/ParseKraken2.py \
                --input ${WD}/results/kraken2/${Name}_1_${i}_report.txt \
                --output ${WD}/results/kraken2/${Name}_1_${i}_summary.tsv \
                --Level S

            python ${WD}/scripts/ParseKraken2.py \
                --input ${WD}/results/kraken2/${Name}_1_${i}_report.txt \
                --output ${WD}/results/kraken2/${Name}_1_${i}_summary_genus.tsv \
                --Level G

        fi

    done <${WD}/data/datasets.csv
done

for i in {1..10}; do
    echo $i

    kraken2 \
        --db ${WD}/data/mt_kraken_db \
        --threads 150 \
        --report ${WD}/results/kraken2/Dmel${i}_report.txt \
        --output ${WD}/results/kraken2/Dmel${i}_output.txt \
        --paired \
        /media/inter/SeqData/raw/Macrogen/Illumina/250212/reads/D_mel_${i}_1.fastq.gz \
        /media/inter/SeqData/raw/Macrogen/Illumina/250212/reads/D_mel_${i}_2.fastq.gz

    python ${WD}/scripts/ParseKraken2.py \
        --input ${WD}/results/kraken2/Dmel${i}_report.txt \
        --output ${WD}/results/kraken2/Dmel${i}_summary.tsv \
        --Level S

    python ${WD}/scripts/ParseKraken2.py \
        --input ${WD}/results/kraken2/Dmel${i}_report.txt \
        --output ${WD}/results/kraken2/Dmel${i}_summary_genus.tsv \
        --Level S
done

## rename filenames that Start with "VCFC_" to "VCF"
for file in ${WD}/results/kraken2/VCF_*.tsv; do
    mv "$file" "${file/VCF_/VCF}"
done
