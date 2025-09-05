# This script trims paired-end reads using fastp.
WD=$1

# Create the results directory if it doesn't exist
mkdir -p ${WD}/results/trimmed_reads
cd ${WD}/results/trimmed_reads

conda activate ${WD}/scripts/conda_env

while IFS=$"," read -r Library Name Age City Country Wolb Type SRA; do
    ## skip the header line
    if [[ "${Library}" != "Library" ]]; then
        echo "Processing library ${Library}"

        ## trim the reads using fastp
        fastp \
            -i ${WD}/data/${Type}/${Library}_1.fastq.gz \
            -I ${WD}/data/${Type}/${Library}_2.fastq.gz \
            -o ${WD}/results/trimmed_reads/${Name}_1_trimmed.fastq.gz \
            -O ${WD}/results/trimmed_reads/${Name}_2_trimmed.fastq.gz \
            --dedup \
            --trim_poly_g \
            --html ${WD}/results/trimmed_reads/${Name}.html \
            --json ${WD}/results/trimmed_reads/${Name}.json \
            --thread 20 \
            --detect_adapter_for_pe
    fi

done <${WD}/data/datasets.csv
