#!/bin/bash
# Test dataset for ECMSD pipeline
#
# This dataset is used to test the ECMSD pipeline. It contains a small set of
# reads from hDNA sequencing data of unknown origin. The reads are in FASTQ format
# and are stored in the TestData directory. The dataset is not intended for
# production use.

# Install ECMSD via conda (conda install -c bioconda ecmsd) or by running install.sh

# Set the working directory (parent folder containing TestData/)
WD="<path_to_your_working_directory>"

# Set the path to the database folder (must be created first with --create-db)
DB="<path_to_your_db_folder>"

if [ -z "$WD" ] || [ "$WD" = "<path_to_your_working_directory>" ]; then
    echo "Error: Please set the WD variable to your working directory."
    exit 1
fi

if [ -z "$DB" ] || [ "$DB" = "<path_to_your_db_folder>" ]; then
    echo "Error: Please set the DB variable to the path of your database folder."
    exit 1
fi

# Step 1: Build the reference database (only needed once)
ECMSD \
    --create-db \
    --db-folder "${DB}"

# Step 2: Run the ECMSD pipeline with test data
ECMSD \
    --fwd "${WD}/TestData/merged.fastq.gz" \
    --out "${WD}/TestOutput" \
    --db-folder "${DB}" \
    --threads 20 \
    --cov-threshold 10 \
    --top-n 25 \
    --mapping_quality 20 \
    --taxonomic-hierarchy genus \
    --force
