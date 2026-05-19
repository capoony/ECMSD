#!/bin/bash
# Test dataset for ECMSD pipeline
#
# This dataset is used to test the ECMSD pipeline. It contains a small set of
# reads from hDNA sequencing data of unknown origin. The reads are in FASTQ format
# and are stored in the TestData directory. The dataset is not intended for
# production use.

# Set the working directory to your ECMSD installation
WD="<path_to_your_ECMSD_directory>"

# Set the path to the database folder (must be created first with --create-db)
DB="<path_to_your_db_folder>"

# Step 1: Build the reference database (only needed once)
# bash "${WD}/shell/ECMSD.sh" \
#     --create-db \
#     --db-folder "${DB}"

# Step 2: Run the ECMSD pipeline with test data
bash "${WD}/shell/ECMSD.sh" \
    --fwd "${WD}/TestData/merged.fastq.gz" \
    --out "${WD}/TestOutput" \
    --db-folder "${DB}" \
    --threads 20 \
    --cov-threshold 10 \
    --top-n 25 \
    --mapping_quality 20 \
    --taxonomic-hierarchy genus \
    --force
