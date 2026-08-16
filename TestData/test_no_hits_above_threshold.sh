#!/bin/bash
# Regression test: no reference passes --cov-threshold
#
# Sets --cov-threshold above what any reference in the test data can reach,
# forcing the "zero hits" path (fixed in fix/no-hits-above-threshold — ECMSD.sh
# used to exit 1 / process_files.R used to stop() here, which broke Snakemake
# rules expecting declared outputs to exist). Asserts the pipeline instead
# exits 0 and still writes header-only summary tables.
#
# Requires the same setup as shell.sh: an existing DB built with --create-db.

set -euo pipefail

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

OUT="${WD}/TestOutput_no_hits"

ECMSD \
    --fwd "${WD}/TestData/merged.fastq.gz" \
    --out "${OUT}" \
    --db-folder "${DB}" \
    --threads 20 \
    --cov-threshold 100 \
    --mapping_quality 20 \
    --taxonomic-hierarchy species \
    --force

# --- Assertions -------------------------------------------------------------
summary="${OUT}/mapping/Mito_summary.txt"
taxon_summary="${OUT}/mapping/Mito_summary_species.txt"
proportions="${OUT}/mapping/Mito_summary_species_proportions.txt"

for f in "$summary" "$taxon_summary" "$proportions"; do
    if [[ ! -f "$f" ]]; then
        echo "FAIL: expected output '$f' was not created."
        exit 1
    fi
    if [[ $(wc -l <"$f") -gt 1 ]]; then
        echo "FAIL: '$f' has data rows; --cov-threshold 100 should pass none."
        exit 1
    fi
done

echo "PASS: zero-hit run exited 0 and left header-only summary tables."
