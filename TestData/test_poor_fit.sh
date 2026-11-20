#!/bin/bash
# Regression test: reads that fit no reference well are labelled poor_fit.
#
# merged_poor_fit.fastq.gz is merged.fastq.gz with 6% of the bases substituted
# (make_poor_fit.py), mimicking a haplotype missing from the database.
# Asserts that Mito_calls.txt has poor_fit rows and no confident row.
#
# Requires the same setup as shell.sh: an existing DB built with --create-db.

set -euo pipefail

WD="<path_to_your_working_directory>"
DB="<path_to_your_db_folder>"

if [ "$WD" = "<path_to_your_working_directory>" ] || [ "$DB" = "<path_to_your_db_folder>" ]; then
    echo "Error: Please set the WD and DB variables."
    exit 1
fi

OUT="${WD}/TestOutput_poor_fit"

ECMSD \
    --fwd "${WD}/TestData/merged_poor_fit.fastq.gz" \
    --out "${OUT}" \
    --db-folder "${DB}" \
    --threads 20 \
    --taxonomic-hierarchy species \
    --force

CALLS="${OUT}/mapping/Mito_calls.txt"
grep -v '^#' "$CALLS" | cut -f2 | grep -qx poor_fit || { echo "FAIL: no poor_fit call"; exit 1; }
! grep -v '^#' "$CALLS" | cut -f2 | grep -qx confident || { echo "FAIL: confident call in poor-fit sample"; exit 1; }
echo "OK: poor_fit sample labelled correctly"
