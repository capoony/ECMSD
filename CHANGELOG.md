# ECMSD Changelog

---

## v1.2.0 — Coverage-based filtering & conda integration (2026-05-20)

### Added

- **`ECMSD.sh`**: `--cov-threshold` / `-u` — per-reference breadth-of-coverage (%) as the reference retention criterion (replaces `--RMUS-threshold`)
- **`ECMSD.sh`**: `--skip-mapping` / `-k` — reuse existing `Mito.paf.gz` and `Mito_coverage.txt` without re-mapping
- **`ECMSD.sh`**: `--top-n` / `-n` — controls how many top references receive per-reference alignment/depth PDFs
- **`ECMSD.sh`**: inline Python coverage calculation block (vectorised difference-array, O(n + Σref_len)); coverage step was absent in the old pipeline
- **`scripts/plot_paf_alignments.R`**: new script generating per-reference alignment-breadth and depth PNGs for the top-N references (was absent in `ECMSD_old.sh`)
- **`scripts/LinkTaxonomy.py`**: `--MapQuality` argument added so MQ filtering is consistent with the awk post-filter
- **`scripts/LinkTaxonomy.py`**: `taxon_trace()` memoised with `functools.lru_cache`
- **`conda/ecmsd_env.yaml`**: added `r-gridextra`, `r-data.table` dependencies
- **`conda/meta.yaml`**: added `r-gridextra`, `r-data.table` to run requirements
- **`shell/install.sh`**: added `r-gridextra`, `r-data.table` to conda install block
- **`.github/workflows/auto-tag.yml`**: CI step that automatically rewrites the `version=` line in `ECMSD.sh` on each new tag and commits the change back to `main`

### Changed

- **`ECMSD.sh` — `run_mapping()`**: minimap2 runs with `--secondary=no`; only primary alignments are retained; awk post-filter simplified to a single MQ check (`$12 >= Q {print}`)
- **`ECMSD.sh` — coverage Python block**: MQ threshold applied uniformly to all alignments (no secondary detection logic)
- **`ECMSD.sh`**: output directory handling refactored — `--force` now clears only `mapping/` and `logs/` subdirectories rather than deleting the entire output tree
- **`ECMSD.sh`**: `--skip-mapping` preserves the existing PAF and coverage file while clearing only the logs directory
- **`scripts/LinkTaxonomy.py`**: `--Bins` / `--RMUS` arguments replaced by `--Coverage`, `--CovThreshold`, and `--MapQuality`
- **`scripts/process_files.R`**: `read.table()` replaced with `data.table::fread()`; proportional abundance computed from unique primary-alignment reads only (`distinct(SeqID)`)
- **`scripts/plot_paf_alignments.R`**: `pafr::read_paf()` replaced with a `gzip -dc | awk | data.table::fread` pipeline (replaces `zcat` for macOS compatibility); depth calculation vectorised with `tabulate()` + `cumsum()`; fixed null-reference guard before `nrow()` check
- **`conda/ecmsd_env.yaml`**: removed `numpy` and `openjdk=17`
- **`shell/install.sh`**: removed `numpy` from conda install block

### Removed

- **`ECMSD.sh`**: `--Binsize` / `-b` parameter (bin-based RMUS analysis removed)
- **`ECMSD.sh`**: `--RMUS-threshold` / `-u` parameter (replaced by `--cov-threshold`)
- **`scripts/LinkTaxonomy.py`**: `--Bins` and `--RMUS` arguments
- **`scripts/LinkTaxonomy_old.py`**: legacy script removed from repository
- **`shell/requirements.sh`**: superseded by `install.sh` (removed)
- **`pafr`** R package dependency (removed entirely)

---

## v1.1.0 — Database integration & installation overhaul (2026-03-24)

### Added

- **`ECMSD.sh`**: `--create-db` / `-z` — triggers `MakeRef.sh` to build the reference database in `--db-folder` and exits; no analysis is run
- **`ECMSD.sh`**: `--db-folder` / `-d` — explicit path to the reference database folder; if required files are absent the pipeline auto-invokes `MakeRef.sh` to generate them before proceeding
- **`ECMSD.sh`**: `--prefix` / `-p` — prefix for output file names (PAF, summary, plots)
- **`ECMSD.sh`**: dynamic script-path resolution — detects conda install (`$CONDA_PREFIX/lib/ecmsd/`) vs direct repo clone and sets `SCRIPT_DIR` / `SHELL_DIR` accordingly
- **`shell/install.sh`**: new script that installs all dependencies (minimap2, bbmap, R, pigz, …) into the active conda environment via `conda install` and copies scripts + `ECMSD` entry-point into `$CONDA_PREFIX`
- **`conda/build.sh`**: conda recipe build script — copies shell scripts and R/Python scripts into `$PREFIX/lib/ecmsd/` and installs the `ECMSD` entry-point in `$PREFIX/bin/`
- **`conda/ecmsd_env.yaml`**: conda environment definition (Python 3.11, NumPy, R 4.3, r-tidyverse, minimap2, pigz, openjdk 17)
- **`LICENSE`**: MIT licence added to the repository

### Changed

- **`shell/MakeRef.sh`**: completely refactored from a linear script into modular, resumable functions — each step (genome download, accession-to-taxid download, taxid extraction, FASTA header renaming, low-complexity masking, NCBI taxonomy dump) checks whether its output already exists and skips if so; signature changed to `MakeRef.sh <DB_FOLDER> <SCRIPT_DIR> <THREADS>`
- **`shell/MakeRef.sh`**: framing corrected from "Kraken2 database" to "minimap2 reference database"
- **`ECMSD.sh`**: `--db-folder` is now required for analysis runs; the old hardcoded fallback path (`data/refseq_mito`) is removed
- **`scripts/LinkTaxonomy.py`**: validation added for required options; exits with a help message if any are missing
- **`scripts/renameFASTA_taxid.py`**: removed unused `--logical` and `--param` arguments; added required-option and file-existence checks; fixed dict-assignment bug (`IDlist[TaxidDict[name]]` had no assignment — now `= True`)

### Removed

- **`ECMSD.sh`**: `--skip_environment` / `-s` (conda activation is now handled entirely by `install.sh` / `build.sh`)
- **`shell/kraken2.sh`**: Kraken2-based workflow retired
- **`shell/requirements.sh`**: superseded by `install.sh`
- **`shell/trim.sh`**: trimming step removed from the pipeline

---

## v1.0.0 — Initial release

- minimap2 short-read alignment (`-x sr`, `--secondary=no`), MQ-filtered with awk
- `LinkTaxonomy.py` taxonomy linking via NCBI nodes/names dump
- `process_files.R` proportional abundance plots
- Parameters: `--fwd`, `--rev`, `--merged`, `--Binsize`, `--RMUS-threshold`, `--mapping_quality`, `--threads`, `--force`, `--taxonomic-hierarchy`
