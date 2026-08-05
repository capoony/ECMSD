# ECMSD Changelog

---

## v1.3.0 — Taxonomic hierarchy fidelity & curated reference database (2026-08-06)

### Added

- **`ECMSD.sh`**: `--force-rebuild` / `-b` — discards an existing database and rebuilds it from scratch. Implies `--create-db`, so `ECMSD --force-rebuild -d /path/to/db` is sufficient. Only files ECMSD is known to create are removed, by name; the database folder itself is never deleted, since `--db-folder` may point at a shared directory. The flag is never propagated from an analysis run, so a job that repairs a missing database cannot discard one that parallel jobs are reading
- **`scripts/LinkTaxonomy.py`**: subspecies is now a reported rank — a `subspecies` column in the per-read summary and a `subspecies_name` column in `ref_summary`. `--TaxonomicHierarchy subspecies` previously fell through to species silently; it now selects the rank it names. References with no subspecies rank fall back to the species name

### Changed

- **`shell/MakeRef.sh`**: the reference database is now restricted to curated `NC_` accessions; `NW_` records are excluded (`extract_taxids()` filters on `^NC_` instead of `^N[CW]`). `NC_` marks a complete, reviewed genomic molecule — for this dataset, a finished mitogenome. `NW_` marks an unplaced WGS scaffold that automated screening flagged as mitochondrial: uncurated, often partial, and a common carrier of NUMTs. Such references made `--cov-threshold` inconsistent across the database and, because only one reference is kept per taxid, an `NW_` scaffold could silently displace a taxon's curated `NC_` mitogenome
- **`README.md`**: new "What goes into the reference database" section documenting the `NC_`-only policy, the RefSeq accession prefixes, and why `NW_` scaffolds are excluded
- **`scripts/renameFASTA_taxid.py`**: FASTA headers are now written as the bare taxid (`>{taxid}`) instead of the legacy `>kraken:taxid|{taxid}` format inherited from the retired Kraken2 workflow (see v1.1.0 removal of `shell/kraken2.sh`)
- **`scripts/LinkTaxonomy.py`**: reads the reference name directly as the taxid instead of extracting it via `ref_name.split("|")[1]`
- **`scripts/plot_paf_alignments.R`**: sequencing-depth panel now uses a log1p y-axis scale with power-of-ten breaks, so low-depth regions stay visible next to sharp high-depth peaks instead of flattening out on a linear axis
- **`scripts/plot_paf_alignments.R`**: alignment-coverage panel subtitle now reports breadth-of-coverage percentage instead of the reference name, for clarity
- **`scripts/plot_paf_alignments.R`**: alignment plots are now labelled with the rank requested via `-x`. `ECMSD.sh` passes the rank through as a fifth argument (default `species` when the script is called by hand) and the title takes the `taxon_name` column, which `LinkTaxonomy.py` already writes for the requested rank with the walk-up to the next populated coarser rank. The `subspecies_name` column is only consulted when subspecies was actually requested, as a fallback for `ref_summary` files written before `taxon_name` became rank-aware
- **`scripts/plot_paf_alignments.R`**: plot file names follow `-x` too — `%04d_taxid<taxid>_<label>.png`, using the same label as the title. The species epithet is no longer appended unconditionally, which duplicated the label at species rank (`..._Homo_sapiens_sapiens.png`) and named a finer rank than requested at genus and above. The rank index and taxid still keep names unique where a coarse rank collapses several references
- **`ECMSD.sh`**: default `--cov-threshold` lowered from 50% to 25%
- **`shell/MakeRef.sh`**: build intermediates are now deleted once the database is verified complete. Only the three files an analysis run actually reads are retained — `mitochondrion_refseq_taxid_masked.fna.gz`, `NCBI_taxdump/nodes.dmp` and `NCBI_taxdump/names.dmp`. Removed: `nucl_gb.accession2taxid.gz`, `nucl_refseq.accession2taxid.tsv`, `mitochondrion_refseq.fa.gz`, the unmasked `mitochondrion_refseq_taxid.fna`, and the seven unused NCBI taxdump files (`citations.dmp`, `delnodes.dmp`, `division.dmp`, `gencode.dmp`, `merged.dmp`, `gc.prt`, `readme.txt`). This shrinks a finished database folder from ~5–10 GB to well under 1 GB
- **`shell/MakeRef.sh`**: the unmasked `mitochondrion_refseq_taxid.fna` is no longer gzipped (it was only ever an input to `bbmask`); it is deleted instead
- **`shell/MakeRef.sh`**: a complete database now short-circuits the build instead of re-running steps. Running `--create-db` against an existing database folder reports it as complete and prunes any leftover intermediates, which also migrates folders built by earlier versions

### Fixed

- **`shell/MakeRef.sh`**: `mask_low_complexity_regions()` and `compress_fasta()` tested for the *uncompressed* `.fna` files, which `pigz` had already consumed. Re-running `--create-db` on a complete database therefore repeated the renaming and the slow `bbmask` step, then skipped compression, leaving stale uncompressed FASTA files behind. The guards now test for the compressed artefacts
- **`shell/MakeRef.sh`**: cleanup only runs after the three required files are verified present, so an interrupted build keeps its intermediates and can be resumed without re-downloading several GB
- **`scripts/process_files.R`**: the requested taxonomic rank is now resolved before summarising. `Mito_summary.txt` stores one raw column per rank with a literal `"NA"` where a reference is not resolved that far, and the script grouped straight on that column — so with `-x subspecies` the bulk of the reads (most references stop at species) collapsed into a single `NA` category in both plots and in the written summary tables. `""`/`"NA"` are now normalised to real `NA` and each read walks up to the next populated coarser rank, falling back to `"Unknown"` only when every rank is empty. The rank column's existence is also checked up front (a bad `-x` previously failed cryptically inside `select()`), and mixed-case values are tolerated, since `ECMSD.sh` passes `-x` through verbatim while `LinkTaxonomy.py` lowercases it
- **`scripts/LinkTaxonomy.py`**: the same walk-up replaces the hardcoded subspecies→species special case, so the `ref_summary` label and the `process_files.R` plots agree and coarser ranks get the fallback too
- **`scripts/plot_paf_alignments.R`**: `ref_name` is now read with an explicit `character` column class. With headers written as bare taxids, `read.table()` type-inferred the column as numeric and the reference names no longer matched the PAF

> NOTE: Existing database folders must be rebuilt — they use the old header format and still contain `NW_` records. Because a complete database now short-circuits the build, a plain `--create-db` will not update them; run `ECMSD --create-db --force-rebuild --db-folder /path/to/db` instead

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
