# ECMSD Changelog

---

## Current release — changes relative to v1.0  (2026-05-19)

### Added

- **`ECMSD.sh`**: `--cov-threshold` / `-u` — per-reference breadth-of-coverage (%) as the reference retention criterion (replaces `--RMUS-threshold`)
- **`ECMSD.sh`**: `--skip-mapping` / `-k` — reuse existing `Mito.paf.gz` and `Mito_coverage.txt` without re-mapping
- **`ECMSD.sh`**: `--top-n` / `-n` — controls how many top references receive per-reference alignment/depth PDFs
- **`ECMSD.sh`**: inline Python coverage calculation block (vectorised difference-array, O(n + Σref_len)); coverage step was absent in the old pipeline
- **`scripts/plot_paf_alignments.R`**: new script generating per-reference alignment-breadth and depth PNGs for the top-N references (was absent in `ECMSD_old.sh`)
- **`scripts/LinkTaxonomy.py`**: `--MapQuality` argument added so MQ filtering is consistent with the awk post-filter
- **`scripts/LinkTaxonomy.py`**: `taxon_trace()` memoised with `functools.lru_cache`

### Changed

- **`ECMSD.sh` — `run_mapping()`**: minimap2 runs with `--secondary=no`; only primary alignments are retained; awk post-filter simplified to a single MQ check (`$12 >= Q {print}`)
- **`ECMSD.sh` — coverage Python block**: MQ threshold applied uniformly to all alignments (no secondary detection logic)
- **`ECMSD.sh`**: output directory handling refactored — `--force` now clears only `mapping/` and `logs/` subdirectories rather than deleting the entire output tree
- **`ECMSD.sh`**: `--skip-mapping` preserves the existing PAF and coverage file while clearing only the logs directory
- **`scripts/LinkTaxonomy.py`**: `--Bins` / `--RMUS` arguments replaced by `--Coverage`, `--CovThreshold`, and `--MapQuality`
- **`scripts/process_files.R`**: `read.table()` replaced with `data.table::fread()`; proportional abundance computed from unique primary-alignment reads only (`distinct(SeqID)`)
- **`scripts/plot_paf_alignments.R`**: `pafr::read_paf()` replaced with a `zcat | awk | data.table::fread` pipeline; depth calculation vectorised with `tabulate()` + `cumsum()`

### Removed

- **`ECMSD.sh`**: `--Binsize` / `-b` parameter (bin-based RMUS analysis removed)
- **`ECMSD.sh`**: `--RMUS-threshold` / `-u` parameter (replaced by `--cov-threshold`)
- **`scripts/LinkTaxonomy.py`**: `--Bins` and `--RMUS` arguments
- **`pafr`** R package dependency (removed entirely)

---

## v1.0 — Initial release (`ECMSD_old.sh`)

- minimap2 short-read alignment (`-x sr`, `--secondary=no`), MQ-filtered with awk
- `LinkTaxonomy.py` taxonomy linking via NCBI nodes/names dump
- `process_files.R` proportional abundance plots
- Parameters: `--fwd`, `--rev`, `--merged`, `--Binsize`, `--RMUS-threshold`, `--mapping_quality`, `--threads`, `--force`, `--taxonomic-hierarchy`
