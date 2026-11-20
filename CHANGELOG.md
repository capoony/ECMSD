# ECMSD Changelog

---

## Unreleased

### Added

- **`scripts/assess_references.py`**, **`ECMSD.sh`**: spatial statistics of the discriminatory reads in `Mito_reference_evidence.txt` (`breadth_ratio`, `discriminatory_windows`, `window_ratio`, `max_1kb_share`) and a new credibility requirement `--min-disc-window-ratio` (default 0.5): a reference whose discriminatory reads hit less than half the 1 kb windows that as many random reads would hit is not credible
- **`shell/install.sh`**: new `-p/--prefix PATH` option creates a new local conda environment at `PATH` and installs the dependencies and the `ECMSD` command there, without an active conda environment. An existing environment at `PATH` is reused; a non-empty directory that is not a conda environment is refused. New `-h/--help`
- **`shell/install.sh`**: unknown options exit with an error
- **`ECMSD.sh`**: two-stage mapping. Stage 1 maps all reads with `-x sr -k 17 -w 7 -n 1 -m 20 -s 20 --secondary=no`, no MAPQ filter and no `-c`, and collects the names of all reads with an alignment. Stage 2 remaps only those reads, extracted from the original FASTQ file(s), against the full database with `-c --secondary=yes -N 20`
- **`ECMSD.sh`**: new options `--seed-k` (default 17), `--seed-w` (default 7), `--min-chain-score` (default 20), `--max-alns` (default 20) and `--keep-all-alignments`
- **`scripts/extract_candidates.py`**: streams the original single-end, merged or paired-end FASTQ (plain or gzip) and writes the Stage 1 hits unchanged; a pair is kept when either mate hit, mate order is verified, only read names are held in memory
- **`scripts/select_best_alignments.py`**: keeps all best-scoring Stage 2 alignments of each read (one per reference, ties included) by `AS`, never by MAPQ, and tags them with `nc`, `dm` and `na`; `--one-per-read` keeps a single alignment (higher `AS`, then primary, then lower `NM`); optional read filters `--min-mapq`, `--min-aln-len`, `--min-aln-frac` and `--min-identity`; reports and handles missing `tp`/`AS`/`NM` tags and stops if no alignment has `AS`; reads from stdin and writes to stdout with `-`; with `--table` it writes a table of all alignments with `tp`, `AS`, `NM`, MAPQ and the `AS` difference to the best one
- **`scripts/assess_references.py`**, **`ECMSD.sh`**: competitive, sample-wide assessment of each reference. Reads are classified as compatible (best `AS`, ties included), discriminatory (the only reference with the best `AS`, margin of at least 10, alignment of at least 20 bp, 50% of the read and 90% identity), weak unique, or shared; per reference the table `Mito_reference_evidence.txt` gives counts, breadth, depth and evenness for compatible and discriminatory reads, how many discriminatory reads were at the `-N` limit, and whether the reference is credible
- **`ECMSD.sh`**: new options `--min-disc-reads` (default 20), `--min-disc-breadth` (default 5) and `--min-disc-evenness` (default 0.5), which define a credible reference
- **`ECMSD.sh`**: new option `--min-aln-len` (default 20): the minimum alignment length of a read that counts as discriminatory evidence. It was fixed at 30 bp, which left samples of very short, adapter-trimmed reads (about 28 bp on average) without any discriminatory read
- **`scripts/assess_references.py`**: new columns `discriminatory_identity` (pooled 1 - sum `NM` / sum alignment length), `discriminatory_exact` (fraction of reads without mismatch) and `edit_distance_decline` (the HOPS -Δ% over 0-3 mismatches; below 0.9 the reads do not fit the reference well) in `Mito_reference_evidence.txt`, computed from the discriminatory reads only. Reported only, not used to decide. They can show that the species of a sample is missing from the database, but cannot separate species that share mitochondrial haplotypes
- **`scripts/assess_references.py`**, **`ECMSD.sh`**: PCR duplicates (same references, start, end and strand) are counted once in `Mito_reference_evidence.txt`; new column `duplicate_fraction`; new option `--no-dedup` to count every read
- **`scripts/LinkTaxonomy.py`**: new `--Evidence` option; only credible references are reported, and a read compatible with several credible references is split between them in the proportion of their discriminatory reads, as one row per reference with its share in the new column `Weight` of `Mito_summary.txt` (no assignment to a higher taxonomic level). `scripts/process_files.R` sums the weights; files without `Weight` count every read once
- **tests**: unit tests for the new scripts and options

### Changed

- **`ECMSD.sh`**: the default `--cov-threshold` is 25% again (it was 5% since v1.3.3). With the discriminatory-evidence step, 5% still reported two spurious taxa with 5-8% breadth on the test sample; 25% reports only *D. melanogaster* and *G. bimaculatus*. Because breadth depends on depth, a rare species can fall below 25% in a shallow sample (*G. bimaculatus* at 10% of the reads: 19%); lower it for such samples
- **`ECMSD.sh`**: the default `-q/--mapping_quality` is now 0 (no MAPQ filtering) and drops reads whose compatible alignments all have a lower MAPQ. It was 20 and applied to a single `--secondary=no` pass, which discarded correct reads with ambiguous (MAPQ 0) alignments between related mitogenomes
- **`ECMSD.sh`**: the `sr` chaining thresholds are relaxed to `-n 1 -m 20 -s 20` in both stages (preset: `-n 2 -m 25 -s 40`); with the preset values reads under ~30 bp are almost never mapped. A first choice of `-m 10` passed 85% of the reads of the test sample to Stage 2, 14% of them with AS < 40; `-m 20` passes 59% and none of them
- **`scripts/process_files.R`**: no longer drops reads with `MappingQuality` 0, which are now valid best alignments
- **`ECMSD.sh`**, **`scripts/coverage.py`**: coverage is calculated from all compatible alignments, so reads that tie between references contribute to each of them; the number of reads in `Mito.paf.gz` is therefore smaller than the number of lines
- **`ECMSD.sh`**: `minimap2.log` is replaced by `minimap2_stage1.log` and `minimap2_stage2.log`
- **`ECMSD.sh`**: the Stage 2 output is piped straight into `select_best_alignments.py` and no longer written to disk unless `--keep-all-alignments` is given. On the 1.09 M read test sample this avoids a 337 MB gzip of 16.6 M alignments; Stage 2 plus best-alignment selection takes about 55 s instead of about 5 min
- **`scripts/select_best_alignments.py`**: keeps only the best alignment of the current read in memory and parses the tags with string searches instead of splitting every line (about 2x faster)

---

## v1.3.3 — Audit fixes (2026-10-25)

### Changed

- **`ECMSD.sh`**: default `--cov-threshold` lowered from 25% to 5%
- **`scripts/LinkTaxonomy.py`**: `load_best_hits` keeps only the first passing alignment per read instead of every line, which cuts memory on large samples
- **`ECMSD.sh`**: the coverage filter is now `scripts/coverage.py`, covered by unit tests
- **`shell/MakeRef.sh`**: at the end of `--create-db`, `nodes.dmp` and `names.dmp` are trimmed to the lineages of the taxids in the FASTA, so each run parses a few thousand taxa instead of the whole NCBI taxonomy. A database built earlier is not trimmed until it is rebuilt

### Fixed

- **`ECMSD.sh`**: log text from `run_mapping` went to stdout and ended up inside `Mito.paf.gz`; it now goes to stderr
- **`ECMSD.sh`**: `-c/--force` ran `rm -rf` on the whole output folder and could delete input reads or a project folder; it now removes only `mapping/` and `logs/`
- **`ECMSD.sh`**: an analysis run no longer starts `MakeRef.sh` when database files are missing (parallel jobs raced on the same `.tmp` files); it prints "Run `ECMSD --create-db` first" and exits 1
- **`README.md`**: removed the "Paired-end with additional merged reads" example, which the script rejects
- **`ECMSD.sh`**: option values (`-x`, `-u`, `-q`, `-n`, `-t`) and missing option values are checked before mapping instead of failing after hours of compute
- **`scripts/LinkTaxonomy.py`**: an unknown `--TaxonomicHierarchy` rank now exits with an error instead of falling back to species
- **`scripts/coverage.py`**: reads with the same name (paired mates) on one reference were deduplicated to one, undercounting breadth; the dedupe key now includes start and end
- **`TestData/shell.sh`**: replaced the removed `--binsize` and `--RMUS-threshold` options with `--cov-threshold 10 --top-n 25`
- **`shell/MakeRef.sh`**: the taxdump is extracted into `NCBI_taxdump.tmp/` and moved into place, so a killed extraction no longer leaves a truncated `names.dmp`
- **`.github/workflows/auto-tag.yml`**: commit author name, email and tag are passed through `env:` instead of being interpolated into shell commands
- **`ECMSD.sh`**: the "no references passed" warning said plots were skipped although both R scripts still run; it now says plots will be empty
- **`ECMSD.sh`**: `LinkTaxonomy.py` is called with `python3`

### Removed

- **`ECMSD.sh`**, **`scripts/LinkTaxonomy.py`**, **`scripts/renameFASTA_taxid.py`**: unreachable `rev`/`merged` checks, the duplicate merged-file check, the unused `get_alignment_score` and the empty `OptionGroup`
- **`conda/meta.yaml`**: stale v1.1.0 recipe; the bioconda recipe is the maintained one

### Docs

- **`README.md`**: new "Known limitations" section: one genome is kept per taxid, so multi-record mitochondrial genomes lose all but the first record

---

## v1.3.2 — Database build fixes (2026-09-14)

### Fixed

- **`shell/MakeRef.sh`**: downloads now use HTTPS and are checked with `gzip -t` before they are kept. The old FTP fetch could deliver a corrupt `mitochondrion_refseq.fa.gz`, which crashed `renameFASTA_taxid.py` with a `zlib.error`
- **`shell/MakeRef.sh`**: each step writes to a `.tmp` file and renames it only on success, so a re-run no longer accepts a partial output left by a crashed step

---

## v1.3.1 — Test coverage & no-hits threshold fix (2026-08-17)

### Added

- **`scripts/tests/`**: unit tests for `LinkTaxonomy.py` and `renameFASTA_taxid.py` (stdlib `unittest`, no new dependency). Covers PAF alignment-score parsing, coverage-threshold filtering, mapping-quality filtering, taxonomy-path tracing, the taxonomic rank walk-up fallback, and FASTA header rename/dedup-by-taxid
- **`TestData/test_no_hits_above_threshold.sh`**: regression script for the zero-references-pass-threshold fix below. Runs the pipeline with an unreachable `--cov-threshold` and asserts a clean exit plus header-only summary tables

### Changed

- **`scripts/LinkTaxonomy.py`**, **`scripts/renameFASTA_taxid.py`**: driver logic moved behind `if __name__ == "__main__"`; both scripts previously ran their full pipeline at import time, which made the underlying functions impossible to unit test. CLI behaviour is unchanged. The rank walk-up logic in `LinkTaxonomy.py` is now a standalone `fallback_label()` function and the FASTA dedup loop in `renameFASTA_taxid.py` is now `rename_fasta()`

### Fixed

- **`ECMSD.sh`**, **`scripts/process_files.R`**: a run where no reference clears `--cov-threshold` no longer fails. `ECMSD.sh` exited 1 on an empty summary file and `process_files.R` called `stop()` on zero taxa; both now log a warning/message, skip plot generation, and exit 0, still writing the declared (header-only) summary tables. A Snakemake rule around ECMSD (e.g. pastforward) requires its declared outputs to exist and treats a nonzero exit as rule failure, so a legitimate "no hits" result was indistinguishable from a crash
- **`shell/install.sh`**: the `conda-forge::python>=3.6` dependency spec could resolve to a GraalPy build instead of CPython, which breaks every `python3` call the pipeline makes (GraalPy's launcher fails outright without its own JVM setup). Pinned to `conda-forge::python>=3.6=*cpython*`

### Chore

- **`scripts/tests/`**: further coverage for both scripts' `load_data()` (plain/gzip/stdin dispatch), `LinkTaxonomy.py`'s `main()` (nodes.dmp/names.dmp parsing, including that non-"scientific name" `names.dmp` entries are ignored) and `taxon_trace()`'s error exit on a broken taxonomy path, and `renameFASTA_taxid.py`'s `run()` (missing input/taxid file handling, plus an end-to-end run producing a renamed FASTA). No functional code changed

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
