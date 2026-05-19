# ECMSD - Efficient Comprehensive Mitochondrial Sequence Detection Pipeline

## Overview

ECMSD (Efficient Comprehensive Mitochondrial Sequence Detection) is a bioinformatics pipeline for the detection and taxonomic classification of mitochondrial sequences from high-throughput Illumina sequencing data. It aligns reads against a comprehensive NCBI RefSeq mitochondrial reference database using minimap2, filters alignments by coverage and quality, assigns taxonomy via the NCBI taxonomy database, and generates quantitative summaries and visualisation plots.

## Features

- **Fast mitochondrial sequence detection** using minimap2 short-read alignment
- **Comprehensive reference database** built from NCBI RefSeq mitochondrial genomes (repeat-masked)
- **Taxonomic classification** at configurable hierarchy levels (species → kingdom)
- **Quality-based filtering**: separate, correct handling of primary (MQ threshold) and secondary alignments (chaining-score threshold)
- **Coverage-based reference filtering**: only references with sufficient breadth of coverage are retained
- **Secondary alignment support**: configurable number of secondary alignments per read for highly similar reference genomes
- **Per-reference alignment and depth plots** for the top-N references
- **Proportional abundance plots** based on unique primary-alignment reads only
- **Skip-mapping mode** to rerun analysis steps on existing PAF output without re-mapping
- **Conda environment management** for reproducible results

---

## Pipeline Workflow

1. **Reference Database Construction** (`MakeRef.sh`): Downloads NCBI RefSeq mitochondrial genomes, repeat-masks them with BBMask, and embeds taxIDs in the FASTA headers
2. **Sequence Mapping**: Aligns reads to the reference with minimap2 (`-x sr`, secondary alignments enabled); primary alignments are filtered by mapping quality (MQ ≥ threshold); secondary alignments are filtered by the s1 chaining score
3. **Coverage Calculation**: Computes per-reference breadth of coverage and mean/std depth using a vectorised difference-array algorithm; references below the coverage threshold are discarded
4. **Taxonomic Assignment** (`LinkTaxonomy.py`): Links each alignment to the full NCBI taxonomic lineage; primary alignments are again MQ-filtered to ensure consistency
5. **Proportional Abundance Analysis** (`process_files.R`): Summarises read counts and proportions per taxon using **unique primary-alignment reads only** (secondary alignments excluded)
6. **Per-Reference Alignment Plots** (`plot_paf_alignments.R`): Generates coverage breadth and sequencing depth PDFs for the top-N references using **all alignments** (including secondaries) read directly from the PAF file

---

## Requirements

- Linux/Unix operating system
- Conda package manager
- Internet connection (for initial database download)

### Software Dependencies (automatically installed via Conda)

| Tool | Purpose |
|---|---|
| minimap2 | Read alignment |
| BBTools (bbmask) | Reference repeat masking |
| pigz | Parallel gzip |
| Python 3 | Coverage calculation and taxonomy linking |
| R + tidyverse | Proportional abundance plots |
| R + data.table | Fast file I/O for large PAF/summary files |
| R + ggplot2 + gridExtra | Alignment and depth plots |

---

## Installation

```bash
git clone <repository-url>
cd ECMSD
```

Dependencies and the reference database are installed automatically on the first run.

---

## Usage

### Basic Usage

```bash
# Single-end reads
./shell/ECMSD.sh -f reads.fastq -o output_directory/

# Paired-end reads
./shell/ECMSD.sh -f reads_R1.fastq -r reads_R2.fastq -o output_directory/

# Paired-end + separately merged reads
./shell/ECMSD.sh -f reads_R1.fastq -r reads_R2.fastq -m merged.fastq -o output_directory/
```

### Command Line Options

#### Required Arguments

| Flag | Description |
|---|---|
| `-f, --fwd` | Path to the forward (or single-end/merged) FASTQ file |
| `-o, --out` | Path to the output folder |

#### Optional Arguments

| Flag | Default | Description |
|---|---|---|
| `-r, --rev` | — | Path to the reverse FASTQ file (paired-end) |
| `-m, --merged` | — | Path to an additional merged FASTQ file |
| `-u, --cov-threshold` | 50 | Minimum % of a reference covered by reads for it to be retained (0–100) |
| `-s, --score-ratio` | 0.9 | Minimum ratio of secondary-to-primary chaining score for a secondary alignment to be retained (0–1) |
| `-a, --min-score` | 20 | Minimum s1 chaining score for secondary alignments |
| `-N, --max-secondary` | 5 | Maximum number of secondary alignments per read emitted by minimap2. **Increase this (e.g. 20) when many highly similar reference genomes are present** |
| `-n, --top-n` | 25 | Number of top references for which alignment/depth PDFs are generated |
| `-q, --mapping_quality` | 20 | Minimum mapping quality for **primary** alignments |
| `-t, --threads` | 10 | Number of CPU threads |
| `-c, --force` | — | Force overwrite of existing output files |
| `-k, --skip-mapping` | — | Skip the mapping and coverage steps; reuse existing `Mito.paf.gz` and `Mito_coverage.txt`. Useful for re-running analysis with different parameters without re-mapping |
| `-x, --taxonomic-hierarchy` | species | Taxonomic level for classification (`species`, `genus`, `family`, `order`, `phylum`, `kingdom`) |
| `-v, --version` | — | Show version and exit |
| `-h, --help` | — | Show help and exit |

### Example — Full Run

```bash
bash shell/ECMSD.sh \
    --fwd reads_R1.fastq.gz \
    --rev reads_R2.fastq.gz \
    --merged merged.fastq.gz \
    --out results/ \
    --threads 20 \
    --cov-threshold 10 \
    --mapping_quality 20 \
    --score-ratio 0.9 \
    --max-secondary 20 \
    --taxonomic-hierarchy genus \
    --top-n 25 \
    --force
```

### Example — Rerun Analysis on Existing Mapping Results

```bash
bash shell/ECMSD.sh \
    --fwd reads_R1.fastq.gz \
    --rev reads_R2.fastq.gz \
    --out results/ \
    --cov-threshold 5 \
    --taxonomic-hierarchy genus \
    --skip-mapping
```

---

## Output Files

```text
output_directory/
├── mapping/
│   ├── Mito.paf.gz                          # Compressed PAF alignment file
│   ├── Mito_coverage.txt                    # Per-reference coverage statistics
│   ├── Mito_summary.txt                     # Per-alignment taxonomic assignments (all alignments)
│   ├── Mito_summary.ref_summary.txt         # Per-reference read counts and taxon labels (ranked)
│   ├── Mito_summary_<taxon>.txt             # Read counts per taxon × read length (primary reads only)
│   ├── Mito_summary_<taxon>_proportions.txt # Proportional abundance per taxon (primary reads only)
│   ├── Mito_summary_<taxon>_ReadLengths.png # Read length histogram for top 10 taxa
│   ├── Mito_summary_<taxon>_Proportions.png # Proportional abundance bar chart
│   └── alignment_plots/
│       └── <rank>_taxid<id>_<taxon>_<species>.pdf  # Coverage + depth PDF per reference
└── logs/
    └── logfile.log                          # minimap2 stderr and pipeline log
```

### Output File Descriptions

| File | Description |
|---|---|
| `Mito.paf.gz` | All filtered alignments in PAF format (primary MQ-filtered; secondaries s1-filtered) |
| `Mito_coverage.txt` | Columns: `reference`, `ref_length`, `mean_coverage`, `std_coverage`, `pct_covered` |
| `Mito_summary.txt` | One row per alignment: `SeqID`, `TaxID`, `Length`, `MappingQuality`, full lineage columns |
| `Mito_summary.ref_summary.txt` | Ranked list of references: `ref_name`, `taxid`, `taxon_name`, `species_name`, `read_count` |
| `Mito_summary_<taxon>_proportions.txt` | Proportions based on unique primary-alignment reads |
| `alignment_plots/*.pdf` | Two-panel PDF per reference: alignment coverage breadth (top) + sequencing depth (bottom) |

---

## Filtering Logic in Detail

### Primary Alignments

Retained when **mapping quality ≥ `--mapping_quality`** (MQ field, PAF column 12). This threshold is applied consistently at three points:

1. In the minimap2 awk post-filter (writing the PAF)
2. In the inline Python coverage calculation
3. In `LinkTaxonomy.py` (via `--MapQuality`)

### Secondary Alignments

minimap2 always sets MQ = 0 for secondary alignments. They are therefore **not** subject to the MQ filter. Instead they are retained when the **s1 chaining score ≥ `--min-score`** and the secondary/primary score ratio ≥ `--score-ratio` (controlled by minimap2 `-p`).

Use **`--max-secondary`** to control how many secondaries per read minimap2 emits. The default is 5; for datasets with many highly similar reference genomes (e.g. congeneric mitogenomes) increasing this to 20 or more ensures that all close relatives receive mapped reads.

### Coverage Threshold

After mapping, per-reference breadth of coverage is calculated. Only references where `pct_covered ≥ --cov-threshold` are passed to the taxonomic assignment step.

### Proportional Abundance

Read proportions (`process_files.R`) are computed from **unique primary-alignment reads only**: secondary alignments (MQ = 0) are excluded, and each read ID is counted once (`distinct(SeqID)`). This avoids inflating counts from multi-mapping reads.

### Alignment/Depth Plots

The per-reference PDFs (`plot_paf_alignments.R`) use **all alignments** from the PAF (including secondaries) to accurately represent true sequencing depth across the reference, which is the biologically relevant quantity for these plots.

---

## Parameters

### `--cov-threshold`

Minimum percentage of reference positions covered by at least one read. References below this threshold are excluded from taxonomy assignment and downstream plots. Lowering this value retains more references but may include spurious hits.

### `--score-ratio` / `--min-score`

Control which secondary alignments are retained. `--score-ratio` sets the minimap2 `-p` parameter (secondary score must be ≥ ratio × primary score). `--min-score` sets the minimum absolute s1 chaining score for secondaries.

### `--max-secondary`

Sets minimap2 `-N`. With a large reference database containing many closely related genomes, the default of 5 may silently discard valid hits to other congeners. Increase to 20–50 for genus- or family-level analyses.

### `--mapping_quality`

Filters primary alignments by reliability. Higher values reduce false positives but may miss divergent true positives. Not applied to secondaries (which always have MQ = 0 in minimap2 PAF output).

### `--taxonomic-hierarchy`

Determines the grouping level for summary tables and plots. Supported values: `species`, `genus`, `family`, `order`, `phylum`, `kingdom`.

### `--skip-mapping`

Reuses existing `Mito.paf.gz` and `Mito_coverage.txt` without re-running minimap2. Only the log directory is cleared. Useful when iterating over different `--cov-threshold`, `--taxonomic-hierarchy`, or `--top-n` values.

---

## Test Data

```bash
cd /path/to/ECMSD
bash shell/ECMSD.sh \
    --fwd TestData/merged.fastq.gz \
    --out TestOutput \
    --threads 20 \
    --cov-threshold 10 \
    --mapping_quality 20 \
    --taxonomic-hierarchy genus \
    --force
```

---

## Troubleshooting

| Symptom | Solution |
|---|---|
| `Conda environment could not be activated` | Run `conda init bash` and restart your shell |
| `Reference data not found` | Pipeline auto-downloads on first run; check internet connectivity and disk space (~5–10 GB) |
| `Output directory already exists` | Use `--force` to overwrite, or `--skip-mapping` to reuse existing mapping results |
| PAF reading very slow or out-of-memory in R | The PAF is pre-filtered with `awk` before entering R; ensure `zcat` is available |
| Too few secondary hits for closely related species | Increase `--max-secondary` (e.g. `--max-secondary 20`) |

---

## Version History

### v1.1 — Performance and correctness overhaul

**New parameters**

- `--max-secondary` (`-N`): exposes minimap2 `-N` (previously hardcoded to 5) to allow capturing all close secondary hits when many similar reference genomes are present
- `--skip-mapping` (`-k`): reuse existing PAF and coverage files, re-run only analysis and plotting steps

**Mapping / filtering fixes**

- Primary alignments are now MQ-filtered consistently at three points: awk post-filter, Python coverage script, and `LinkTaxonomy.py`
- Secondary alignments (MQ always 0 in minimap2) are correctly exempted from the MQ filter and instead filtered by s1 chaining score only
- `LinkTaxonomy.py` gained a `--MapQuality` argument to enforce the same MQ threshold as the awk filter

**Performance improvements**

- **Coverage calculation** (inline Python): replaced O(n × read_length) per-position dict updates with a vectorised difference-array algorithm — O(n + Σref_len)
- **`LinkTaxonomy.py`**: `taxon_trace()` is now memoised with `functools.lru_cache`, reducing repeated taxonomy-tree traversals from O(reads) to O(unique taxIds)
- **`process_files.R`**: `read.table()` replaced with `data.table::fread()` for fast I/O on large summary files
- **`plot_paf_alignments.R`**: `pafr::read_paf()` (slow pure-R parser that loads the entire file) replaced with a `zcat | awk | data.table::fread` pipeline that pre-filters to the top-N references before any data enters R, avoiding the 2^31-byte R string limit on large PAF files
- **`plot_paf_alignments.R`**: `compute_depth()` row-loop replaced with a vectorised difference-array using `tabulate()` + `cumsum()`
- **`pafr`** package removed as a dependency entirely

**Analysis changes**

- Proportional abundance (`process_files.R`): now uses **unique primary-alignment reads only** — secondary alignments (MQ = 0) are excluded and each read ID is deduplicated (`distinct(SeqID)`) before counting
- Alignment/depth PDFs (`plot_paf_alignments.R`): continue to use **all alignments** including secondaries for accurate depth representation
- PDF filenames now include both the taxon-rank label and the species binomial: `<rank>_taxid<id>_<taxon>_<species>.pdf`

### v1.0 — Initial release

---

## License

MIT License

## Author

Martin Kapun

## Contact

For questions, issues, or feature requests, please open a GitHub issue.

## Acknowledgments

- [minimap2](https://github.com/lh3/minimap2) — sequence alignment
- [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) — sequence masking
- [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) — mitochondrial reference genomes
- [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) — taxonomic classification
