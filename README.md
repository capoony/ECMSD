# ECMSD - Efficient Comprehensive Mitochondrial Sequence Detection Pipeline

## Overview

ECMSD (Efficient Comprehensive Mitochondrial Sequence Detection) is a bioinformatics pipeline designed for the detection and taxonomic classification of mitochondrial sequences from high-throughput sequencing data. The pipeline uses minimap2 for fast sequence alignment against a comprehensive mitochondrial reference database, filters references by breadth of coverage, and provides quantitative summaries and visualisation plots.

## Features

- **Fast mitochondrial sequence detection** using minimap2 alignment
- **Comprehensive reference database** built from NCBI RefSeq mitochondrial genomes
- **Taxonomic classification** with customizable hierarchy levels
- **Quality-based filtering** with configurable mapping quality thresholds
- **Coverage-based reference filtering**: only references with sufficient breadth of coverage are retained
- **Automated visualization** of results with R-based plotting
- **Support for multiple input formats**: single-end, paired-end, and merged reads
- **Decoupled database management**: the reference database is built once and shared across runs

## Pipeline Workflow

The workflow has two distinct phases:

### Phase 1 — One-time database setup

Build the reference database before running any analysis:

```
ECMSD --create-db --db-folder /path/to/db_folder
```

This downloads and processes NCBI RefSeq mitochondrial genomes and the NCBI taxonomy database. It only needs to be run once (or whenever you want to update the reference data). The database folder can then be shared across all subsequent analysis runs and across multiple parallel jobs.

### Phase 2 — Per-sample analysis

1. **Sequence Mapping**: Aligns input reads to the mitochondrial reference using minimap2
2. **Quality Filtering**: Filters alignments based on mapping quality threshold
3. **Coverage Calculation**: Computes per-reference breadth of coverage; references below the threshold are discarded
4. **Taxonomic Assignment**: Links aligned sequences to taxonomic information using NCBI taxonomy
5. **Visualization**: Generates summary plots and per-reference alignment coverage plots

## Requirements

- Linux/Unix operating system
- Conda package manager
- Internet connection (for database download)

### Software Dependencies (installed automatically via Bioconda or `install.sh`)

- minimap2
- BBTools (bbmask)
- Python 3
- R with tidyverse
- pigz
- wget

## Installation

### Recommended — Install via Bioconda

```bash
conda create -n ecmsd_env -c bioconda -c conda-forge ecmsd
conda activate ecmsd_env
```

This installs ECMSD and all dependencies in a single step.

### Alternative — Manual installation via `install.sh`

<details>
<summary>Click to expand manual installation steps</summary>

**Step 1 — Clone the repository**

```bash
git clone https://github.com/capoony/ECMSD.git
cd ECMSD
```

**Step 2 — Create and activate a conda environment**

You can use the provided environment file:

```bash
conda env create -f conda/ecmsd_env.yaml
conda activate ecmsd_env
```

Or create a minimal environment manually:

```bash
conda create -n ecmsd_env
conda activate ecmsd_env
```

**Step 3 — Install dependencies and scripts**

With your conda environment active, run the installation script:

```bash
bash shell/install.sh
```

This installs all required dependencies (minimap2, bbmap, R, Python packages, etc.) into the active conda environment and registers the `ECMSD` command.

</details>

### Build the reference database

Regardless of the installation method, the reference database must be created before running any analysis. This is a one-time step:

```bash
ECMSD --create-db --db-folder /path/to/db_folder
```

The database folder can be placed anywhere and reused across all future runs. This step downloads several GB of data from NCBI, so ensure sufficient disk space (~5–10 GB) and a stable internet connection.

> **Why is the database a separate step?**
> Separating database creation from analysis prevents race conditions when multiple ECMSD instances run in parallel (e.g. in a Snakemake workflow). In the old design, all parallel jobs would compete to create the same database simultaneously. By building the database once upfront, each analysis job simply reads from the shared, pre-built database. This also enables pipelines such as the **[pastForward Snakemake pipeline](https://github.com/SarahSaadain/PastForward)** to run contamination checks from scratch, downloading all required reference files in a single dedicated setup step before launching parallel analysis jobs.

## Test Data

ECMSD includes a test dataset to verify the installation and demonstrate pipeline functionality. The test data consists of historic DNA (hDNA) sequencing reads from an unknown origin.

### Test Dataset Details

- **File**: `TestData/merged.fastq.gz` (~37 MB)
- **Content**: Merged FASTQ reads from hDNA sequencing
- **Purpose**: Pipeline validation and demonstration
- **Origin**: Unknown sample (for testing purposes only)

### Running the Test

```bash
# Navigate to ECMSD directory
cd /path/to/ECMSD

# Step 1: create the database (first time only)
ECMSD --create-db --db-folder /path/to/db_folder

# Step 2: run the pipeline on the test data
ECMSD --fwd /path/to/TestData/merged.fastq.gz \
      --out TestOutput \
      --db-folder /path/to/db_folder \
      --threads 20 \
      --cov-threshold 10 \
      --top-n 25 \
      --mapping_quality 20 \
      --taxonomic-hierarchy genus \
      --force
```

### Expected Test Results

The test run should complete successfully and generate:

- Alignment results in `TestOutput/mapping/Mito.paf.gz`
- Taxonomic summary in `TestOutput/mapping/Mito_summary.txt`
- Visualization plots in `TestOutput/mapping/`
- Log files in `TestOutput/logs/`

### Test Data Validation

A successful test run indicates that:

- All dependencies are properly installed
- The reference database is correctly built
- The pipeline can process FASTQ data
- Output files are generated in the expected format

## Usage

### Basic Usage

```bash
# First: create the database (once)
ECMSD --create-db --db-folder /path/to/db

# Single-end reads
ECMSD -f reads.fastq -o output_directory/ -d /path/to/db

# Paired-end reads
ECMSD -f reads_R1.fastq -r reads_R2.fastq -o output_directory/ -d /path/to/db

# Merged reads only
ECMSD -m merged_reads.fastq -o output_directory/ -d /path/to/db

# Paired-end with additional merged reads
ECMSD -f reads_R1.fastq -r reads_R2.fastq -m merged_reads.fastq -o output_directory/ -d /path/to/db
```

### Command Line Options

#### Required Arguments for Analysis

- `-f, --fwd`: Path to the forward or single-end FASTQ file (required unless `--merged` is used)
- `-m, --merged`: Path to a merged FASTQ file (required unless `--fwd` is used; cannot be combined with `--fwd`/`--rev`)
- `-o, --out`: Path to the output folder
- `-d, --db-folder`: Path to the database folder (must be created first with `--create-db`)

#### Required Arguments for Database Creation

- `-z, --create-db`: Create a new reference database
- `-d, --db-folder`: Path to the folder where the database will be stored

#### Optional Arguments

- `-r, --rev`: Path to the reverse FASTQ file (for paired-end data; requires `--fwd`)
- `-u, --cov-threshold`: Minimum percentage of a reference covered by reads to retain it (default: 50)
- `-n, --top-n`: Number of top references to generate per-reference alignment plots for (default: 25)
- `-q, --mapping_quality`: Mapping quality threshold (default: 20)
- `-t, --threads`: Number of threads to use (default: 10)
- `-c, --force`: Force overwrite of existing output files
- `-k, --skip-mapping`: Skip mapping and coverage steps; reuse existing PAF and coverage files (implies `--force`)
- `-x, --taxonomic-hierarchy`: Taxonomic hierarchy level (default: species)
- `-p, --prefix`: Prefix for output files (default: None)
- `-v, --version`: Show version and exit
- `-h, --help`: Show help message and exit

### Advanced Examples

```bash
# High-stringency analysis with custom thresholds
ECMSD -f reads_R1.fastq -r reads_R2.fastq -o results/ -d /path/to/db \
  -q 30 -u 75 -n 50 -t 16

# Genus-level taxonomic classification
ECMSD -f reads.fastq -o results/ -d /path/to/db -x genus

# Force overwrite existing results
ECMSD -f reads.fastq -o existing_results/ -d /path/to/db -c

# Add a prefix to output files
ECMSD -f reads.fastq -o results/ -d /path/to/db -p sample01
```

### Running Multiple Instances in Parallel

Because the database is maintained separately, multiple ECMSD jobs can safely run in parallel against the same database without any risk of file conflicts:

```bash
# Build the database once
ECMSD --create-db --db-folder /shared/db

# Then run many samples concurrently (e.g. via a scheduler or Snakemake)
ECMSD -f sample1.fastq -o out_sample1/ -d /shared/db
ECMSD -f sample2.fastq -o out_sample2/ -d /shared/db
ECMSD -f sample3.fastq -o out_sample3/ -d /shared/db
```

## Output Files

The pipeline generates the following output structure:

```text
output_directory/
├── mapping/
│   ├── Mito.paf.gz                              # Compressed alignment results
│   ├── Mito_coverage.txt                        # Per-reference coverage statistics
│   ├── Mito_summary.txt                         # Taxonomic assignment per read
│   ├── Mito_summary.ref_summary.txt             # Ranked reference summary for plotting
│   ├── Mito_summary_<taxon>.txt                 # Per-taxon read count summary
│   ├── Mito_summary_<taxon>_ReadLengths.png     # Read length distribution plot
│   ├── Mito_summary_<taxon>_Proportions.png     # Proportions bar plot
│   └── alignment_plots/                         # Per-reference coverage and depth plots
└── logs/
    └── minimap2.log                             # Mapping log
```

### Output File Descriptions

- **Mito.paf.gz**: Compressed PAF format alignment file containing all mitochondrial sequence alignments
- **Mito_coverage.txt**: Per-reference coverage statistics — columns: `reference`, `ref_length`, `mean_coverage`, `std_coverage`, `pct_covered`
- **Mito_summary.txt**: Tab-separated file with taxonomic assignments and read counts per read
- **Mito_summary.ref_summary.txt**: Ranked reference summary used for generating alignment plots — columns: `ref_name`, `taxid`, `taxon_name`, `species_name`, `subspecies_name`, `read_count` (`subspecies_name` is `NA` for references not resolved below species level)
- **Mito_summary_\<taxon\>.txt**: Aggregated read counts per taxon at the chosen taxonomic level
- **Mito_summary_\<taxon\>_ReadLengths.png**: Faceted histogram of read length distributions for the top 10 taxa
- **Mito_summary_\<taxon\>_Proportions.png**: Bar plot of relative read proportions per taxon
- **alignment_plots/**: Two-panel PNG per reference showing alignment coverage breadth (top) and sequencing depth (bottom) for the top-N references
- **minimap2.log**: Detailed minimap2 alignment log

## Parameters Explanation

### Coverage-based Reference Filtering (`--cov-threshold`)

After mapping, ECMSD calculates per-reference breadth of coverage from the PAF file. For each reference, it computes the percentage of reference positions covered by at least one read (`pct_covered`). Only references where `pct_covered ≥ --cov-threshold` are passed to the taxonomic assignment step.

This replaces the previous RMUS-based filtering. Coverage breadth is a more direct and interpretable measure for distinguishing authentic mitochondrial hits (which should show broad, even coverage across the genome) from mapping artifacts or NUMTs (which tend to pile up on a small region of the reference).

The `Mito_coverage.txt` output file contains the full per-reference statistics for inspection.

### `--top-n`

Controls how many of the top-ranked references (by read count) receive a detailed per-reference alignment plot. Plots are written to `mapping/alignment_plots/` and show coverage breadth and sequencing depth across the reference. Increase this value to inspect more candidate references.

### `--skip-mapping`

Skips the minimap2 mapping and coverage calculation steps and reuses the existing `Mito.paf.gz` and `Mito_coverage.txt` files. Only the log directory is cleared. Useful when iterating over different `--cov-threshold`, `--taxonomic-hierarchy`, or `--top-n` values without re-running the (slow) mapping step.

### Mapping Quality

The mapping quality threshold filters alignments based on their reliability. Higher values increase stringency but may reduce sensitivity for divergent sequences.

### Taxonomic Hierarchy

Determines the taxonomic level for classification (e.g., species, genus, family, order, phylum or kingdom). Lower levels provide more specific identification but may have reduced sensitivity.

## Troubleshooting

### Common Issues

1. **"No conda environment is active"**
   - Activate a conda environment before running `install.sh`: `conda activate <your-env>`
   - If no environment exists yet, create one: `conda create -n ecmsd_env && conda activate ecmsd_env`

2. **"No database folder provided"** or **database files missing**
   - Build the database before running analysis: `ECMSD --create-db --db-folder /path/to/db`
   - Ensure sufficient disk space (~5–10 GB) and internet connectivity
   - If a previous database build was interrupted, re-run `--create-db` to complete it

3. **"Output directory already exists"**
   - Use the `-c` or `--force` flag to overwrite existing results
   - Or specify a different output directory

4. **Memory or disk space issues**
   - Reduce the number of threads with `-t`
   - Ensure sufficient disk space for the reference database (~5–10 GB)

5. **Parallel jobs failing or conflicting**
   - Ensure the database is fully built before launching parallel jobs
   - Each job must write to its own `--out` directory

### Performance Optimization

- Use more threads (`-t`) on systems with multiple CPU cores
- Use `--skip-mapping` (`-k`) to re-run analysis with different thresholds without repeating the mapping step
- Adjust `--cov-threshold` (`-u`) to control how strictly references are filtered by breadth of coverage
- Consider mapping quality threshold (`-q`) based on your sequencing platform
- Build the database on fast local storage to reduce I/O latency during mapping

## Citation

If you use ECMSD in your research, please cite:

[Citation information to be added]

## License

MIT License (MIT)

## Authors

- Martin Kapun (idea, design, algorithm development, implementation, testing, documentation)
- Sarah Saadain (conda environment setup, implementation, bioconda installation, testing, documentation)

## Contact

For questions, issues, or feature requests, please open a GitHub issue.

## Version History

- **v1.2.0**: Replaced RMUS-based filtering with coverage-based reference filtering (`--cov-threshold`). Added per-reference alignment plots (`--top-n`), `--skip-mapping` flag for faster iteration, and `plot_paf_alignments.R` for visualising coverage breadth and sequencing depth.
- **v1.1.0**: Decoupled database creation from analysis runs. The database is now built once with `--create-db --db-folder` and provided as a required parameter to all analysis runs. This prevents race conditions in parallel workflows and enables pipeline-level database management (e.g. [pastForward Snakemake pipeline](https://github.com/SarahSaadain/PastForward)). Added `install.sh` for explicit dependency installation.
- **v1.0.0**: Initial release with core mitochondrial detection functionality

## Acknowledgments

This pipeline uses several open-source tools:

- minimap2 for sequence alignment
- BBTools for sequence masking
- NCBI RefSeq database for reference genomes
- NCBI taxonomy database for taxonomic classification
