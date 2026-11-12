# ECMSD - Efficient Comprehensive Mitochondrial Sequence Detection Pipeline

## Overview

ECMSD (Efficient Comprehensive Mitochondrial Sequence Detection) is a bioinformatics pipeline designed for the detection and taxonomic classification of mitochondrial sequences from high-throughput sequencing data. The pipeline uses minimap2 for fast sequence alignment against a comprehensive mitochondrial reference database, filters references by breadth of coverage, and provides quantitative summaries and visualisation plots.

## Features

- **Fast mitochondrial sequence detection** using minimap2 alignment
- **Curated reference database** built from complete, reviewed NCBI RefSeq mitochondrial genomes (`NC_` accessions only)
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

This downloads and processes NCBI RefSeq mitochondrial genomes and the NCBI taxonomy database. It only needs to be run once (add `--force-rebuild` when you want to refresh the reference data). The database folder can then be shared across all subsequent analysis runs and across multiple parallel jobs.

### Phase 2 — Per-sample analysis

1. **Stage 1 — liberal screen**: Maps all reads to the mitochondrial reference with minimap2 using sensitive seeds (`-k 17 -w 7`) and relaxed chaining (`-n 1 -m 20 -s 20`), without MAPQ filtering (MAPQ 0 is kept) and without base-level alignment. Every read with at least one alignment is a candidate; the candidates are extracted unchanged from the original FASTQ file(s) (pairs stay together)
2. **Stage 2 — accurate mapping**: Remaps only the candidates to the complete database with base-level alignment (`-c`) and secondary alignments (`--secondary=yes -N 20`). The output is streamed straight into the best-alignment selection, so the very large all-alignments file is not stored. For each read all alignments with the best alignment score (`AS`) are kept (ties included), never ranked by MAPQ
3. **Reference Evidence**: For each reference, separates reads that are *compatible* with it (it has the best `AS`, possibly tied with other references) from reads that *discriminate* it from all other references. A reference is only reported when enough discriminatory reads, spread over enough of its genome, support it (see [Taxonomic assignment](#taxonomic-assignment))
4. **Coverage Calculation**: Computes per-reference breadth of coverage from all compatible alignments; references below the threshold are discarded
5. **Taxonomic Assignment**: Links reads to taxonomic information using NCBI taxonomy; only credible references are reported, always at the level of the reference (never at a higher taxonomic level); a read compatible with several credible references is assigned to the one with the most discriminatory reads
6. **Visualization**: Generates summary plots and per-reference alignment coverage plots

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

**Alternative — let `install.sh` create a new local environment**

Instead of steps 2 and 3, you can give `install.sh` a path with `--prefix` (`-p`). It creates a new conda environment at that location (e.g. on a shared project drive), installs all dependencies and registers `ECMSD` there. No active conda environment is needed:

```bash
bash shell/install.sh --prefix /path/to/ecmsd_env
conda activate /path/to/ecmsd_env
```

If the path is already a conda environment, ECMSD is (re)installed into it. A non-empty directory that is not a conda environment is refused.

</details>

### Build the reference database

Regardless of the installation method, the reference database must be created before running any analysis. This is a one-time step:

```bash
ECMSD --create-db --db-folder /path/to/db_folder
```

The database folder can be placed anywhere and reused across all future runs. This step downloads several GB of data from NCBI, so ensure sufficient disk space (~5–10 GB during the build) and a stable internet connection. Build intermediates are deleted once the database is complete, leaving a finished folder of well under 1 GB.

Running `--create-db` against an existing, complete database is a no-op — it reports the database as ready and exits without re-downloading anything. This makes it safe to call unconditionally at the start of a workflow. To refresh the reference data against the current NCBI release, add `--force-rebuild`:

```bash
ECMSD --create-db --force-rebuild --db-folder /path/to/db_folder
```

This discards the existing database files and downloads everything again. Only files ECMSD created are removed, so a database folder shared with other content is left otherwise intact. Do not use `--force-rebuild` while other ECMSD jobs are reading the same database folder.

> **Why is the database a separate step?**
> Separating database creation from analysis prevents race conditions when multiple ECMSD instances run in parallel (e.g. in a Snakemake workflow). In the old design, all parallel jobs would compete to create the same database simultaneously. By building the database once upfront, each analysis job simply reads from the shared, pre-built database. This also enables pipelines such as the **[pastForward Snakemake pipeline](https://github.com/SarahSaadain/PastForward)** to run contamination checks from scratch, downloading all required reference files in a single dedicated setup step before launching parallel analysis jobs.

### What goes into the reference database

The database is built from the NCBI RefSeq mitochondrion release, but **only `NC_` accessions are retained**. `NW_` records are discarded.

RefSeq accession prefixes encode how a sequence was produced and how much curation it received:

| Prefix | Meaning | Kept? |
| --- | --- | --- |
| `NC_` | Complete genomic molecule, curated RefSeq — a finished, reviewed mitogenome that someone deliberately sequenced, assembled and checked | ✅ |
| `NW_` | Unplaced/unlocalised genomic scaffold, auto-generated from a whole-genome shotgun assembly — a leftover contig that automated screening flagged as mitochondrial, with no curation and no guarantee of completeness | ❌ |

For reference, other prefixes you may encounter in RefSeq: `NM_`/`NP_` (curated mRNA/protein), `NG_` (curated gene region), `NT_` (assembly scaffold, like `NW_`), `XM_`/`XP_` (automatically predicted, not curated at all).

Including `NW_` scaffolds actively degrades results:

- **Partial and low-quality references.** An `NW_` entry may cover a fraction of the mitogenome, or include nuclear sequence. Reads mapping to such a reference still count towards a taxon, but its breadth of coverage is not comparable to that of a complete mitogenome, which makes `--cov-threshold` behave inconsistently across references.
- **NUMT contamination.** WGS scaffolds are exactly where nuclear mitochondrial insertions (NUMTs) end up. Mapping to a NUMT-bearing scaffold produces confident-looking hits that are not mitochondrial in origin.
- **They can displace the curated genome.** The database keeps one reference per taxid, so whichever record for a taxon is encountered first wins. An `NW_` scaffold appearing before the taxon's `NC_` mitogenome would silently replace it.

The filter is applied in `shell/MakeRef.sh` (`extract_taxids()`), which restricts the accession-to-taxid mapping to `NC_`. Since `scripts/renameFASTA_taxid.py` only writes sequences whose accession appears in that mapping, `NW_` records never reach the FASTA.

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
```

### Known limitations

- The database keeps one genome per taxid (the first sequence encountered wins), so multi-record mitogenomes lose their other records.

### Command Line Options

#### Required Arguments for Analysis

- `-f, --fwd`: Path to the forward or single-end FASTQ file (required unless `--merged` is used)
- `-m, --merged`: Path to a merged FASTQ file (required unless `--fwd` is used; cannot be combined with `--fwd`/`--rev`)
- `-o, --out`: Path to the output folder
- `-d, --db-folder`: Path to the database folder (must be created first with `--create-db`)

#### Required Arguments for Database Creation

- `-z, --create-db`: Create a new reference database
- `-d, --db-folder`: Path to the folder where the database will be stored
- `-b, --force-rebuild`: Optional — discard an existing database and rebuild it from scratch (implies `--create-db`)

#### Optional Arguments

- `-r, --rev`: Path to the reverse FASTQ file (for paired-end data; requires `--fwd`)
- `-u, --cov-threshold`: Minimum percentage of a reference covered by reads to retain it (default: 25; see [Coverage-based Reference Filtering](#coverage-based-reference-filtering---cov-threshold) for low-depth samples)
- `-n, --top-n`: Number of top references to generate per-reference alignment plots for (default: 25)
- `-q, --mapping_quality`: Drop reads whose best alignments all have a lower MAPQ (default: 0, i.e. no MAPQ filtering; see [Mapping Quality](#mapping-quality))
- `--seed-k`: minimap2 k-mer size, both stages (default: 17)
- `--seed-w`: minimap2 minimizer window size, both stages (default: 7)
- `--min-chain-score`: minimap2 minimal chaining score `-m`, both stages (default: 20). Higher values make Stage 1 stricter and Stage 2 faster; see [Two-stage mapping parameters](#two-stage-mapping-parameters)
- `--max-alns`: Stage 2 `minimap2 -N`, maximum number of secondary alignments per read (default: 20)
- `--min-disc-reads`: A reference is only reported when at least this many reads discriminate it from all other references (default: 20)
- `--min-disc-breadth`: ... and these reads cover at least this % of the reference (default: 5)
- `--min-disc-evenness`: ... and are spread over the reference with at least this evenness, 0-1 (default: 0.5)
- `--min-aln-len`: Minimum alignment length in bp of a read that counts as discriminatory evidence (default: 20). Lower it for very short reads, raise it to be stricter
- `--keep-all-alignments`: Also keep every Stage 2 alignment as `Mito_stage2_all_alignments.paf.gz` and `Mito_stage2_alignments.tsv.gz` (large, and adds several minutes on deep samples)
- `-t, --threads`: Number of threads to use (default: 10)
- `-c, --force`: Force overwrite of existing results (removes only `<out>/mapping` and `<out>/logs`)
- `-k, --skip-mapping`: Skip mapping and coverage steps; reuse existing PAF and coverage files (implies `--force`)
- `-x, --taxonomic-hierarchy`: Taxonomic hierarchy level (default: species)
- `-p, --prefix`: Prefix for output files (default: None)
- `-v, --version`: Show version and exit
- `-h, --help`: Show help message and exit

### Advanced Examples

```bash
# High-stringency analysis with custom thresholds
ECMSD -f reads_R1.fastq -r reads_R2.fastq -o results/ -d /path/to/db \
  -q 30 -u 25 -n 50 -t 16

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
│   ├── Mito.paf.gz                              # Best Stage 2 alignment per read
│   ├── Mito_stage1_hits.txt                     # Names of reads with a Stage 1 alignment
│   ├── Mito_stage1_candidates[_R1|_R2].fastq.gz # Candidate reads, input of Stage 2
│   ├── Mito_stage2_all_alignments.paf.gz        # Only with --keep-all-alignments
│   ├── Mito_stage2_alignments.tsv.gz            # Only with --keep-all-alignments
│   ├── Mito_reference_evidence.txt              # Compatible and discriminatory evidence per reference
│   ├── Mito_coverage.txt                        # Per-reference coverage statistics
│   ├── Mito_summary.txt                         # Taxonomic assignment per read
│   ├── Mito_summary.ref_summary.txt             # Ranked reference summary for plotting
│   ├── Mito_summary_<taxon>.txt                 # Per-taxon read count summary
│   ├── Mito_summary_<taxon>_ReadLengths.png     # Read length distribution plot
│   ├── Mito_summary_<taxon>_Proportions.png     # Proportions bar plot
│   └── alignment_plots/                         # Per-reference coverage and depth plots
└── logs/
    ├── minimap2_stage1.log                      # Stage 1 mapping log
    └── minimap2_stage2.log                      # Stage 2 mapping log
```

### Output File Descriptions

- **Mito.paf.gz**: Compressed PAF with, for each read, the alignments to all references that tie for its best `AS` (the references the read is compatible with). Each line carries the extra tags `nc` (number of compatible references), `dm` (`AS` margin to the best lower-scoring reference, if one was reported) and `na` (number of alignments minimap2 reported for the read). Evidence, coverage, taxonomy and plots are based on this file
- **Mito_reference_evidence.txt**: One row per reference with read counts (`compatible_reads`, `discriminatory_reads`, `weak_unique_reads`, `shared_reads`, `discriminatory_at_limit`), breadth, depth and evenness of the compatible and of the discriminatory reads, whether the reference is `credible` and, if not, the `reason`
- **Mito_stage1_hits.txt**, **Mito_stage1_candidates*.fastq.gz**: Stage 1 output. Hit names (any MAPQ) and the candidate reads extracted from the original FASTQ file(s)
- **Mito_stage2_all_alignments.paf.gz**, **Mito_stage2_alignments.tsv.gz** (only with `--keep-all-alignments`): all Stage 2 alignments, including secondary ones (`tp:A:S`). The table lists reference, coordinates, aligned length, MAPQ, `tp`, `AS`, `NM`, whether the alignment was chosen as best and its `AS` difference to the best alignment of the read, so ambiguous reads can be inspected
- **Mito_coverage.txt**: Per-reference coverage statistics — columns: `reference`, `ref_length`, `mean_coverage`, `std_coverage`, `pct_covered`
- **Mito_summary.txt**: Tab-separated file with taxonomic assignments and read counts per read
- **Mito_summary.ref_summary.txt**: Ranked reference summary used for generating alignment plots — columns: `ref_name`, `taxid`, `taxon_name`, `species_name`, `subspecies_name`, `read_count` (`subspecies_name` is `NA` for references not resolved below species level)
- **Mito_summary_\<taxon\>.txt**: Aggregated read counts per taxon at the chosen taxonomic level
- **Mito_summary_\<taxon\>_ReadLengths.png**: Faceted histogram of read length distributions for the top 10 taxa
- **Mito_summary_\<taxon\>_Proportions.png**: Bar plot of relative read proportions per taxon
- **alignment_plots/**: Two-panel PNG per reference showing alignment coverage breadth (top) and sequencing depth (bottom) for the top-N references. Titles and file names use the label for the rank requested with `-x`, so a run at species level never shows a subspecies name
- **minimap2_stage1.log**, **minimap2_stage2.log**: Detailed minimap2 logs of the two stages

## Parameters Explanation

### Coverage-based Reference Filtering (`--cov-threshold`)

After mapping, ECMSD calculates per-reference breadth of coverage from the PAF file. For each reference, it computes the percentage of reference positions covered by at least one read (`pct_covered`). Only references where `pct_covered ≥ --cov-threshold` are passed to the taxonomic assignment step.

This replaces the previous RMUS-based filtering. Coverage breadth is a more direct and interpretable measure for distinguishing authentic mitochondrial hits (which should show broad, even coverage across the genome) from mapping artifacts or NUMTs (which tend to pile up on a small region of the reference).

The default is 25%. This threshold is on absolute breadth, and breadth grows with sequencing depth, so a species that is present at low abundance can fall below it in a shallow sample. On the 1.09 M read test sample (*D. melanogaster* plus a rare *G. bimaculatus*), 25% reports exactly these two species. Randomly subsampled to 10% of the reads, *G. bimaculatus* is still credible by its discriminatory reads (27 reads, 19% breadth) but is dropped by the 25% threshold; at 3% and 1% it has too few discriminatory reads to be reported at any threshold. In that range the discriminatory-evidence criteria (`--min-disc-*`) already reject the weak candidates that 25% is needed for at full depth, so for shallow samples a lower `--cov-threshold` (for example 5-10) is reasonable. `Mito_reference_evidence.txt` shows both values for every reference, so a borderline result can be checked.

The `Mito_coverage.txt` output file contains the full per-reference statistics for inspection.

### `--top-n`

Controls how many of the top-ranked references (by read count) receive a detailed per-reference alignment plot. Plots are written to `mapping/alignment_plots/` and show coverage breadth and sequencing depth across the reference. Increase this value to inspect more candidate references.

### `--skip-mapping`

Skips the minimap2 mapping and coverage calculation steps and reuses the existing `Mito.paf.gz` and `Mito_coverage.txt` files. Only the log directory is cleared. Useful when iterating over different `--cov-threshold`, `--taxonomic-hierarchy`, or `--top-n` values without re-running the (slow) mapping step.

### Mapping Quality

MAPQ is low (often 0) whenever a read fits several closely related mitogenomes equally well, even if the alignment is correct. ECMSD therefore does not filter on MAPQ by default (`-q 0`) and chooses between references by alignment score (`AS`). `-q` can still be used to drop reads whose best alignments all have a lower MAPQ, which discards ambiguous reads.

### Taxonomic assignment

Closely related mitogenomes share many reads. A species that is not in the sample but is almost identical to one that is can therefore reach a high breadth of coverage on reads that do not tell the two apart, so breadth alone cannot separate it from a species that is present. ECMSD thus counts, for each reference, two kinds of reads:

- **compatible**: the reference has the best alignment score (`AS`) for the read, including reads that tie with other references;
- **discriminatory**: the reference is the only one with the best `AS`, by at least 10 `AS` points over every competing reference that was reported, and the alignment is adequate (at least `--min-aln-len` bp, 20 by default, at least half of the read aligned, at least 90% identity). Short chance matches are `weak_unique_reads`, not evidence.

A reference is reported (`credible`) when it has at least `--min-disc-reads` discriminatory reads that cover at least `--min-disc-breadth` % of it and are spread over it with an evenness (normalised Shannon entropy over 20 equal windows) of at least `--min-disc-evenness`. No threshold depends on how abundant other species are, so a rare species with well-spread discriminatory reads is kept. Reads shared between references still count towards the coverage of every credible reference they are compatible with, and a read that is compatible with several credible references is assigned to the one with the most discriminatory reads (ties: the first one reported). Only credible references are reported, always as the reference itself; nothing is reported at a higher taxonomic level. If a sample only contains shared reads, no reference is credible and nothing is reported.

Limits: minimap2 reports at most `--max-alns` (+1) alignments per read, chosen before base-level alignment, so a competing reference may not be reported for a read. A unique best alignment is therefore only provisional; `discriminatory_at_limit` in `Mito_reference_evidence.txt` counts discriminatory reads for which the limit was reached. Species that are credible only because of reads at the limit, or because of a few reads, should be treated with caution. The thresholds are defaults chosen on one test data set and may need adjusting.

### Two-stage mapping parameters

The `sr` preset requires two minimizers per chain (`-n 2`), a chaining score of 25 (`-m 25`) and a peak score of 40 (`-s 40`). Reads shorter than about 30 bp are almost never mapped with these. Both stages use `-n 1 -m 20 -s 20` together with `-k 17 -w 7` instead.

The chaining score `-m` is what controls how many reads Stage 1 passes on. On the 1.09 M read test sample (Drosophila melanogaster and Gryllus bimaculatus reads against the full RefSeq mitochondrial database), `-m 10` passed 85% of all reads, of which 14% could only align with a very low alignment score (AS < 40), which made Stage 2 slow (16.6 M alignments). `-m 20` passes 59% of the reads, none of them with AS < 40, and still keeps 99% of the *D. melanogaster* reads, 99% of the reads with AS >= 60 and 94% of the *G. bimaculatus* reads (98.7% of those with an alignment of AS >= 40). `-n` and `-s` have no effect once `-m` is 20 or higher. Use `--min-chain-score 18` to recover a few percent more of the weakest hits (about 97% of *G. bimaculatus*) at the cost of 5 percentage points more candidates, or 25 (the preset value) for a stricter screen that loses about 7% of the *G. bimaculatus* reads.

### Taxonomic Hierarchy

Determines the taxonomic level for classification (e.g., species, genus, family, order, phylum or kingdom). Lower levels provide more specific identification but may have reduced sensitivity.

## Troubleshooting

### Common Issues

1. **"No conda environment is active"**
   - Activate a conda environment before running `install.sh`: `conda activate <your-env>`
   - If no environment exists yet, create one: `conda create -n ecmsd_env && conda activate ecmsd_env`
   - Or let the script create one at a path of your choice: `bash shell/install.sh --prefix /path/to/ecmsd_env`

2. **"No database folder provided"** or **database files missing**
   - Analysis never builds the database; it exits with "Database not found. Run `ECMSD --create-db` first". Build it first: `ECMSD --create-db --db-folder /path/to/db`
   - Ensure sufficient disk space (~5–10 GB) and internet connectivity
   - If a previous database build was interrupted, re-run `--create-db` to complete it

3. **"Output directory already exists"**
   - Use the `-c` or `--force` flag to overwrite existing results
   - Or specify a different output directory

4. **Memory or disk space issues**
   - Reduce the number of threads with `-t`
   - Ensure sufficient disk space for the database build (~5–10 GB peak; the finished database is well under 1 GB)

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

> Saadain, S., Kapun, M., & Kofler, R. (2026). pastForward: a Snakemake pipeline for ancient and historical DNA with eukaryote-wide taxonomic screening and tracking of copy-number variation. *bioRxiv*. https://doi.org/10.64898/2026.08.07.743613


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
