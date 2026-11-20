# ECMSD - Efficient Comprehensive Mitochondrial Sequence Detection Pipeline

## Overview

ECMSD (Efficient Comprehensive Mitochondrial Sequence Detection) is a bioinformatics pipeline designed for the detection and taxonomic classification of mitochondrial sequences from high-throughput sequencing data. The pipeline uses minimap2 for fast sequence alignment against a comprehensive mitochondrial reference database, filters references by breadth of coverage, and provides quantitative summaries and visualisation plots.

## Features

- Fast mitochondrial sequence detection with minimap2
- A reference database of complete, reviewed NCBI RefSeq mitochondrial genomes (`NC_` accessions only)
- Taxonomic classification at a hierarchy level you choose
- Configurable mapping quality threshold
- Only references with sufficient breadth of coverage are retained
- Plots of the results, made with R
- Single-end, paired-end and merged reads as input
- The reference database is built once and shared across runs

## Pipeline Workflow

The workflow has two distinct phases:

### Phase 1: one-time database setup

Build the reference database before running any analysis:

```
ECMSD --create-db --db-folder /path/to/db_folder
```

This downloads and processes NCBI RefSeq mitochondrial genomes and the NCBI taxonomy database. It only needs to be run once (add `--force-rebuild` when you want to refresh the reference data). The database folder can then be shared across all subsequent analysis runs and across multiple parallel jobs.

### Phase 2: per-sample analysis

1. **Mapping**: Maps all reads once to the complete mitochondrial database with minimap2, using sensitive seeds (`-k 17 -w 7`), relaxed chaining (`-n 1 -m 20 -s 20`), base-level alignment (`-c`) and secondary alignments (`--secondary=yes -N 20`), without MAPQ filtering (MAPQ 0 is kept). The output is streamed straight into the best-alignment selection, so the very large all-alignments file is not stored. For each read all alignments with the best alignment score (`AS`) are kept (ties included), never ranked by MAPQ
2. **Reference Evidence**: For each reference, separates reads that are *compatible* with it (it has the best `AS`, possibly tied with other references) from reads that *discriminate* it from all other references. A reference is only reported when enough discriminatory reads, spread over enough of its genome, support it (see [Taxonomic assignment](#taxonomic-assignment))
3. **Coverage Calculation**: Computes per-reference breadth of coverage from all compatible alignments; references below the threshold are discarded
4. **Taxonomic Assignment**: Links reads to taxonomic information using NCBI taxonomy; only credible references are reported, always at the level of the reference (never at a higher taxonomic level); a read compatible with several credible references is split between them in the proportion of their discriminatory reads
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

### Recommended: install via Bioconda

```bash
conda create -n ecmsd_env -c bioconda -c conda-forge ecmsd
conda activate ecmsd_env
```

This installs ECMSD and all dependencies in a single step.

### Alternative: manual installation via `install.sh`

<details>
<summary>Click to expand manual installation steps</summary>

**Step 1: clone the repository**

```bash
git clone https://github.com/capoony/ECMSD.git
cd ECMSD
```

**Step 2: create and activate a conda environment**

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

**Step 3: install dependencies and scripts**

With your conda environment active, run the installation script:

```bash
bash shell/install.sh
```

This installs all required dependencies (minimap2, bbmap, R, Python packages, etc.) into the active conda environment and registers the `ECMSD` command.

**Alternative: let `install.sh` create a new local environment**

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

The database folder can be placed anywhere and reused across all future runs. This step downloads several GB of data from NCBI, so ensure sufficient disk space (about 5 to 10 GB during the build) and a stable internet connection. Build intermediates are deleted once the database is complete, leaving a finished folder of well under 1 GB.

Running `--create-db` against an existing, complete database is a no-op: it reports the database as ready and exits without re-downloading anything. This makes it safe to call unconditionally at the start of a workflow. To refresh the reference data against the current NCBI release, add `--force-rebuild`:

```bash
ECMSD --create-db --force-rebuild --db-folder /path/to/db_folder
```

This discards the existing database files and downloads everything again. Only files ECMSD created are removed, so a database folder shared with other content is left otherwise intact. Do not use `--force-rebuild` while other ECMSD jobs are reading the same database folder.

> **Why is the database a separate step?**
> Separating database creation from analysis prevents race conditions when multiple ECMSD instances run in parallel (e.g. in a Snakemake workflow). If each job built its own database, they would compete to create the same files. Because the database is built once upfront, each analysis job only reads from the shared, pre-built database. This also enables pipelines such as the **[pastForward Snakemake pipeline](https://github.com/SarahSaadain/PastForward)** to run contamination checks from scratch, downloading all required reference files in a single dedicated setup step before launching parallel analysis jobs.

### What goes into the reference database

The database is built from the NCBI RefSeq mitochondrion release, but **only `NC_` accessions are retained**. `NW_` records are discarded.

RefSeq accession prefixes encode how a sequence was produced and how much curation it received:


| Prefix | Meaning                                                                                                                                                                                                              | Kept? |
| ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----- |
| `NC_`  | Complete genomic molecule, curated RefSeq: a finished, reviewed mitogenome that someone deliberately sequenced, assembled and checked                                                                               | Yes   |
| `NW_`  | Unplaced/unlocalised genomic scaffold, auto-generated from a whole-genome shotgun assembly: a leftover contig that automated screening flagged as mitochondrial, with no curation and no guarantee of completeness  | No    |

For reference, other prefixes you may encounter in RefSeq: `NM_`/`NP_` (curated mRNA/protein), `NG_` (curated gene region), `NT_` (assembly scaffold, like `NW_`), `XM_`/`XP_` (automatically predicted, not curated at all).

Including `NW_` scaffolds actively degrades results:

- **Partial and low-quality references.** An `NW_` entry may cover a fraction of the mitogenome, or include nuclear sequence. Reads mapping to such a reference still count towards a taxon, but its breadth of coverage is not comparable to that of a complete mitogenome, which makes `--cov-threshold` behave inconsistently across references.
- **NUMT contamination.** WGS scaffolds are exactly where nuclear mitochondrial insertions (NUMTs) end up. Mapping to a NUMT-bearing scaffold produces confident-looking hits that are not mitochondrial in origin.
- **They can displace the curated genome.** The database keeps one reference per taxid, so whichever record for a taxon is encountered first wins. An `NW_` scaffold appearing before the taxon's `NC_` mitogenome would silently replace it.

The filter is applied in `shell/MakeRef.sh` (`extract_taxids()`), which restricts the accession-to-taxid mapping to `NC_`. Since `scripts/renameFASTA_taxid.py` only writes sequences whose accession appears in that mapping, `NW_` records never reach the FASTA.

### Database version, updates and extensions

Every build writes `db_info.txt` into the database folder: the ECMSD version, the build date, the RefSeq release, the date of the NCBI taxonomy, the number of references and of extra references, and the date of the last extension. Each analysis run prints it, and the second header line of `Mito_calls.txt` summarises it (`# db: RefSeq release 238, built ...`), so every result can be traced to its database. The full NCBI taxonomy (`nodes.dmp`, `names.dmp`, `merged.dmp`, gzipped, about 70 MB) is kept in `NCBI_taxdump/full/`; the analysis reads only the trimmed copies.

**Updating.** `ECMSD --update-db --db-folder DB` compares the RefSeq release of the database with the current one and does nothing when they are the same. Otherwise it builds the new release in `DB/update.new/` while the old database stays usable, then moves the finished files into place and removes the staging folder. Extra references are kept: their taxids are checked against the new taxonomy (merged taxids are replaced). The new taxonomy also covers the old references, so a job that is running during the swap still finds its lineages. Still, results of one study should come from one database version: do not update while a workflow is running. Updating needs the disk space of a full build next to the existing database.

**Extending.** The RefSeq database has one mitogenome per species. When a sample carries a haplotype that is not in it (for example *D. simulans* siIII, while RefSeq has siI), its reads fit the species' reference poorly, and the call is labelled `poor_fit` instead of `confident`. Known haplotypes or species can be added:

```bash
# NCBI accessions, one per line; sequences and taxids are fetched from NCBI
ECMSD --extend-db --accessions accessions.txt --db-folder DB
# or your own sequences, with headers >taxid|id (e.g. >7240|AF200848.1)
ECMSD --extend-db --fasta my_mitogenomes.fa --db-folder DB
```

- The new sequences go into `extra_refs.fna.gz`; RefSeq is not changed. `extensions.log` lists every added sequence with date, source, taxid, length and warnings. Delete `extra_refs.fna.gz` to undo all extensions (the extra taxonomy lines are harmless).
- Each record is checked against the stored taxonomy: unknown taxids are rejected, merged taxids are replaced, ids already in the database are skipped, and sequences under 5 kb and RefSeq `NC_` accessions are added with a warning. Low-complexity regions are masked as in the main database.
- With several references of one taxid, the analysis turns on `--group-by-taxid` and sets `--max-alns` to 20 plus the largest number of references of one taxid, unless you set it.
- Only add what you have a reason to expect. A mitochondrial haplotype that introgressed from another species gives the species of the mitochondrion, not of the nuclear genome: a simulated *D. mauritiana* sample with the maI haplotype, which is almost identical to *D. simulans* siIII, is called *D. simulans* with confidence once siIII is in the database.
- Every file is replaced, never edited in place, so a running analysis sees either the old or the new set. To keep the original database unchanged, extend a hard-linked copy: `cp -al DB DB_ext`.

One database job (build, update or extension) runs at a time per folder; a second one waits for the first.

## Test Data

ECMSD includes a test dataset to verify the installation and demonstrate pipeline functionality. The test data consists of historic DNA (hDNA) sequencing reads from an unknown origin.

### Test Dataset Details

- **File**: `TestData/merged.fastq.gz` (~37 MB)
- **Content**: Merged FASTQ reads from hDNA sequencing
- **Purpose**: Pipeline validation and demonstration
- **Origin**: Unknown sample (for testing purposes only)
- **Poor-fit sample**: `TestData/merged_poor_fit.fastq.gz` (~7 MB) holds the first 200,000 reads of `merged.fastq.gz` with 6% of the bases substituted (`TestData/make_poor_fit.py`, seeded). Its reads fit every reference poorly, so `TestData/test_poor_fit.sh` expects `poor_fit` calls and no `confident` call.

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

Only `confident` rows of `Mito_calls.txt` are species calls. The limits below were found in tests on simulated and real reads with known species (25 scenarios, plus depth tests). Each is marked as addressed, partly addressed or open, with the proposed remedy where there is one.

**Addressed**

- **A species missing from the database is reported as its closest relatives.** The relatives used to be reported as credible hits. They now get the label `poor_fit` or `ambiguous`, and the sample verdict names the genus. In the tests no wrong species was labelled `confident`; the old version did this in 12 of 25 scenarios. The relatives are still listed, and they are not species calls. See [Identity](#identity-right-species-or-closest-relative).
- **A haplotype missing from the database.** RefSeq holds one mitogenome per species, so the right species can fit worse than a sister species (the *D. simulans* sample called *D. mauritiana*). Without added references it is labelled `poor_fit`. `--extend-db` adds known haplotypes, which are counted as one unit per species; with 27 *D. simulans* complex haplotypes added, that sample is called *D. simulans* with confidence.
- **Contaminated reference sequences.** A stretch of foreign sequence in a reference (*Exorista sorbillans*, taxid 191189, has a Drosophila segment of about 4.5 kb; *D. formosana*, 157056, a COI-COII stretch 99% identical to *D. melanogaster*) gave credible hits. `--min-disc-window-ratio` removes references whose discriminatory reads pile up in a few windows.
- **Shared reads.** A read shared by several credible references is split in proportion to their discriminatory reads, not given to one winner.
- **PCR duplicates** no longer inflate the discriminatory read counts and identity columns (`--no-dedup` switches this off).
- **Multi-record mitogenomes.** The database keeps one genome per taxid (the first sequence encountered wins). Further mitogenomes of a taxid can be added with `--extend-db`.
- **Competing alignments cut off by minimap2** (`--max-alns`). Raising the limit hardly changes the calls and costs 4-5 times the mapping time, so it is kept. Calls that depend on reads at the limit are labelled `low_evidence`. With extra references, `--max-alns` is raised automatically.
- **Aligner damage.** C to T and G to A damage at read ends did not change the calls in tests with simulated reads at 47 and 70 bp. Wrong aDNA calls came from short reads and missing haplotypes. A damage-aware score was tried and rejected (it creates reads that look unique by mistake).

**Open, with a proposed solution**

- **Mitochondrial introgression.** A mitochondrion that moved between species gives the species of the mitochondrion, not of the nuclear genome. *D. mauritiana* maI is almost identical to *D. simulans* siIII, so maI reads are called *D. simulans* with confidence once siIII is in the database. No read-level rule can detect this. Use nuclear markers to decide, and add haplotypes only when you expect them.
- **A species or haplotype that is not in the database is never `confident`.** Absent species are still named, as closest relatives with `poor_fit` or `ambiguous`; the sample verdict says "an unsampled genus species or haplotype". Proposed: extend the database with the species or haplotype you expect (`--extend-db`), and run a Metazoa-only screen for references with long near-identical stretches to another genus (the all-against-all screen of the full database did not finish).
- **Species complexes cannot be separated.** Species that share most of their reads are grouped as `ambiguous` and reported as "one of the group". The label shows the problem; it does not solve it. Proposed: report such groups at a higher level, or add a nuclear marker step.
- **Thresholds were chosen on the test data.** `--min-decline` (0.9), `--min-breadth-ratio` (0.5) and `--max-overlap` (0.5) were set after looking at the benchmark, and only one real in-database species (*D. melanogaster*) was tested; its breadth ratio is 0.05-0.15 above the cutoff. A species whose reachable fraction of the mitogenome is below 0.5 is always labelled `poor_fit` (this errs on the safe side). The decline reacts to batches of poor reads (0.76 in one batch of the *D. simulans* sample). Proposed: calibrate on more real samples of known species, and on a second simulation seed.
- **Breadth depends on depth.** `--cov-threshold` (default 25) is an absolute breadth, which grows with sequencing depth, so a rare species in a shallow sample can fall below it and is not reported. In the depth test the *D. melanogaster* hit was hidden by the threshold at 10% of the reads and had too few discriminatory reads at 5%. The label rules use discriminatory reads and do not depend on this threshold. Proposed: lower `--cov-threshold` for low-depth samples, or replace it by the expected breadth given the depth.
- **Closest relative ranks first when reads fit equally.** Equal identities on several references say nothing about which is right. Proposed: EM reassignment over the credible references only. In tests it ranked the real *D. simulans* sample correctly (*D. simulans* 0.37 against *D. mauritiana* 0.33) but did not reduce wrong calls, and it has not been tested together with the labels.
- **NUMTs are not detected directly.** ECMSD only shows that reads do not fit a reference the way reads fit their own species (`poor_fit`, window rule). This was not tested on a NUMT benchmark. Proposed: a benchmark with known NUMT-rich species.
- **Amplicon and capture data.** ECMSD expects shotgun data. For amplicon or capture data the breadth and evenness thresholds do not apply.
- **Low-quality references.** The database cannot be fully checked: a RefSeq mitogenome may itself be chimeric or contaminated. Proposed: flag or remove the two known entries above and screen the rest.
- **Reads shorter than about 30 bp** are rarely mapped, and at 30 bp a single mismatch is 3.3%, so identity columns are noisy. At 47 bp the decline is 1.000 for right and wrong references alike.

**Not tested**

- Samples with very few mitochondrial reads (about 10-50).
- Large mitogenomes of plants and fungi.

### Command Line Options

#### Required Arguments for Analysis

- `-f, --fwd`: Path to the forward or single-end FASTQ file (required unless `--merged` is used)
- `-m, --merged`: Path to a merged FASTQ file (required unless `--fwd` is used; cannot be combined with `--fwd`/`--rev`)
- `-o, --out`: Path to the output folder
- `-d, --db-folder`: Path to the database folder (must be created first with `--create-db`)

#### Required Arguments for Database Creation

- `-z, --create-db`: Create a new reference database
- `-d, --db-folder`: Path to the folder where the database will be stored
- `-b, --force-rebuild`: Optional. Discard an existing database and rebuild it from scratch (implies `--create-db`). Extra references are kept

#### Managing the Database

- `--update-db`: update the database to the current NCBI RefSeq release, if it is newer (see [Database version, updates and extensions](#database-version-updates-and-extensions))
- `--extend-db --accessions FILE`: add mitogenomes by NCBI nucleotide accession, one per line
- `--extend-db --fasta FASTA`: add mitogenomes from a FASTA file with headers `>taxid|id`

#### Optional Arguments

- `-r, --rev`: Path to the reverse FASTQ file (for paired-end data; requires `--fwd`)
- `-u, --cov-threshold`: Minimum percentage of a reference covered by reads to retain it (default: 25; see [Coverage-based Reference Filtering](#coverage-based-reference-filtering---cov-threshold) for low-depth samples)
- `-n, --top-n`: Number of top references to generate per-reference alignment plots for (default: 25)
- `-q, --mapping_quality`: Drop reads whose best alignments all have a lower MAPQ (default: 0, i.e. no MAPQ filtering; see [Mapping Quality](#mapping-quality))
- `--seed-k`: minimap2 k-mer size (default: 17)
- `--seed-w`: minimap2 minimizer window size (default: 7)
- `--min-chain-score`: minimap2 minimal chaining score `-m` (default: 20). Higher values align fewer weak reads and make mapping faster; see [Mapping parameters](#mapping-parameters)
- `--max-alns`: `minimap2 -N`, maximum number of secondary alignments per read (default: 20)
- `--min-disc-reads`: A reference is only reported when at least this many reads discriminate it from all other references (default: 20)
- `--min-disc-breadth`: ... and these reads cover at least this % of the reference (default: 5)
- `--min-disc-evenness`: ... and are spread over the reference with at least this evenness, 0-1 (default: 0.5)
- `--min-disc-window-ratio`: ... and hit at least this ratio of the 1 kb windows that as many reads placed at random would hit; 0 disables (default: 0.5)
- `--min-aln-len`: Minimum alignment length in bp of a read that counts as discriminatory evidence (default: 20). Lower it for very short reads, raise it to be stricter
- `--no-dedup`: Count PCR duplicates as separate reads in `Mito_reference_evidence.txt` (default: a read with the same references, start, end and strand as an earlier read is counted once, see [Identity](#identity-right-species-or-closest-relative))
- `--keep-all-alignments`: Also keep every alignment as `Mito_all_alignments.paf.gz` and `Mito_alignments.tsv.gz` (large, and adds several minutes on deep samples)
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
│   ├── Mito_calls.txt                           # Main result: one row per call with its status
│   ├── Mito.paf.gz                              # Best alignments per read
│   ├── Mito_all_alignments.paf.gz               # Only with --keep-all-alignments
│   ├── Mito_alignments.tsv.gz                   # Only with --keep-all-alignments
│   ├── Mito_reference_evidence.txt              # Compatible and discriminatory evidence per reference
│   ├── Mito_coverage.txt                        # Per-reference coverage statistics
│   ├── Mito_summary.txt                         # Taxonomic assignment per read
│   ├── Mito_summary.ref_summary.txt             # Ranked reference summary for plotting
│   ├── Mito_summary_<taxon>.txt                 # Per-taxon read count summary
│   ├── Mito_summary_<taxon>_ReadLengths.png     # Read length distribution plot
│   ├── Mito_summary_<taxon>_proportions.txt     # Proportions per taxon with call status
│   ├── Mito_summary_<taxon>_Proportions.png     # Proportions bar plot, coloured by call status
│   └── alignment_plots/                         # Per-reference coverage and depth plots
└── logs/
    └── minimap2.log                             # Mapping log
```

### Output File Descriptions

- **Mito_calls.txt**: the main result. The first line gives the verdict for the sample (`# sample: all calls confident`, `some calls are not confident, see status`, or that the reads fit no reference well and come from an unsampled species or haplotype of one genus); the second line gives the database version. Then one row per call: `call`, `status`, `proportion`, `weighted_reads`, `taxids`, `members` and `reason`. **Only `confident` rows are species calls.** The status is:

  - `confident`: the reads fit this reference as well as reads fit their own species, and they separate it from the other reported references;
  - `ambiguous`: references that share most of their reads (one row per group, e.g. `Drosophila group: D. mauritiana / D. simulans`); the reads cannot tell which of them is present;
  - `poor_fit`: the reads fit the reference worse than reads fit their own species, or only in one region: a relative of the source, not the source (the species or its haplotype is probably missing from the database, see `--extend-db`);
  - `low_evidence`: reported only through reads at the minimap2 `-N` limit.

  The same summary is printed at the end of the run.
- **Mito.paf.gz**: Compressed PAF with, for each read, the alignments to all references that tie for its best `AS` (the references the read is compatible with). Each line carries the extra tags `nc` (number of compatible references), `dm` (`AS` margin to the best lower-scoring reference, if one was reported) and `na` (number of alignments minimap2 reported for the read). Evidence, coverage, taxonomy and plots are based on this file
- **Mito_reference_evidence.txt**: One row per reference with read counts (`compatible_reads`, `discriminatory_reads`, `weak_unique_reads`, `shared_reads`, `discriminatory_at_limit`), breadth, depth, evenness of the compatible and of the discriminatory reads, spatial statistics of the discriminatory reads (`breadth_ratio`: observed breadth / breadth expected for reads placed at random, 1 - exp(-depth); `discriminatory_windows`: fraction of 1 kb windows with a discriminatory read; `window_ratio`: that fraction / the fraction expected for as many random reads; `max_1kb_share`: share of discriminatory reads in the densest 1 kb window), the identity of the discriminatory reads (`discriminatory_identity`, `discriminatory_exact`, `edit_distance_decline`, see [Identity](#identity-right-species-or-closest-relative)), the fraction of compatible reads skipped as PCR duplicates (`duplicate_fraction`), whether the reference is `credible` and, if not, the `reason`
- **Mito_all_alignments.paf.gz**, **Mito_alignments.tsv.gz** (only with `--keep-all-alignments`): all alignments, including secondary ones (`tp:A:S`). The table lists reference, coordinates, aligned length, MAPQ, `tp`, `AS`, `NM`, whether the alignment was chosen as best and its `AS` difference to the best alignment of the read, so ambiguous reads can be inspected
- **Mito_coverage.txt**: Per-reference coverage statistics. Columns: `reference`, `ref_length`, `mean_coverage`, `std_coverage`, `pct_covered`
- **Mito_summary.txt**: Tab-separated file with taxonomic assignments and read counts per read
- **Mito_summary.ref_summary.txt**: Ranked reference summary used for generating alignment plots. Columns: `ref_name`, `taxid`, `taxon_name`, `species_name`, `subspecies_name`, `read_count` (`subspecies_name` is `NA` for references not resolved below species level)
- **Mito_summary_\<taxon\>.txt**: Aggregated read counts per taxon at the chosen taxonomic level
- **Mito_summary_\<taxon\>_ReadLengths.png**: Faceted histogram of read length distributions for the top 10 taxa, with the call status in each facet title
- **Mito_summary_\<taxon\>_proportions.txt**: Reads and proportion per taxon at the chosen level, with the call `status` from `Mito_calls.txt`; a taxon of a coarser level whose calls have different statuses is `mixed`
- **Mito_summary_\<taxon\>_Proportions.png**: Bar plot of relative read proportions per taxon, coloured by call status, with the sample verdict as subtitle
- **alignment_plots/**: Two-panel PNG per reference showing alignment coverage breadth (top) and sequencing depth (bottom) for the top-N references; the subtitle starts with the call status. Titles and file names use the label for the rank requested with `-x`, so a run at species level never shows a subspecies name
- **minimap2.log**: Detailed minimap2 log

## Parameters Explanation

### Coverage-based Reference Filtering (`--cov-threshold`)

After mapping, ECMSD calculates per-reference breadth of coverage from the PAF file. For each reference, it computes the percentage of reference positions covered by at least one read (`pct_covered`). Only references where `pct_covered ≥ --cov-threshold` are passed to the taxonomic assignment step.

Coverage breadth is a more direct and interpretable measure for distinguishing authentic mitochondrial hits (which should show broad, even coverage across the genome) from mapping artifacts or NUMTs (which tend to pile up on a small region of the reference).

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

A reference is reported (`credible`) when it has at least `--min-disc-reads` discriminatory reads that cover at least `--min-disc-breadth` % of it and are spread over it with an evenness (normalised Shannon entropy over 20 equal windows) of at least `--min-disc-evenness`, and hit at least `--min-disc-window-ratio` (0.5) times the number of 1 kb windows that as many reads placed at random would hit (`window_ratio`). For 20 or more reads placed at random the ratio stays above 0.6, so the default only removes references whose discriminatory reads pile up in a few regions, such as a conserved region that is almost identical to the species present, or a stretch of foreign sequence in the reference. No threshold depends on how abundant other species are, so a rare species with well-spread discriminatory reads is kept. Reads shared between references still count towards the coverage of every credible reference they are compatible with, and a read that is compatible with several credible references is split between them in the proportion of their discriminatory reads: with 300 and 100 discriminatory reads, a shared read counts 0.75 and 0.25. `Mito_summary.txt` has one row per read and reference, with this share in the column `Weight`, and the species tables sum the weights. Only credible references are reported, always as the reference itself; nothing is reported at a higher taxonomic level. If a sample only contains shared reads, no reference is credible and nothing is reported.

Limits: minimap2 reports at most `--max-alns` (+1) alignments per read, chosen before base-level alignment, so a competing reference may not be reported for a read. A unique best alignment is therefore only provisional; `discriminatory_at_limit` in `Mito_reference_evidence.txt` counts discriminatory reads for which the limit was reached. Species that are credible only because of reads at the limit, or because of a few reads, should be treated with caution. The thresholds are defaults chosen on one test data set and may need adjusting.

### Identity: right species or closest relative

The database only holds species with a curated RefSeq mitogenome. When the species in a sample is not in it, its reads go to the closest relative that is. That relative then gets many discriminatory reads, spread over its whole genome, and is reported as `credible`, exactly like a correct hit. Counts, breadth and evenness cannot tell "this is species X" from "this is something close to X".

How well the reads match the reference can sometimes help. For every reference, `Mito_reference_evidence.txt` gives three columns, computed from the alignments of its discriminatory reads only:

- `discriminatory_identity`: pooled identity, 1 - (sum of `NM`) / (sum of alignment lengths). Unlike a median, it does not move in steps of one mismatch per read.
- `discriminatory_exact`: fraction of the discriminatory reads without any mismatch.
- `edit_distance_decline`: the edit distance decline (-Δ%) of HOPS (Hübler et al. 2019, doi:10.1186/s13059-019-1903-0), computed as in its post-processing script: over the numbers of reads with 0, 1, 2 and 3 mismatches, the sum of the decreases between neighbouring numbers divided by the sum of all changes. 1 means the numbers fall at every step, as expected for reads of the reference species; reads of a relative peak at 1 or 2 mismatches and give a lower value. HOPS accepts 0.9 or more. `NA` when there is no change (no reads with 0-3 mismatches).

Shared reads are left out because they come from conserved regions, which match distant references as well as the right one. For the identity columns, `NA` means the reference has no discriminatory read with an `NM` tag. The columns are reported only; they do not change which references are credible or how reads are assigned.

What the columns can and cannot tell:

- **They can show a poor match.** A credible reference with `edit_distance_decline` below 0.9 does not fit the reads well: the species in the sample, or its haplotype, is probably not in the database. In the test below this held for every relative of a missing species but one (next points). `discriminatory_identity` and `discriminatory_exact` show how poor the fit is; compare them only within a sample or with a sample of known species of the same read length, platform and processing.
- **Compare identities within a sample, never with a fixed number.** Sequencing errors, damage and the haplotype of the sample shift all references of a sample together: the *D. simulans* sample below matches its own species at 0.979, the *D. melanogaster* sample at 0.996.
- **A good decline is not proof.** A reference whose discriminatory reads all come from one region nearly identical to the sample species can still reach 1 (*D. formosana* below); check `discriminatory_breadth` as well. With few reads (about 10-20) the value is noisy.
- **They cannot tell close relatives apart.** Species that share haplotypes give the same identity. Equal identities on several references say nothing about which one is right.

Measured on 4M read pairs of *D. simulans* (76 bp) and *D. melanogaster* (75 bp), against the full database and against a database without the sample's own species. Credible references, `discriminatory_identity` / `discriminatory_exact` / `edit_distance_decline`:


| Sample                               | Database                 | Credible references                                                                                                                                           |
| ------------------------------------ | ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| *D. melanogaster*                    | full                     | *D. melanogaster* 0.9956 / 0.764 / 1.000                                                                                                                      |
| *D. melanogaster*                    | without*D. melanogaster* | *D. formosana* 0.9941 / 0.686 / 1.000, *D. mauritiana* 0.9774 / 0.274 / 0.667, *D. simulans* 0.9747 / 0.205 / 0.357                                           |
| *D. simulans*                        | full                     | *D. melanogaster* 0.9822 / 0.299 / 0.744, *D. mauritiana* 0.9813 / 0.288 / 0.738, *D. sechellia* 0.9802 / 0.261 / 0.813, *D. simulans* 0.9789 / 0.244 / 0.666 |
| *D. simulans*                        | without*D. simulans*     | *D. sechellia* 0.9805 / 0.260 / 0.758, *D. mauritiana* 0.9797 / 0.263 / 0.715, *D. melanogaster* 0.9786 / 0.208 / 0.631                                       |
| 3*D. melanogaster* : 1 *D. simulans* | full                     | *D. melanogaster* 0.9855 / 0.571 / 1.000, *D. mauritiana*, *D. sechellia*, *D. simulans* 0.969-0.971 / 0.19-0.21 / 0.58-0.69                                  |

*D. melanogaster* fits its own reference (decline 1.000); with it removed, its closest relatives fall to 0.36-0.67. *D. formosana* stays at 1.000 because its discriminatory reads come from one region (17% breadth) that is nearly identical to *D. melanogaster*. The *D. simulans* sample fits no reference (0.67-0.81, also with the full database), because its haplogroup is not in the database (see below); the same holds for the *D. simulans* part of the mix.

Why the expected value is below 1, and what affects it:

- **Variation within a species.** The reference is one individual. The haplotype of a sample can differ from it by several percent even for the right species: the *D. simulans* haplogroups siI and siII differ by about 3%, more than *D. simulans* siI differs from *D. sechellia* (about 1.3%) (Ballard 2000, doi:10.1007/s002390010067).
- **Sequencing errors and ancient DNA damage.** Both add mismatches to every read. C to T and G to A changes at read ends lower the identity of aDNA samples, so compare a sample with a sample of known species processed the same way, not with 1.
- **NUMTs.** Nuclear copies of mitochondrial DNA have changed since they were inserted and map with lower identity. They give a tail of low-identity reads even on the correct reference. Mitochondria have hundreds to thousands of copies per cell and a nuclear region two, and *Drosophila* has few NUMTs, so the identity is rarely affected; samples with very few mitochondrial reads are the exception.
- **Short reads.** At 30 bp one mismatch is 3.3%, so `discriminatory_exact` is higher for short reads at the same divergence. Compare samples of similar read length.
- **Divergent reads are not mapped.** Reads that differ too much from every reference do not pass the minimap2 seeding and chaining thresholds, and discriminatory reads need at least `--min-identity`. The identity therefore underestimates the divergence of a missing species, and a reference can look close from one conserved region alone; check `discriminatory_breadth` as well.

#### PCR duplicates

Copies of one DNA molecule add no evidence, but would count towards `--min-disc-reads` and repeat the same errors in the identity columns. By default, a read whose compatible alignments have the same references, start, end and strand as an earlier read is skipped in `Mito_reference_evidence.txt`; `duplicate_fraction` gives the share of compatible reads skipped per reference. Breadth is not affected. The read assignment and the summary proportions still count every read.

The key is exact for merged reads, whose start and end are the ends of the molecule. Unmerged mates of the same trimmed length are keyed on one end only, so at high depth unrelated reads collide: in the *D. simulans* test sample (about 130x on the mitogenome) 47-55% of the reads were skipped, against 2-3% at about 2x in the *D. melanogaster* sample. There, `duplicate_fraction` is an upper bound. In these tests no reference changed its `credible` status. Use `--no-dedup` to count every read, for example when duplicates were already removed.

#### Species that share mitochondrial haplotypes

*D. simulans* carries three divergent haplogroups (siI, siII, siIII), and *D. mauritiana* haplogroup maI is nearly identical to siIII after introgression (Ballard 2000, doi:10.1093/oxfordjournals.molbev.a026394). The database has one mitogenome per species (the *D. simulans* RefSeq is probably siI, a rare Pacific lineage), so the reference of the right species can match a sample worse than the reference of its sister species: in the test above, the *D. simulans* sample fits the *D. mauritiana* reference better than its own. Read counts and identity cannot separate such species. Treat credible hits within such a group as "one of the group", and use nuclear markers to decide between them. The same applies to other species complexes with shared or introgressed mitochondria.

### Mapping parameters

The `sr` preset requires two minimizers per chain (`-n 2`), a chaining score of 25 (`-m 25`) and a peak score of 40 (`-s 40`). Reads shorter than about 30 bp are almost never mapped with these. ECMSD uses `-n 1 -m 20 -s 20` together with `-k 17 -w 7` instead.

The chaining score `-m` is what controls how many reads are aligned. On the 1.09 M read test sample (Drosophila melanogaster and Gryllus bimaculatus reads against the full RefSeq mitochondrial database), `-m 10` aligned 85% of all reads, of which 14% could only align with a very low alignment score (AS < 40), which made mapping slow (16.6 M alignments). `-m 20` aligns 59% of the reads, none of them with AS < 40, and still keeps 99% of the *D. melanogaster* reads, 99% of the reads with AS >= 60 and 94% of the *G. bimaculatus* reads (98.7% of those with an alignment of AS >= 40). `-n` and `-s` have no effect once `-m` is 20 or higher. Use `--min-chain-score 18` to recover a few percent more of the weakest hits (about 97% of *G. bimaculatus*) at the cost of 5 percentage points more aligned reads, or 25 (the preset value) for a stricter screen that loses about 7% of the *G. bimaculatus* reads.

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
   - Ensure sufficient disk space (about 5 to 10 GB) and internet connectivity
   - If a previous database build was interrupted, re-run `--create-db` to complete it
3. **"Output directory already exists"**

   - Use the `-c` or `--force` flag to overwrite existing results
   - Or specify a different output directory
4. **Memory or disk space issues**

   - Reduce the number of threads with `-t`
   - Ensure sufficient disk space for the database build (about 5 to 10 GB peak; the finished database is well under 1 GB)
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
