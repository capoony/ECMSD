# ECMSD - Efficient Comprehensive Mitochondrial Sequence Detection Pipeline

## Overview

ECMSD (Efficient Comprehensive Mitochondrial Sequence Detection) is a bioinformatics pipeline designed for the detection and taxonomic classification of mitochondrial sequences from high-throughput sequencing data. The pipeline uses minimap2 for fast sequence alignment against a comprehensive mitochondrial reference database and provides quantitative analysis through Read Mapping Uniformity Score (RMUS) calculations.

## Features

- **Fast mitochondrial sequence detection** using minimap2 alignment
- **Comprehensive reference database** built from NCBI RefSeq mitochondrial genomes
- **Taxonomic classification** with customizable hierarchy levels
- **Quality-based filtering** with configurable mapping quality thresholds
- **RMUS (Read Mapping Uniformity Score)** calculation for quantitative analysis
- **Automated visualization** of results with R-based plotting
- **Support for multiple input formats**: single-end, paired-end, and merged reads
- **Conda environment management** for reproducible results

## Pipeline Workflow

1. **Reference Database Construction**: Downloads and processes NCBI RefSeq mitochondrial genomes
2. **Sequence Mapping**: Aligns input reads to mitochondrial reference using minimap2
3. **Quality Filtering**: Filters alignments based on mapping quality threshold
4. **Taxonomic Assignment**: Links aligned sequences to taxonomic information using NCBI taxonomy
5. **RMUS Calculation**: Computes Read Mapping Uniformity Score for quantitative assessment
6. **Visualization**: Generates plots and summary statistics

## Requirements

- Linux/Unix operating system
- Conda package manager
- Internet connection (for initial database download)

### Software Dependencies (automatically installed)

- minimap2
- BBTools (bbmask)
- Python 3 with numpy
- R with tidyverse
- pigz
- wget

## Installation

1. Clone or download the ECMSD pipeline:

```bash
git clone <repository-url>
cd ECMSD
```

2. The pipeline will automatically install dependencies and build the reference database on first run.

## Test Data

ECMSD includes a test dataset to verify the installation and demonstrate pipeline functionality. The test data consists of historic DNA (hDNA) sequencing reads from an unknown origin, providing a realistic example for mitochondrial sequence detection.

### Test Dataset Details

- **File**: `TestData/merged.fastq.gz` (~37 MB)
- **Content**: Merged FASTQ reads from hDNA sequencing
- **Purpose**: Pipeline validation and demonstration
- **Origin**: Unknown sample (for testing purposes only)

### Running the Test

To test the pipeline with the provided dataset:

```bash
# Navigate to ECMSD directory
cd /path/to/ECMSD

# Replace <path_to_your_ECMSD_directory> with your actual ECMSD installation path
bash shell/ECMSD.sh \
    --fwd TestData/merged.fastq.gz \
    --out TestOutput \
    --threads 20 \
    --Binsize 1000 \
    --RMUS-threshold 0.15 \
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

**Note**: Before running the test script, make sure to edit the `WD` variable in `TestData/shell.sh` to point to your actual ECMSD installation directory.

### Test Data Validation

A successful test run indicates that:

- All dependencies are properly installed
- The reference database is correctly built
- The pipeline can process FASTQ data
- Output files are generated in the expected format

## Usage

### Basic Usage

```bash
# Single-end reads
./shell/ECMSD.sh -f reads.fastq -o output_directory/

# Paired-end reads
./shell/ECMSD.sh -f reads_R1.fastq -r reads_R2.fastq -o output_directory/

# With merged reads
./shell/ECMSD.sh -f reads_R1.fastq -r reads_R2.fastq -m merged_reads.fastq -o output_directory/
```

### Command Line Options

#### Required Arguments

- `-f, --fwd`: Path to the forward or single end (e.g. merged) FASTQ file
- `-o, --out`: Path to the output folder

#### Optional Arguments

- `-r, --rev`: Path to the reverse FASTQ file (for paired-end data)
- `-m, --merged`: Path to merged FASTQ file that should be mapped in addition to the forward and reverse reads
- `-b, --Binsize`: Bin size for RMUS analysis (default: 1000)
- `-u, --RMUS-threshold`: RMUS threshold for analysis (default: 0.15)
- `-q, --mapping_quality`: Mapping quality threshold (default: 20)
- `-t, --threads`: Number of threads to use (default: 10)
- `-c, --force`: Force overwrite of existing output files
- `-x, --taxonomic-hierarchy`: Taxonomic hierarchy level (default: species)
- `-v, --version`: Show version and exit
- `-h, --help`: Show help message and exit

### Advanced Examples

```bash
# High-stringency analysis with custom thresholds
./shell/ECMSD.sh -f reads_R1.fastq -r reads_R2.fastq -o results/ \
  -q 30 -u 0.25 -b 500 -t 16

# Genus-level taxonomic classification
./shell/ECMSD.sh -f reads.fastq -o results/ -x genus

# Force overwrite existing results
./shell/ECMSD.sh -f reads.fastq -o existing_results/ -c
```

## Output Files

The pipeline generates the following output structure:

```text
output_directory/
├── mapping/
│   ├── Mito.paf.gz           # Compressed alignment results
│   ├── Mito_summary.txt      # Taxonomic summary with RMUS scores
│   └── Mito_summary_plot.*   # Visualization plots
└── logs/
    └── logfile.log           # Pipeline execution logs
```

### Output File Descriptions

- **Mito.paf.gz**: Compressed PAF format alignment file containing all mitochondrial sequence alignments
- **Mito_summary.txt**: Tab-separated file with taxonomic assignments, read counts, and RMUS scores
- **Mito_summary_plot.***: Various plot formats showing taxonomic composition and abundance
- **logfile.log**: Detailed log of pipeline execution including any errors or warnings

## Parameters Explanation

### RMUS (Read Mapping Uniformity Score)

RMUS is a novel metric based on [Shannon's entropy](https://arxiv.org/pdf/1405.2061) that assesses the uniformity of read distribution across genomes. Higher RMUS values indicate more uniform coverage, which can help distinguish authentic mitochondrial content from mapping artifacts, contamination or overall low-quality alignments.

RMUS is calculated as follows:

```math
RMUS = \frac{H}{H_{max}}
```

where:

- \( H \) is the Shannon entropy of the read distribution across taxonomic categories

```math
H = -\sum_{i=1}^{n} p_i \log_2(p_i)
```

- \( p_i \) is the proportion of reads assigned to category \( i \)
- \( n \) is the total number of categories
  
- \( H_{max} \) is the maximum possible entropy for the given number of categories  
- If \( H_{max} \) is zero, RMUS is set to zero to avoid division by zero errors.

RMUS values range from 0 to 1, where 1 indicates perfect uniformity of reads across all categories (i.e. genomic bins). This metric is particularly useful for assessing if the mapping of reads is an artifact (mapping only to a single specific region), as it provides a quantitative measure of how evenly reads are distributed across the reference genomes. We assume that only a high RMUS value indicates that the original reads are likely from a specific mitochondrial genome, while low RMUS values may suggest mapping artifacts or numts. It helps in identifying potential contamination or uneven coverage that may affect downstream analyses. RMUS is calculated for each taxonomic category during the taxonomic assignment step and is included in the output summary file for easy interpretation.

### Bin Size

The bin size parameter determines the resolution for RMUS calculation. Smaller bins provide higher resolution but may be more sensitive to sequencing artifacts.

### Mapping Quality

The mapping quality threshold filters alignments based on their reliability. Higher values increase stringency but may reduce sensitivity for divergent sequences.

### Taxonomic Hierarchy

Determines the taxonomic level for classification (e.g., species, genus, family, order, phylum or kingdom). Lower levels provide more specific identification but may have reduced sensitivity.

## Troubleshooting

### Common Issues

1. **"Conda environment could not be activated"**
   - Ensure conda is properly installed and initialized
   - Run `conda init bash` and restart your shell

2. **"Reference data not found"**
   - The pipeline will automatically download reference data on first run
   - Ensure internet connectivity and sufficient disk space

3. **"Output directory already exists"**
   - Use the `-c` or `--force` flag to overwrite existing results
   - Or specify a different output directory

4. **Memory or disk space issues**
   - Reduce the number of threads with `-t`
   - Ensure sufficient disk space for reference database (~5-10 GB)

### Performance Optimization

- Use more threads (`-t`) on systems with multiple CPU cores
- Adjust bin size (`-b`) based on your data characteristics
- Consider mapping quality threshold (`-q`) based on your sequencing platform

## Citation

If you use ECMSD in your research, please cite:

[Citation information to be added]

## License

MIT License (MIT)

## Author

Martin Kapun

## Contact

For questions, issues, or feature requests, please write an issue.

## Version History

- **v1.0**: Initial release with core mitochondrial detection functionality

## Acknowledgments

This pipeline uses several open-source tools:

- minimap2 for sequence alignment
- BBTools for sequence masking
- NCBI RefSeq database for reference genomes
- NCBI taxonomy database for taxonomic classification
