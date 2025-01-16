# README

## Overview
The script processes NGS reads. It is primarily meant for Biocenter Oulu Sequencing Center ONT customers but can be easily modified for other purposes. The script has been tested on macOS (Sequoia) on both Intel and ARM systems, as well as on the CSC Puhti server.

The 12nt-barcodes are custom indexes used in BCO to index amplicons in ONT projects.

## Steps
The script executes the following steps:
1. Demultiplexing forward and reverse reads according to linked index combinations.
2. Reverse complementing the reverse reads.
3. Merging corresponding forward and reverse reads.
4. Trimming primers from the merged reads.
5. Generating a summary file with statistics for each sequence file.
6. Cleaning up intermediate directories.

## Requirements
- **Python**: Version 3.6 or higher.
- **Libraries**:
  - `os`
  - `subprocess`
  - `argparse`
  - `glob`
  - `gzip`
  - `logging`
  - `yaml`
- **External Tools**:
  - `cutadapt`: Version 3.0 or higher.
  - `seqkit`: Version 2.6.0 or higher (supports `AvgQual`).

## Installation
The script expects external tools to be available in the PATH variable.

To install the required Python libraries, you can use `pip`:
```bash
pip install pyyaml
```

To install the external tools, you can use `conda` or `pip`:
```bash
conda install -c bioconda cutadapt seqkit
```

## Usage
The script can be executed from the command line with the following arguments:
- `--config`: Path to the YAML configuration file (optional). Default is `config.yaml` in the same directory as the script.
- `--project`: Project type (e.g., 16S, ITS) (required).
- `--input`: Input sequence file (FASTQ.gz format) (required).
- `--output`: Output directory for processed files (required).
- `--log`: Log file to store the output (optional, default: `process.log`).

Example command:
```bash
./process_reads.py --project 16S --input reads.fastq.gz --output output_dir
```

## Configuration File
The configuration file should be in YAML format and contain the following parameters for each project type:
- `min_len`: Minimum length of reads.
- `max_len`: Maximum length of reads.
- `forward_primer`: Forward primer sequence.
- `reverse_primer`: Reverse primer sequence.
- `forward_barcodes`: Path to the file containing linked barcodes in the forward direction.
- `reverse_barcodes`: Path to the file containing linked barcodes in the reverse direction.
- `cores`: Number of CPU cores to use (optional, default: 1).
- `error_rate`: Allowed error rate for the primer trimming step (optional, default: 0.1).
