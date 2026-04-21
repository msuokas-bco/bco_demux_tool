# BCO Demux Tool

## Overview

`bco_demux.py` processes Oxford Nanopore amplicon reads sequenced at Biocenter Oulu Sequencing Center (BCO). It demultiplexes reads using a samplesheet that maps sample IDs to BCO's 10×12 custom barcode plate (F1–F10 × R1–R12), extracts Oxford Nanopore DNA CS control reads for QC (if they have been included in the sequencing library), generates reverse-complements from reverse-oriented reads, merges sample sequences together, and trims PCR primers from the merged reads, generating analysis-ready amplicon sequences.

The tool has currently been tested on macOS version 26 and Ubuntu 22.04 LTS.

---

## Pipeline steps

1. **Demultiplex** — A combined barcode FASTA (both strand orientations) is built from the samplesheet and passed to a single `cutadapt` run. Reads are routed to per-sample files in `demux_temp/`.
2. **Write barcode record** — A `run_barcodes.tsv` summarising the barcode combinations used is written to the output directory.
3. **DNA CS extraction** *(optional)* — If `minimap2` and `samtools` are available, all input reads are aligned to the ONT DNA CS reference sequence (3,584 bp). Mapped reads are saved to `dna_cs/dna_cs.fastq.gz` for use with nanoporeQC. This step runs on the original input because DNA CS reads carry no sample barcodes and are discarded by the demux step. Skipped automatically if either tool is missing; if the dataset contains no DNA CS reads the output file will be empty.
4. **Reverse complement** — Reads matched on the reverse strand (`_REV` files) are reverse-complemented with `seqkit` so all reads are arranged 5'→3'.
5. **Merge** — Forward and RC-corrected reads for each sample are concatenated into `merged_reads/merged_{sample_id}.fastq.gz`.
6. **Trim primers** — PCR primers are trimmed with `cutadapt` (`--trimmed-only`). Samples with fewer reads than `min_reads` after trimming are removed and excluded from the summary.
7. **Summarise** — `seqkit stats` is run on all trimmed files and selected columns are written to `trimmed_reads/summary.tsv`.
8. **Clean up** — Intermediate directories (`demux_temp/`, `merged_reads/`) are removed.

---

## Requirements

- **Python** 3.8 or higher
- **Python libraries**: `pyyaml` (all other imports are from the standard library)
- **External tools (required)**:
  - `cutadapt` ≥ 3.0
  - `seqkit` ≥ 2.6.0 (required for `AvgQual` column in stats)
- **External tools (optional)**:
  - `minimap2` — required for DNA CS extraction
  - `samtools` — required for DNA CS extraction

If `minimap2` or `samtools` is not found, the DNA CS extraction step is skipped automatically and the rest of the pipeline runs unchanged.

### Installation

```bash
# Python dependency
pip install pyyaml

# Required external tools (via conda/mamba recommended)
conda install -c bioconda cutadapt seqkit

# Optional: enable DNA CS extraction
conda install -c bioconda minimap2 samtools
```

Required tools must be available on `PATH`. The script checks this at startup and exits with a clear error if either is missing. Optional tools are checked separately; a warning is logged if they are absent.

---

## Samplesheet format

The samplesheet is a **TSV** (`.tsv`) or **CSV** (`.csv`) file with three required columns:

| SampleID | ForwardBarcode | ReverseBarcode |
|----------|---------------|----------------|
| Sample_A | F1            | R1             |
| Sample_B | F3            | R7             |
| Sample_C | F10           | R12            |

- `ForwardBarcode`: one of `F1`–`F10` (row barcodes)
- `ReverseBarcode`: one of `R1`–`R12` (column barcodes)
- Each SampleID and each barcode combination must be unique
- Characters unsafe for filenames (anything other than letters, digits, `_`, `-`) are replaced with `_` in output file names

Additional columns beyond the three required ones are permitted and will be carried through to `run_barcodes.tsv` unchanged.

The 10×12 barcode plate supports up to **120 samples per run**.

---

## Usage

```
python bco_demux.py --project <TYPE> --samplesheet <FILE> --input <READS.fastq.gz> --output <DIR> [options]
```

### Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--project` | Yes | — | Project type key in the config file (e.g. `16S`, `ITS`) |
| `--samplesheet` | Yes | — | Samplesheet file (`.tsv` or `.csv`) |
| `--input` | Yes | — | Input reads in FASTQ.gz format |
| `--output` | Yes | — | Output directory (created if absent) |
| `--config` | No | `config.yaml` next to script | Path to YAML configuration file |
| `--log` | No | `process.log` | Path to log file |
| `--dry-run` | No | off | Log all commands without executing them |
| `--version` | No | — | Print version and exit |

### Example

```bash
python bco_demux.py \
  --project 16S \
  --samplesheet my_plate.tsv \
  --input run_reads.fastq.gz \
  --output results/
```

Dry-run (validate samplesheet and config without processing any data):

```bash
python bco_demux.py \
  --project 16S \
  --samplesheet my_plate.tsv \
  --input run_reads.fastq.gz \
  --output results/ \
  --dry-run
```

---

## Configuration file

The configuration file (`config.yaml`) contains one block per project type. All parameters except `cores`, `error_rate`, and `min_reads` are required.

```yaml
16S:
  forward_primer: "RGTTYGATYMTGGCTCAG"   # 5'→3' sequence
  reverse_primer: "RGYTACCTTGTTACGACTT"  # 5'→3' sequence (RC computed automatically)
  min_len: 1200        # Minimum read length (barcode demux and primer trimming)
  max_len: 1700        # Maximum read length (barcode demux only)
  cores: 4             # CPU threads for cutadapt and minimap2 (default: 1)
  error_rate: 0.1      # Allowed mismatch rate for primer trimming (default: 0.1)
  min_reads: 20        # Minimum reads per sample after primer trimming (default: 20)
```

> **Note:** `forward_primer` and `reverse_primer` are stored 5'→3' as on the primer. The reverse primer is automatically reverse-complemented before being passed to `cutadapt -a`.

---

## Output structure

```
output_dir/
├── run_barcodes.tsv              # Barcode combinations used in this run
├── dna_cs/                       # Only present when minimap2 and samtools are available
│   └── dna_cs.fastq.gz           # ONT DNA CS control reads (for nanoporeQC)
└── trimmed_reads/
    ├── summary.tsv               # Per-sample stats (num_seqs, lengths, Q20, AvgQual)
    ├── SampleA_trimmed.fastq.gz
    ├── SampleB_trimmed.fastq.gz
    └── ...
```

Intermediate directories (`demux_temp/`, `merged_reads/`) are removed automatically on successful completion.

---

## Barcode reference

BCO uses a plate of **10 forward** (row) and **12 reverse** (column) 12-nt barcodes, giving 120 unique combinations. The sequences are hardcoded in `bco_demux.py` and do not require external FASTA files.

| Label | Sequence (5'→3') | | Label | Sequence (5'→3') |
|-------|------------------|-|-------|------------------|
| F1  | TACAACACTCGT | | R1  | GTCCACGAGTGA |
| F2  | ACGCGTACCATA | | R2  | TGTCTGCTGTGA |
| F3  | GGTCACCTCCAT | | R3  | ACGCTCTTGTGA |
| F4  | GGAGAAGAAGAA | | R4  | CTGACATTCTGA |
| F5  | ACCACAGAAGAA | | R5  | TGAACCGATTGA |
| F6  | CTCTAGGAAGAA | | R6  | TTCGTGAGTTGA |
| F7  | CATGCCGAAGAA | | R7  | GATTGCAGTTGA |
| F8  | CCGTTACAAGAA | | R8  | AAGAAGCCTTGA |
| F9  | AACCTGTTCAGA | | R9  | TCCTGCGATTGG |
| F10 | TTGTCCTTAGGA | | R10 | GACTCTGGTTGG |
|     |              | | R11 | TAAGTGTCTTGG |
|     |              | | R12 | AACCTAACTGCG |

---

## Troubleshooting

### "Too many open files" during demultiplexing

`cutadapt` holds all output files open simultaneously during demultiplexing — one file per sample per strand orientation. With a full 120-sample run this requires approximately 240 file descriptors plus overhead.

The script attempts to raise the file descriptor limit automatically at startup. On most Linux systems this succeeds without any manual intervention. On macOS or in restricted environments (e.g. HPC cluster login nodes) the system hard limit may be too low for the automatic raise to succeed, in which case the log will include a warning.

**Fix:** increase the limit in your shell before running the script:

```bash
ulimit -n 8192
python bco_demux.py ...
```

If you do not have permission to raise the limit (common on shared HPC systems), contact your system administrator and ask them to increase `nofile` for your account, or run the script in a job submission environment where higher limits are typically set by default.

> **Technical note:** This is an architectural property of how `cutadapt` handles demultiplexing into named output files, not a limitation of `bco_demux.py` itself. The current 10×12 barcode matrix (max 120 samples × 2 orientations = 240 open files) fits within default limits on most modern Linux systems. If the barcode matrix were expanded substantially, a two-stage demultiplexing strategy — first by forward barcode (10 files), then by reverse barcode within each group (12 files) — would eliminate the file descriptor concern entirely regardless of matrix size.

---

## Acknowledgements

This tool relies on the following open-source software. If you use `bco_demux.py` in work leading to a publication, please cite these tools accordingly.

**cutadapt**
> Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1), 10–12. https://doi.org/10.14806/ej.17.1.200

**seqkit**
> Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLOS ONE*, 11(10), e0163962. https://doi.org/10.1371/journal.pone.0163962

**minimap2**
> Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

**samtools**
> Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

---

## License

GPL-3.0 — see [LICENSE](LICENSE).
