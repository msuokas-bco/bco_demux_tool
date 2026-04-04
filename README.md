# BCO Demux Tool

## Overview

`bco_demux.py` processes Oxford Nanopore amplicon reads sequenced at Biocenter Oulu Sequencing Center (BCO). It demultiplexes reads using a samplesheet that maps sample IDs to BCO's 10×12 custom barcode plate (F1–F10 × R1–R12), then generates reverse-complements from reverse-oriented reads, merges sample sequences together, and finally trims PCR primers from the merged reads, generating analysis-ready amplicon sequences.

The tool has currently been tested on macOS version 26.

---

## Pipeline steps

1. **Demultiplex** — A combined barcode FASTA (both strand orientations) is built from the samplesheet and passed to a single `cutadapt` run. Reads are routed to per-sample files in `demux_temp/`.
2. **Reverse complement** — Reads matched on the reverse strand (`_REV` files) are reverse-complemented with `seqkit` so all reads are arranged 5'→3'.
3. **Merge** — Forward and RC-corrected reads for each sample are concatenated into `merged_reads/merged_{sample_id}.fastq.gz`.
4. **Trim primers** — Samples meeting the minimum read count threshold have forward and reverse primers trimmed with `cutadapt`. Output goes to `trimmed_reads/{sample_id}_trimmed.fastq.gz`.
5. **Summarise** — `seqkit stats` is run on all trimmed files and selected columns are written to `trimmed_reads/summary.tsv`.
6. **Clean up** — Intermediate directories (`demux_temp/`, `merged_reads/`) are removed.

A `run_barcodes.tsv` file recording the barcode combinations used is written to the output directory at the end of step 1.

---

## Requirements

- **Python** 3.6 or higher
- **Python libraries**: `pyyaml` (all other imports are from the standard library)
- **External tools**:
  - `cutadapt` ≥ 3.0
  - `seqkit` ≥ 2.6.0 (required for `AvgQual` column in stats)

### Installation

```bash
# Python dependency
pip install pyyaml

# External tools (via conda/mamba recommended)
conda install -c bioconda cutadapt seqkit
```

Both tools must be available on `PATH`. The script checks this at startup and exits with a clear error if either is missing.

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
  cores: 4             # CPU threads for cutadapt (default: 1)
  error_rate: 0.1      # Allowed mismatch rate for primer trimming (default: 0.1)
  min_reads: 20        # Minimum reads per sample to proceed to primer trimming (default: 20)
```

> **Note:** `forward_primer` and `reverse_primer` are stored 5'→3' as on the primer. The reverse primer is automatically reverse-complemented before being passed to `cutadapt -a`.

---

## Output structure

```
output_dir/
├── run_barcodes.tsv              # Barcode combinations used in this run
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

## Acknowledgements

This tool relies on the following open-source software. If you use `bco_demux.py` in work leading to a publication, please cite these tools accordingly.

**cutadapt**
> Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1), 10–12. https://doi.org/10.14806/ej.17.1.200

**seqkit**
> Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLOS ONE*, 11(10), e0163962. https://doi.org/10.1371/journal.pone.0163962

---

## License

GPL-3.0 — see [LICENSE](LICENSE).
