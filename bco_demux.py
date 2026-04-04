#!/usr/bin/env python3

import os
import re
import csv
import subprocess
import argparse
import glob
import gzip
import logging
import yaml
import sys
import shutil

# Version information
__version__ = "0.6.1"

# Barcode Definitions
# 10 Forward Barcodes: F1–F10 (rows)
FORWARD_BARCODES = {
    "F1":  "TACAACACTCGT", "F2":  "ACGCGTACCATA", "F3":  "GGTCACCTCCAT",
    "F4":  "GGAGAAGAAGAA", "F5":  "ACCACAGAAGAA", "F6":  "CTCTAGGAAGAA",
    "F7":  "CATGCCGAAGAA", "F8":  "CCGTTACAAGAA", "F9":  "AACCTGTTCAGA",
    "F10": "TTGTCCTTAGGA"
}

# 12 Reverse Barcodes: R1–R12 (columns) - stored as primer sequences (5'-3')
REVERSE_BARCODES = {
    "R1":  "GTCCACGAGTGA", "R2":  "TGTCTGCTGTGA", "R3":  "ACGCTCTTGTGA",
    "R4":  "CTGACATTCTGA", "R5":  "TGAACCGATTGA", "R6":  "TTCGTGAGTTGA",
    "R7":  "GATTGCAGTTGA", "R8":  "AAGAAGCCTTGA", "R9":  "TCCTGCGATTGG",
    "R10": "GACTCTGGTTGG", "R11": "TAAGTGTCTTGG", "R12": "AACCTAACTGCG"
}

# The '...' spacer is required by cutadapt to restrict barcode matching to the
# ends of reads (non-internal matching). Do not modify.
BARCODE_SPACER = "..."

# Characters not safe for use in filenames
_UNSAFE_CHARS = re.compile(r"[^\w\-]")  # keep alphanumerics, underscore, hyphen


def reverse_complement(seq):
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M', 'S': 'S', 'W': 'W',
        'H': 'D', 'D': 'H', 'B': 'V', 'V': 'B', 'N': 'N',
        '.': '.'
    }
    return "".join(complement[base] for base in reversed(seq))


def sanitize_sample_id(sample_id):
    """Replace characters unsafe for filenames with underscores."""
    sanitized = _UNSAFE_CHARS.sub("_", sample_id)
    if sanitized != sample_id:
        logging.warning(
            f"SampleID '{sample_id}' contains unsafe characters — "
            f"renamed to '{sanitized}' for output files."
        )
    return sanitized


def load_samplesheet(path):
    """Parse a samplesheet (TSV or CSV, auto-detected by extension).

    Expected columns: SampleID, ForwardBarcode, ReverseBarcode
    ForwardBarcode values must be F1–F10; ReverseBarcode values must be R1–R12.

    Returns a list of dicts:
        [{'sample_id': str, 'fwd_label': str, 'rev_label': str,
          'fwd_seq': str, 'rev_seq': str}, ...]
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".tsv":
        delimiter = "\t"
    elif ext == ".csv":
        delimiter = ","
    else:
        raise ValueError(
            f"Unrecognised samplesheet extension '{ext}'. "
            "Use .tsv (tab-separated) or .csv (comma-separated)."
        )

    required_cols = {"SampleID", "ForwardBarcode", "ReverseBarcode"}
    samples = []
    seen_ids = set()
    seen_combos = set()

    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)

        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Samplesheet is missing required column(s): {', '.join(sorted(missing))}"
            )

        # Preserve any extra columns (beyond the three required ones) so they
        # can be passed through to run_barcodes.tsv unchanged.
        extra_cols = [f for f in reader.fieldnames if f not in required_cols]

        for lineno, row in enumerate(reader, start=2):  # 2 = first data row
            raw_id = row["SampleID"].strip()
            fwd_label = row["ForwardBarcode"].strip().upper()
            rev_label = row["ReverseBarcode"].strip().upper()

            if not raw_id:
                raise ValueError(f"Samplesheet line {lineno}: SampleID is empty.")

            sample_id = sanitize_sample_id(raw_id)

            if sample_id in seen_ids:
                raise ValueError(
                    f"Samplesheet line {lineno}: duplicate SampleID '{sample_id}'."
                )
            seen_ids.add(sample_id)

            if fwd_label not in FORWARD_BARCODES:
                raise ValueError(
                    f"Samplesheet line {lineno} ('{sample_id}'): "
                    f"unknown ForwardBarcode '{fwd_label}'. "
                    f"Valid values: {', '.join(sorted(FORWARD_BARCODES))}."
                )
            if rev_label not in REVERSE_BARCODES:
                raise ValueError(
                    f"Samplesheet line {lineno} ('{sample_id}'): "
                    f"unknown ReverseBarcode '{rev_label}'. "
                    f"Valid values: {', '.join(sorted(REVERSE_BARCODES))}."
                )

            combo = (fwd_label, rev_label)
            if combo in seen_combos:
                raise ValueError(
                    f"Samplesheet line {lineno} ('{sample_id}'): "
                    f"duplicate barcode combination {fwd_label}/{rev_label}."
                )
            seen_combos.add(combo)

            samples.append({
                "sample_id": sample_id,
                "fwd_label":  fwd_label,
                "rev_label":  rev_label,
                "fwd_seq":    FORWARD_BARCODES[fwd_label],
                "rev_seq":    REVERSE_BARCODES[rev_label],
                "extra":      {k: row[k] for k in extra_cols},
            })

    if not samples:
        raise ValueError("Samplesheet contains no data rows.")

    logging.info(f"Loaded {len(samples)} sample(s) from samplesheet: {path}")
    return samples


def validate_input_file(file_path):
    """Check that the input file exists and is a readable gzip file."""
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")
    try:
        with gzip.open(file_path, 'rt') as f:
            f.read(1024)
    except gzip.BadGzipFile:
        raise ValueError(f"Input file is not a valid gzip file: {file_path}")
    except Exception as e:
        raise RuntimeError(f"Could not read input file {file_path}: {e}")


def count_sequences(file_path):
    """Count the number of sequences in a gzipped FASTQ file."""
    count = 0
    try:
        with gzip.open(file_path, 'rt') as f:
            for i, _ in enumerate(f):
                if i % 4 == 0:  # FASTQ records are 4 lines: header, seq, '+', quality
                    count += 1
    except gzip.BadGzipFile:
        raise ValueError(f"File is not a valid gzip file: {file_path}")
    except Exception as e:
        raise RuntimeError(f"Failed to count sequences in {file_path}: {e}")
    return count


def check_dependencies():
    """Verify that required external tools are available."""
    required_tools = ["cutadapt", "seqkit"]
    missing_tools = [t for t in required_tools if shutil.which(t) is None]
    if missing_tools:
        raise EnvironmentError(
            f"Missing required tools: {', '.join(missing_tools)}. "
            "Please install them and ensure they are in your PATH."
        )


def run_command(cmd, dry_run=False):
    """Run a command (as an argument list), logging stderr on success. Skips execution in dry-run mode."""
    logging.info(f"Running command: {' '.join(cmd)}")
    if dry_run:
        logging.info("[DRY RUN] Command not executed.")
        return
    try:
        result = subprocess.run(
            cmd, shell=False, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if result.stderr:
            logging.info(f"STDERR: {result.stderr.strip()}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {' '.join(cmd)}")
        logging.error(f"Return code: {e.returncode}")
        if e.stdout:
            logging.error(f"STDOUT: {e.stdout.strip()}")
        if e.stderr:
            logging.error(f"STDERR: {e.stderr.strip()}")
        raise


def merge_files(files_to_merge, output_file, dry_run=False):
    """Concatenate gzip files using Python file I/O (avoids shell redirection issues)."""
    logging.info(f"Merging {len(files_to_merge)} file(s) into {output_file}")
    if dry_run:
        logging.info("[DRY RUN] Merge not executed.")
        return
    with open(output_file, 'wb') as out:
        for f in files_to_merge:
            with open(f, 'rb') as inp:
                shutil.copyfileobj(inp, out)


# Step 1: Demultiplex reads using samplesheet-derived barcodes
def demux_all_reads(samples, reads_file, demux_temp_dir, output_dir,
                    min_len, max_len, cores, dry_run=False):
    """Builds a barcode FASTA from samplesheet entries only (both orientations)
    and runs a single cutadapt pass to demultiplex reads by SampleID."""
    combined_barcode_temp = os.path.join(output_dir, "generated_barcodes.fasta")
    logging.info(f"Generating barcode FASTA for {len(samples)} sample(s): {combined_barcode_temp}")

    try:
        with open(combined_barcode_temp, 'w') as outfile:
            for s in samples:
                # Reverse barcodes are stored 5'→3' as on the primer, but they
                # appear at the 3' end of a forward-oriented read as their RC.
                r_seq_rc = reverse_complement(s["rev_seq"])
                # Forward-oriented read: fwd_barcode...RC(rev_barcode)
                seq_fwd = f"{s['fwd_seq']}{BARCODE_SPACER}{r_seq_rc}"
                # Reverse-oriented read (opposite strand): full RC of seq_fwd.
                # Both entries are required because nanopore reads can arrive in
                # either orientation; cutadapt routes each match to a named output
                # file using the FASTA entry name ({name} in the -o pattern).
                seq_rev = reverse_complement(seq_fwd)
                outfile.write(f">{s['sample_id']}\n{seq_fwd}\n")
                outfile.write(f">{s['sample_id']}_REV\n{seq_rev}\n")

        # -e 0: zero mismatches required (exact barcode matching)
        # -O 24: minimum overlap = full barcode length (12 fwd + 12 rev)
        # -g: anchored 5' adapter (barcode at read start, spacer handles 3' end)
        # --discard-untrimmed: discard reads that did not match any barcode
        # {name}: cutadapt placeholder replaced with the matched FASTA entry name
        cmd = [
            "cutadapt",
            "-e", "0",
            "-O", "24",
            "-g", f"file:{combined_barcode_temp}",
            "--discard-untrimmed",
            "-m", str(min_len),
            "-M", str(max_len),
            "-j", str(cores),
            "-o", f"{demux_temp_dir}/{{name}}.fastq.gz",
            reads_file,
        ]
        run_command(cmd, dry_run=dry_run)

    finally:
        if os.path.exists(combined_barcode_temp):
            os.remove(combined_barcode_temp)


# Step 2: Write a record of which barcodes were used in this run
def write_run_barcodes(samples, output_dir):
    """Saves a TSV summarising the barcode combinations used in this run."""
    barcode_list_file = os.path.join(output_dir, "run_barcodes.tsv")
    logging.info(f"Writing run barcode summary: {barcode_list_file}")
    # Derive extra column names from the first sample (all samples share the same keys).
    extra_cols = list(samples[0]["extra"].keys()) if samples else []

    with open(barcode_list_file, 'w') as fh:
        # Column order: SampleID | [extra cols] | ForwardLabel | ForwardSequence | ReverseLabel | ReverseSequence
        # Barcode columns are kept last for consistent readability across runs.
        extra_header = ("\t" + "\t".join(extra_cols)) if extra_cols else ""
        fh.write("SampleID" + extra_header + "\tForwardLabel\tForwardSequence\tReverseLabel\tReverseSequence\n")
        for s in samples:
            extra_vals = ("\t" + "\t".join(s["extra"].get(k, "") for k in extra_cols)) if extra_cols else ""
            fh.write(f"{s['sample_id']}{extra_vals}\t{s['fwd_label']}\t{s['fwd_seq']}\t{s['rev_label']}\t{s['rev_seq']}\n")


# Step 3: Reverse complement _REV reads and merge with forward-oriented reads
def process_and_merge_reads(samples, demux_temp_dir, merged_dir, dry_run=False):
    """For each sample, reverse-complements any _REV-matched reads and merges
    them with their forward-oriented counterparts."""
    for s in samples:
        sample_id = s["sample_id"]
        fwd_file = os.path.join(demux_temp_dir, f"{sample_id}.fastq.gz")
        rev_file = os.path.join(demux_temp_dir, f"{sample_id}_REV.fastq.gz")
        merged_file = os.path.join(merged_dir, f"merged_{sample_id}.fastq.gz")
        rev_comp_file = os.path.join(demux_temp_dir, f"{sample_id}_REV_RC.fastq.gz")

        files_to_merge = []
        try:
            if os.path.exists(fwd_file):
                files_to_merge.append(fwd_file)
            if os.path.exists(rev_file):
                # -r: reverse, -p: complement → together = reverse complement.
                # Brings _REV-matched reads into the same 5'→3' orientation as
                # forward-matched reads before merging.
                cmd = ["seqkit", "seq", "-rp", "--seq-type", "DNA", "-o", rev_comp_file, rev_file]
                run_command(cmd, dry_run=dry_run)
                if not dry_run:
                    files_to_merge.append(rev_comp_file)

            if files_to_merge:
                merge_files(files_to_merge, merged_file, dry_run=dry_run)
            elif not dry_run:
                logging.warning(
                    f"No demultiplexed files found for sample '{sample_id}'. "
                    "It may have had zero reads assigned."
                )
        finally:
            if os.path.exists(rev_comp_file):
                os.remove(rev_comp_file)


# Step 4: Trim forward and reverse primers from merged reads
def trim_primers(merged_dir, trimmed_dir, forward_primer, reverse_primer_rc,
                 min_len, cores, erate, min_reads, dry_run=False):
    """Filters merged read files by a minimum read count, then trims primers."""
    merged_files = glob.glob(f"{merged_dir}/merged_*.fastq.gz")
    for merged_file in merged_files:
        num_sequences = count_sequences(merged_file)
        if num_sequences >= min_reads:
            base_name = os.path.basename(merged_file)
            sample_id = base_name.replace("merged_", "").replace(".fastq.gz", "")
            output_file = os.path.join(trimmed_dir, f"{sample_id}_trimmed.fastq.gz")
            # -g: trim forward primer from the 5' end of each read
            # -a: trim reverse primer (passed as its RC) from the 3' end
            # --trimmed-only: discard reads where neither primer was found
            cmd = [
                "cutadapt",
                "-g", forward_primer,
                "-a", reverse_primer_rc,
                "--trimmed-only",
                "-j", str(cores),
                "--error-rate", str(erate),
                "-m", str(min_len),
                "-o", output_file,
                merged_file,
            ]
            try:
                run_command(cmd, dry_run=dry_run)
            except subprocess.CalledProcessError as e:
                logging.error(f"Error during trimming {merged_file}: {e}")
        else:
            logging.warning(
                f"'{os.path.basename(merged_file)}' has only {num_sequences} reads "
                f"(threshold: {min_reads}). Skipping primer trimming."
            )


# Step 5: Generate summary file
def generate_summary(trimmed_dir, dry_run=False):
    """Runs 'seqkit stats' on all trimmed files and writes a summary TSV."""
    if dry_run:
        logging.info("[DRY RUN] Summary generation skipped.")
        return

    summary_file = os.path.join(trimmed_dir, "summary.tsv")
    trimmed_files = sorted(glob.glob(f"{trimmed_dir}/*_trimmed.fastq.gz"))

    if not trimmed_files:
        logging.warning("No trimmed files found — summary will be empty.")
        return

    with open(summary_file, 'w') as summary:
        for idx, trimmed_file in enumerate(trimmed_files):
            result = subprocess.run(
                # -a: compute all statistics (including Q20, AvgQual); -T: TSV output
                ["seqkit", "stats", "-aT", trimmed_file],
                shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                universal_newlines=True
            )
            if result.returncode != 0:
                logging.error(
                    f"Failed to generate stats for {trimmed_file}. Error: {result.stderr}"
                )
                continue

            lines = result.stdout.strip().split('\n')
            if len(lines) < 2:
                logging.warning(
                    f"Unexpected seqkit stats output for {trimmed_file}: {result.stdout}"
                )
                continue

            raw_header = lines[0].strip().split('\t')
            raw_stats = lines[1].strip().split('\t')

            desired_cols = ['file', 'num_seqs', 'min_len', 'avg_len', 'max_len', 'Q20(%)', 'AvgQual']
            indices = [i for i, col in enumerate(raw_header) if col in desired_cols]

            header = [raw_header[i] for i in indices]
            stats = [raw_stats[i] for i in indices]

            if idx == 0:
                summary.write("\t".join(header) + "\n")  # write header once, from the first file

            if 'file' in header:
                # seqkit stats writes the full path; replace with basename for readability
                stats[header.index('file')] = os.path.basename(trimmed_file)

            summary.write("\t".join(stats) + "\n")
            logging.info(f"Stats written for {os.path.basename(trimmed_file)}.")


def cleanup_intermediate_dirs(demux_temp_dir, merged_dir):
    """Removes temporary directories created during processing."""
    for directory in [demux_temp_dir, merged_dir]:
        if os.path.exists(directory):
            logging.info(f"Removing intermediate directory: {directory}")
            shutil.rmtree(directory)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Demultiplex amplicon reads using a samplesheet, then reverse complement, "
            "merge, and trim primers from matched reads."
        )
    )
    parser.add_argument('--version', action='store_true',
                        help="Show version information and exit.")
    parser.add_argument('--config', type=str,
                        help="Path to the YAML configuration file "
                             "(default: config.yaml alongside this script).")
    parser.add_argument('--project', type=str, required=True,
                        help="Project type matching a key in the config file (e.g. 16S, ITS).")
    parser.add_argument('--samplesheet', type=str, required=True,
                        help="Samplesheet file (.tsv or .csv) with columns: "
                             "SampleID, ForwardBarcode (F1–F10), ReverseBarcode (R1–R12).")
    parser.add_argument('--input', type=str, required=True,
                        help="Input reads file in FASTQ.gz format.")
    parser.add_argument('--output', type=str, required=True,
                        help="Output directory for all processed files.")
    parser.add_argument('--log', type=str, default='process.log',
                        help="Path to the log file (default: process.log).")
    parser.add_argument('--dry-run', action='store_true',
                        help="Log all commands without executing them. "
                             "Useful for verifying samplesheet and configuration.")
    return parser.parse_args()


def main():
    # Handle --version before parse_args() because argparse would otherwise error
    # on the required positional arguments (--project, --samplesheet, etc.).
    if '--version' in sys.argv:
        print(f"BCO Demux Tool v{__version__}")
        print("Samplesheet-driven demultiplexing, RC merging, and primer trimming for amplicon reads.")
        sys.exit(0)

    args = parse_args()
    dry_run = args.dry_run

    # Setup logging
    logging.basicConfig(
        filename=args.log, level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s', filemode='w'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logging.getLogger('').addHandler(console)

    logging.info(f"BCO Demux Tool v{__version__} | Execution started")
    logging.info(f"Command: {' '.join(sys.argv)}")
    if dry_run:
        logging.info("DRY RUN mode enabled — no commands will be executed.")

    # Load configuration
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = args.config if args.config else os.path.join(script_dir, "config.yaml")

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with open(config_path, 'r') as fh:
        config = yaml.safe_load(fh)

    if args.project not in config:
        raise ValueError(f"Project type '{args.project}' not found in configuration file.")
    pipeline_config = config[args.project]

    required_params = ['min_len', 'max_len', 'forward_primer', 'reverse_primer']
    for param in required_params:
        if param not in pipeline_config:
            raise ValueError(
                f"Missing required config parameter for project '{args.project}': {param}"
            )

    forward_primer = pipeline_config['forward_primer']
    reverse_primer = pipeline_config['reverse_primer']
    # Config stores the reverse primer 5'→3'; cutadapt -a expects the sequence as
    # it appears on the read (i.e. the RC), so we pre-compute it here.
    reverse_primer_rc = reverse_complement(reverse_primer)
    min_len = pipeline_config['min_len']
    max_len = pipeline_config['max_len']
    cores = pipeline_config.get('cores', 1)
    erate = pipeline_config.get('error_rate', 0.1)
    min_reads = pipeline_config.get('min_reads', 20)

    output_dir = args.output
    demux_temp_dir = os.path.join(output_dir, "demux_temp")
    merged_dir = os.path.join(output_dir, "merged_reads")
    trimmed_dir = os.path.join(output_dir, "trimmed_reads")

    os.makedirs(demux_temp_dir, exist_ok=True)
    os.makedirs(merged_dir, exist_ok=True)
    os.makedirs(trimmed_dir, exist_ok=True)

    try:
        check_dependencies()
        validate_input_file(args.input)

        # Load and validate samplesheet — fail early before any processing
        samples = load_samplesheet(args.samplesheet)

        demux_all_reads(samples, args.input, demux_temp_dir, output_dir,
                        min_len, max_len, cores, dry_run=dry_run)
        write_run_barcodes(samples, output_dir)
        process_and_merge_reads(samples, demux_temp_dir, merged_dir, dry_run=dry_run)
        trim_primers(merged_dir, trimmed_dir, forward_primer, reverse_primer_rc,
                     min_len, cores, erate, min_reads, dry_run=dry_run)
        generate_summary(trimmed_dir, dry_run=dry_run)
        cleanup_intermediate_dirs(demux_temp_dir, merged_dir)

        logging.info(f"BCO Demux Tool v{__version__} | Pipeline completed successfully")

    except Exception as e:
        logging.critical(f"Pipeline execution failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
