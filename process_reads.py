#!/usr/bin/env python3

import os
import subprocess
import argparse
import glob
import gzip
import logging
import yaml

# Function to reverse complement a primer sequence
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M', 'S': 'S', 'W': 'W', 'H': 'D', 'D': 'H', 'B': 'V', 'V': 'B', 'N': 'N'}
    return "".join(complement[base] for base in reversed(seq))

# Function to count the number of sequences in a gzipped FASTQ file
def count_sequences(file_path):
    count = 0
    try:
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if line.startswith('@'):
                    count += 1
    except Exception as e:
        logging.error(f"Failed to count sequences in {file_path}: {e}")
    return count

# Function to run shell commands and log output
def run_command(cmd):
    logging.info(f"Running command: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}\nError: {e}")

# Command-line arguments parser
parser = argparse.ArgumentParser(description="Process reads by demultiplexing, reverse complementing, merging, and trimming primers.")
parser.add_argument('--config', type=str, help="Path to the YAML configuration file")
parser.add_argument('--project', type=str, required=True, help="Project type (e.g., 16S, ITS)")
parser.add_argument('--input', type=str, required=True, help="Input reads file (FASTQ.gz format)")
parser.add_argument('--output', type=str, required=True, help="Output directory for processed files")
parser.add_argument('--log', type=str, default='process.log', help="Log file to store the output")
args = parser.parse_args()

# Determine script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))
default_config_path = os.path.join(script_dir, "config.yaml")

# Check if --config argument is provided, else use the default config path
config_path = args.config if args.config else default_config_path

# Load configuration from YAML file
if not os.path.exists(config_path):
    raise FileNotFoundError(f"Configuration file not found: {config_path}")
with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

# Validate project type
if args.project not in config:
    raise ValueError(f"Project type '{args.project}' not found in configuration file")
pipeline_config = config[args.project]

# Setup logging
logging.basicConfig(filename=args.log, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Validate required parameters
required_params = ['min_len', 'max_len', 'forward_primer', 'reverse_primer']
for param in required_params:
    if param not in pipeline_config:
        raise ValueError(f"Missing required parameter in configuration for project '{args.project}': {param}")

forward_primer = pipeline_config['forward_primer']
reverse_primer = pipeline_config['reverse_primer']
reverse_primer_rc = reverse_complement(reverse_primer)
min_len = pipeline_config['min_len']
max_len = pipeline_config['max_len']
cores = pipeline_config.get('cores', 1)
erate = pipeline_config.get('error_rate', 0.1)

reads_file = args.input
output_dir = args.output

fwd_demux_dir = os.path.join(output_dir, "fdemuxed")
rev_demux_dir = os.path.join(output_dir, "rdemuxed")
merged_dir = os.path.join(output_dir, "merged_reads")
trimmed_dir = os.path.join(output_dir, "trimmed_reads")

# Create necessary directories if they don't exist
os.makedirs(fwd_demux_dir, exist_ok=True)
os.makedirs(rev_demux_dir, exist_ok=True)
os.makedirs(merged_dir, exist_ok=True)
os.makedirs(trimmed_dir, exist_ok=True)

# Step 1: Demultiplex forward reads
def demux_forward():
    cmd = (
        f"cutadapt -e 0 -O 24 -g file:{pipeline_config.get('forward_barcodes')} --trimmed-only -m {min_len} -M {max_len} "
        f"-j {cores} -o {fwd_demux_dir}/{{name}}.fastq.gz {reads_file}"
    )
    run_command(cmd)

# Step 2: Demultiplex reverse reads
def demux_reverse():
    cmd = (
        f"cutadapt -e 0 -O 24 -g file:{pipeline_config.get('reverse_barcodes')} --trimmed-only -m {min_len} -M {max_len} "
        f"-j {cores} -o {rev_demux_dir}/{{name}}.fastq.gz {reads_file}"
    )
    run_command(cmd)

# Step 3: Reverse complement the reverse reads
def reverse_complement_reads():
    reverse_files = glob.glob(f"{rev_demux_dir}/*.fastq.gz")
    for rev_file in reverse_files:
        base_name = os.path.basename(rev_file)
        rev_comp_file = os.path.join(rev_demux_dir, f"rev_comp_{base_name}")
        cmd = f"seqkit seq -rp --seq-type DNA -o {rev_comp_file} {rev_file}"
        run_command(cmd)

# Step 4: Merge corresponding forward and reverse reads
def merge_reads():
    forward_files = glob.glob(f"{fwd_demux_dir}/*.fastq.gz")
    for fwd_file in forward_files:
        base_name = os.path.basename(fwd_file)
        rev_comp_file = os.path.join(rev_demux_dir, f"rev_comp_{base_name}")
        merged_file = os.path.join(merged_dir, f"merged_{base_name}")
        if os.path.exists(rev_comp_file):
            cmd = f"cat {fwd_file} {rev_comp_file} > {merged_file}"
            run_command(cmd)

# Step 5: Trim forward and reverse primers from merged reads
def trim_primers():
    merged_files = glob.glob(f"{merged_dir}/*.fastq.gz")
    trimmed_files = []
    for merged_file in merged_files:
        num_sequences = count_sequences(merged_file)
        if num_sequences >= 5:
            base_name = os.path.basename(merged_file)
            barcode = base_name.replace("merged_", "").replace(".fastq.gz", "")
            output_file = os.path.join(trimmed_dir, f"{barcode}_trimmed.fastq.gz")
            output_directory = os.path.dirname(output_file)
            os.makedirs(output_directory, exist_ok=True)
            cmd = (
                f"cutadapt -g {forward_primer} -a {reverse_primer_rc} --trimmed-only "
                f"-j {cores} --error-rate {erate} -m {min_len} -o {output_file} {merged_file}"
            )
            try:
                run_command(cmd)
                trimmed_files.append(output_file)
            except subprocess.CalledProcessError as e:
                logging.error(f"Error during trimming {merged_file}: {e}")
        else:
            logging.warning(f"{merged_file} contains only {num_sequences} sequences. Deleting the file.")
            os.remove(merged_file)

# Step 6: Generate summary file
def generate_summary():
    summary_file = os.path.join(trimmed_dir, "summary.tsv")
    with open(summary_file, 'w') as summary:
        trimmed_files = sorted(glob.glob(f"{trimmed_dir}/*.fastq.gz"))
        
        for idx, trimmed_file in enumerate(trimmed_files):
            cmd = f"seqkit stats -aT {trimmed_file}"
            result = subprocess.run(
                cmd, 
                shell=True, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                universal_newlines=True
            )
            if result.returncode != 0:
                logging.error(f"Failed to generate stats for {trimmed_file}. Error: {result.stderr}")
                continue

            # Parse seqkit output
            lines = result.stdout.strip().split('\n')
            if len(lines) < 2:
                logging.warning(f"Unexpected output from seqkit stats for {trimmed_file}: {result.stdout}")
                continue

            header = lines[0].split()  # Column names
            stats = lines[1].split()  # Data row
            
            # Remove unimportant columns
            header = [col for col in header if col in ('file', 'num_seqs', 'min_len',  'avg_len', 'max_len', 'Q20(%)',  'AvgQual')]
            stats = [stat for i, stat in enumerate(stats) if i in (1, 3, 5, 6, 7, 14, 16)]

            # Write header to the file only once
            if idx == 0:
                summary.write("\t".join(header) + "\n")
            
            # Use just the base name of the file for the first column
            stats[0] = os.path.basename(trimmed_file)
            summary.write("\t".join(stats) + "\n")
            logging.info(f"Stats for {trimmed_file} written to summary.")

# Cleanup intermediate directories
def cleanup_intermediate_dirs():
    intermediate_dirs = [fwd_demux_dir, rev_demux_dir, merged_dir]
    for directory in intermediate_dirs:
        if os.path.exists(directory):
            logging.info(f"Deleting intermediate directory: {directory}")
            for file in glob.glob(f"{directory}/*"):
                os.remove(file)
            os.rmdir(directory)

if __name__ == "__main__":
    # Clear the log file
    open(args.log, 'w').close()
    
    # Execute the steps
    demux_forward()
    demux_reverse()
    reverse_complement_reads()
    merge_reads()
    trim_primers()
    generate_summary()
    cleanup_intermediate_dirs()
