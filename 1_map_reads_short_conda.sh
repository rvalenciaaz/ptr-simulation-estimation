#!/bin/bash -l
# Usage:
#   ./map_reads.sh <reads_1.fastq> <reads_2.fastq> <reference> <output.sam> [threads] <conda_binary> <conda_env>
#
# Example:
#   ./map_reads.sh R1.fq R2.fq ref.fasta out.sam 8 /path/to/miniconda3/condabin/conda bowtie2_env
#
# Arguments:
#   1) Paired-end reads #1          (e.g. R1.fq)
#   2) Paired-end reads #2          (e.g. R2.fq)
#   3) Reference file               (any extension, e.g. ref.fasta, ref.fa, etc.)
#   4) Output SAM file              (e.g. out.sam)
#   5) [OPTIONAL] Number of threads (default: 1)
#   6) Path to conda binary         (e.g. /path/to/miniconda3/condabin/conda)
#   7) Conda environment name       (e.g. bowtie2_env or /path/to/envs/bowtie2_env)
#
# Description:
#   This script maps paired-end reads to a reference genome using Bowtie2.
#   It will build the Bowtie2 index automatically if it doesn't already exist.
#   Instead of activating the Conda environment in the current shell, it uses
#   `conda run -n <env>` for each Bowtie2 command.

READS_1=$1
READS_2=$2
REFERENCE=$3
OUTPUT=$4
THREADS=$5
CONDA_BIN=$6
CONDA_ENV=$7

# Basic argument checks
if [ -z "$READS_1" ] || [ -z "$READS_2" ] || [ -z "$REFERENCE" ] || [ -z "$OUTPUT" ] || [ -z "$CONDA_BIN" ] || [ -z "$CONDA_ENV" ]; then
  echo "Error: Missing required arguments."
  echo "Usage: $0 <reads_1.fastq> <reads_2.fastq> <reference> <output.sam> [threads] <conda_binary> <conda_env>"
  exit 1
fi

# Default THREADS to 1 if not provided
if [ -z "$THREADS" ]; then
  THREADS=1
fi

# ------------------------------------------------------------------------------
# Derive reference basename (strips the final extension)
# ------------------------------------------------------------------------------
BASENAME=$(basename "$REFERENCE")
REFERENCE_BASENAME="${BASENAME%.*}"

# ------------------------------------------------------------------------------
# Build Bowtie2 index only if not present
# ------------------------------------------------------------------------------
if ! ls "${REFERENCE_BASENAME}".*.bt2* >/dev/null 2>&1; then
    echo "Bowtie2 index for '$REFERENCE' not found. Building index..."
    "$CONDA_BIN" run -n "$CONDA_ENV" bowtie2-build "$REFERENCE" "$REFERENCE_BASENAME"
else
    echo "Bowtie2 index for '$REFERENCE' already exists. Skipping index building."
fi

# ------------------------------------------------------------------------------
# Map the paired-end reads with Bowtie2 (using conda run)
# ------------------------------------------------------------------------------
echo "Running Bowtie2 alignment..."
"$CONDA_BIN" run -n "$CONDA_ENV" bowtie2 \
    -x "$REFERENCE_BASENAME" \
    -1 "$READS_1" \
    -2 "$READS_2" \
    -S "$OUTPUT" \
    -p "$THREADS" \
    --fast

echo "Alignment finished. Output written to $OUTPUT"
