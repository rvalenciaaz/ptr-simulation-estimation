#!/bin/bash
# Usage:
#   ./sort_bam.sh <input.sam> <sorted.bam> <conda_binary> <conda_env>
#
# Example:
#   ./sort_bam.sh input.sam output.bam /path/to/miniconda3/condabin/conda samtools_env
#
# Description:
#   Converts a SAM file to a BAM file, then sorts it, using `samtools`.
#   The commands are run through `conda run -n <env>` so that `samtools`
#   is used from the specified Conda environment.

INPUT_SAM=$1
OUTPUT_BAM=$2
CONDA_BIN=$3
CONDA_ENV=$4

# Check for required arguments
if [ -z "$INPUT_SAM" ] || [ -z "$OUTPUT_BAM" ] || [ -z "$CONDA_BIN" ] || [ -z "$CONDA_ENV" ]; then
  echo "Error: Missing required arguments."
  echo "Usage: $0 <input.sam> <sorted.bam> <conda_binary> <conda_env>"
  exit 1
fi

# Convert SAM to BAM and sort, using conda run
# `samtools view -bS` = convert SAM to BAM
# `samtools sort -o` = sort and write to output file
"$CONDA_BIN" run -n "$CONDA_ENV" samtools view -bS "$INPUT_SAM" | \
"$CONDA_BIN" run -n "$CONDA_ENV" samtools sort -o "$OUTPUT_BAM"

echo "Sorted BAM file created: $OUTPUT_BAM"
