#!/bin/bash
# Usage:
#   ./compute_coverage.sh <sorted.bam> <coverage.txt> <conda_binary> <conda_env>
#
# Example:
#   ./compute_coverage.sh sorted.bam coverage.txt /path/to/miniconda3/condabin/conda samtools_env
#
# Description:
#   Computes per-base coverage from a sorted BAM file using `samtools depth`.
#   The command is run through `conda run -n <env>` so that `samtools`
#   is used from the specified Conda environment.

SORTED_BAM=$1
COVERAGE_FILE=$2
CONDA_BIN=$3
CONDA_ENV=$4

# Check for required arguments
if [ -z "$SORTED_BAM" ] || [ -z "$COVERAGE_FILE" ] || [ -z "$CONDA_BIN" ] || [ -z "$CONDA_ENV" ]; then
  echo "Error: Missing required arguments."
  echo "Usage: $0 <sorted.bam> <coverage.txt> <conda_binary> <conda_env>"
  exit 1
fi

# Compute per-base coverage using conda run
"$CONDA_BIN" run -n "$CONDA_ENV" samtools depth "$SORTED_BAM" > "$COVERAGE_FILE"

echo "Coverage information written to: $COVERAGE_FILE"
