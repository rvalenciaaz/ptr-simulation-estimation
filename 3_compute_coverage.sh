#!/bin/bash
# Usage: ./compute_coverage.sh <sorted.bam> <coverage.txt>

SORTED_BAM=$1
COVERAGE_FILE=$2

# Compute per-base coverage
samtools depth "$SORTED_BAM" > "$COVERAGE_FILE"
