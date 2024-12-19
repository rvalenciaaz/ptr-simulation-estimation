#!/bin/bash
# Usage: ./map_reads.sh <reads.fastq> <reference.fasta> <output.sam>

READS=$1
REFERENCE=$2
OUTPUT=$3

# Map reads using minimap2
minimap2 -ax map-ont "$REFERENCE" "$READS" > "$OUTPUT"
