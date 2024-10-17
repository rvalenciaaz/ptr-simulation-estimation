#!/bin/bash
# Usage: ./sort_bam.sh <input.sam> <sorted.bam>

INPUT_SAM=$1
OUTPUT_BAM=$2

# Convert SAM to BAM and sort
samtools view -bS "$INPUT_SAM" | samtools sort -o "$OUTPUT_BAM"
