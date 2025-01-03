#!/bin/bash
# Usage: ./map_reads.sh <reads_1.fastq> <reads_2.fastq> <reference.fasta> <output.sam> [threads]
# Example: ./map_reads.sh R1.fq R2.fq ref.fasta out.sam 8
#
# Notes:
# - If the index has been previously built, you can skip rebuilding it every time.
#   Just ensure the `reference_index.*` files are present.
# - Adjust the Bowtie2 parameters (like --fast) for your accuracy/performance needs.

READS_1=$1
READS_2=$2
REFERENCE=$3
OUTPUT=$4
THREADS=$5

# If THREADS is not specified, default to 1
if [ -z "$THREADS" ]; then
  THREADS=1
fi

# Build the Bowtie2 index only if not already present
# If you run this script multiple times on the same reference, 
# consider commenting out the index-building step after the first run.
if [ ! -f "reference_index.1.bt2" ] && [ ! -f "reference_index.1.bt2l" ]; then
    bowtie2-build "$REFERENCE" reference_index
fi

# Map the paired-end short reads using Bowtie2
# -x reference_index: specifies the prefix of the index files
# -1 and -2: paired-end reads
# -S: output SAM format
# -p: number of threads for parallel alignment
# --fast: a preset for faster alignment (optional)
bowtie2 -x reference_index -1 "$READS_1" -2 "$READS_2" -S "$OUTPUT" -p "$THREADS" 

#--fast

