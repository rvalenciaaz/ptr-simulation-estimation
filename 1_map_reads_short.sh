#!/bin/bash
# Usage: ./map_reads.sh <reads_1.fastq> <reads_2.fastq> <reference> <output.sam> [threads]
# Example: ./map_reads.sh R1.fq R2.fq ref.fasta out.sam 8
#
# Notes:
# - The index name is derived from the reference file name (minus the extension).
# - If the index has already been built, the script will skip rebuilding it.
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

# Derive the base name from the reference (strip any path and final extension)
BASENAME=$(basename "$REFERENCE")
REFERENCE_BASENAME="${BASENAME%.*}"

# Check for existing Bowtie2 index (could be *.bt2 or *.bt2l)
if ! ls "${REFERENCE_BASENAME}".*.bt2* >/dev/null 2>&1; then
    echo "Bowtie2 index for $REFERENCE not found. Building index..."
    bowtie2-build "$REFERENCE" "$REFERENCE_BASENAME"
else
    echo "Bowtie2 index for $REFERENCE already exists. Skipping index building."
fi

# Map the paired-end short reads using Bowtie2
# -x : prefix of the index files
# -1 and -2: paired-end reads
# -S : output SAM format
# -p : number of threads
# --fast : preset for faster alignment (optional)
echo "Running Bowtie2 alignment..."
bowtie2 -x "$REFERENCE_BASENAME" \
        -1 "$READS_1" \
        -2 "$READS_2" \
        -S "$OUTPUT" \
        -p "$THREADS" \
        --fast

echo "Alignment finished. Output written to $OUTPUT"
