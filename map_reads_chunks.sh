#!/bin/bash
# Usage: ./map_reads_chunks.sh <reads_1.fastq> <reads_2.fastq> <chunks_folder> <output_sams_folder> [threads]
# Example: ./map_reads_chunks.sh R1.fq R2.fq output_chunks sams_folder 8
#
# This script will:
# 1. Iterate over each FASTA file in the chunks_folder.
# 2. Build a Bowtie2 index for that chunk if not already present.
# 3. Map the paired-end reads against that chunk reference.
# 4. Save the output SAM file in the output_sams_folder.
# 5. Remove the Bowtie2 index files to save space after mapping is done.

READS_1=$1
READS_2=$2
CHUNKS_FOLDER=$3
OUTPUT_SAMS_FOLDER=$4
THREADS=$5

# If THREADS is not specified, default to 1
if [ -z "$THREADS" ]; then
  THREADS=1
fi

# Create output directory if it does not exist
mkdir -p "$OUTPUT_SAMS_FOLDER"

# Loop over each FASTA file (chunk) in the chunks_folder
for CHUNK_REF in "$CHUNKS_FOLDER"/*.fa; do
    # Extract the base name of the chunk (without path and extension)
    CHUNK_BASENAME=$(basename "$CHUNK_REF" .fa)
    
    # Index prefix will be based on the chunk basename
    INDEX_PREFIX="$CHUNKS_FOLDER/${CHUNK_BASENAME}_index"
    
    # Check if index files exist; if not, build them
    if [ ! -f "${INDEX_PREFIX}.1.bt2" ] && [ ! -f "${INDEX_PREFIX}.1.bt2l" ]; then
        bowtie2-build "$CHUNK_REF" "$INDEX_PREFIX"
    fi
    
    # Output SAM file will be named after the chunk
    OUTPUT_SAM="$OUTPUT_SAMS_FOLDER/${CHUNK_BASENAME}.sam"
    
    # Map the reads to this chunk
    bowtie2 -x "$INDEX_PREFIX" -1 "$READS_1" -2 "$READS_2" -S "$OUTPUT_SAM" -p "$THREADS"
    
    echo "Mapping completed for $CHUNK_REF -> $OUTPUT_SAM"

    # Remove the index files to save space
    rm -f ${INDEX_PREFIX}.*bt2*
done
