#!/usr/bin/env bash
#
# Usage:
#   ./run_one_folder.sh <mother_directory> <folder_name> <reference_genome> [threads]
#
# Example:
#   ./run_one_folder.sh /path/to/mother E4T11 GCF_014131755.1_ASM1413175v1_genomic.fna 16
#
# Description:
#   1) Receives four arguments (with the fourth optional):
#      - mother_directory (e.g., /path/to/mother)
#      - folder_name      (e.g., E4T11)
#      - reference_genome (e.g., GCF_014131755.1_ASM1413175v1_genomic.fna)
#      - threads          (optional, default=16)
#   2) It does NOT cd into mother_directory.
#      Instead, it references the FASTQ files as:
#         mother_directory/folder_name/folder_name_R1_trimmed.fastq
#         mother_directory/folder_name/folder_name_R2_trimmed.fastq
#   3) Creates three output folders in the current directory:
#         SAM_files/, BAM_files/, Coverage_files/
#   4) Runs:
#       - 1_map_reads_short.sh  -> produce .sam -> in SAM_files
#       - 2_sort_bam.sh        -> produce .bam -> in BAM_files (removes .sam)
#       - 3_compute_coverage.sh -> produce coverage.txt -> in Coverage_files
#   5) Each step can use **its own** Conda environment by passing the Conda binary and
#      environment name to the sub-scripts.
#   6) ALL output (stdout/stderr) and logs are written to a timestamped .log file,
#      with each line timestamped.

# --------------------------------------------------------------------
# 1. Parse arguments
# --------------------------------------------------------------------
MOTHER_DIR="$1"    # e.g., /path/to/mother
FOLDER_NAME="$2"   # e.g., E4T11
REF_GENOME="$3"
THREADS="$4"       # optional

if [ -z "$MOTHER_DIR" ] || [ -z "$FOLDER_NAME" ] || [ -z "$REF_GENOME" ]; then
  echo "ERROR: Missing required arguments."
  echo "Usage: $0 <mother_directory> <folder_name> <reference_genome> [threads]"
  exit 1
fi

# Default threads to 16 if not provided
if [ -z "$THREADS" ]; then
  THREADS=16
fi

# --------------------------------------------------------------------
# 2. (Optional) Define separate Conda info for each step
# --------------------------------------------------------------------
# Bowtie2 environment
BOWTIE2_CONDA_BIN="/home/rvalenciaaz/miniconda3/condabin/mamba"
BOWTIE2_CONDA_ENV="biotools"

# Samtools environment for sorting
SAMTOOLS_SORT_CONDA_BIN="/home/rvalenciaaz/miniconda3/condabin/mamba"
SAMTOOLS_SORT_CONDA_ENV="biotools"

# Samtools environment for coverage
SAMTOOLS_COV_CONDA_BIN="/home/rvalenciaaz/miniconda3/condabin/mamba"
SAMTOOLS_COV_CONDA_ENV="biotools"

# --------------------------------------------------------------------
# 3. Record the directory from which we're running & set up logging
# --------------------------------------------------------------------
TOP_DIR="$(pwd)"  # The directory where this script is invoked
LOGFILE="run_one_folder_$(date +%Y%m%d_%H%M%S).log"

# Redirect ALL stdout & stderr into the log file (nothing on console)
exec >"$LOGFILE" 2>&1

# --------------------------------------------------------------------
# 4. A function to timestamp log messages
# --------------------------------------------------------------------
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# --------------------------------------------------------------------
# 5. Create output directories in the current directory
# --------------------------------------------------------------------
mkdir -p "$TOP_DIR/SAM_files"
mkdir -p "$TOP_DIR/BAM_files"
mkdir -p "$TOP_DIR/Coverage_files"

# --------------------------------------------------------------------
# 6. Print initial info
# --------------------------------------------------------------------
log "-----------------------------------------------------------------"
log "Starting run_one_folder.sh"
log "Mother directory:          $MOTHER_DIR"
log "Folder name:               $FOLDER_NAME"
log "Reference genome:          $REF_GENOME"
log "Threads:                   $THREADS"
log "Current dir (TOP_DIR):     $TOP_DIR"
log "Log file:                  $TOP_DIR/$LOGFILE"
log "SAM directory:             $TOP_DIR/SAM_files/"
log "BAM directory:             $TOP_DIR/BAM_files/"
log "Coverage directory:        $TOP_DIR/Coverage_files/"
log "-----------------------------------------------------------------"
log ""

# --------------------------------------------------------------------
# 7. Check for FASTQ files
#    mother_directory/folder_name/folder_name_R1_trimmed.fastq
# --------------------------------------------------------------------
R1_FASTQ="${MOTHER_DIR}/${FOLDER_NAME}/${FOLDER_NAME}_R1_trimmed.fastq"
R2_FASTQ="${MOTHER_DIR}/${FOLDER_NAME}/${FOLDER_NAME}_R2_trimmed.fastq"

if [[ ! -f "$R1_FASTQ" || ! -f "$R2_FASTQ" ]]; then
  log "ERROR: FASTQ files not found at:"
  log "       $R1_FASTQ"
  log "       $R2_FASTQ"
  log "Exiting..."
  exit 1
fi

# --------------------------------------------------------------------
# 8. Define output filenames
# --------------------------------------------------------------------
SAM_OUTPUT="$TOP_DIR/SAM_files/${FOLDER_NAME}_mapped.sam"
BAM_OUTPUT="$TOP_DIR/BAM_files/${FOLDER_NAME}_sorted.bam"
COVERAGE_OUTPUT="$TOP_DIR/Coverage_files/${FOLDER_NAME}_coverage.txt"

# --------------------------------------------------------------------
# 9. Map Reads => .sam
# --------------------------------------------------------------------
log "[${FOLDER_NAME}] Starting: 1_map_reads_short_conda.sh"
bash "$TOP_DIR/1_map_reads_short_conda.sh" \
  "$R1_FASTQ" \
  "$R2_FASTQ" \
  "$REF_GENOME" \
  "$SAM_OUTPUT" \
  "$THREADS" \
  "$BOWTIE2_CONDA_BIN" \
  "$BOWTIE2_CONDA_ENV"
log "[${FOLDER_NAME}] Finished: 1_map_reads_short_conda.sh"

if [ ! -f "$SAM_OUTPUT" ]; then
  log "ERROR: SAM file not created -> $SAM_OUTPUT"
  log "Exiting..."
  exit 1
fi

# --------------------------------------------------------------------
# 10. Sort => .bam
# --------------------------------------------------------------------
log "[${FOLDER_NAME}] Starting: 2_sort_bam_conda.sh"
bash "$TOP_DIR/2_sort_bam_conda.sh" \
  "$SAM_OUTPUT" \
  "$BAM_OUTPUT" \
  "$SAMTOOLS_SORT_CONDA_BIN" \
  "$SAMTOOLS_SORT_CONDA_ENV"
log "[${FOLDER_NAME}] Finished: 2_sort_bam_conda.sh"

if [ ! -f "$BAM_OUTPUT" ]; then
  log "ERROR: Sorted BAM not created -> $BAM_OUTPUT"
  log "Exiting..."
  exit 1
fi

# Remove the SAM file to free space
rm -f "$SAM_OUTPUT"
log "[${FOLDER_NAME}] Removed SAM file: $SAM_OUTPUT"

# --------------------------------------------------------------------
# 11. Coverage => coverage.txt
# --------------------------------------------------------------------
log "[${FOLDER_NAME}] Starting: 3_compute_coverage_conda.sh"
bash "$TOP_DIR/3_compute_coverage_conda.sh" \
  "$BAM_OUTPUT" \
  "$COVERAGE_OUTPUT" \
  "$SAMTOOLS_COV_CONDA_BIN" \
  "$SAMTOOLS_COV_CONDA_ENV"
log "[${FOLDER_NAME}] Finished: 3_compute_coverage_conda.sh"

if [ ! -f "$COVERAGE_OUTPUT" ]; then
  log "ERROR: Coverage file not created -> $COVERAGE_OUTPUT"
  log "Exiting..."
  exit 1
fi

# --------------------------------------------------------------------
# 12. Wrap up
# --------------------------------------------------------------------
log "Completed pipeline for folder: $FOLDER_NAME"
log "BAM file:          $BAM_OUTPUT"
log "Coverage file:     $COVERAGE_OUTPUT"
log ""
log "Script finished at: $(date '+%Y-%m-%d %H:%M:%S')"
log "-----------------------------------------------------------------"
