#!/usr/bin/env bash
#
# Usage:
#   ./run_all_folders.sh <mother_directory> <reference_genome> [threads]
#
# Example:
#   ./run_all_folders.sh /path/to/mother GCF_014131755.1_ASM1413175v1_genomic.fna 16
#
# Description:
#   1) Receives three arguments (with the third optional):
#      - mother_directory (e.g., /path/to/mother)
#      - reference_genome (e.g., GCF_014131755.1_ASM1413175v1_genomic.fna)
#      - threads          (optional, default=16)
#   2) Finds all subfolders in mother_directory whose names match the pattern:
#         stringNumberStringNumber
#         stringNumberStringNumberelse
#   3) For each matching folder, it references the FASTQ files as:
#         mother_directory/folder_name/folder_name_R1_trimmed.fastq
#         mother_directory/folder_name/folder_name_R2_trimmed.fastq
#   4) Creates three output folders in the current directory:
#         SAM_files/, BAM_files/, Coverage_files/
#   5) Runs:
#       - 1_map_reads_short.sh  -> produce .sam -> in SAM_files
#       - 2_sort_bam.sh        -> produce .bam -> in BAM_files (removes .sam)
#       - 3_compute_coverage.sh -> produce coverage.txt -> in Coverage_files
#   6) Each step can use **its own** Conda environment by passing the Conda binary and
#      environment name to the sub-scripts.
#   7) ALL output (stdout/stderr) and logs are written to a timestamped .log file
#      (per folder), with each line timestamped.

# --------------------------------------------------------------------
# 1. Parse arguments
# --------------------------------------------------------------------
MOTHER_DIR="$1"    # e.g., /path/to/mother
REF_GENOME="$2"
THREADS="$3"       # optional

if [ -z "$MOTHER_DIR" ] || [ -z "$REF_GENOME" ]; then
  echo "ERROR: Missing required arguments."
  echo "Usage: $0 <mother_directory> <reference_genome> [threads]"
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
# 3. Record the directory from which we're running
# --------------------------------------------------------------------
TOP_DIR="$(pwd)"  # The directory where this script is invoked

# --------------------------------------------------------------------
# 4. Create output directories in the current directory (if not exist)
# --------------------------------------------------------------------
mkdir -p "$TOP_DIR/SAM_files"
mkdir -p "$TOP_DIR/BAM_files"
mkdir -p "$TOP_DIR/Coverage_files"

# --------------------------------------------------------------------
# 5. A helper function to timestamp log messages
# --------------------------------------------------------------------
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# --------------------------------------------------------------------
# 6. A function to process *one* folder
# --------------------------------------------------------------------
process_one_folder() {
  local folder_name="$1"

  # We create a unique log file for each folder
  local logfile="run_${folder_name}_$(date +%Y%m%d_%H%M%S).log"
  
  # Redirect ALL stdout & stderr into the log file (nothing on console)
  exec >"$logfile" 2>&1

  # Print header info
  log "-----------------------------------------------------------------"
  log "Processing folder:         $folder_name"
  log "Mother directory:          $MOTHER_DIR"
  log "Reference genome:          $REF_GENOME"
  log "Threads:                   $THREADS"
  log "Current dir (TOP_DIR):     $TOP_DIR"
  log "Log file:                  $TOP_DIR/$logfile"
  log "SAM directory:             $TOP_DIR/SAM_files/"
  log "BAM directory:             $TOP_DIR/BAM_files/"
  log "Coverage directory:        $TOP_DIR/Coverage_files/"
  log "-----------------------------------------------------------------"
  log ""

  # ------------------------------------------------------------------
  # Define FASTQ input paths
  # ------------------------------------------------------------------
  R1_FASTQ="${MOTHER_DIR}/${folder_name}/${folder_name}_R1_trimmed.fastq"
  R2_FASTQ="${MOTHER_DIR}/${folder_name}/${folder_name}_R2_trimmed.fastq"

  if [[ ! -f "$R1_FASTQ" || ! -f "$R2_FASTQ" ]]; then
    log "ERROR: FASTQ files not found at:"
    log "       $R1_FASTQ"
    log "       $R2_FASTQ"
    log "Skipping folder: $folder_name"
    log "-----------------------------------------------------------------"
    return 1
  fi

  # ------------------------------------------------------------------
  # Define output filenames
  # ------------------------------------------------------------------
  SAM_OUTPUT="$TOP_DIR/SAM_files/${folder_name}_mapped.sam"
  BAM_OUTPUT="$TOP_DIR/BAM_files/${folder_name}_sorted.bam"
  COVERAGE_OUTPUT="$TOP_DIR/Coverage_files/${folder_name}_coverage.txt"

  # ------------------------------------------------------------------
  # 1. Map Reads => .sam
  # ------------------------------------------------------------------
  log "[${folder_name}] Starting: 1_map_reads_short_conda.sh"
  bash "$TOP_DIR/1_map_reads_short_conda.sh" \
    "$R1_FASTQ" \
    "$R2_FASTQ" \
    "$REF_GENOME" \
    "$SAM_OUTPUT" \
    "$THREADS" \
    "$BOWTIE2_CONDA_BIN" \
    "$BOWTIE2_CONDA_ENV"
  log "[${folder_name}] Finished: 1_map_reads_short_conda.sh"

  if [ ! -f "$SAM_OUTPUT" ]; then
    log "ERROR: SAM file not created -> $SAM_OUTPUT"
    log "Skipping folder: $folder_name"
    log "-----------------------------------------------------------------"
    return 1
  fi

  # ------------------------------------------------------------------
  # 2. Sort => .bam
  # ------------------------------------------------------------------
  log "[${folder_name}] Starting: 2_sort_bam_conda.sh"
  bash "$TOP_DIR/2_sort_bam_conda.sh" \
    "$SAM_OUTPUT" \
    "$BAM_OUTPUT" \
    "$SAMTOOLS_SORT_CONDA_BIN" \
    "$SAMTOOLS_SORT_CONDA_ENV"
  log "[${folder_name}] Finished: 2_sort_bam_conda.sh"

  if [ ! -f "$BAM_OUTPUT" ]; then
    log "ERROR: Sorted BAM not created -> $BAM_OUTPUT"
    log "Skipping folder: $folder_name"
    log "-----------------------------------------------------------------"
    return 1
  fi

  # Remove the SAM file to free space
  rm -f "$SAM_OUTPUT"
  log "[${folder_name}] Removed SAM file: $SAM_OUTPUT"

  # ------------------------------------------------------------------
  # 3. Coverage => coverage.txt
  # ------------------------------------------------------------------
  log "[${folder_name}] Starting: 3_compute_coverage_conda.sh"
  bash "$TOP_DIR/3_compute_coverage_conda.sh" \
    "$BAM_OUTPUT" \
    "$COVERAGE_OUTPUT" \
    "$SAMTOOLS_COV_CONDA_BIN" \
    "$SAMTOOLS_COV_CONDA_ENV"
  log "[${folder_name}] Finished: 3_compute_coverage_conda.sh"

  if [ ! -f "$COVERAGE_OUTPUT" ]; then
    log "ERROR: Coverage file not created -> $COVERAGE_OUTPUT"
    log "Skipping folder: $folder_name"
    log "-----------------------------------------------------------------"
    return 1
  fi

  # ------------------------------------------------------------------
  # Finished with this folder
  # ------------------------------------------------------------------
  log "Completed pipeline for folder: $folder_name"
  log "BAM file:          $BAM_OUTPUT"
  log "Coverage file:     $COVERAGE_OUTPUT"
  log ""
  log "Script finished at: $(date '+%Y-%m-%d %H:%M:%S')"
  log "-----------------------------------------------------------------"

  # Return success
  return 0
}

# --------------------------------------------------------------------
# 7. Main Loop: find all folders matching the pattern
# --------------------------------------------------------------------
PATTERN='^[A-Za-z]+[0-9]+[A-Za-z]+[0-9]+(else)?$'

# List all directories in $MOTHER_DIR that match $PATTERN
for folder_name in "$(ls -1 "$MOTHER_DIR" | grep -E "$PATTERN")"; do
  # Verify it's actually a directory (just in case)
  if [ -d "${MOTHER_DIR}/${folder_name}" ]; then
    # Process the folder
    process_one_folder "$folder_name"
    rc=$?
    if [ $rc -ne 0 ]; then
      echo "[Main Script] Pipeline failed for folder: $folder_name"
      # Uncomment the line below if you want to stop as soon as a folder fails
      # break
    fi
  fi
done

# echo "All matching folders processed!"
exit 0
