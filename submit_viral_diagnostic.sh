#!/bin/bash
#SBATCH --job-name=viral_diagnostic
#SBATCH --output=viral_diagnostic_%j.out
#SBATCH --error=viral_diagnostic_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general

# Submit script for viral contamination diagnostic module
# Usage: sbatch submit_viral_diagnostic.sh <R1> <R2> <accession> <sample_name> [threads]

# Check arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> <sample_name> [threads]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 SMS_14 4"
    echo ""
    echo "This diagnostic script runs in parallel with the main pipeline to:"
    echo "  1. Check mapping statistics to expected reference"
    echo "  2. Perform de novo assembly with MEGAHIT"
    echo "  3. BLAST contigs against viral database for contamination detection"
    echo "  4. Generate comprehensive diagnostic report"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
SAMPLE_NAME=$4
THREADS=${5:-4}

# Change to the directory where the script was submitted from
cd ${SLURM_SUBMIT_DIR}

echo "========================================="
echo "VIRAL DIAGNOSTIC MODULE"
echo "========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Job started at: $(date)"
echo "Working directory: $(pwd)"
echo "========================================="

# Extract pipeline directory from the original sbatch command
# SLURM preserves the original command in the job details
ORIGINAL_SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | grep -oP 'Command=\K[^ ]+')
PIPELINE_DIR="$(dirname "$ORIGINAL_SCRIPT_PATH")"
DIAGNOSTIC_SCRIPT="${PIPELINE_DIR}/viral_diagnostic.sh"

echo "Original script path: $ORIGINAL_SCRIPT_PATH"
echo "Pipeline directory: $PIPELINE_DIR"
echo "Diagnostic script: $DIAGNOSTIC_SCRIPT"

# Verify the script exists
if [ ! -f "$DIAGNOSTIC_SCRIPT" ]; then
    echo "ERROR: Cannot find viral_diagnostic.sh at: $DIAGNOSTIC_SCRIPT"
    echo "This usually means the script wasn't submitted with full path"
    echo "Please use: sbatch /full/path/to/submit_viral_diagnostic.sh"
    exit 1
fi

# Run the diagnostic script
echo "Running viral diagnostic with:"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Reference: $ACCESSION"
echo "  Sample: $SAMPLE_NAME"
echo "  Threads: $THREADS"
echo "  Working directory: $(pwd)"

# Execute the diagnostic script
"$DIAGNOSTIC_SCRIPT" "$R1" "$R2" "$ACCESSION" "$SAMPLE_NAME" "$THREADS"
DIAGNOSTIC_EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job completed at: $(date)"
if [ $DIAGNOSTIC_EXIT_CODE -eq 0 ]; then
    echo "Diagnostic analysis completed successfully"
else
    echo "Diagnostic analysis failed with exit code: $DIAGNOSTIC_EXIT_CODE"
fi
echo "========================================="

exit $DIAGNOSTIC_EXIT_CODE