#!/bin/bash
#SBATCH --job-name=viral_diagnostic_extreme
#SBATCH --output=viral_diagnostic_%j.out
#SBATCH --error=viral_diagnostic_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=8
#SBATCH --partition=general

# Submit script for viral contamination diagnostic module (EXTREME MEMORY)
# Usage: sbatch submit_viral_diagnostic_extreme.sh <R1> <R2> <accession> <sample_name> [threads]

# Check arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> <sample_name> [threads]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 SMS_14 8"
    echo ""
    echo "This EXTREME MEMORY diagnostic script runs in parallel with the main pipeline to:"
    echo "  1. Check mapping statistics to expected reference"
    echo "  2. Perform de novo assembly with MEGAHIT (256GB memory)"
    echo "  3. BLAST contigs against viral database for contamination detection"
    echo "  4. Generate comprehensive diagnostic report"
    echo ""
    echo "Use this version for very large files (>50GB) that fail with regular diagnostic."
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
SAMPLE_NAME=$4
THREADS=${5:-8}

# Change to the directory where the script was submitted from
cd ${SLURM_SUBMIT_DIR}

echo "========================================="
echo "VIRAL DIAGNOSTIC MODULE (EXTREME MEMORY)"
echo "========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Job started at: $(date)"
echo "Working directory: $(pwd)"
echo "Memory allocated: 256G"
echo "Partition: general"
echo "Time limit: 6:00:00"
echo "========================================="
echo "Running viral diagnostic with:"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Reference: $ACCESSION"
echo "  Sample: $SAMPLE_NAME"
echo "  Threads: $THREADS"
echo "  Memory mode: extreme (256GB)"
echo "  Working directory: $(pwd)"

# Extract pipeline directory from the original sbatch command
ORIGINAL_SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | grep -oP 'Command=\K[^ ]+')
PIPELINE_DIR="$(dirname "$ORIGINAL_SCRIPT_PATH")"

echo "Original script path: $ORIGINAL_SCRIPT_PATH"
echo "Pipeline directory: $PIPELINE_DIR"

# Call the main diagnostic script with extreme memory flag
"${PIPELINE_DIR}/viral_diagnostic.sh" "$R1" "$R2" "$ACCESSION" "$SAMPLE_NAME" "$THREADS" --extremely-large-files

echo ""
echo "========================================="
echo "Job completed at: $(date)"
if [ $? -eq 0 ]; then
    echo "Diagnostic analysis completed successfully"
else
    echo "Diagnostic analysis failed with exit code: $?"
fi
echo "========================================="