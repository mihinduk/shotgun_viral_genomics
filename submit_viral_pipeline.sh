#!/bin/bash
#SBATCH --job-name=viral_pipeline
#SBATCH --output=viral_pipeline_%j.out
#SBATCH --error=viral_pipeline_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general

# Submit script for viral pipeline with proper --large-files handling
# Usage: sbatch submit_viral_pipeline.sh <R1> <R2> <accession> <threads> [--large-files]

# Check arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> <threads> [--large-files]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4 --large-files"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=$4
LARGE_FILES_FLAG=""

# Check if --large-files flag is provided
if [ $# -eq 5 ] && [ "$5" == "--large-files" ]; then
    LARGE_FILES_FLAG="--large-files"
    echo "Using --large-files flag for increased memory allocation"
    # Update SLURM resources for large files
    export SLURM_MEM_PER_NODE=64G
fi

# Change to the directory where the script was submitted from
cd ${SLURM_SUBMIT_DIR}

echo "========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Job started at: $(date)"
echo "Working directory: $(pwd)"
echo "========================================="

# Set up environment
echo "Setting up HTCF environment..."
# Note: The wrapper script handles conda/mamba activation internally

# Set paths
PIPELINE_DIR="/scratch/sahlab/kathie/Diamond_test/shotgun_viral_genomics"
WRAPPER_SCRIPT="${PIPELINE_DIR}/run_pipeline_htcf_consolidated.sh"

# Validate input files exist
if [ ! -f "$R1" ]; then
    echo "Error: R1 file not found: $R1"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "Error: R2 file not found: $R2"
    exit 1
fi

echo "Running viral pipeline with:"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Accession: $ACCESSION"
echo "  Threads: $THREADS"
echo "  Large files flag: $LARGE_FILES_FLAG"
echo "  Working directory: $(pwd)"
echo ""

# Run the pipeline wrapper script
if [ -n "$LARGE_FILES_FLAG" ]; then
    # Run with large files support
    bash "$WRAPPER_SCRIPT" "$R1" "$R2" "$ACCESSION" "$THREADS" "$LARGE_FILES_FLAG"
else
    # Run with standard settings
    bash "$WRAPPER_SCRIPT" "$R1" "$R2" "$ACCESSION" "$THREADS"
fi

PIPELINE_EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job completed at: $(date)"
if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo "Pipeline completed successfully!"
else
    echo "Pipeline failed with exit code: $PIPELINE_EXIT_CODE"
fi
echo "========================================="

exit $PIPELINE_EXIT_CODE