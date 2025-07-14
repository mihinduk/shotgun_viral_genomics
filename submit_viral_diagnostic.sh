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

# Set paths
PIPELINE_DIR="/scratch/sahlab/kathie/Diamond_test/shotgun_viral_genomics"
DIAGNOSTIC_SCRIPT="${PIPELINE_DIR}/viral_diagnostic.sh"

# Validate input files exist
if [ ! -f "$R1" ]; then
    echo "Error: R1 file not found: $R1"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "Error: R2 file not found: $R2"
    exit 1
fi

echo "Running viral diagnostic with:"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Reference: $ACCESSION"
echo "  Sample: $SAMPLE_NAME"
echo "  Threads: $THREADS"
echo "  Working directory: $(pwd)"
echo ""

# Run the diagnostic script
bash "$DIAGNOSTIC_SCRIPT" "$R1" "$R2" "$ACCESSION" "$SAMPLE_NAME" "$THREADS"

DIAGNOSTIC_EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job completed at: $(date)"
if [ $DIAGNOSTIC_EXIT_CODE -eq 0 ]; then
    echo "Diagnostic analysis completed successfully!"
    echo ""
    echo "Key outputs:"
    if [ -d "./diagnostic_${SAMPLE_NAME}" ]; then
        echo "  Report: ./diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_diagnostic_report.txt"
        echo "  BLAST results: ./diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_viral_blast.tsv"
        echo "  Assembly: ./diagnostic_${SAMPLE_NAME}/assembly_${SAMPLE_NAME}/final.contigs.fa"
        echo ""
        
        # Show quick summary if report exists
        if [ -f "./diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_diagnostic_report.txt" ]; then
            echo "Quick Summary:"
            grep -E "(Mapping Percentage|Total Contigs|Contigs >1000bp)" "./diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_diagnostic_report.txt" || true
        fi
    fi
else
    echo "Diagnostic analysis failed with exit code: $DIAGNOSTIC_EXIT_CODE"
fi
echo "========================================="

exit $DIAGNOSTIC_EXIT_CODE