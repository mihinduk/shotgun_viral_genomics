#!/bin/bash
# Launch viral pipeline and diagnostic module in parallel
# Usage: ./run_parallel_analysis.sh <R1> <R2> <accession> <sample_name> [threads] [memory_flag]

# Check arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> <sample_name> [threads] [memory_flag]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 SMS_14"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 SMS_14 4"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 SMS_14 4 --large-files"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 SMS_14 4 --extremely-large-files"
    echo ""
    echo "This script launches two jobs in parallel:"
    echo "  1. Main viral genomics pipeline (variant calling)"
    echo "  2. Contamination diagnostic module (assembly + BLAST)"
    echo ""
    echo "Both jobs will run simultaneously and provide complementary analysis."
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
SAMPLE_NAME=$4
THREADS=${5:-4}
MEMORY_FLAG=$6

# Set paths
PIPELINE_DIR="/scratch/sahlab/kathie/Diamond_test/shotgun_viral_genomics"

# Validate input files exist
if [ ! -f "$R1" ]; then
    echo "Error: R1 file not found: $R1"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "Error: R2 file not found: $R2"
    exit 1
fi

echo "========================================="
echo "PARALLEL VIRAL ANALYSIS LAUNCHER"
echo "========================================="
echo "Sample: $SAMPLE_NAME"
echo "R1: $R1"
echo "R2: $R2"
echo "Reference: $ACCESSION"
echo "Threads: $THREADS"
echo "Memory flag: ${MEMORY_FLAG:-none}"
echo "Working directory: $(pwd)"
echo "Launch time: $(date)"
echo "========================================="

# Launch diagnostic module first
echo "Launching diagnostic module..."
if [ -n "$MEMORY_FLAG" ]; then
    DIAGNOSTIC_JOB=$(sbatch --parsable "${PIPELINE_DIR}/submit_viral_diagnostic.sh" "$R1" "$R2" "$ACCESSION" "$SAMPLE_NAME" "$THREADS")
else
    DIAGNOSTIC_JOB=$(sbatch --parsable "${PIPELINE_DIR}/submit_viral_diagnostic.sh" "$R1" "$R2" "$ACCESSION" "$SAMPLE_NAME" "$THREADS")
fi

echo "  Diagnostic job submitted: $DIAGNOSTIC_JOB"

# Launch main pipeline
echo "Launching main viral pipeline..."
if [ -n "$MEMORY_FLAG" ]; then
    # Use the extreme submit script for memory flags
    PIPELINE_JOB=$(sbatch --parsable "${PIPELINE_DIR}/submit_viral_pipeline_extreme.sh" "$R1" "$R2" "$ACCESSION" "$THREADS" "$MEMORY_FLAG")
else
    # Use standard submit script
    PIPELINE_JOB=$(sbatch --parsable "${PIPELINE_DIR}/submit_viral_pipeline.sh" "$R1" "$R2" "$ACCESSION" "$THREADS")
fi

echo "  Pipeline job submitted: $PIPELINE_JOB"

echo ""
echo "========================================="
echo "JOBS LAUNCHED SUCCESSFULLY"
echo "========================================="
echo "Diagnostic Job ID: $DIAGNOSTIC_JOB"
echo "Pipeline Job ID: $PIPELINE_JOB"
echo ""
echo "Monitor job status with:"
echo "  squeue -u \$USER"
echo "  squeue -j $DIAGNOSTIC_JOB,$PIPELINE_JOB"
echo ""
echo "Check job outputs:"
echo "  # Diagnostic module"
echo "  tail -f viral_diagnostic_${DIAGNOSTIC_JOB}.out"
echo "  # Main pipeline"
echo "  tail -f viral_pipeline_${PIPELINE_JOB}.out"
echo ""
echo "Expected outputs:"
echo "  # Diagnostic results"
echo "  ./diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_diagnostic_report.txt"
echo "  ./diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_viral_blast.tsv"
echo "  # Pipeline results"
echo "  ./cleaned_seqs/variants/${SAMPLE_NAME}*.snpEFF.ann.tsv"
echo ""
echo "Both jobs will complete independently and can be analyzed together."
echo "========================================="

# Create a simple status check script
cat > "check_${SAMPLE_NAME}_status.sh" << EOF
#!/bin/bash
echo "========================================="
echo "STATUS CHECK FOR SAMPLE: $SAMPLE_NAME"
echo "========================================="
echo "Job Status:"
squeue -j $DIAGNOSTIC_JOB,$PIPELINE_JOB 2>/dev/null || echo "Jobs completed or not found"
echo ""
echo "Diagnostic Results:"
if [ -f "./diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_diagnostic_report.txt" ]; then
    echo "  ✓ Diagnostic report available"
    echo "  Preview:"
    grep -E "(Mapping Percentage|Total Contigs)" "./diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_diagnostic_report.txt" 2>/dev/null || echo "    Report parsing in progress..."
else
    echo "  ⏳ Diagnostic analysis in progress..."
fi
echo ""
echo "Pipeline Results:"
if ls ./cleaned_seqs/variants/*${SAMPLE_NAME}*.snpEFF.ann.tsv 2>/dev/null; then
    echo "  ✓ Pipeline variants available"
    VARIANT_COUNT=\$(wc -l < \$(ls ./cleaned_seqs/variants/*${SAMPLE_NAME}*.snpEFF.ann.tsv | head -1) 2>/dev/null || echo "0")
    echo "  Variants found: \$VARIANT_COUNT"
else
    echo "  ⏳ Pipeline analysis in progress..."
fi
echo ""
echo "Run this script again to check status: ./check_${SAMPLE_NAME}_status.sh"
echo "========================================="
EOF

chmod +x "check_${SAMPLE_NAME}_status.sh"

echo "Status check script created: ./check_${SAMPLE_NAME}_status.sh"
echo ""