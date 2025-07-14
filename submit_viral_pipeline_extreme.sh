#!/bin/bash
# Submit script for viral pipeline with dynamic resource allocation
# Usage: sbatch submit_viral_pipeline_extreme.sh <R1> <R2> <accession> <threads> [--large-files|--extremely-large-files]

# Check arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> <threads> [--large-files|--extremely-large-files]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4 --large-files"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4 --extremely-large-files"
    echo ""
    echo "Memory allocation:"
    echo "  Standard: 32GB (default)"
    echo "  --large-files: 64GB"
    echo "  --extremely-large-files: 256GB on high-memory partition"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=$4
LARGE_FILES_FLAG=""

# Set default resources
MEMORY="32G"
PARTITION="general"
TIME="4:00:00"

# Check if memory flags are provided
if [ $# -eq 5 ]; then
    if [ "$5" == "--large-files" ]; then
        LARGE_FILES_FLAG="--large-files"
        MEMORY="64G"
        TIME="8:00:00"
        echo "Using --large-files flag: 64GB memory allocation"
    elif [ "$5" == "--extremely-large-files" ]; then
        LARGE_FILES_FLAG="--extremely-large-files"
        MEMORY="256G"
        PARTITION="general"
        TIME="12:00:00"
        echo "Using --extremely-large-files flag: 256GB memory on general partition"
    fi
fi

# Submit the job with dynamic resources
sbatch --job-name=viral_pipeline \
       --output=viral_pipeline_%j.out \
       --error=viral_pipeline_%j.err \
       --time=$TIME \
       --mem=$MEMORY \
       --cpus-per-task=$THREADS \
       --partition=$PARTITION \
       --wrap="
# Change to the directory where the script was submitted from
cd ${PWD}

echo '========================================='
echo 'SLURM Job ID: '\$SLURM_JOB_ID
echo 'Job started at: '\$(date)
echo 'Working directory: '\$(pwd)
echo 'Memory allocated: $MEMORY'
echo 'Partition: $PARTITION'
echo 'Time limit: $TIME'
echo '========================================='

# Set up environment
echo 'Setting up HTCF environment...'
# Note: The wrapper script handles conda/mamba activation internally

# Set paths
PIPELINE_DIR='/scratch/sahlab/kathie/Diamond_test/shotgun_viral_genomics'
WRAPPER_SCRIPT=\${PIPELINE_DIR}/run_pipeline_htcf_consolidated.sh

# Validate input files exist
if [ ! -f '$R1' ]; then
    echo 'Error: R1 file not found: $R1'
    exit 1
fi

if [ ! -f '$R2' ]; then
    echo 'Error: R2 file not found: $R2'
    exit 1
fi

echo 'Running viral pipeline with:'
echo '  R1: $R1'
echo '  R2: $R2'
echo '  Accession: $ACCESSION'
echo '  Threads: $THREADS'
echo '  Memory flag: $LARGE_FILES_FLAG'
echo '  Working directory: '\$(pwd)
echo ''

# Run the pipeline wrapper script
if [ -n '$LARGE_FILES_FLAG' ]; then
    # Run with memory flag
    bash \$WRAPPER_SCRIPT '$R1' '$R2' '$ACCESSION' '$THREADS' '$LARGE_FILES_FLAG'
else
    # Run with standard settings
    bash \$WRAPPER_SCRIPT '$R1' '$R2' '$ACCESSION' '$THREADS'
fi

PIPELINE_EXIT_CODE=\$?

echo ''
echo '========================================='
echo 'Job completed at: '\$(date)
if [ \$PIPELINE_EXIT_CODE -eq 0 ]; then
    echo 'Pipeline completed successfully!'
else
    echo 'Pipeline failed with exit code: '\$PIPELINE_EXIT_CODE
fi
echo '========================================='

exit \$PIPELINE_EXIT_CODE
"