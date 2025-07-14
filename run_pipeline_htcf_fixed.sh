#\!/bin/bash
# Generic run script for viral pipeline on HTCF using mamba

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=${4:-4}  # Default to 4 threads if not specified

# Set paths
PIPELINE_DIR="$(cd "$(dirname "$0")"; pwd)"
SNPEFF_JAR="/home/mihindu/software/snpEff/snpEff.jar"
JAVA_PATH="java"  # Will use conda environment java

# Validate input files exist
if [ \! -f "$R1" ]; then
    echo "Error: R1 file not found: $R1"
    exit 1
fi

if [ \! -f "$R2" ]; then
    echo "Error: R2 file not found: $R2"
    exit 1
fi

echo "Running viral pipeline with:"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Accession: $ACCESSION"
echo "  Threads: $THREADS"
echo "  Working directory: $(pwd)"

# Run the pipeline using mamba run
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics ${PIPELINE_DIR}/viral_pipeline.py \
    --r1 "$R1" \
    --r2 "$R2" \
    --accession "$ACCESSION" \
    --threads $THREADS \
    --snpeff-jar "$SNPEFF_JAR" \
    --java-path "$JAVA_PATH" \
    --add-to-snpeff

echo "Pipeline completed\!"
