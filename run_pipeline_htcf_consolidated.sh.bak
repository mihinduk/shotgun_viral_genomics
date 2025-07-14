#!/bin/bash
# Consolidated viral pipeline for HTCF with automatic snpEff database building

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    echo ""
    echo "This script will:"
    echo "  1. Check if snpEff database exists for the accession"
    echo "  2. Build the database if needed (viral-friendly settings)"
    echo "  3. Run the complete viral genomics pipeline"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=${4:-4}  # Default to 4 threads if not specified

# Set up environment - use mamba run instead of activate
echo "Setting up HTCF environment..."
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"

# Set paths
PIPELINE_DIR="$(cd "$(dirname "$0")"; pwd)"
SNPEFF_JAR="/home/mihindu/software/snpEff/snpEff.jar"
SNPEFF_DIR="/home/mihindu/software/snpEff"
JAVA_PATH="java"

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
echo "  Working directory: $(pwd)"

# Function to build snpEff database
build_snpeff_database() {
    local acc=$1
    echo "Building snpEff database for $acc..."
    
    # Create database directory
    mkdir -p "$SNPEFF_DIR/data/$acc"
    
    # Download GenBank file to variants folder (for collaborators)
    mkdir -p ./cleaned_seqs/variants
    if [ ! -f "./cleaned_seqs/variants/$acc.gb" ]; then
        echo "Downloading GenBank file for $acc..."
        $MAMBA_CMD efetch -db nucleotide -id "$acc" -format gb > "./cleaned_seqs/variants/$acc.gb"
        if [ $? -ne 0 ] || [ ! -s "./cleaned_seqs/variants/$acc.gb" ]; then
            echo "Error: Failed to download GenBank file for $acc"
            return 1
        fi
        echo "GenBank file downloaded successfully ($(wc -l < "./cleaned_seqs/variants/$acc.gb") lines)"
    fi
    
    # Copy GenBank file to snpEff database directory
    cp "./cleaned_seqs/variants/$acc.gb" "$SNPEFF_DIR/data/$acc/genes.gbk"
    
    # Download FASTA file if not exists
    if [ ! -f "$SNPEFF_DIR/data/$acc/sequences.fa" ]; then
        echo "Downloading FASTA file for $acc..."
        $MAMBA_CMD efetch -db nucleotide -id "$acc" -format fasta > "$SNPEFF_DIR/data/$acc/sequences.fa"
        if [ $? -ne 0 ] || [ ! -s "$SNPEFF_DIR/data/$acc/sequences.fa" ]; then
            echo "Error: Failed to download FASTA file for $acc"
            return 1
        fi
        echo "FASTA file downloaded successfully"
    fi
    
    # Add to snpEff config BEFORE building - proper syntax
    if ! grep -q "$acc.genome" "$SNPEFF_DIR/snpEff.config"; then
        echo "Adding $acc to snpEff configuration..."
        echo "" >> "$SNPEFF_DIR/snpEff.config"
        echo "# $acc" >> "$SNPEFF_DIR/snpEff.config"
        echo "$acc.genome : $acc" >> "$SNPEFF_DIR/snpEff.config"
    fi
    
    # Verify config was updated
    echo "Verifying config file contains $acc..."
    if grep -q "$acc.genome" "$SNPEFF_DIR/snpEff.config"; then
        echo "Config file updated successfully"
    else
        echo "Error: Failed to update config file"
        return 1
    fi
    
    # Build database with viral-friendly settings using mamba run
    echo "Building snpEff database with viral-friendly settings..."
    cd "$SNPEFF_DIR"
    $MAMBA_CMD java -jar snpEff.jar build -v -noCheckProtein -noCheckCds -noLog "$acc"
    local build_result=$?
    
    # Return to original directory
    cd - > /dev/null
    
    if [ $build_result -eq 0 ]; then
        echo "Successfully built snpEff database for $acc"
        return 0
    else
        echo "Error: Failed to build snpEff database for $acc"
        return 1
    fi
}

# Check if snpEff database exists using mamba run
echo "Checking snpEff database for $ACCESSION..."
if $MAMBA_CMD java -jar "$SNPEFF_JAR" databases | grep -q "$ACCESSION"; then
    echo "snpEff database found for $ACCESSION"
    
    # Check if the database is actually built (has snpEffectPredictor.bin)
    if [ ! -f "$SNPEFF_DIR/data/$ACCESSION/snpEffectPredictor.bin" ]; then
        echo "Database exists but not built. Building now..."
        build_snpeff_database "$ACCESSION"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to build snpEff database"
            exit 1
        fi
    fi
else
    echo "snpEff database not found for $ACCESSION. Building now..."
    build_snpeff_database "$ACCESSION"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to build snpEff database"
        exit 1
    fi
fi

echo "Running the viral pipeline..."

# Run the pipeline without --add-to-snpeff since we handled it above
$MAMBA_CMD python ${PIPELINE_DIR}/viral_pipeline.py \
    --r1 "$R1" \
    --r2 "$R2" \
    --accession "$ACCESSION" \
    --threads $THREADS \
    --snpeff-jar "$SNPEFF_JAR" \
    --java-path "$JAVA_PATH"

PIPELINE_EXIT_CODE=$?

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "Pipeline completed successfully!"
    echo ""
    echo "Output files:"
    find ./cleaned_seqs/variants -name "*.tsv" -o -name "*.vcf" 2>/dev/null | sort
    echo ""
    echo "Key results:"
    if ls ./cleaned_seqs/variants/*snpEFF.ann.tsv >/dev/null 2>&1; then
        echo "  Annotation TSV: $(ls ./cleaned_seqs/variants/*snpEFF.ann.tsv | head -1)"
        echo "  Total annotated variants: $(wc -l < $(ls ./cleaned_seqs/variants/*snpEFF.ann.tsv | head -1))"
    fi
else
    echo "Pipeline failed with exit code $PIPELINE_EXIT_CODE"
    exit $PIPELINE_EXIT_CODE
fi