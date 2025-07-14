#!/bin/bash
# Build local BLAST database for microbial contamination detection (no mammals)
# This version excludes large mammalian genomes for faster builds and smaller database size
# Version 2: Fixed file combination and error handling

set -euo pipefail

# Configuration
DB_NAME="microbial_contaminants_v2"
WORK_DIR="/ref/sahlab/data/microbes_db_v2_build"
FINAL_DB_DIR="/ref/sahlab/data/microbes_db_v2"
THREADS=4

# Create working directories
echo "========================================="
echo "Microbial Contamination Database Build v2"
echo "========================================="
echo "Started at: $(date)"
echo "Creating directories..."
mkdir -p "$WORK_DIR"
mkdir -p "$FINAL_DB_DIR"
cd "$WORK_DIR"

# Function to download with retry
download_with_retry() {
    local url=$1
    local output=$2
    local max_attempts=3
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        echo "  Downloading $output (attempt $attempt/$max_attempts)..."
        if wget -q --show-progress "$url" -O "$output"; then
            echo "  ✓ Downloaded successfully"
            return 0
        else
            echo "  ✗ Download failed"
            rm -f "$output"
            attempt=$((attempt + 1))
            [ $attempt -le $max_attempts ] && sleep 5
        fi
    done
    return 1
}

# Step 1: Download microbial genomes (no mammals)
echo ""
echo "Step 1: Downloading microbial genomes"
echo "====================================="

# Create accession list with direct FTP URLs - MICROBES ONLY
declare -A GENOMES=(
    ["Mycoplasma_hyorhinis"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/575/GCF_000195575.1_ASM19557v1/GCF_000195575.1_ASM19557v1_genomic.fna.gz"
    ["Pseudomonas_aeruginosa_PAO1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz"
    ["Candida_albicans_SC5314"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz"
    ["Saccharomyces_cerevisiae_S288C"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
    ["Cryptococcus_neoformans_JEC21"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/245/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.fna.gz"
    ["Escherichia_coli_K12_MG1655"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    ["Staphylococcus_aureus_NCTC8325"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz"
)

# Download each genome
DOWNLOADED_FILES=()
for organism in "${!GENOMES[@]}"; do
    echo "Downloading $organism..."
    if download_with_retry "${GENOMES[$organism]}" "${organism}.fna.gz"; then
        gunzip -f "${organism}.fna.gz"
        if [ -f "${organism}.fna" ]; then
            echo "  ✓ Extracted ${organism}.fna"
            DOWNLOADED_FILES+=("${organism}.fna")
        else
            echo "  ✗ Failed to extract ${organism}.fna"
        fi
    else
        echo "  ✗ Failed to download $organism - continuing with remaining genomes"
    fi
done

# Step 2: Download viral sequences
echo ""
echo "Step 2: Downloading viral RefSeq genomes"
echo "========================================"
VIRAL_URL="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
if download_with_retry "$VIRAL_URL" "viral.1.1.genomic.fna.gz"; then
    echo "Extracting viral sequences..."
    gunzip -f viral.1.1.genomic.fna.gz
    if [ -f "viral.1.1.genomic.fna" ]; then
        echo "  ✓ Extracted viral sequences"
        DOWNLOADED_FILES+=("viral.1.1.genomic.fna")
    else
        echo "  ✗ Failed to extract viral sequences"
    fi
else
    echo "  ✗ Failed to download viral sequences - continuing without them"
fi

# Step 3: Combine sequences and create BLAST database
echo ""
echo "Step 3: Building BLAST database"
echo "==============================="

# Check if we have files to combine
if [ ${#DOWNLOADED_FILES[@]} -eq 0 ]; then
    echo "Error: No FASTA files were successfully downloaded"
    exit 1
fi

# Combine all successfully downloaded FASTA files
echo "Combining ${#DOWNLOADED_FILES[@]} FASTA files..."
echo "Files to combine: ${DOWNLOADED_FILES[*]}"

# Use explicit file list to avoid globbing issues
cat "${DOWNLOADED_FILES[@]}" > "${DB_NAME}_combined.fna"

# Verify combined file was created and has content
if [ ! -f "${DB_NAME}_combined.fna" ] || [ ! -s "${DB_NAME}_combined.fna" ]; then
    echo "Error: Failed to create combined FASTA file"
    exit 1
fi

# Count sequences
TOTAL_SEQS=$(grep -c "^>" "${DB_NAME}_combined.fna" || echo "0")
echo "Total sequences in database: $TOTAL_SEQS"

if [ "$TOTAL_SEQS" -eq 0 ]; then
    echo "Error: No sequences found in combined file"
    exit 1
fi

# Build BLAST database
echo "Building BLAST database..."
# Activate BLAST environment
eval "$(conda shell.bash hook)"
conda activate Blast

makeblastdb -in "${DB_NAME}_combined.fna" \
            -dbtype nucl \
            -out "${FINAL_DB_DIR}/${DB_NAME}" \
            -title "Microbial Contamination Database v2" \
            -parse_seqids

# Check if database was created successfully
if [ -f "${FINAL_DB_DIR}/${DB_NAME}.nhr" ] || [ -f "${FINAL_DB_DIR}/${DB_NAME}.00.nhr" ]; then
    echo "  ✓ BLAST database created successfully"
    
    # Create database info file
    echo "Creating database info file..."
    cat > "${FINAL_DB_DIR}/${DB_NAME}_info.txt" << EOF
Microbial Contamination BLAST Database v2 (No Mammals)
=====================================================
Created: $(date)
Total sequences: $TOTAL_SEQS

Organisms successfully included:
$(for file in "${DOWNLOADED_FILES[@]}"; do
    if [ -f "$file" ]; then
        count=$(grep -c "^>" "$file" 2>/dev/null || echo "0")
        basename_file=$(basename "$file" .fna)
        echo "  - $basename_file: $count sequences"
    fi
done)

Database location: ${FINAL_DB_DIR}/${DB_NAME}
Database name: ${DB_NAME}

Note: This database excludes mammalian genomes to reduce size and build time.

Usage:
blastn -query your_contigs.fa -db ${FINAL_DB_DIR}/${DB_NAME} -outfmt 6 -max_target_seqs 5 -evalue 1e-10 -num_threads 4
EOF

else
    echo "  ✗ Failed to create BLAST database"
    exit 1
fi

# Clean up working directory
echo ""
echo "Cleaning up temporary files..."
rm -rf "$WORK_DIR"
echo "  ✓ Temporary files removed"

echo ""
echo "========================================="
echo "Microbial database v2 build completed!"
echo "========================================="
echo "Database location: ${FINAL_DB_DIR}/${DB_NAME}"
echo "Info file: ${FINAL_DB_DIR}/${DB_NAME}_info.txt"
echo "Total sequences: $TOTAL_SEQS"
echo ""
echo "To use this database:"
echo "blastn -query contigs.fa -db ${FINAL_DB_DIR}/${DB_NAME} -outfmt 6 -max_target_seqs 5 -evalue 1e-10 -num_threads 4"
echo "========================================="