#!/bin/bash
# Build local BLAST database for cell culture contamination detection
# This script downloads representative genomes and creates a custom BLAST database

set -euo pipefail

# Configuration
DB_NAME="cell_culture_contaminants"
WORK_DIR="/ref/sahlab/data/contamination_db_build"
FINAL_DB_DIR="/ref/sahlab/data/contamination_db"
THREADS=4

# Create working directories
echo "========================================="
echo "Cell Culture Contamination Database Build"
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

# Step 1: Download representative genomes
echo ""
echo "Step 1: Downloading representative genomes"
echo "=========================================="

# Create accession list with direct FTP URLs
declare -A GENOMES=(
    ["Mycoplasma_hyorhinis"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/575/GCF_000195575.1_ASM19557v1/GCF_000195575.1_ASM19557v1_genomic.fna.gz"
    ["Pseudomonas_aeruginosa_PAO1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz"
    ["Candida_albicans_SC5314"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz"
    ["Saccharomyces_cerevisiae_S288C"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
    ["Cryptococcus_neoformans_JEC21"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/245/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.fna.gz"
    ["Escherichia_coli_K12_MG1655"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    ["Staphylococcus_aureus_NCTC8325"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz"
    ["Homo_sapiens_rRNA"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz"
    ["Chlorocebus_sabaeus"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/252/025/GCF_015252025.1_Vero_WHO_p1.0/GCF_015252025.1_Vero_WHO_p1.0_genomic.fna.gz"
    ["Mesocricetus_auratus"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/785/GCF_017639785.1_BCM_Maur_2.0/GCF_017639785.1_BCM_Maur_2.0_genomic.fna.gz"
)

# Download each genome
for organism in "${!GENOMES[@]}"; do
    echo "Downloading $organism..."
    if download_with_retry "${GENOMES[$organism]}" "${organism}.fna.gz"; then
        gunzip -f "${organism}.fna.gz"
        echo "  ✓ Extracted ${organism}.fna"
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
    echo "  ✓ Extracted viral sequences"
else
    echo "  ✗ Failed to download viral sequences - continuing without them"
fi

# Step 3: Combine sequences and create BLAST database
echo ""
echo "Step 3: Building BLAST database"
echo "==============================="

# Combine all FASTA files
echo "Combining all sequences..."
cat *.fna > "${DB_NAME}_combined.fna" 2>/dev/null || {
    echo "Error: No FASTA files found to combine"
    exit 1
}

# Count sequences
TOTAL_SEQS=$(grep -c "^>" "${DB_NAME}_combined.fna" || echo "0")
echo "Total sequences in database: $TOTAL_SEQS"

# Build BLAST database
echo "Building BLAST database..."
makeblastdb -in "${DB_NAME}_combined.fna" \
            -dbtype nucl \
            -out "${FINAL_DB_DIR}/${DB_NAME}" \
            -title "Cell Culture Contamination Database" \
            -parse_seqids

# Check if database was created successfully
if [ -f "${FINAL_DB_DIR}/${DB_NAME}.nhr" ]; then
    echo "  ✓ BLAST database created successfully"
else
    echo "  ✗ Failed to create BLAST database"
    exit 1
fi

# Create database info file
echo "Creating database info file..."
cat > "${FINAL_DB_DIR}/${DB_NAME}_info.txt" << EOF
Cell Culture Contamination BLAST Database
=========================================
Created: $(date)
Total sequences: $TOTAL_SEQS

Organisms included:
$(for org in "${!GENOMES[@]}"; do
    if [ -f "${org}.fna" ]; then
        count=$(grep -c "^>" "${org}.fna" 2>/dev/null || echo "0")
        echo "  - $org: $count sequences"
    fi
done)

Viral sequences: $(grep -c "^>" viral.1.1.genomic.fna 2>/dev/null || echo "0")

Database location: ${FINAL_DB_DIR}/${DB_NAME}
Database name: ${DB_NAME}

Usage:
blastn -query your_contigs.fa -db ${FINAL_DB_DIR}/${DB_NAME} -outfmt 6 -max_target_seqs 5 -evalue 1e-10 -num_threads 4
EOF

# Clean up working directory (optional)
echo ""
echo "Cleaning up temporary files..."
# In batch mode, always clean up
rm -rf "$WORK_DIR"
echo "  ✓ Temporary files removed"

echo ""
echo "========================================="
echo "Database build completed!"
echo "========================================="
echo "Database location: ${FINAL_DB_DIR}/${DB_NAME}"
echo "Info file: ${FINAL_DB_DIR}/${DB_NAME}_info.txt"
echo ""
echo "To use this database:"
echo "blastn -query contigs.fa -db ${FINAL_DB_DIR}/${DB_NAME} -outfmt 6 -max_target_seqs 5 -evalue 1e-10 -num_threads 4"
echo "========================================="