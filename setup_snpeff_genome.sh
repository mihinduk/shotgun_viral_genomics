#!/bin/bash
# This script helps set up a genome in snpEff when the automatic approach fails

# Check arguments
if [ $# -lt 2 ]; then
  echo "Usage: $0 <accession> <snpeff_jar_path> [genbank_file]"
  echo "Example: $0 MK573243.1 /home/kathiem/snpEff/snpEff.jar ./MK573243.1.gb"
  exit 1
fi

ACCESSION=$1
SNPEFF_JAR=$2
GB_FILE=$3

# Directory where snpEff is installed
SNPEFF_DIR=$(dirname "$SNPEFF_JAR")
CONFIG_FILE="$SNPEFF_DIR/snpEff.config"

# If no GenBank file was provided, download it
if [ -z "$GB_FILE" ]; then
  echo "No GenBank file provided. Downloading from NCBI..."
  GB_FILE="./${ACCESSION}.gb"
  efetch -db nucleotide -id "$ACCESSION" -format gb > "$GB_FILE"
  if [ $? -ne 0 ]; then
    echo "Failed to download GenBank file. Trying alternative approach..."
    wget -O "$GB_FILE" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACCESSION}&rettype=gb&retmode=text"
  fi
fi

# Extract organism information if available
ORGANISM=$(grep -A1 "ORGANISM" "$GB_FILE" | tail -1 | tr -d '[:space:]' || echo "$ACCESSION")
STRAIN=$(grep -o 'strain="[^"]*"' "$GB_FILE" | head -1 | sed 's/strain="\([^"]*\)"/\1/' || echo "")

# Create an abbreviated name
if [ -n "$STRAIN" ]; then
  # Extract first letter of each word in organism
  ABBR=$(echo "$ORGANISM" | tr ' ' '\n' | grep -v '^$' | grep -o '^.' | tr -d '\n')
  GENOME_NAME="${ABBR}-${STRAIN}"
else
  GENOME_NAME="$ACCESSION"
fi

echo "Using genome name: $GENOME_NAME"

# Create local config file
LOCAL_CONFIG="./snpEff_${ACCESSION}.config"
cat > "$LOCAL_CONFIG" << EOF
# ${ACCESSION} genome configuration
data.dir = ./data/
${ACCESSION}.genome: ${GENOME_NAME}
${ACCESSION}.chromosomes: ${ACCESSION}
${ACCESSION}.codonTable: Standard
EOF

echo "Created local config: $LOCAL_CONFIG"

# Download the genome sequence if needed
echo "Downloading genome sequence..."
FASTA_FILE="./${ACCESSION}.fasta"
efetch -db nucleotide -id "$ACCESSION" -format fasta > "$FASTA_FILE"
if [ $? -ne 0 ]; then
  echo "Failed to download FASTA file. Trying alternative approach..."
  wget -O "$FASTA_FILE" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACCESSION}&rettype=fasta&retmode=text"
fi

# Create data directory structure for snpEff
echo "Setting up data directory..."
mkdir -p "./data/${ACCESSION}"
cp "$FASTA_FILE" "./data/${ACCESSION}/sequences.fa"

# Convert GenBank to GFF3 if available
if command -v biopython &> /dev/null || python -c "import Bio" &> /dev/null; then
  echo "Biopython found, generating GFF3 file..."
  if [ -f "convert_gb_to_gff_advanced.py" ]; then
    python3 convert_gb_to_gff_advanced.py "$GB_FILE" "./data/${ACCESSION}/genes.gff" --verbose
  else
    echo "convert_gb_to_gff_advanced.py not found, using simple approach..."
    # Create a very simple GFF3 file
    cat > "./data/${ACCESSION}/genes.gff" << EOF
##gff-version 3
##sequence-region ${ACCESSION} 1 10000
${ACCESSION}\tGenBank\tregion\t1\t10000\t.\t+\t.\tID=region_1
${ACCESSION}\tGenBank\tCDS\t1\t10000\t.\t+\t0\tID=CDS_1;Name=viral_protein;product=viral protein
EOF
  fi
else
  echo "Biopython not found, creating simple GFF3 file..."
  # Create a very simple GFF3 file
  cat > "./data/${ACCESSION}/genes.gff" << EOF
##gff-version 3
##sequence-region ${ACCESSION} 1 10000
${ACCESSION}\tGenBank\tregion\t1\t10000\t.\t+\t.\tID=region_1
${ACCESSION}\tGenBank\tCDS\t1\t10000\t.\t+\t0\tID=CDS_1;Name=viral_protein;product=viral protein
EOF
fi

# Build the database
echo "Building snpEff database..."
export SNPEFF_CONFIG="$LOCAL_CONFIG"
java -jar "$SNPEFF_JAR" build -gff3 -v -noCheckProtein -noCheckCds -noLog -treatAllAsProteinCoding "$ACCESSION"

# Verify the database was built
if [ $? -eq 0 ]; then
  echo "Success! The snpEff database for $ACCESSION was built successfully."
  echo ""
  echo "To use this database with viral_pipeline.py, you need to:"
  echo "1. Run your command with SNPEFF_CONFIG set:"
  echo "   export SNPEFF_CONFIG=\"$LOCAL_CONFIG\""
  echo "   ./shotgun_viral_genomics/viral_pipeline.py --r1 \"NovaSeq_N917*_R1*.fastq.gz\" --r2 \"NovaSeq_N917*_R2*.fastq.gz\" \\"
  echo "       --accession $ACCESSION \\"
  echo "       --threads 4 \\"
  echo "       --snpeff-jar $SNPEFF_JAR \\"
  echo "       --add-to-snpeff"
  echo ""
  echo "You can also use this database directly with snpEff:"
  echo "   java -jar $SNPEFF_JAR -c $LOCAL_CONFIG $ACCESSION your_variants.vcf"
else
  echo "Failed to build the snpEff database. Check the error messages above."
fi