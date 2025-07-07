#!/bin/bash
# This script directly builds a snpEff database without any fancy config file handling

# Check arguments
if [ $# -lt 2 ]; then
  echo "Usage: $0 <accession> <snpeff_jar_path>"
  echo "Example: $0 MK573243.1 /home/kathiem/snpEff/snpEff.jar"
  exit 1
fi

ACCESSION=$1
SNPEFF_JAR=$2
SNPEFF_DIR=$(dirname "$SNPEFF_JAR")
CONFIG_FILE="$SNPEFF_DIR/snpEff.config"

echo "Creating directories and files for snpEff database..."

# Create sequences directory for your genome
mkdir -p "$SNPEFF_DIR/data/$ACCESSION"

# Download genome
echo "Downloading genome..."
efetch -db nucleotide -id "$ACCESSION" -format fasta > "$SNPEFF_DIR/data/$ACCESSION/sequences.fa"
efetch -db nucleotide -id "$ACCESSION" -format gb > "$SNPEFF_DIR/data/$ACCESSION/$ACCESSION.gb"

# Create backup of original config
cp "$CONFIG_FILE" "${CONFIG_FILE}.bak"

# Add to config file - trying a format that works with your version
echo "Adding genome to snpEff configuration..."
printf "\n# $ACCESSION\n$ACCESSION.genome : $ACCESSION\n" >> "$CONFIG_FILE"

# Verify it was added properly
echo "Checking config file update:"
tail -3 "$CONFIG_FILE"

# Run the build command directly
echo "Building database using GenBank file directly (most reliable approach)..."
cd "$SNPEFF_DIR"
java -jar snpEff.jar build -genbank -v "$ACCESSION"

if [ $? -eq 0 ]; then
  echo "Success! The snpEff database for $ACCESSION was built successfully."
  echo ""
  echo "Now you can run your pipeline without --add-to-snpeff flag:"
  echo "./shotgun_viral_genomics/viral_pipeline.py --r1 \"NovaSeq_N917*_R1*.fastq.gz\" --r2 \"NovaSeq_N917*_R2*.fastq.gz\" \\"
  echo "    --accession $ACCESSION \\"
  echo "    --threads 4 \\"
  echo "    --snpeff-jar $SNPEFF_JAR \\"
  echo "    --java-path \"/usr/lib/jvm/java-11-openjdk-11.0.20.0.8-1.el7_9.x86_64/bin/java\""
else
  echo "Failed to build database. Let's try an alternative approach..."
  
  # If the first method fails, try with minimal config 
  echo "Trying alternative method..."
  
  # Create a minimal config file
  MINIMAL_CONFIG="/tmp/snpEff_minimal.config"
  cat > "$MINIMAL_CONFIG" << EOF
# Minimal snpEff config file

# Data directories
data.dir = $SNPEFF_DIR/data/

# Genome
$ACCESSION.genome : $ACCESSION
EOF

  # Try building with minimal config
  echo "Building with minimal config..."
  java -jar "$SNPEFF_JAR" -c "$MINIMAL_CONFIG" build -genbank -v "$ACCESSION"
  
  if [ $? -eq 0 ]; then
    echo "Success with alternative method!"
    echo "To use this database, run your pipeline with this environment variable:"
    echo "export SNPEFF_CONFIG=\"$MINIMAL_CONFIG\""
    echo "./shotgun_viral_genomics/viral_pipeline.py --r1 \"NovaSeq_N917*_R1*.fastq.gz\" --r2 \"NovaSeq_N917*_R2*.fastq.gz\" \\"
    echo "    --accession $ACCESSION \\"
    echo "    --threads 4 \\"
    echo "    --snpeff-jar $SNPEFF_JAR \\"
    echo "    --java-path \"/usr/lib/jvm/java-11-openjdk-11.0.20.0.8-1.el7_9.x86_64/bin/java\""
  else
    echo "Both approaches failed."
    echo "You might need to manually configure snpEff following these instructions:"
    echo "1. Create directory: $SNPEFF_DIR/data/$ACCESSION/"
    echo "2. Download files: "
    echo "   - Genome: efetch -db nucleotide -id $ACCESSION -format fasta > $SNPEFF_DIR/data/$ACCESSION/sequences.fa"
    echo "   - GenBank: efetch -db nucleotide -id $ACCESSION -format gb > $SNPEFF_DIR/data/$ACCESSION/genes.gbk"
    echo "3. Manually edit $CONFIG_FILE to add: $ACCESSION.genome : $ACCESSION"
    echo "4. Run build: java -jar $SNPEFF_JAR build -genbank $ACCESSION"
  fi
fi