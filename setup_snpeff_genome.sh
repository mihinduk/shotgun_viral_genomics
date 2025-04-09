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

# Genome ${ACCESSION}
${ACCESSION}.genome = ${GENOME_NAME}
${ACCESSION}.chromosomes = ${ACCESSION}
${ACCESSION}.codonTable = Standard
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
python3 -c "import Bio" 2>/dev/null
if [ $? -eq 0 ]; then
  echo "Biopython found, generating GFF3 file..."
  
  # Look for converter script in shotgun_viral_genomics directory
  CONVERTER=""
  for script in "./convert_gb_to_gff_advanced.py" "./shotgun_viral_genomics/convert_gb_to_gff_advanced.py" "./convert_gb_to_gff_simple.py" "./shotgun_viral_genomics/convert_gb_to_gff_simple.py"; do
    if [ -f "$script" ]; then
      CONVERTER="$script"
      echo "Found converter script: $CONVERTER"
      break
    fi
  done
  
  if [ -n "$CONVERTER" ]; then
    python3 "$CONVERTER" "$GB_FILE" "./data/${ACCESSION}/genes.gff" --verbose
  else
    echo "No converter script found, creating one..."
    # Create a simple converter script
    cat > "./convert_gb_to_gff_simple.py" << 'EOF'
#!/usr/bin/env python3
import sys
import re
from Bio import SeqIO

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 convert_gb_to_gff_simple.py input.gb output.gff3")
        sys.exit(1)
    
    gb_file = sys.argv[1]
    gff_file = sys.argv[2]
    
    # Read GenBank file
    record = SeqIO.read(gb_file, "genbank")
    
    # Create GFF3 file
    with open(gff_file, 'w') as out:
        out.write('##gff-version 3\n')
        out.write(f'##sequence-region {record.id} 1 {len(record.seq)}\n')
        
        # Add source feature
        out.write(f'{record.id}\tGenBank\tregion\t1\t{len(record.seq)}\t.\t+\t.\tID=region_1\n')
        
        # Process features
        feature_counter = {}
        for feature in record.features:
            if feature.type not in feature_counter:
                feature_counter[feature.type] = 0
            feature_counter[feature.type] += 1
            
            # Skip source feature (already added)
            if feature.type == 'source':
                continue
                
            # Start and end positions (1-based for GFF)
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            
            # Strand
            strand = '+' if feature.location.strand == 1 else '-'
            
            # Phase for CDS features
            phase = '0'
            if feature.type == 'CDS' and 'codon_start' in feature.qualifiers:
                phase = str(int(feature.qualifiers['codon_start'][0]) - 1)
            else:
                phase = '.'
            
            # Create ID
            feature_id = f"{feature.type.lower()}_{feature_counter[feature.type]}"
            if feature.type == 'CDS' and 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0]
                clean_product = re.sub(r'[^\w\s]', '', product).replace(' ', '_')
                feature_id = f"{clean_product}_{feature_counter[feature.type]}"
            
            # Create attributes
            attributes = [f"ID={feature_id}"]
            
            # Add Name and product if available
            if 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0]
                attributes.append(f"Name={product}")
                attributes.append(f"product={product}")
            
            # Add other qualifiers
            for key, values in feature.qualifiers.items():
                if key not in ['product'] and key != 'translation':
                    for value in values:
                        value = value.replace(';', '%3B').replace('=', '%3D')
                        attributes.append(f"{key}={value}")
            
            # Write the feature
            out.write(f'{record.id}\tGenBank\t{feature.type}\t{start}\t{end}\t.\t{strand}\t{phase}\t{";".join(attributes)}\n')
        
        # Add sequence
        out.write('##FASTA\n')
        out.write(f'>{record.id}\n')
        seq_str = str(record.seq)
        for i in range(0, len(seq_str), 60):
            out.write(seq_str[i:i+60] + '\n')

if __name__ == "__main__":
    main()
EOF
    
    chmod +x "./convert_gb_to_gff_simple.py"
    python3 "./convert_gb_to_gff_simple.py" "$GB_FILE" "./data/${ACCESSION}/genes.gff"
  fi
else
  echo "Biopython not found, creating simple GFF3 file..."
  # Get sequence length from FASTA file
  SEQ_LENGTH=$(grep -v "^>" "$FASTA_FILE" | tr -d '\n' | wc -c)
  if [ -z "$SEQ_LENGTH" ] || [ "$SEQ_LENGTH" -eq 0 ]; then
    SEQ_LENGTH=10000
  fi
  
  # Create a very simple GFF3 file
  cat > "./data/${ACCESSION}/genes.gff" << EOF
##gff-version 3
##sequence-region ${ACCESSION} 1 ${SEQ_LENGTH}
${ACCESSION}\tGenBank\tregion\t1\t${SEQ_LENGTH}\t.\t+\t.\tID=region_1
${ACCESSION}\tGenBank\tCDS\t1\t${SEQ_LENGTH}\t.\t+\t0\tID=CDS_1;Name=viral_protein;product=viral protein
EOF
fi

# Build the database
echo "Building snpEff database..."
export SNPEFF_CONFIG="$LOCAL_CONFIG"
java -jar "$SNPEFF_JAR" build -gff3 -v -noCheckProtein -noCheckCds -noLog "$ACCESSION"

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