#\!/bin/bash
# Setup script for HTCF environment

# Activate conda
source /ref/sahlab/software/anaconda3/bin/activate

# Create viral_genomics environment if it doesn't exist
if \! conda env list  < /dev/null |  grep -q viral_genomics; then
    echo "Creating viral_genomics conda environment..."
    conda env create -f environment.yml
else
    echo "viral_genomics environment already exists"
fi

# Download snpEff if not present
SNPEFF_DIR="/scratch/sahlab/kathie/Diamond_test/snpEff"
if [ ! -f "$SNPEFF_DIR/snpEff.jar" ]; then
    echo "Downloading snpEff..."
    mkdir -p $SNPEFF_DIR
    cd $SNPEFF_DIR
    wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    rm snpEff_latest_core.zip
    cd -
else
    echo "snpEff already installed at $SNPEFF_DIR"
fi

echo "Setup complete!"
echo "To run the pipeline, use:"
echo "  conda activate viral_genomics"
echo "  Java path: /usr/bin/java"
echo "  snpEff path: $SNPEFF_DIR/snpEff.jar"
