#!/bin/bash
# Setup script for viral assembly environment

echo "Creating viral assembly environment..."

# Load anaconda
source /ref/sahlab/software/anaconda3/bin/activate

# Create new environment with assemblers
conda create -n viral_assembly -y python=3.9

# Activate it
conda activate viral_assembly

# Install assemblers and tools
echo "Installing SPAdes, MEGAHIT, and supporting tools..."
conda install -c bioconda -c conda-forge -y \
    spades \
    megahit \
    seqkit \
    blast \
    biopython \
    pandas

echo "Environment setup complete!"
echo "To use: conda activate viral_assembly"