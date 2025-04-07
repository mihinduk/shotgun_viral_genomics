#!/bin/bash
set -e

# Installation script for Shotgun Viral Genomics Pipeline
echo "Installing Shotgun Viral Genomics Pipeline..."

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH."
    echo "Please install Miniconda or Anaconda first: https://docs.conda.io/projects/conda/en/latest/user-guide/install/"
    exit 1
fi

# Create conda environment
echo "Creating conda environment 'viral_genomics'..."
conda env create -f environment.yml

# Activate the environment
echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate viral_genomics

# Make scripts executable
echo "Making scripts executable..."
chmod +x *.sh *.pl viral_pipeline.py

# Check for snpEff
echo "Checking for snpEff..."
if ! command -v snpEff &> /dev/null && ! [ -f snpEff.jar ]; then
    echo "Warning: snpEff not found in PATH or current directory."
    echo "Please install snpEff or make sure snpEff.jar is in the current directory."
    echo "You can download it from: http://pcingola.github.io/SnpEff/"
fi

echo ""
echo "Installation complete!"
echo ""
echo "To use the pipeline:"
echo "1. Activate the environment: conda activate viral_genomics"
echo "2. Run the pipeline: ./viral_pipeline.py --r1 <READ1_FILES> --r2 <READ2_FILES> --accession <ACCESSION_NUMBER> --threads <THREADS>"
echo ""
echo "For more information, see the README.md file."