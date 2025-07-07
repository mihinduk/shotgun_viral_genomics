#!/bin/bash
# This script installs Biopython and other dependencies

# Check if conda is active
if [ -z "$CONDA_PREFIX" ]; then
  echo "WARNING: No active conda environment detected. It's recommended to install within your conda environment."
  echo "If you want to proceed anyway, press Enter. Otherwise, exit and activate your conda environment."
  read -r
fi

# Install Biopython
echo "Installing Biopython..."
pip install biopython

# Verify installation
python -c "from Bio import SeqIO; print('Biopython installed successfully!')"

if [ $? -eq 0 ]; then
  echo "Installation successful!"
else
  echo "Installation failed. Try using conda instead:"
  echo "conda install -c conda-forge biopython"
fi