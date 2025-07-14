#!/bin/bash
#SBATCH --job-name=build_microbes_db
#SBATCH --output=build_microbes_db_%j.out
#SBATCH --error=build_microbes_db_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general

echo "========================================="
echo "MICROBIAL CONTAMINATION DATABASE BUILD"
echo "========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Job started at: $(date)"
echo "Working directory: $PWD"
echo "========================================="

# Set up environment
echo "Setting up HTCF environment..."
source /ref/sahlab/software/anaconda3/bin/activate

# Run the database build script
echo "Running database build script..."
bash build_microbes_db.sh

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "Job completed at: $(date)"
    echo "Database build completed successfully!"
    echo "========================================="
else
    echo ""
    echo "========================================="
    echo "Job completed at: $(date)"
    echo "Database build failed with exit code: $?"
    echo "========================================="
fi