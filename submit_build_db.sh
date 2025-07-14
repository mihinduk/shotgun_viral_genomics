#!/bin/bash
#SBATCH --job-name=build_contamination_db
#SBATCH --output=build_contamination_db_%j.out
#SBATCH --error=build_contamination_db_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general

# Submit script for building contamination database
# Usage: sbatch submit_build_db.sh

# Change to the directory where the script was submitted from
cd ${SLURM_SUBMIT_DIR}

echo "========================================="
echo "CONTAMINATION DATABASE BUILD"
echo "========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Job started at: $(date)"
echo "Working directory: $(pwd)"
echo "========================================="

# Set up environment for BLAST
echo "Setting up HTCF environment..."
source /ref/sahlab/software/anaconda3/bin/activate
eval "$(conda shell.bash hook)"
conda activate Blast

# Check if makeblastdb is available
if ! command -v makeblastdb &> /dev/null; then
    echo "Error: makeblastdb not found in Blast environment"
    exit 1
fi

# Run the database build script
echo "Running database build script..."
bash build_contamination_db.sh

BUILD_EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job completed at: $(date)"
if [ $BUILD_EXIT_CODE -eq 0 ]; then
    echo "Database build completed successfully!"
    echo ""
    # Check if database was created
    if [ -f "/ref/sahlab/databases/contamination_db/cell_culture_contaminants.nhr" ]; then
        echo "Database created at: /ref/sahlab/databases/contamination_db/cell_culture_contaminants"
        echo "Database is ready to use in diagnostic script!"
    fi
else
    echo "Database build failed with exit code: $BUILD_EXIT_CODE"
fi
echo "========================================="

exit $BUILD_EXIT_CODE