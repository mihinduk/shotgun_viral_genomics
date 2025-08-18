# Viral Diagnostic Module - Quick Start Guide

## Prerequisites

### 1. Clone the Repository
```bash
git clone https://github.com/your-repo/shotgun_viral_genomics.git
cd shotgun_viral_genomics
```

### 2. Set Up HTCF Environment

#### Step 1: Set up the main viral genomics environment
```bash
./setup_htcf_env.sh
```

#### Step 2: Set up the viral assembly environment (REQUIRED for diagnostic module)
```bash
./setup_assembly_env.sh
```

#### Step 3: Install MEGAHIT in viral_genomics environment (if not already present)
```bash
source /ref/sahlab/software/anaconda3/bin/activate
/home/mihindu/miniforge3/bin/mamba install -n viral_genomics -c bioconda megahit -y
```

### 3. Verify Environments
```bash
# Check viral_genomics environment
conda activate viral_genomics
which bwa samtools megahit
conda deactivate

# Check viral_assembly environment  
conda activate viral_assembly
which megahit
conda deactivate
```

## Running the Diagnostic Module

### Basic Usage
```bash
sbatch submit_viral_diagnostic.sh <R1_fastq> <R2_fastq> <accession> <sample_name> [threads]
```

### Example Command
```bash
# Navigate to your data directory
cd /path/to/your/data

# Submit diagnostic job
sbatch /full/path/to/submit_viral_diagnostic.sh \
  "./sample_R1.fastq.gz" \
  "./sample_R2.fastq.gz" \
  "NC_001477.1" \
  "my_sample" \
  4
```

### Input Requirements
- **Cleaned reads (preferred)**: Place QC-cleaned reads in `../cleaned_seqs/` directory
  - Named as: `<original_name>.qc.fastq.gz`
- **Raw reads (fallback)**: Will be used if cleaned reads not found

### Expected Directory Structure
```
your_data_directory/
├── sample_R1.fastq.gz          # Raw reads
├── sample_R2.fastq.gz          # Raw reads
└── cleaned_seqs/
    ├── sample_R1.qc.fastq.gz   # Cleaned reads (preferred)
    └── sample_R2.qc.fastq.gz   # Cleaned reads (preferred)
```

### Output Files
The diagnostic will create a directory `diagnostic_<sample_name>/` containing:
- `<sample_name>_quick.bam` - Mapping to reference
- `assembly_<sample_name>/` - MEGAHIT assembly results
- `<sample_name>_contigs_filtered.fa` - Filtered contigs >1kb
- `<sample_name>_blast_all.tsv` - BLAST results
- `<sample_name>_diagnostic_report.txt` - Final contamination report

### Key Success Indicators
- Assembly generates contigs (check `assembly_<sample_name>/final.contigs.fa`)
- BLAST analysis completes (check `<sample_name>_blast_all.tsv`)
- Diagnostic report generated with contamination assessment

## Troubleshooting

### Mamba Lock Issues
If you see "Cannot lock" warnings, the script handles these automatically. No action needed.

### Assembly Fails
1. Check if viral_assembly environment exists: `conda env list | grep viral_assembly`
2. If missing, run: `./setup_assembly_env.sh`
3. Verify MEGAHIT is installed: `conda activate viral_assembly && which megahit`

### No Cleaned Reads Found
The script will use raw reads as fallback. For better results, run QC first:
```bash
fastp -i R1.fastq.gz -I R2.fastq.gz \
      -o cleaned_seqs/R1.qc.fastq.gz \
      -O cleaned_seqs/R2.qc.fastq.gz \
      --html qc_report.html
```

## Tested Configuration
- Successfully tested on HTCF with WNV sample (August 18, 2025)
- Job 26029658 completed with contamination detection working as expected
- Environments: viral_genomics (BWA/samtools), viral_assembly (MEGAHIT)
