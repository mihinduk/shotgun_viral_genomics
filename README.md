# Shotgun Viral Genomics Pipeline

A comprehensive pipeline for processing shotgun sequencing data from viral samples, performing quality control, mapping to reference genomes, variant calling, and annotation.

## Overview

This pipeline takes Illumina shotgun sequencing data and an accession number for a viral genome, then:
1. Calculates basic statistics on input reads (seqkit stats)
2. Cleans the input sequences (fastp)
3. Maps reads to the reference genome (bwa)
4. Calls variants (samtools, lofreq)
5. Filters variants (lofreq)
6. Annotates variants (snpEff)
7. Parses the annotated variants (custom Perl script)

## Installation

### Prerequisites

The pipeline requires the following tools:
- seqkit (v2.0.0+)
- fastp (v0.22.0+)
- bwa (v0.7.17+)
- samtools (v1.18)
- lofreq (v2.1.5+)
- snpEff (v5.0+)
- perl (v5.32.1+)
- python (v3.9+)
- pandas, numpy (for parsing outputs)
- entrez-direct (for downloading reference genomes)

### Method 1: Clone the Repository

```bash
git clone https://github.com/mihinduk/shotgun_viral_genomics
cd shotgun_viral_genomics
```

### Method 2: Install via Conda (Recommended)

```bash
# Create a conda environment
conda create -n viral_genomics -c bioconda -c conda-forge \
  seqkit=2.0.0 fastp=0.22.0 bwa=0.7.17 samtools=1.18 lofreq=2.1.5 \
  perl=5.32.1 entrez-direct python=3.9 numpy pandas
conda activate viral_genomics

# Make scripts executable
chmod +x *.sh *.pl viral_pipeline.py
```

## Quick Start

```bash
# Activate the environment (if using conda)
conda activate viral_genomics

# Run the pipeline with default parameters
./viral_pipeline.py --r1 "*_R1.fastq.gz" --r2 "*_R2.fastq.gz" --accession NC_045512.2 --threads 12

# Or run with explicit file paths
./viral_pipeline.py --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --accession NC_045512.2 --threads 12
```

## Usage

### Basic Usage

```bash
./viral_pipeline.py --r1 <READ1_FILES> --r2 <READ2_FILES> --accession <ACCESSION_NUMBER> --threads <THREADS>
```

### Parameters

- `--r1`: Forward reads (R1) FASTQ files (accepts wildcards like *_R1.fastq.gz or *_R1_001.fastq.gz)
- `--r2`: Reverse reads (R2) FASTQ files (accepts wildcards like *_R2.fastq.gz or *_R2_001.fastq.gz)
- `--accession`: Reference genome accession number (e.g., NC_045512.2)
- `--reference`: Path to local reference genome (when not using accession)
- `--threads`: Number of CPU threads to use (default: 1)
- `--outdir`: Output directory (default: current directory)
- `--min-depth`: Minimum read depth for final variant reporting (default: 200)

### Advanced Parameters

- `--skip-download`: Skip genome download (use local reference)
- `--force-download`: Force download even if reference exists locally
- `--skip-qc`: Skip quality control step (useful for pre-processed data)
- `--skip-stats`: Skip read statistics step
- `--skip-mapping`: Skip read mapping step
- `--skip-variants`: Skip variant calling step
- `--skip-annotation`: Skip variant annotation step
- `--add-to-snpeff`: Attempt to add genome to snpEff if not already present
- `--snpeff-jar`: Path to snpEff.jar (default: "snpEff.jar")
- `--keep-tmp`: Keep temporary files (useful for debugging)

## Output Structure

- `./cleaned_seqs/` - BWA indexed genome, cleaned fastqs, fastp reports
- `./cleaned_seqs/mapping/` - Alignment files (.sam, .bam)
- `./cleaned_seqs/variants/` - Variant calls and annotations

## Pipeline Steps

### 1. Reference Genome Preparation

- Downloads reference genome from NCBI based on accession number
- Checks if the reference genome is in the snpEff database
- Optionally adds the genome to snpEff database if missing

### 2. Quality Control & Statistics

- Generates statistics on input reads using seqkit
- Cleans reads with fastp (adapter trimming, quality filtering)
- Produces QC reports with key metrics:
  - Number of reads passing filter
  - Number of reads failing due to quality, N content, or length
  - Number of reads with adapter trimming
  - Duplication rate
  - Insert size peak

### 3. Alignment & Variant Calling

- Indexes reference genome with BWA
- Maps reads to reference with BWA MEM
- Processes alignments:
  - Fixing mate information
  - Sorting BAM files
  - Marking duplicate reads
- Improves alignment quality with LoFreq:
  - Viterbi realignment
  - Indel quality calibration
  - Alignment quality calibration
- Performs variant calling with LoFreq

### 4. Variant Filtering & Annotation

- Filters variants based on quality criteria
- Annotates variants using snpEff
- Parses annotated variants into user-friendly TSV format
- Generates summary reports of variant effects

## Examples

### Example 1: SARS-CoV-2 Analysis

```bash
./viral_pipeline.py --r1 "*_R1.fastq.gz" --r2 "*_R2.fastq.gz" --accession NC_045512.2 --threads 12
```

### Example 2: Custom Reference with Minimum Depth

```bash
./viral_pipeline.py --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --reference custom.fasta --min-depth 100
```

### Example 3: Add Missing Genome to snpEff

```bash
./viral_pipeline.py --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --accession AF253419.1 --add-to-snpeff --threads 8
```

## Common Issues

### The needed genome is NOT in the snpEff database

**Detection:**
```bash
java -jar snpEff.jar databases | grep "AF253419.1"
```

**Fix:** Follow the protocol at [SnpEff-Viral-Genomes](https://github.com/mihinduk/Bioinformatics_protocols/blob/main/SnpEff/snpEff-Viral-Genomes.md)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- Based on the viral_genome_sop pipeline written by: Ana Jung
- Protocol for adding viral genomes to snpEff: [SnpEff-Viral-Genomes](https://github.com/mihinduk/Bioinformatics_protocols/blob/main/SnpEff/snpEff-Viral-Genomes.md)
- [LoFreq](https://csb5.github.io/lofreq/) for variant calling
- [snpEff](http://pcingola.github.io/SnpEff/) for variant annotation