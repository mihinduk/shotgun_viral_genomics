# Annotate from BLASTN

This tool takes a FASTA file, finds the most closely related genome using BLASTN, and applies that annotation to the query sequence to generate a GenBank file.

## Requirements

- Python 3.6+
- Required Python packages:
  - Biopython
  - requests
- BLAST+ (installed and in PATH)

## Installation

1. Make sure BLAST+ is installed and in your PATH
2. Install required Python packages:

```bash
pip install biopython requests
```

Or use conda:

```bash
conda install -c conda-forge -c bioconda biopython requests blast
```

3. Make the script executable:

```bash
chmod +x annotate_from_blastn.py
chmod +x create_viral_db.py
```

## Usage

### Annotate a FASTA file

```bash
python annotate_from_blastn.py [input.fasta] [output.gbk] [options]
```

#### Arguments

- `input.fasta`: Path to the input FASTA file containing the sequence to annotate
- `output.gbk`: Path where the output GenBank file will be written

#### Options

- `--blast_db`: BLAST database to use (default: nt)
- `--email`: Email for NCBI Entrez queries

### Create a Custom Viral BLAST Database

```bash
python create_viral_db.py [options]
```

#### Options for creating a database from local files

```bash
python create_viral_db.py --input viral_genomes/*.fasta --output viral_db
```

#### Options for downloading and creating a database from NCBI

```bash
python create_viral_db.py --download --email your.email@example.com --output viral_db
```

Optional parameters:
- `--taxid`: NCBI Taxonomy ID to restrict download (e.g., 11118 for coronaviridae)
- `--max_sequences`: Maximum number of sequences to download (0 for all available)

## How It Works

1. Reads the input FASTA file
2. Runs BLASTN against the specified database to find the closest related genome
3. Downloads the GenBank record for the top hit
4. Extracts features from the reference GenBank record
5. Maps features from the reference to the query sequence, adjusting coordinates
6. Generates a new GenBank file with the mapped annotations

## Example

```bash
# Create a custom viral database
python create_viral_db.py --download --email your.email@example.com --output rsv_db --taxid 11250

# Annotate a sequence using the custom database
python annotate_from_blastn.py sequence.fa sequence.gbk --blast_db rsv_db --email your.email@example.com
```

## Notes

- If multiple sequences are found in the input FASTA, only the first one will be used
- Features are mapped based on BLAST alignment coordinates
- Features that don't overlap with the aligned region are excluded
- The script handles both forward and reverse orientations

## Integration with Shotgun Viral Genomics Pipeline

This script can be used as part of the shotgun viral genomics pipeline or as a standalone tool for quickly annotating viral sequences when a closely related reference exists in the NCBI database.