#!/usr/bin/env python3

"""
create_viral_db.py

This script creates a custom BLAST database from viral sequences.
It can either use a local collection of FASTA files or download 
viral sequences from NCBI.

Dependencies:
- Python packages: biopython, requests
- External tools: BLAST+ (makeblastdb)

Usage:
python create_viral_db.py [options]

Examples:
# Create from local FASTA files
python create_viral_db.py --input viral_genomes/*.fasta --output viral_db

# Download and create from NCBI viral genomes
python create_viral_db.py --download --email your.email@example.com --output viral_db
"""

import argparse
import os
import sys
import glob
import subprocess
import tempfile
from typing import List, Optional

try:
    from Bio import Entrez, SeqIO
except ImportError:
    sys.stderr.write("Error: This script requires Biopython. "
                     "Please install with 'pip install biopython'.\n")
    sys.exit(1)


def check_dependencies():
    """Check if BLAST+ is installed."""
    try:
        subprocess.run(['makeblastdb', '-version'], 
                       stdout=subprocess.PIPE, 
                       stderr=subprocess.PIPE,
                       check=True)
        print("BLAST+ is available.")
    except (subprocess.SubprocessError, FileNotFoundError):
        sys.stderr.write("Error: BLAST+ (makeblastdb) not found. Please install BLAST+ "
                         "and ensure it's in your PATH.\n")
        sys.exit(1)


def create_db_from_files(input_files: List[str], output_db: str):
    """Create a BLAST database from a collection of FASTA files."""
    if not input_files:
        sys.stderr.write("Error: No input files found.\n")
        sys.exit(1)
    
    # Concatenate all input files into a single FASTA
    temp_fasta = tempfile.NamedTemporaryFile(delete=False, suffix='.fasta')
    temp_fasta.close()
    
    sequence_count = 0
    
    with open(temp_fasta.name, 'w') as out_handle:
        for fasta_file in input_files:
            try:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    SeqIO.write(record, out_handle, "fasta")
                    sequence_count += 1
            except Exception as e:
                sys.stderr.write(f"Warning: Error processing {fasta_file}: {str(e)}\n")
    
    if sequence_count == 0:
        os.unlink(temp_fasta.name)
        sys.stderr.write("Error: No valid sequences found in input files.\n")
        sys.exit(1)
    
    print(f"Processing {sequence_count} sequences from {len(input_files)} files.")
    
    # Create the BLAST database
    cmd = [
        'makeblastdb',
        '-in', temp_fasta.name,
        '-dbtype', 'nucl',
        '-out', output_db,
        '-title', 'Custom viral BLAST database'
    ]
    
    print("Creating BLAST database...")
    try:
        subprocess.run(cmd, check=True)
        print(f"Successfully created BLAST database: {output_db}")
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error creating BLAST database: {str(e)}\n")
        os.unlink(temp_fasta.name)
        sys.exit(1)
    
    # Clean up
    os.unlink(temp_fasta.name)


def download_viral_sequences(email: str, taxid: Optional[str] = None, max_sequences: int = 0):
    """Download viral sequences from NCBI."""
    # Set email for Entrez
    Entrez.email = email
    
    # Prepare the query
    if taxid:
        query = f"txid{taxid}[Organism:exp] AND complete[Title] AND genome[Title]"
    else:
        query = "viruses[Organism] AND complete[Title] AND genome[Title]"
    
    print(f"Searching NCBI for viral sequences with query: {query}")
    
    # First, get the list of IDs
    search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_sequences if max_sequences > 0 else 10000)
    search_record = Entrez.read(search_handle)
    search_handle.close()
    
    id_list = search_record["IdList"]
    
    if not id_list:
        sys.stderr.write("Error: No sequences found with the given query.\n")
        sys.exit(1)
    
    print(f"Found {len(id_list)} sequences. Downloading...")
    
    # Create a temporary file to store the downloaded sequences
    temp_fasta = tempfile.NamedTemporaryFile(delete=False, suffix='.fasta')
    temp_fasta.close()
    
    # Download in batches to avoid timeout issues
    batch_size = 100
    sequence_count = 0
    
    with open(temp_fasta.name, 'w') as out_handle:
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i+batch_size]
            try:
                # Download batch
                print(f"Downloading batch {i//batch_size + 1} of {(len(id_list)-1)//batch_size + 1}...")
                fetch_handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="fasta", retmode="text")
                data = fetch_handle.read()
                fetch_handle.close()
                
                # Write to file
                out_handle.write(data)
                sequence_count += len(batch_ids)
                
            except Exception as e:
                sys.stderr.write(f"Warning: Error downloading batch {i//batch_size + 1}: {str(e)}\n")
    
    print(f"Downloaded {sequence_count} sequences to {temp_fasta.name}")
    return temp_fasta.name


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Create a custom BLAST database from viral sequences.')
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--input', nargs='+', help='Input FASTA files or wildcards')
    input_group.add_argument('--download', action='store_true', help='Download viral sequences from NCBI')
    
    # Download options
    parser.add_argument('--email', help='Email address for NCBI Entrez queries (required with --download)')
    parser.add_argument('--taxid', help='NCBI Taxonomy ID to restrict download (e.g., 11118 for coronaviridae)')
    parser.add_argument('--max_sequences', type=int, default=0, 
                        help='Maximum number of sequences to download (0 for all available)')
    
    # Output options
    parser.add_argument('--output', required=True, help='Output BLAST database name')
    
    args = parser.parse_args()
    
    # Check dependencies
    check_dependencies()
    
    # Process based on input mode
    try:
        if args.download:
            if not args.email:
                sys.stderr.write("Error: --email is required with --download option.\n")
                sys.exit(1)
            
            # Download viral sequences
            temp_fasta = download_viral_sequences(args.email, args.taxid, args.max_sequences)
            
            # Create DB from the downloaded file
            create_db_from_files([temp_fasta], args.output)
            
            # Clean up
            os.unlink(temp_fasta)
            
        else:
            # Expand any wildcards in input files
            expanded_files = []
            for pattern in args.input:
                files = glob.glob(pattern)
                if not files:
                    sys.stderr.write(f"Warning: No files match pattern '{pattern}'\n")
                expanded_files.extend(files)
            
            # Create DB from input files
            create_db_from_files(expanded_files, args.output)
        
        print("\nDatabase creation completed successfully!")
        print(f"To use this database with annotate_from_blastn.py:")
        print(f"python annotate_from_blastn.py input.fasta output.gbk --blast_db {args.output}")
        
    except Exception as e:
        sys.stderr.write(f"Error: {str(e)}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()