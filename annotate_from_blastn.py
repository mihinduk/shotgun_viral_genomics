#!/usr/bin/env python3

"""
annotate_from_blastn.py

This script takes a FASTA file, finds the most closely related genome using blastn,
and applies that annotation to the query sequence to generate a GenBank file.

Dependencies:
- Python packages: biopython, requests
- External tools: BLAST+ (blastn)

Usage:
python annotate_from_blastn.py <input.fasta> <output.gbk>
"""

import argparse
import os
import sys
import tempfile
import subprocess
import re
import datetime
from io import StringIO
from typing import List, Dict, Tuple, Optional, Any

try:
    from Bio import SeqIO, Entrez
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    from Bio.Blast import NCBIXML
    import requests
except ImportError:
    sys.stderr.write("Error: This script requires Biopython and requests. "
                     "Please install with 'pip install biopython requests'.\n")
    sys.exit(1)

# Set your email for Entrez queries
Entrez.email = "your_email@example.com"  # Please replace with a valid email


def check_dependencies():
    """Check if all required dependencies are installed."""
    # Check if BLAST+ is installed
    try:
        subprocess.run(['blastn', '-version'], 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE,
                        check=True)
        print("All dependencies satisfied.")
    except (subprocess.SubprocessError, FileNotFoundError):
        sys.stderr.write("Error: BLAST+ not found. Please install BLAST+ and ensure it's in your PATH.\n")
        sys.exit(1)


def validate_input(args):
    """Validate command line arguments."""
    if not os.path.exists(args.input_file):
        sys.stderr.write(f"Error: Input file '{args.input_file}' not found.\n")
        sys.exit(1)
    
    if not re.search(r'\.fa$|\.fasta$', args.input_file, re.IGNORECASE):
        sys.stderr.write(f"Warning: Input file does not have a .fa or .fasta extension.\n")
    
    if not re.search(r'\.gb$|\.gbk$|\.genbank$', args.output_file, re.IGNORECASE):
        sys.stderr.write(f"Warning: Output file does not have a .gb, .gbk, or .genbank extension.\n")


def read_fasta(file_path):
    """Read a FASTA file and return a list of sequence records."""
    try:
        sequences = list(SeqIO.parse(file_path, "fasta"))
        if not sequences:
            sys.stderr.write(f"Error: No sequences found in FASTA file.\n")
            sys.exit(1)
        print(f"Read {len(sequences)} sequence(s) from {file_path}")
        return sequences
    except Exception as e:
        sys.stderr.write(f"Error reading FASTA file: {str(e)}\n")
        sys.exit(1)


def find_closest_genome(sequence, tmp_dir=None, max_targets=5, blast_db='nt'):
    """Run BLAST to find the closest genome and return the top hit."""
    if tmp_dir is None:
        tmp_dir = tempfile.mkdtemp()
    
    tmp_fasta = os.path.join(tmp_dir, "query.fasta")
    tmp_blast_out = os.path.join(tmp_dir, "blast_results.xml")
    
    # Write sequence to temporary FASTA file
    with open(tmp_fasta, 'w') as f:
        SeqIO.write(sequence, f, "fasta")
    
    # Run blastn against specified database
    print(f"Running BLAST to find closest genomes using database: {blast_db}...")
    blast_cmd = [
        'blastn',
        '-query', tmp_fasta,
        '-db', blast_db,
        '-outfmt', '5',  # XML format
        '-max_target_seqs', str(max_targets),
        '-out', tmp_blast_out
    ]
    
    try:
        subprocess.run(blast_cmd, check=True)
    except subprocess.CalledProcessError:
        sys.stderr.write("Error: BLAST search failed.\n")
        sys.exit(1)
    
    if not os.path.exists(tmp_blast_out) or os.path.getsize(tmp_blast_out) == 0:
        sys.stderr.write("Error: BLAST search returned no hits.\n")
        sys.exit(1)
    
    # Parse BLAST results
    with open(tmp_blast_out) as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))
    
    if not blast_records or not blast_records[0].alignments:
        sys.stderr.write("Error: No BLAST hits found.\n")
        sys.exit(1)
    
    record = blast_records[0]
    hits = []
    
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            # Extract accession number
            # First check if we have a proper accession, if not use hit_id or title
            accession = alignment.accession
            hit_id = alignment.hit_id.split('|')
            title = alignment.title
            
            # First try: Look for accession in hit_id
            if accession.isdigit() and len(hit_id) > 1:
                # Look for accession format in hit_id
                for part in hit_id:
                    # Look for typical accession format (letters followed by numbers, possibly with dots/underscores)
                    if re.match(r'^[A-Za-z]+[0-9._]+$', part):
                        accession = part
                        break
            
            # Second try: Look for accession in title
            if accession.isdigit() and title:
                # Extract accession patterns like XX123456.1 from the title
                accession_match = re.search(r'([A-Z]{2}\d{6,}\.\d+)', title)
                if accession_match:
                    accession = accession_match.group(1)
                    print(f"Extracted accession {accession} from title")
            
            # Store hit information
            hits.append({
                'accession': accession,
                'hit_id': alignment.hit_id,
                'title': alignment.title,
                'pident': (hsp.identities / hsp.align_length) * 100,
                'length': hsp.align_length,
                'mismatch': hsp.align_length - hsp.identities,
                'gapopen': hsp.gaps,
                'qstart': hsp.query_start,
                'qend': hsp.query_end,
                'sstart': hsp.sbjct_start,
                'send': hsp.sbjct_end,
                'evalue': hsp.expect,
                'bitscore': hsp.bits
            })
    
    # Sort hits by bitscore (highest first)
    hits.sort(key=lambda x: x['bitscore'], reverse=True)
    
    # Print top hits
    print(f"Found {len(hits)} BLAST hits.")
    if hits:
        print(f"Top {min(5, len(hits))} hits:")
        for i, hit in enumerate(hits[:5], 1):
            print(f"{i}. {hit['accession']} ({hit['pident']:.2f}% identity, E-value: {hit['evalue']:.2e}): {hit['title']}")
        
        # Return top hit
        return hits[0]
    else:
        sys.stderr.write("Error: No BLAST hits found after parsing.\n")
        sys.exit(1)


def fetch_genbank(accession, hit_id=None, tmp_dir=None):
    """Fetch GenBank record for the closest match."""
    if tmp_dir is None:
        tmp_dir = tempfile.mkdtemp()
    
    tmp_gb = os.path.join(tmp_dir, "reference.gb")
    
    print(f"Fetching GenBank record for accession {accession}...")
    
    # If accession is numeric only (likely a GI number), try to get the corresponding accession
    if accession.isdigit() and hit_id:
        print(f"Detected numeric-only accession (GI number). Attempting to find proper accession from hit_id: {hit_id}")
        
        # First try: extract from hit_id
        parts = hit_id.split('|')
        for part in parts:
            if re.match(r'^[A-Za-z]+[0-9._]+$', part):
                accession = part
                print(f"Using alternative accession: {accession}")
                break
        
        # Second try: Look for accession in hit_id string
        if accession.isdigit():  # Still numeric, try another approach
            accession_match = re.search(r'([A-Z]{2}\d{6,}\.\d+)', hit_id)
            if accession_match:
                accession = accession_match.group(1)
                print(f"Extracted accession {accession} from hit_id")
    
    # Try to fetch using Entrez
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        gb_record = handle.read()
        handle.close()
        
        with open(tmp_gb, 'w') as f:
            f.write(gb_record)
        
        print(f"Downloaded GenBank record to {tmp_gb}")
        return tmp_gb
    
    except Exception as e:
        print(f"First attempt failed: {str(e)}")
        # If Entrez fails, try using direct URL
        try:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=gb&retmode=text"
            response = requests.get(url)
            response.raise_for_status()
            
            with open(tmp_gb, 'w') as f:
                f.write(response.text)
            
            print(f"Downloaded GenBank record to {tmp_gb}")
            return tmp_gb
        
        except Exception as e2:
            # If both methods fail and we have a hit_id, try searching by hit description
            if hit_id:
                try:
                    print(f"Trying to search for GenBank record using hit description...")
                    # Extract organism name from hit title if possible
                    parts = hit_id.split()
                    if len(parts) > 2:
                        search_term = " ".join(parts[1:3]) + " complete genome"
                        print(f"Searching for: {search_term}")
                        
                        # Search for similar sequences
                        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1)
                        record = Entrez.read(handle)
                        handle.close()
                        
                        if record["IdList"]:
                            new_id = record["IdList"][0]
                            print(f"Found alternative ID: {new_id}")
                            
                            # Fetch the record
                            handle = Entrez.efetch(db="nucleotide", id=new_id, rettype="gb", retmode="text")
                            gb_record = handle.read()
                            handle.close()
                            
                            with open(tmp_gb, 'w') as f:
                                f.write(gb_record)
                            
                            print(f"Downloaded alternative GenBank record to {tmp_gb}")
                            return tmp_gb
                    
                except Exception as e3:
                    print(f"Alternative search failed: {str(e3)}")
            
            sys.stderr.write(f"Error: Failed to fetch GenBank record for accession {accession}: {str(e2)}\n")
            sys.exit(1)


def parse_genbank_features(gb_file):
    """Parse GenBank features from a file."""
    print("Parsing GenBank features...")
    
    try:
        record = next(SeqIO.parse(gb_file, "genbank"))
        features = [feat for feat in record.features if feat.type != "source"]
        
        if not features:
            sys.stderr.write("Error: No features found in GenBank file (excluding source).\n")
            sys.exit(1)
        
        print(f"Extracted {len(features)} features from GenBank file.")
        return features, record
    
    except Exception as e:
        sys.stderr.write(f"Error parsing GenBank file: {str(e)}\n")
        sys.exit(1)


def map_annotations(query_seq, reference_features, blast_hit):
    """Map annotations from reference to query sequence."""
    print("Mapping annotations to query sequence...")
    
    # Calculate query length
    query_length = len(query_seq)
    
    # Get coordinates from BLAST hit
    query_start = blast_hit['qstart']
    reference_start = blast_hit['sstart']
    reference_end = blast_hit['send']
    
    # Determine if reference sequence is in reverse orientation
    is_reverse = reference_start > reference_end
    
    if is_reverse:
        # Swap start and end for reverse orientation
        reference_start, reference_end = reference_end, reference_start
    
    # Map features
    mapped_features = []
    
    for feature in reference_features:
        # Skip features without location
        if not feature.location:
            continue
        
        # Get feature location
        if isinstance(feature.location, CompoundLocation):
            # Skip complex features for now (potential future enhancement)
            continue
        
        feat_start = int(feature.location.start)
        feat_end = int(feature.location.end)
        
        # Skip features outside of the reference alignment region
        if feat_end < reference_start or feat_start > reference_end:
            continue
        
        # Calculate new positions on query sequence
        if is_reverse:
            # For reverse orientation
            query_feat_start = query_start + (reference_end - feat_end)
            query_feat_end = query_start + (reference_end - feat_start)
            strand = -feature.location.strand if feature.location.strand else 0
        else:
            # For forward orientation
            query_feat_start = query_start + (feat_start - reference_start)
            query_feat_end = query_start + (feat_end - reference_start)
            strand = feature.location.strand
        
        # Ensure feature is within query bounds
        if query_feat_start < 0:
            query_feat_start = 0
        if query_feat_end > query_length:
            query_feat_end = query_length
        
        # Skip features that would be entirely outside query bounds
        if query_feat_start >= query_length or query_feat_end <= 0:
            continue
        
        # Create new feature
        new_location = FeatureLocation(query_feat_start, query_feat_end, strand)
        new_feature = SeqFeature(
            location=new_location,
            type=feature.type,
            qualifiers=feature.qualifiers.copy()
        )
        
        mapped_features.append(new_feature)
    
    print(f"Mapped {len(mapped_features)} features to query sequence.")
    return mapped_features


def generate_genbank(query_seq, query_name, mapped_features, output_file, blast_hit):
    """Generate a GenBank file with the mapped annotations."""
    print("Generating GenBank file...")
    
    # Create sequence record
    if isinstance(query_seq, SeqRecord):
        seq_record = query_seq
    else:
        seq_record = SeqRecord(query_seq)
    
    # Update record attributes
    seq_record.id = query_name
    seq_record.name = query_name
    seq_record.description = "Automatically annotated using annotate_from_blastn.py"
    
    # Add annotations
    seq_record.annotations["molecule_type"] = "DNA"
    seq_record.annotations["topology"] = "linear"
    seq_record.annotations["date"] = datetime.datetime.now().strftime("%d-%b-%Y")
    
    # Clear existing features and add source feature
    seq_record.features = []
    
    source_feature = SeqFeature(
        FeatureLocation(0, len(seq_record.seq)),
        type="source",
        qualifiers={
            "organism": ["Unknown"],
            "mol_type": ["genomic DNA"],
            "note": [f"Annotated using closest BLAST hit: {blast_hit['accession']} "
                    f"({blast_hit['pident']:.2f}% identity, E-value: {blast_hit['evalue']:.2e})"]
        }
    )
    seq_record.features.append(source_feature)
    
    # Add mapped features
    for feature in mapped_features:
        seq_record.features.append(feature)
    
    # Write GenBank file
    SeqIO.write(seq_record, output_file, "genbank")
    print(f"GenBank file written to {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Annotate a FASTA file using BLAST and GenBank.')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('output_file', help='Output GenBank file')
    parser.add_argument('--blast_db', help='BLAST database to use (default: nt)', default='nt')
    parser.add_argument('--email', help='Email for NCBI Entrez queries', default=Entrez.email)
    
    args = parser.parse_args()
    
    # Set email for Entrez
    if args.email:
        Entrez.email = args.email
    
    # Check dependencies
    check_dependencies()
    
    # Validate input
    validate_input(args)
    
    # Create temporary directory
    tmp_dir = tempfile.mkdtemp()
    
    try:
        # Read FASTA file
        sequences = read_fasta(args.input_file)
        
        # For simplicity, use only the first sequence if there are multiple
        if len(sequences) > 1:
            sys.stderr.write("Warning: Multiple sequences found in input file. Using only the first sequence.\n")
            query_seq = sequences[0]
        else:
            query_seq = sequences[0]
        
        # Extract query name
        query_name = query_seq.id
        
        # Run BLAST to find closest genome
        blast_hit = find_closest_genome(query_seq, tmp_dir, blast_db=args.blast_db)
        
        # Download GenBank record for closest match
        gb_file = fetch_genbank(blast_hit['accession'], blast_hit.get('hit_id', ''), tmp_dir)
        
        # Parse GenBank features
        reference_features, reference_record = parse_genbank_features(gb_file)
        
        # Map annotations to query sequence
        mapped_features = map_annotations(query_seq, reference_features, blast_hit)
        
        # Generate GenBank file
        generate_genbank(query_seq, query_name, mapped_features, args.output_file, blast_hit)
        
        print("\nAnnotation completed successfully!")
        
    finally:
        # Clean up temporary files
        import shutil
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()