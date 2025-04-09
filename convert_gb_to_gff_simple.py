#!/usr/bin/env python3
"""
Convert GenBank to GFF3 with advanced handling for viral polyproteins.
Simplified version that doesn't require Biopython.

This script specifically handles the complex annotation structures in viral genomes,
particularly those with polyproteins and mature peptides (mat_peptide features).
It creates proper parent-child relationships in the GFF3 output to ensure the
annotation is correctly interpreted by tools like snpEff.

Example usage:
    python3 convert_gb_to_gff_simple.py input.gb output.gff3
"""

import argparse
import re
import sys
import os
from datetime import datetime

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Convert GenBank to GFF3 with advanced handling for viral polyproteins.')
    parser.add_argument('genbank_file', help='Input GenBank file')
    parser.add_argument('gff_output', help='Output GFF3 file')
    parser.add_argument('--source', default='GenBank', help='Source field for GFF (default: GenBank)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    return parser.parse_args()

def clean_product_name(product):
    """Clean product name to make it suitable for an ID."""
    # Remove " protein" suffix if present
    product = re.sub(r'\s+protein$', '', product)
    # Remove special characters, replace spaces with underscores
    product = re.sub(r'[^\w\s]', '', product)
    product = re.sub(r'\s+', '_', product)
    return product

def extract_features_from_genbank(gb_file):
    """Extract features from GenBank file without using Biopython."""
    with open(gb_file, 'r') as f:
        content = f.read()
    
    # Extract sequence ID
    accession_match = re.search(r'ACCESSION\s+(\S+)', content)
    if accession_match:
        sequence_id = accession_match.group(1)
    else:
        # Try alternative approach
        locus_match = re.search(r'LOCUS\s+(\S+)', content)
        if locus_match:
            sequence_id = locus_match.group(1)
        else:
            # Last resort - use the file name
            sequence_id = os.path.basename(gb_file).split('.')[0]
    
    # Extract sequence length
    length_match = re.search(r'LOCUS\s+\S+\s+(\d+)\s+bp', content)
    if length_match:
        sequence_length = length_match.group(1)
    else:
        # Default fallback
        sequence_length = "10000"
    
    # Extract sequence
    sequence = ""
    in_sequence = False
    for line in content.split('\n'):
        if line.startswith('ORIGIN'):
            in_sequence = True
            continue
        elif line.startswith('//'):
            in_sequence = False
            continue
        
        if in_sequence:
            # Remove numbers and whitespace from sequence lines
            sequence_part = re.sub(r'\d+', '', line).strip().replace(' ', '')
            sequence += sequence_part
    
    # Parse features
    features = []
    feature_blocks = re.findall(r'     (\S+)\s+([^/]+)((?:/[^=]+="[^"]*")+)?', content)
    
    for feature_type, location, qualifiers_text in feature_blocks:
        # Parse location
        location = location.strip()
        start = 0
        end = 0
        strand = '+'
        
        # Handle complement
        if 'complement' in location:
            strand = '-'
            location = re.sub(r'complement\(|\)', '', location)
        
        # Handle join - just take the overall range for simplicity
        if 'join' in location:
            location = re.sub(r'join\(|\)', '', location)
            parts = location.split(',')
            starts = []
            ends = []
            for part in parts:
                part_start, part_end = part.split('..')
                starts.append(int(part_start))
                ends.append(int(part_end))
            start = min(starts)
            end = max(ends)
        else:
            # Regular location
            if '..' in location:
                start, end = location.split('..')
            else:
                # Single position
                start = end = location
        
        # Convert to integers, making start 0-based for internal use
        try:
            start = int(start) - 1  # GFF is 1-based, but we'll convert when writing
            end = int(end)
        except ValueError:
            # Try to handle expressions like '<1' or '>10000'
            start = int(re.sub(r'[<>]', '', start)) - 1
            end = int(re.sub(r'[<>]', '', end))
        
        # Parse qualifiers
        qualifiers = {}
        if qualifiers_text:
            qualifier_matches = re.findall(r'/([^=]+)="([^"]*)"', qualifiers_text)
            for key, value in qualifier_matches:
                if key not in qualifiers:
                    qualifiers[key] = []
                qualifiers[key].append(value)
        
        features.append({
            'type': feature_type,
            'start': start,
            'end': end,
            'strand': strand,
            'qualifiers': qualifiers
        })
    
    return sequence_id, sequence_length, sequence, features

def convert_genbank_to_gff3(genbank_file, gff_output, source="GenBank", verbose=False):
    """
    Convert GenBank file to GFF3 with special handling for viral polyproteins.
    """
    # Extract features from GenBank
    if verbose:
        print(f"Reading GenBank file: {genbank_file}")
    
    try:
        record_id, sequence_length, sequence, features = extract_features_from_genbank(genbank_file)
    except Exception as e:
        print(f"Error extracting features from GenBank: {str(e)}")
        sys.exit(1)
    
    if verbose:
        print(f"Found {len(features)} features")
    
    # Open output file
    with open(gff_output, 'w') as gff_handle:
        # Write GFF3 header
        gff_handle.write('##gff-version 3\n')
        gff_handle.write(f'##date {datetime.now().strftime("%Y-%m-%d")}\n')
        gff_handle.write(f'##source-version convert_gb_to_gff_simple.py\n')
        gff_handle.write(f'##sequence-region {record_id} 1 {sequence_length}\n')
        
        # Track feature counters and IDs
        feature_counter = {}
        feature_ids = {}
        polyproteins = {}
        
        # First pass: count features by type
        for feature in features:
            feature_type = feature['type']
            if feature_type not in feature_counter:
                feature_counter[feature_type] = 0
            feature_counter[feature_type] += 1
        
        # Reset counters for actual processing
        actual_counter = {}
        
        # Process source feature first
        gff_handle.write(f"{record_id}\t{source}\tregion\t1\t{sequence_length}\t.\t+\t.\tID=region_1\n")
        
        # Store all feature lines to write later
        gff_lines = []
        
        # First process gene features to establish parent relationships
        for feature in features:
            if feature['type'] != 'gene':
                continue
                
            feature_type = feature['type']
            if feature_type not in actual_counter:
                actual_counter[feature_type] = 0
            actual_counter[feature_type] += 1
            
            # Generate ID for gene
            gene_id = None
            if 'gene' in feature['qualifiers']:
                gene_name = feature['qualifiers']['gene'][0]
                gene_id = gene_name
                feature_ids[gene_name] = gene_id
            elif 'locus_tag' in feature['qualifiers']:
                locus_tag = feature['qualifiers']['locus_tag'][0]
                gene_id = locus_tag
                feature_ids[locus_tag] = gene_id
            else:
                gene_id = f"gene_{actual_counter[feature_type]}"
            
            # Format attributes
            attributes = [f"ID={gene_id}"]
            if 'gene' in feature['qualifiers']:
                attributes.append(f"Name={feature['qualifiers']['gene'][0]}")
            if 'locus_tag' in feature['qualifiers']:
                attributes.append(f"locus_tag={feature['qualifiers']['locus_tag'][0]}")
            
            # Add other qualifiers
            for key, values in feature['qualifiers'].items():
                if key not in ['gene', 'locus_tag']:
                    for value in values:
                        value = value.replace(';', '%3B').replace('=', '%3D').replace('&', '%26')
                        attributes.append(f"{key}={value}")
            
            # Format GFF line
            strand = feature['strand']
            gff_line = [
                record_id,
                source,
                feature_type,
                str(feature['start'] + 1),  # Convert to 1-based
                str(feature['end']),
                '.',
                strand,
                '.',
                ";".join(attributes)
            ]
            gff_lines.append("\t".join(gff_line))
        
        # Process CDS features
        for feature in features:
            if feature['type'] != 'CDS':
                continue
                
            feature_type = feature['type']
            if feature_type not in actual_counter:
                actual_counter[feature_type] = 0
            actual_counter[feature_type] += 1
            
            # Generate ID for CDS - prioritize product name
            feature_id = None
            if 'product' in feature['qualifiers']:
                product_name = feature['qualifiers']['product'][0]
                clean_name = clean_product_name(product_name)
                feature_id = f"{clean_name}_{actual_counter[feature_type]}"
            elif 'protein_id' in feature['qualifiers']:
                feature_id = feature['qualifiers']['protein_id'][0]
            else:
                feature_id = f"cds_{actual_counter[feature_type]}"
            
            # Find parent gene
            parent_id = None
            if 'gene' in feature['qualifiers'] and feature['qualifiers']['gene'][0] in feature_ids:
                parent_id = feature_ids[feature['qualifiers']['gene'][0]]
            elif 'locus_tag' in feature['qualifiers'] and feature['qualifiers']['locus_tag'][0] in feature_ids:
                parent_id = feature_ids[feature['qualifiers']['locus_tag'][0]]
            
            # Format attributes
            attributes = [f"ID={feature_id}"]
            if parent_id:
                attributes.append(f"Parent={parent_id}")
            
            # Add name attribute based on product for better display in snpEff
            if 'product' in feature['qualifiers']:
                product_name = feature['qualifiers']['product'][0]
                attributes.append(f"product={product_name}")
                attributes.append(f"Name={product_name}")
                
                # Check if polyprotein for later reference by mat_peptide
                if 'polyprotein' in product_name.lower():
                    polyproteins[feature_id] = feature
            
            # Always include protein_id if available
            if 'protein_id' in feature['qualifiers']:
                attributes.append(f"protein_id={feature['qualifiers']['protein_id'][0]}")
            
            # Add other qualifiers
            for key, values in feature['qualifiers'].items():
                if key not in ['protein_id', 'product', 'gene', 'locus_tag', 'translation']:
                    for value in values:
                        value = value.replace(';', '%3B').replace('=', '%3D').replace('&', '%26')
                        attributes.append(f"{key}={value}")
            
            # Determine phase
            phase = '0'
            if 'codon_start' in feature['qualifiers']:
                try:
                    phase = str(int(feature['qualifiers']['codon_start'][0]) - 1)
                except:
                    phase = '0'
            
            # Format GFF line
            strand = feature['strand']
            gff_line = [
                record_id,
                source,
                feature_type,
                str(feature['start'] + 1),  # Convert to 1-based
                str(feature['end']),
                '.',
                strand,
                phase,
                ";".join(attributes)
            ]
            gff_lines.append("\t".join(gff_line))
        
        # Process mat_peptide features
        for feature in features:
            if feature['type'] != 'mat_peptide':
                continue
                
            feature_type = feature['type']
            if feature_type not in actual_counter:
                actual_counter[feature_type] = 0
            actual_counter[feature_type] += 1
            
            # Generate ID for mature peptide
            feature_id = None
            if 'product' in feature['qualifiers']:
                product_name = feature['qualifiers']['product'][0]
                clean_name = clean_product_name(product_name)
                feature_id = f"{clean_name}_{actual_counter[feature_type]}"
            else:
                feature_id = f"mat_peptide_{actual_counter[feature_type]}"
            
            # Find parent gene or polyprotein
            parent_id = None
            
            # First try gene parent
            if 'gene' in feature['qualifiers'] and feature['qualifiers']['gene'][0] in feature_ids:
                parent_id = feature_ids[feature['qualifiers']['gene'][0]]
            elif 'locus_tag' in feature['qualifiers'] and feature['qualifiers']['locus_tag'][0] in feature_ids:
                parent_id = feature_ids[feature['qualifiers']['locus_tag'][0]]
            
            # If no gene parent, try to find polyprotein parent
            if not parent_id:
                for poly_id, poly_feature in polyproteins.items():
                    if (feature['start'] >= poly_feature['start'] and 
                        feature['end'] <= poly_feature['end']):
                        parent_id = poly_id
                        break
            
            # Format attributes
            attributes = [f"ID={feature_id}"]
            if parent_id:
                attributes.append(f"Parent={parent_id}")
            
            # Add name and product attributes
            if 'product' in feature['qualifiers']:
                product_name = feature['qualifiers']['product'][0]
                attributes.append(f"product={product_name}")
                attributes.append(f"Name={product_name}")
            
            # Add other qualifiers
            for key, values in feature['qualifiers'].items():
                if key not in ['product', 'gene', 'locus_tag', 'translation']:
                    for value in values:
                        value = value.replace(';', '%3B').replace('=', '%3D').replace('&', '%26')
                        attributes.append(f"{key}={value}")
            
            # Format GFF line
            strand = feature['strand']
            gff_line = [
                record_id,
                source,
                feature_type,
                str(feature['start'] + 1),  # Convert to 1-based
                str(feature['end']),
                '.',
                strand,
                '.',
                ";".join(attributes)
            ]
            gff_lines.append("\t".join(gff_line))
        
        # Process all other features
        for feature in features:
            if feature['type'] in ['gene', 'CDS', 'mat_peptide', 'source']:
                continue
                
            feature_type = feature['type']
            if feature_type not in actual_counter:
                actual_counter[feature_type] = 0
            actual_counter[feature_type] += 1
            
            # Generate simple ID
            feature_id = f"{feature_type.lower()}_{actual_counter[feature_type]}"
            
            # Find parent gene
            parent_id = None
            if 'gene' in feature['qualifiers'] and feature['qualifiers']['gene'][0] in feature_ids:
                parent_id = feature_ids[feature['qualifiers']['gene'][0]]
            elif 'locus_tag' in feature['qualifiers'] and feature['qualifiers']['locus_tag'][0] in feature_ids:
                parent_id = feature_ids[feature['qualifiers']['locus_tag'][0]]
            
            # Format attributes
            attributes = [f"ID={feature_id}"]
            if parent_id:
                attributes.append(f"Parent={parent_id}")
            
            # Add product if available for better display
            if 'product' in feature['qualifiers']:
                product_name = feature['qualifiers']['product'][0]
                attributes.append(f"product={product_name}")
                attributes.append(f"Name={product_name}")
            
            # Add other qualifiers
            for key, values in feature['qualifiers'].items():
                if key not in ['product']:
                    for value in values:
                        value = value.replace(';', '%3B').replace('=', '%3D').replace('&', '%26')
                        attributes.append(f"{key}={value}")
            
            # Format GFF line
            strand = feature['strand']
            gff_line = [
                record_id,
                source,
                feature_type,
                str(feature['start'] + 1),  # Convert to 1-based
                str(feature['end']),
                '.',
                strand,
                '.',
                ";".join(attributes)
            ]
            gff_lines.append("\t".join(gff_line))
        
        # Write all GFF lines
        for line in gff_lines:
            gff_handle.write(line + '\n')
        
        # Write FASTA sequence
        gff_handle.write('##FASTA\n')
        gff_handle.write(f'>{record_id}\n')
        
        # Write sequence in chunks of 60 characters
        for i in range(0, len(sequence), 60):
            gff_handle.write(sequence[i:i+60] + '\n')
    
    if verbose:
        print(f"Conversion complete. GFF3 file written to {gff_output}")
        print(f"Processed {len(features)} features")

def main():
    """Main function."""
    args = parse_args()
    convert_genbank_to_gff3(args.genbank_file, args.gff_output, args.source, args.verbose)

if __name__ == "__main__":
    main()