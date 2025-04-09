#!/usr/bin/env python3
"""
Convert GenBank to GFF3 with advanced handling for viral polyproteins.

This script specifically handles the complex annotation structures in viral genomes,
particularly those with polyproteins and mature peptides (mat_peptide features).
It creates proper parent-child relationships in the GFF3 output to ensure the
annotation is correctly interpreted by tools like snpEff.

Example usage:
    python3 convert_gb_to_gff_advanced.py input.gb output.gff3
"""

import argparse
import re
import sys
from collections import defaultdict
from Bio import SeqIO
from datetime import datetime

def parse_args():
    """Parse command-line arguments."""
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

def get_gene_name(feature):
    """Extract gene name from feature."""
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    elif 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    elif 'product' in feature.qualifiers:
        # Create a gene name from the product if needed
        return clean_product_name(feature.qualifiers['product'][0])
    return None

def get_product_name(feature):
    """Extract product name from feature."""
    if 'product' in feature.qualifiers:
        return feature.qualifiers['product'][0]
    return None

def get_feature_id(feature, feature_counter, feature_type=None):
    """
    Generate a unique ID for a feature.
    For CDS, use the protein_id if available.
    For other features, combine feature type and a counter.
    
    This enhanced version prioritizes product names for better annotation in snpEff.
    """
    if feature_type is None:
        feature_type = feature.type
    
    # First priority: use product name for CDS and mat_peptide features
    if feature_type in ['CDS', 'mat_peptide'] and 'product' in feature.qualifiers:
        product_name = feature.qualifiers['product'][0]
        clean_name = clean_product_name(product_name)
        # Add a unique identifier to ensure uniqueness
        return f"{clean_name}_{feature_counter[feature_type]}"
    
    # Second priority: use protein_id for CDS features    
    if feature_type == 'CDS' and 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    
    # Third priority: use gene name    
    gene_name = get_gene_name(feature)
    product_name = get_product_name(feature)
    
    # Generate IDs based on available information
    if gene_name and product_name:
        clean_product = clean_product_name(product_name)
        if gene_name.lower() != clean_product.lower():
            # If gene name is not just a cleaned version of the product, include both
            return f"{gene_name}_{clean_product}"
        else:
            return gene_name
    elif gene_name:
        return gene_name
    elif product_name:
        return clean_product_name(product_name)
    else:
        # Fallback to feature type and number
        return f"{feature_type.lower()}_{feature_counter[feature_type]}"

def format_qualifiers(qualifiers):
    """Format qualifiers as GFF3 attributes."""
    attributes = []
    
    # Filter out certain qualifiers that we handle separately or don't want in the output
    skip_qualifiers = set(['translation'])
    
    for key, values in qualifiers.items():
        if key in skip_qualifiers:
            continue
            
        for value in values:
            # Escape special characters
            value = value.replace(';', '%3B')
            value = value.replace('=', '%3D')
            value = value.replace('&', '%26')
            value = value.replace(',', '%2C')
            value = value.replace('\n', ' ')  # Replace newlines with spaces
            
            attributes.append(f"{key}={value}")
    
    return ';'.join(attributes)

def convert_genbank_to_gff3(genbank_file, gff_output, source="GenBank", verbose=False):
    """
    Convert GenBank file to GFF3 with special handling for viral polyproteins.
    """
    # Open input and output files
    with open(gff_output, 'w') as gff_handle:
        # Write GFF3 header
        gff_handle.write('##gff-version 3\n')
        gff_handle.write(f'##date {datetime.now().strftime("%Y-%m-%d")}\n')
        gff_handle.write(f'##source-version convert_gb_to_gff_advanced.py\n')
        
        # Read GenBank record
        record = SeqIO.read(genbank_file, "genbank")
        
        # Write sequence-region header
        gff_handle.write(f'##sequence-region {record.id} 1 {len(record.seq)}\n')
        
        # Track feature counters for generating unique IDs
        feature_counter = defaultdict(int)
        
        # First pass: count features for ID generation
        for feature in record.features:
            feature_counter[feature.type] += 1
        
        # Reset counters for the actual processing
        actual_counter = defaultdict(int)
        
        # Store feature IDs for parent-child relationships
        feature_ids = {}
        mat_peptide_parents = {}
        
        # Dictionary to store polyproteins for looking up parent-child relationships later
        polyproteins = {}
        
        # List to track all lines to be written
        gff_lines = []
        
        # Process all features
        for feature in record.features:
            # Skip source feature, we'll add a custom one
            if feature.type == 'source':
                # Add a source feature line
                source_id = "region_1"
                source_line = [
                    record.id,  # seqid
                    source,     # source
                    "region",   # type
                    "1",        # start
                    str(len(record.seq)),  # end
                    ".",        # score
                    "+",        # strand
                    ".",        # phase
                    f"ID={source_id}"  # attributes
                ]
                gff_lines.append("\t".join(source_line))
                continue
                
            # Initialize GFF fields
            actual_counter[feature.type] += 1
            seqid = record.id
            type_source = source
            feature_type = feature.type
            start = str(feature.location.start + 1)  # 1-based coordinates in GFF
            end = str(feature.location.end)
            score = "."
            
            # Determine strand
            if feature.location.strand == 1:
                strand = "+"
            elif feature.location.strand == -1:
                strand = "-"
            else:
                strand = "."
                
            # Determine phase for CDS features
            if feature_type == "CDS":
                phase = str(0) if 'codon_start' not in feature.qualifiers else str(int(feature.qualifiers['codon_start'][0]) - 1)
            else:
                phase = "."
                
            # Generate a unique ID for this feature
            feature_id = get_feature_id(feature, actual_counter)
            
            # Process gene features specially
            if feature_type == 'gene':
                # Store the gene ID for linking to child features
                gene_name = get_gene_name(feature)
                if gene_name:
                    feature_ids[gene_name] = feature_id
                
                # Set up attributes
                attributes = [f"ID={feature_id}"]
                if 'gene' in feature.qualifiers:
                    attributes.append(f"Name={feature.qualifiers['gene'][0]}")
                if 'locus_tag' in feature.qualifiers:
                    attributes.append(f"locus_tag={feature.qualifiers['locus_tag'][0]}")
                
                # Add additional qualifiers
                for key, values in feature.qualifiers.items():
                    if key not in ['gene', 'locus_tag']:
                        for value in values:
                            value = value.replace(';', '%3B').replace('=', '%3D').replace('&', '%26')
                            attributes.append(f"{key}={value}")
                
                # Assemble and store the line
                gff_line = [
                    seqid, type_source, feature_type, start, end, score, strand, phase, 
                    ";".join(attributes)
                ]
                gff_lines.append("\t".join(gff_line))
                
            # Process CDS features specially
            elif feature_type == 'CDS':
                # Create a unique ID for this CDS - prioritize product name
                if 'product' in feature.qualifiers:
                    product_name = feature.qualifiers['product'][0]
                    clean_name = clean_product_name(product_name)
                    feature_id = f"{clean_name}_{actual_counter[feature_type]}"
                elif 'protein_id' in feature.qualifiers:
                    feature_id = feature.qualifiers['protein_id'][0]
                else:
                    feature_id = f"cds_{actual_counter[feature_type]}"
                
                # Find parent gene if possible
                parent_id = None
                if 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] in feature_ids:
                    parent_id = feature_ids[feature.qualifiers['gene'][0]]
                elif 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] in feature_ids:
                    parent_id = feature_ids[feature.qualifiers['locus_tag'][0]]
                
                # Set up attributes
                attributes = [f"ID={feature_id}"]
                if parent_id:
                    attributes.append(f"Parent={parent_id}")
                
                # Add name attribute based on product for better display in snpEff
                if 'product' in feature.qualifiers:
                    product_name = feature.qualifiers['product'][0]
                    attributes.append(f"product={product_name}")
                    attributes.append(f"Name={product_name}")
                    
                    # Check if this is a polyprotein
                    product = product_name.lower()
                    if 'polyprotein' in product:
                        # Store this polyprotein for reference by mat_peptide features
                        polyproteins[feature_id] = {'feature': feature}
                        if verbose:
                            print(f"Found polyprotein: {feature_id} - {product}")
                
                # Always include protein_id if available
                if 'protein_id' in feature.qualifiers:
                    attributes.append(f"protein_id={feature.qualifiers['protein_id'][0]}")
                
                # Add additional qualifiers
                for key, values in feature.qualifiers.items():
                    if key not in ['protein_id', 'product', 'gene', 'locus_tag', 'translation']:
                        for value in values:
                            value = value.replace(';', '%3B').replace('=', '%3D').replace('&', '%26')
                            attributes.append(f"{key}={value}")
                
                # Assemble and store the line
                gff_line = [
                    seqid, type_source, feature_type, start, end, score, strand, phase, 
                    ";".join(attributes)
                ]
                gff_lines.append("\t".join(gff_line))
                
            # Handle mature peptides specially for viral polyproteins
            elif feature_type == 'mat_peptide':
                # Create a unique ID for this mature peptide
                if 'product' in feature.qualifiers:
                    product_name = clean_product_name(feature.qualifiers['product'][0])
                    feature_id = f"{product_name}_{actual_counter[feature_type]}"
                else:
                    feature_id = f"mat_peptide_{actual_counter[feature_type]}"
                
                # Find parent gene or polyprotein
                parent_id = None
                
                # First check for gene parent
                if 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] in feature_ids:
                    parent_id = feature_ids[feature.qualifiers['gene'][0]]
                elif 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] in feature_ids:
                    parent_id = feature_ids[feature.qualifiers['locus_tag'][0]]
                
                # If we don't have a gene parent, try to identify a polyprotein parent
                if not parent_id:
                    # Find a polyprotein CDS that contains this mature peptide
                    for poly_id, poly_data in polyproteins.items():
                        poly_feature = poly_data['feature']
                        if (feature.location.start >= poly_feature.location.start and 
                            feature.location.end <= poly_feature.location.end):
                            parent_id = poly_id
                            break
                
                # Set up attributes
                attributes = [f"ID={feature_id}"]
                if parent_id:
                    attributes.append(f"Parent={parent_id}")
                    # Register this mat_peptide with its parent for lookup
                    if parent_id not in mat_peptide_parents:
                        mat_peptide_parents[parent_id] = []
                    mat_peptide_parents[parent_id].append(feature_id)
                
                if 'product' in feature.qualifiers:
                    product_name = feature.qualifiers['product'][0]
                    attributes.append(f"product={product_name}")
                    attributes.append(f"Name={product_name}")
                
                # Add additional qualifiers
                for key, values in feature.qualifiers.items():
                    if key not in ['product', 'gene', 'locus_tag', 'translation']:
                        for value in values:
                            value = value.replace(';', '%3B').replace('=', '%3D').replace('&', '%26')
                            attributes.append(f"{key}={value}")
                
                # Assemble and store the line
                gff_line = [
                    seqid, type_source, feature_type, start, end, score, strand, phase, 
                    ";".join(attributes)
                ]
                gff_lines.append("\t".join(gff_line))
                
            # Process other features generically
            else:
                # Generate a unique ID
                if feature_type == 'source':
                    feature_id = f"region_{actual_counter[feature_type]}"
                else:
                    feature_id = f"{feature_type.lower()}_{actual_counter[feature_type]}"
                
                # Find parent gene if possible
                parent_id = None
                if 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] in feature_ids:
                    parent_id = feature_ids[feature.qualifiers['gene'][0]]
                elif 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] in feature_ids:
                    parent_id = feature_ids[feature.qualifiers['locus_tag'][0]]
                
                # Set up attributes
                attributes = [f"ID={feature_id}"]
                if parent_id:
                    attributes.append(f"Parent={parent_id}")
                
                # Add remaining qualifiers
                other_attributes = format_qualifiers(feature.qualifiers)
                if other_attributes:
                    attributes.append(other_attributes)
                
                # Assemble and store the line
                gff_line = [
                    seqid, type_source, feature_type, start, end, score, strand, phase, 
                    ";".join(attributes)
                ]
                gff_lines.append("\t".join(gff_line))
        
        # Write all GFF lines
        for line in gff_lines:
            gff_handle.write(line + '\n')
        
        # Write FASTA sequence (optional for tools that need it)
        gff_handle.write('##FASTA\n')
        gff_handle.write(f'>{record.id}\n')
        
        # Write sequence in chunks of 60 characters
        seq_str = str(record.seq)
        for i in range(0, len(seq_str), 60):
            gff_handle.write(seq_str[i:i+60] + '\n')
            
    if verbose:
        print(f"Conversion complete. GFF3 file written to {gff_output}")
        print(f"Processed {sum(actual_counter.values())} features:")
        for feature_type, count in actual_counter.items():
            print(f"  {feature_type}: {count}")

def main():
    """Main function."""
    args = parse_args()
    convert_genbank_to_gff3(args.genbank_file, args.gff_output, args.source, args.verbose)

if __name__ == "__main__":
    main()