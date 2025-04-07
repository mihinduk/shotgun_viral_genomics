#!/usr/bin/env python3
"""
Fix VCF files for snpEff annotation by removing problematic IUB ambiguity code lines.
"""

import os
import re
import sys
import argparse
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def fix_vcf_for_snpeff(vcf_path):
    """
    Fix VCF file to remove problematic IUB ambiguity code lines that cause snpEff to fail.
    
    Args:
        vcf_path: Path to VCF file to fix
        
    Returns:
        Path to fixed VCF file
    """
    logger.info(f"Fixing VCF file for snpEff: {vcf_path}")
    
    # Create a temporary file for the fixed VCF
    fixed_vcf = f"{vcf_path}.fixed"
    
    # Patterns for problematic variant lines - find lines with empty ALT field or non-standard nucleotides
    # Standard nucleotides are A, C, G, T, N
    # IUB ambiguity codes include R, Y, S, W, K, M, B, D, H, V
    problematic_pattern = re.compile(r'\t[RYSWKMBDHV]\t|\t\t')
    
    with open(vcf_path, 'r') as infile, open(fixed_vcf, 'w') as outfile:
        line_num = 0
        skipped_lines = 0
        
        for line in infile:
            line_num += 1
            
            # Always keep header lines (starting with #)
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            # Check if line has problematic variants
            if problematic_pattern.search(line):
                skipped_lines += 1
                logger.debug(f"Skipping problematic variant at line {line_num}: {line.strip()}")
                continue
            
            # Write good lines to output
            outfile.write(line)
    
    if skipped_lines > 0:
        logger.info(f"Removed {skipped_lines} problematic variant lines from VCF file")
    
    # Replace original with fixed version
    os.rename(fixed_vcf, vcf_path)
    
    return vcf_path

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Fix VCF files for snpEff annotation by removing problematic IUB ambiguity code lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("vcf", help="VCF file to fix")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Fix VCF file
    fix_vcf_for_snpeff(args.vcf)
    
    logger.info("VCF fix completed successfully")

if __name__ == "__main__":
    main()