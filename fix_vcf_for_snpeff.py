#!/usr/bin/env python3
"""
Fix VCF files for snpEff annotation by removing problematic IUB ambiguity code lines.
"""

import os
import re
import sys
import argparse
import logging
import subprocess

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
    Handles binary files and non-UTF-8 encodings that can occur in VCF files.
    
    Args:
        vcf_path: Path to VCF file to fix
        
    Returns:
        Path to fixed VCF file
    """
    logger.info(f"Fixing VCF file for snpEff: {vcf_path}")
    
    # Create a temporary file for the fixed VCF
    fixed_vcf = f"{vcf_path}.fixed"
    
    # First approach: Use grep directly (most reliable for binary files)
    try:
        logger.info("Using grep-based filtering to remove problematic lines")
        
        # Extract header lines (starting with #)
        grep_header_cmd = f"grep '^#' {vcf_path} > {fixed_vcf}"
        subprocess.run(grep_header_cmd, shell=True, check=True)
        
        # Extract content lines, filtering out problematic ones
        # - Skip lines with IUB ambiguity codes in the ALT field (tab + IUB code + tab)
        # - Skip lines with empty ALT field (tab + tab)
        grep_content_cmd = f"grep -v '^#' {vcf_path} | grep -v -P '\\t[RYSWKMBDHV]\\t' | grep -v '\\t\\t' >> {fixed_vcf}"
        subprocess.run(grep_content_cmd, shell=True, check=True)
        
        # Count filtered lines by comparing original with fixed
        orig_count = int(subprocess.run(f"grep -v '^#' {vcf_path} | wc -l", shell=True, capture_output=True, text=True).stdout.strip())
        fixed_count = int(subprocess.run(f"grep -v '^#' {fixed_vcf} | wc -l", shell=True, capture_output=True, text=True).stdout.strip())
        filtered_count = orig_count - fixed_count
        
        logger.info(f"Removed {filtered_count} problematic variant lines from VCF file")
        
        # Replace original with fixed version
        os.rename(fixed_vcf, vcf_path)
        return vcf_path
        
    except Exception as e:
        logger.warning(f"Grep-based filtering failed: {str(e)}")
        
        # If the first approach fails, try a second approach with text processing
        try:
            logger.warning("Trying text-based processing with latin-1 encoding")
            
            # Use latin-1 encoding which can handle any byte value
            with open(vcf_path, 'r', encoding='latin-1', errors='replace') as infile, \
                 open(fixed_vcf, 'w', encoding='latin-1') as outfile:
                
                line_num = 0
                skipped_lines = 0
                
                # Pattern for IUB codes in ALT field or empty ALT field
                pattern = re.compile(r'\t[RYSWKMBDHV]\t|\t\t')
                
                for line in infile:
                    line_num += 1
                    
                    # Always keep header lines
                    if line.startswith('#'):
                        outfile.write(line)
                        continue
                    
                    # Skip problematic lines
                    if pattern.search(line):
                        skipped_lines += 1
                        continue
                    
                    # Keep good lines
                    outfile.write(line)
            
            logger.info(f"Removed {skipped_lines} problematic variant lines using text processing")
            
            # Replace original with fixed version
            os.rename(fixed_vcf, vcf_path)
            return vcf_path
            
        except Exception as e:
            logger.error(f"Text-based processing also failed: {str(e)}")
            
            # Clean up any partial files
            if os.path.exists(fixed_vcf):
                os.remove(fixed_vcf)
                
            # Rethrow the exception
            raise RuntimeError(f"Failed to process VCF file: {str(e)}")
    
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