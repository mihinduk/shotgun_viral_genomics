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
    Handles binary files and non-UTF-8 encodings that can occur in VCF files.
    
    Args:
        vcf_path: Path to VCF file to fix
        
    Returns:
        Path to fixed VCF file
    """
    logger.info(f"Fixing VCF file for snpEff: {vcf_path}")
    
    # Create a temporary file for the fixed VCF
    fixed_vcf = f"{vcf_path}.fixed"
    
    # Try to use most permissive reading mode
    try:
        # Method 1: First try reading line by line in binary mode, with robust pattern matching
        with open(vcf_path, 'rb') as infile, open(fixed_vcf, 'wb') as outfile:
            line_num = 0
            skipped_lines = 0
            
            # Patterns for problematic variant lines in binary mode
            # Look for tab + single char + tab or empty alt field
            iub_pattern = re.compile(b'\t[RYSWKMBDHV]\t')
            empty_alt_pattern = re.compile(b'\t\t')
            header_pattern = re.compile(b'^#')
            
            for line in infile:
                line_num += 1
                
                # Always keep header lines (starting with #)
                if header_pattern.match(line):
                    outfile.write(line)
                    continue
                
                # Skip lines with IUB codes or empty ALT fields
                if iub_pattern.search(line) or empty_alt_pattern.search(line):
                    skipped_lines += 1
                    logger.debug(f"Skipping problematic variant at line {line_num}")
                    continue
                
                # Write good lines to output
                outfile.write(line)
        
        if skipped_lines > 0:
            logger.info(f"Removed {skipped_lines} problematic variant lines from VCF file")
        
        # Replace original with fixed version
        os.rename(fixed_vcf, vcf_path)
        return vcf_path
        
    except Exception as e:
        logger.warning(f"Binary processing failed: {str(e)}")
        logger.warning("Trying alternative method...")
        
        # Method 2: If binary processing fails, try using the 'latin-1' encoding which accepts any byte value
        try:
            with open(vcf_path, 'r', encoding='latin-1', errors='replace') as infile, \
                 open(fixed_vcf, 'w', encoding='latin-1') as outfile:
                
                line_num = 0
                skipped_lines = 0
                
                # Compile patterns for problematic lines
                problematic_pattern = re.compile(r'\t[RYSWKMBDHV]\t|\t\t')
                
                for line in infile:
                    line_num += 1
                    
                    # Always keep header lines
                    if line.startswith('#'):
                        outfile.write(line)
                        continue
                    
                    # Skip problematic lines
                    if problematic_pattern.search(line):
                        skipped_lines += 1
                        logger.debug(f"Skipping problematic variant at line {line_num}")
                        continue
                    
                    # Write good lines
                    outfile.write(line)
            
            if skipped_lines > 0:
                logger.info(f"Removed {skipped_lines} problematic variant lines using latin-1 encoding")
            
            # Replace original with fixed version
            os.rename(fixed_vcf, vcf_path)
            return vcf_path
            
        except Exception as e:
            logger.error(f"Failed to process VCF file with both methods: {str(e)}")
            
            # If all else fails, try to use a direct regex-based approach on raw bytes
            try:
                import re
                import subprocess
                
                logger.warning("Attempting direct grep-based filtering...")
                
                # Use grep to remove problematic lines - this works on binary files
                grep_cmd = f"grep -v -P '\\t[RYSWKMBDHV]\\t|\\t\\t' {vcf_path} > {fixed_vcf}.tmp"
                grep_header = f"grep '^#' {vcf_path} > {fixed_vcf}"
                grep_content = f"grep -v '^#' {fixed_vcf}.tmp >> {fixed_vcf}"
                
                subprocess.run(grep_header, shell=True, check=True)
                subprocess.run(grep_cmd, shell=True, check=True)
                subprocess.run(grep_content, shell=True, check=True)
                
                # Clean up temporary file
                if os.path.exists(f"{fixed_vcf}.tmp"):
                    os.remove(f"{fixed_vcf}.tmp")
                
                logger.info("Used grep-based filtering to remove problematic lines")
                
                # Replace original with fixed version
                os.rename(fixed_vcf, vcf_path)
                return vcf_path
                
            except Exception as e:
                logger.error(f"All VCF fixing methods failed: {str(e)}")
                if os.path.exists(fixed_vcf):
                    os.remove(fixed_vcf)
                raise
    
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