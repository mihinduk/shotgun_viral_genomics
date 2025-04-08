#!/usr/bin/env python3
"""
Fix empty snpEff annotation VCF files by running snpEff directly with additional debug options.
"""

import os
import sys
import argparse
import subprocess
import shutil
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def run_command(cmd, shell=False, check=True):
    """
    Run a shell command and log the output.
    
    Args:
        cmd: Command to run (string or list of arguments)
        shell: Whether to run the command in a shell
        check: Whether to check the return code
        
    Returns:
        CompletedProcess instance with stdout and stderr
    """
    logger.debug(f"Running command: {cmd}")
    
    result = subprocess.run(
        cmd, 
        shell=shell, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True,
        check=False
    )
    
    if result.returncode != 0:
        logger.error(f"Command failed with exit code {result.returncode}")
        logger.error(f"STDOUT: {result.stdout}")
        logger.error(f"STDERR: {result.stderr}")
        if check:
            raise subprocess.SubprocessError(f"Command failed: {cmd}")
    
    return result

def fix_annotation(input_vcf, genome_id, output_vcf=None, snpeff_jar="snpEff.jar", java_path="java"):
    """
    Fix annotation by running snpEff directly with additional parameters
    """
    if not output_vcf:
        output_vcf = f"{os.path.splitext(input_vcf)[0]}.snpEFF.ann.vcf"
    
    # Check if input VCF exists and has variants
    if not os.path.exists(input_vcf):
        logger.error(f"Input VCF file not found: {input_vcf}")
        return False
    
    logger.info(f"Checking if VCF file has variants: {input_vcf}")
    variant_count_cmd = f"grep -v '^#' {input_vcf} | wc -l"
    variant_result = run_command(variant_count_cmd, shell=True, check=False)
    variant_count = 0
    
    try:
        variant_count = int(variant_result.stdout.strip())
    except (ValueError, AttributeError):
        pass
    
    logger.info(f"Found {variant_count} variants in input VCF")
    
    if variant_count == 0:
        logger.warning("No variants found in input VCF. Creating empty annotation file.")
        with open(output_vcf, 'w') as f:
            f.write("")
        return False
    
    # Verify genome exists in snpEff database
    verify_cmd = f"{java_path} -jar {snpeff_jar} databases | grep -i {genome_id}"
    verify_result = run_command(verify_cmd, shell=True, check=False)
    
    if verify_result.returncode != 0:
        logger.error(f"Genome {genome_id} not found in snpEff database")
        # List available corona genomes
        logger.info("Available coronavirus genomes:")
        run_command(f"{java_path} -jar {snpeff_jar} databases | grep -i corona", shell=True, check=False)
        return False
    else:
        logger.info(f"Found genome {genome_id} in snpEff database")
    
    # Create a safe copy of the VCF file
    safe_copy = f"{input_vcf}.safe"
    shutil.copy2(input_vcf, safe_copy)
    
    # Check chromosome names
    get_chrom_cmd = f"grep -v '^#' {safe_copy} | head -1 | cut -f1"
    chrom_result = run_command(get_chrom_cmd, shell=True, check=False)
    chrom_name = chrom_result.stdout.strip() if chrom_result.returncode == 0 else "UNKNOWN"
    logger.info(f"Chromosome name in VCF: {chrom_name}")
    
    # Run snpEff with verbose output
    summary_html = f"{os.path.splitext(output_vcf)[0]}_summary.html"
    cmd = f"{java_path} -jar -Xmx4g {snpeff_jar} -v {genome_id} {safe_copy} -s {summary_html} > {output_vcf}"
    
    logger.info(f"Running snpEff command: {cmd}")
    try:
        result = run_command(cmd, shell=True)
        logger.info("snpEff annotation completed successfully")
    except subprocess.SubprocessError as e:
        logger.error(f"snpEff annotation failed: {str(e)}")
        
        # Try with more aggressive filtering
        logger.info("Trying with more aggressive filtering...")
        filtered_vcf = f"{safe_copy}.filtered"
        
        # Extract header
        run_command(f"grep '^#' {safe_copy} > {filtered_vcf}", shell=True)
        
        # Filter content lines to only include standard SNPs
        filter_cmd = (
            f"grep -v '^#' {safe_copy} | "
            f"grep -E '^[^[:space:]]+[[:space:]]+[0-9]+[[:space:]]+\\.[[:space:]]+[ACGT][[:space:]]+[ACGT][[:space:]]' >> {filtered_vcf}"
        )
        run_command(filter_cmd, shell=True, check=False)
        
        # Count filtered variants
        filter_count_cmd = f"grep -v '^#' {filtered_vcf} | wc -l"
        filter_result = run_command(filter_count_cmd, shell=True, check=False)
        filter_count = 0
        try:
            filter_count = int(filter_result.stdout.strip())
        except (ValueError, AttributeError):
            pass
        
        logger.info(f"After filtering: {filter_count} variants (originally {variant_count})")
        
        if filter_count > 0:
            # Try again with filtered VCF
            retry_cmd = f"{java_path} -jar -Xmx4g {snpeff_jar} -v {genome_id} {filtered_vcf} -s {summary_html} > {output_vcf}"
            try:
                logger.info(f"Retrying snpEff with filtered VCF: {retry_cmd}")
                run_command(retry_cmd, shell=True)
                logger.info("Retry successful")
            except subprocess.SubprocessError as retry_error:
                logger.error(f"Retry also failed: {str(retry_error)}")
                # Create empty output file
                with open(output_vcf, 'w') as f:
                    f.write("")
                if os.path.exists(filtered_vcf):
                    os.remove(filtered_vcf)
                return False
            
            # Clean up filtered VCF
            if os.path.exists(filtered_vcf):
                os.remove(filtered_vcf)
        else:
            logger.warning("No variants left after filtering")
            # Create empty output file
            with open(output_vcf, 'w') as f:
                f.write("")
            if os.path.exists(filtered_vcf):
                os.remove(filtered_vcf)
            return False
    
    # Check if the output contains variants
    check_cmd = f"grep -v '^#' {output_vcf} | wc -l"
    check_result = run_command(check_cmd, shell=True, check=False)
    output_count = 0
    
    try:
        output_count = int(check_result.stdout.strip())
    except (ValueError, AttributeError):
        pass
    
    logger.info(f"Annotation contains {output_count} annotated variants")
    
    # Clean up
    if os.path.exists(safe_copy):
        os.remove(safe_copy)
    
    return output_count > 0

def main():
    parser = argparse.ArgumentParser(description="Fix empty snpEff annotation VCF files")
    parser.add_argument("input_vcf", help="Input VCF file to annotate")
    parser.add_argument("genome_id", help="Genome ID/accession for snpEff")
    parser.add_argument("--output", "-o", help="Output VCF file (default: input_file_base.snpEFF.ann.vcf)")
    parser.add_argument("--snpeff-jar", default="snpEff.jar", help="Path to snpEff.jar")
    parser.add_argument("--java-path", default="java", help="Path to Java executable")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    success = fix_annotation(
        args.input_vcf,
        args.genome_id,
        args.output,
        args.snpeff_jar,
        args.java_path
    )
    
    if success:
        logger.info("Annotation completed successfully")
        sys.exit(0)
    else:
        logger.error("Annotation failed")
        sys.exit(1)

if __name__ == "__main__":
    main()