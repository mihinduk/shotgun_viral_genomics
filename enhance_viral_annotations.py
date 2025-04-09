#!/usr/bin/env python3
"""
Enhance Viral Annotations Pipeline

A comprehensive script to handle complex viral annotations in shotgun sequencing data.
This script integrates all the tools needed to properly process viral genomes with
complex features like polyproteins and mature peptides.

It can:
1. Convert GenBank to GFF3 with advanced handling of polyproteins
2. Add genomes to snpEff with the appropriate approach (GenBank or GFF3)
3. Run full annotation on VCF files

Copyright (c) 2024
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Union

__version__ = "0.1.0"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Enhance Viral Annotations Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    
    # Command subparsers
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Convert GenBank to GFF3
    convert_parser = subparsers.add_parser("convert", help="Convert GenBank to GFF3")
    convert_parser.add_argument("--genbank", required=True, help="Input GenBank file")
    convert_parser.add_argument("--output", required=True, help="Output GFF3 file")
    
    # Add genome to snpEff
    add_parser = subparsers.add_parser("add-genome", help="Add genome to snpEff")
    add_parser.add_argument("--accession", required=True, help="Genome accession number")
    add_parser.add_argument("--reference", help="Path to local reference FASTA (optional)")
    add_parser.add_argument("--snpeff-jar", default="snpEff.jar", help="Path to snpEff.jar")
    add_parser.add_argument("--java-path", default="java", help="Path to Java executable")
    add_parser.add_argument("--force-gff3", action="store_true", help="Force using GFF3 approach")
    
    # Annotate VCF
    annotate_parser = subparsers.add_parser("annotate", help="Annotate VCF file")
    annotate_parser.add_argument("--vcf", required=True, help="Input VCF file")
    annotate_parser.add_argument("--accession", required=True, help="Genome accession number")
    annotate_parser.add_argument("--output-dir", help="Output directory (defaults to VCF directory)")
    annotate_parser.add_argument("--snpeff-jar", default="snpEff.jar", help="Path to snpEff.jar")
    annotate_parser.add_argument("--java-path", default="java", help="Path to Java executable")
    annotate_parser.add_argument("--add-genome", action="store_true", help="Add genome to snpEff if not present")
    annotate_parser.add_argument("--min-depth", type=int, default=200, help="Minimum read depth for final variant reporting")
    annotate_parser.add_argument("--skip-parse", action="store_true", help="Skip parsing the annotation")
    
    return parser.parse_args()

def run_command(cmd: Union[str, List[str]], shell: bool = True, check: bool = True) -> subprocess.CompletedProcess:
    """
    Run a shell command and log the output.
    
    Args:
        cmd: Command to run (string or list of arguments)
        shell: Whether to run the command in a shell
        check: Whether to check the return code
        
    Returns:
        CompletedProcess instance with stdout and stderr
    """
    if isinstance(cmd, list):
        cmd_str = " ".join(cmd)
    else:
        cmd_str = cmd
        
    logger.debug(f"Running command: {cmd_str}")
    
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
            raise subprocess.SubprocessError(f"Command failed: {cmd_str}")
    
    return result

def check_genome_complexity(accession: str, fasta_path: str) -> bool:
    """
    Check if a genome requires custom GFF3 processing due to complex annotation structure.
    
    Args:
        accession: Genome accession number
        fasta_path: Path to the genome FASTA file
        
    Returns:
        True if genome needs custom GFF3, False if standard GenBank processing is adequate
    """
    logger.info(f"Analyzing annotation complexity for genome: {accession}")
    
    # Download GenBank file to analyze annotation structure
    gb_path = fasta_path.replace(".fasta", ".gb")
    if not os.path.exists(gb_path):
        try:
            cmd = f"efetch -db nucleotide -id {accession} -format gb > {gb_path}"
            run_command(cmd, shell=True, check=False)
        except:
            logger.warning(f"Could not download GenBank for complexity analysis. Using standard processing.")
            return False
    
    try:
        # Check for specific complex annotation patterns
        complex_patterns = [
            # Check for polyproteins
            "codon_start=1.*product=\".*polyprotein.*\"",
            # Check for mat_peptide features
            "mat_peptide[ ]+[0-9]+\.\.[0-9]+",
            # Check for CDS with gene name in product
            "CDS[ ]+[0-9]+\.\.[0-9]+.*product=\".*protein.*\"",
            # Check for multiple genes with same product
            "product=\"nonstructural polyprotein.*\".*product=\"nonstructural polyprotein",
            # Check for POLY gene
            "gene=\"POLY\"",
            # Check for Alphavirus-like patterns
            "product=\".*nsP[1-4].*\"",
            # Check specifically for common virus families that need custom handling
            "Togaviridae|Picornaviridae|Flaviviridae|Coronaviridae"
        ]
        
        # Check each pattern
        needs_custom_gff = False
        for pattern in complex_patterns:
            cmd = f"grep -E '{pattern}' {gb_path}"
            result = run_command(cmd, shell=True, check=False)
            if result.returncode == 0 and result.stdout.strip():
                logger.info(f"Complex annotation pattern detected: {pattern}")
                needs_custom_gff = True
                break
        
        if needs_custom_gff:
            # Get virus family for additional context
            cmd = f"grep -A2 'ORGANISM' {gb_path}"
            family_result = run_command(cmd, shell=True, check=False)
            if family_result.returncode == 0:
                logger.info(f"Taxonomy info: {family_result.stdout.strip()}")
                
            logger.info(f"Genome {accession} requires custom GFF3 processing due to complex annotation structure")
        else:
            logger.info(f"Genome {accession} has standard annotation structure - using regular GenBank processing")
            
        return needs_custom_gff
    
    except Exception as e:
        logger.warning(f"Error analyzing genome complexity: {str(e)}. Using standard processing.")
        return False

def find_script(script_name: str) -> Optional[str]:
    """
    Find a script by searching in common locations.
    
    Args:
        script_name: Name of the script to find
        
    Returns:
        Path to the script if found, None otherwise
    """
    # Get the directory where this Python script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Potential locations for the script
    script_locations = [
        # Current directory
        f"./{script_name}",
        # Same directory as this Python script
        os.path.join(script_dir, script_name),
        # Parent directory of this Python script
        os.path.join(os.path.dirname(script_dir), script_name),
        # Specific locations for the project
        f"/Users/handley_lab/Handley Lab Dropbox/virome/KM_algorithm_dev/shotgun_viral_genomics/{script_name}"
    ]
    
    for location in script_locations:
        if os.path.exists(location):
            logger.info(f"Found script at: {location}")
            return location
    
    return None

def download_reference_genome(accession: str, output_dir: str) -> str:
    """
    Download a reference genome from NCBI using Entrez Direct utilities.
    
    Args:
        accession: The genome accession number (e.g., NC_045512.2)
        output_dir: Directory to save the downloaded genome
        
    Returns:
        Path to the downloaded genome FASTA file
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{accession}.fasta")
    
    # Check if file already exists
    if os.path.exists(output_file):
        logger.info(f"Reference genome file already exists: {output_file}")
        return output_file
    
    logger.info(f"Downloading reference genome: {accession}")
    
    # Try using Entrez Direct utilities
    try:
        cmd = f"efetch -db nucleotide -id {accession} -format fasta > {output_file}"
        logger.debug(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=True)
    except (subprocess.SubprocessError, FileNotFoundError):
        # Fall back to wget if efetch is not available
        logger.info("efetch not found, falling back to NCBI E-utilities via wget")
        try:
            ncbi_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=fasta&retmode=text"
            cmd = f"wget -O {output_file} \"{ncbi_url}\""
            logger.debug(f"Running command: {cmd}")
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.SubprocessError:
            raise RuntimeError(f"Failed to download reference genome: {accession}")
    
    # Verify the downloaded file
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        raise RuntimeError(f"Downloaded reference genome file is empty or missing: {output_file}")
    
    logger.info(f"Reference genome downloaded to: {output_file}")
    return output_file

def generate_custom_gff3(accession: str, gb_path: str, output_dir: str) -> Optional[str]:
    """
    Generate custom GFF3 for complex viral genomes.
    
    Args:
        accession: Genome accession number
        gb_path: Path to the GenBank file
        output_dir: Directory to save the GFF3 file
        
    Returns:
        Path to the generated GFF3 file or None if failed
    """
    logger.info(f"Generating custom GFF3 for genome: {accession}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find the advanced GFF3 converter script
    converter_script = find_script("convert_gb_to_gff_advanced.py")
    
    if not converter_script:
        logger.error("Could not find convert_gb_to_gff_advanced.py script")
        return None
    
    # Generate the GFF3 file
    gff3_path = os.path.join(output_dir, f"{accession}.gff3")
    
    try:
        cmd = f"python3 {converter_script} {gb_path} {gff3_path} --verbose"
        run_command(cmd, shell=True)
        
        # Verify the GFF3 file was created and is not empty
        if os.path.exists(gff3_path) and os.path.getsize(gff3_path) > 0:
            logger.info(f"Successfully generated custom GFF3 file: {gff3_path}")
            return gff3_path
        else:
            logger.error(f"GFF3 generation failed or produced empty file")
            return None
    except Exception as e:
        logger.error(f"Error generating custom GFF3: {str(e)}")
        return None

def check_snpeff_database(accession: str, snpeff_jar: str, java_path: str = "java") -> bool:
    """
    Check if a genome is in the snpEff database.
    
    Args:
        accession: Genome accession number
        snpeff_jar: Path to snpEff.jar
        java_path: Path to Java executable
        
    Returns:
        True if the genome is in the database, False otherwise
    """
    logger.info(f"Checking if {accession} is in snpEff database")
    
    cmd = f"{java_path} -jar {snpeff_jar} databases | grep {accession}"
    result = run_command(cmd, shell=True, check=False)
    
    return result.returncode == 0

def add_genome_to_snpeff(accession: str, fasta_path: str, snpeff_jar: str, java_path: str = "java", force_gff3: bool = False) -> bool:
    """
    Add a genome to the snpEff database using either GenBank or GFF3 approach.
    
    Args:
        accession: Genome accession number
        fasta_path: Path to the genome FASTA file
        snpeff_jar: Path to snpEff.jar
        java_path: Path to Java executable
        force_gff3: Force using GFF3 approach even if genome is not complex
        
    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Adding genome {accession} to snpEff database")
    
    # Extract snpEff directory
    snpeff_dir = os.path.dirname(os.path.abspath(snpeff_jar))
    
    # Check for genome complexity
    needs_custom_gff = force_gff3 or check_genome_complexity(accession, fasta_path)
    
    if needs_custom_gff:
        logger.info(f"Using custom GFF3 approach for genome {accession}")
        
        # Generate GenBank file path
        gb_path = fasta_path.replace(".fasta", ".gb")
        if not os.path.exists(gb_path):
            cmd = f"efetch -db nucleotide -id {accession} -format gb > {gb_path}"
            run_command(cmd, shell=True)
        
        # Generate custom GFF3
        temp_dir = os.path.dirname(fasta_path)
        gff3_path = generate_custom_gff3(accession, gb_path, temp_dir)
        
        if not gff3_path:
            logger.warning("Custom GFF3 generation failed, falling back to standard GenBank approach")
            cmd = f"{java_path} -jar {snpeff_jar} build -genbank -v -noCheckProtein {accession}"
            result = run_command(cmd, shell=True, check=False)
        else:
            # Set up snpEff with GFF3 and FASTA
            data_dir = os.path.join(snpeff_dir, "data", accession)
            os.makedirs(data_dir, exist_ok=True)
            
            # Copy FASTA file to snpEff data directory
            sequences_path = os.path.join(data_dir, "sequences.fa")
            shutil.copy2(fasta_path, sequences_path)
            
            # Copy GFF3 file to snpEff data directory
            genes_path = os.path.join(data_dir, "genes.gff")
            shutil.copy2(gff3_path, genes_path)
            
            # Modify snpEff.config
            config_path = os.path.join(snpeff_dir, "snpEff.config")
            
            # Check if accession already exists in config
            check_cmd = f"grep -E '^{accession}\\.genome:' {config_path}"
            check_result = run_command(check_cmd, shell=True, check=False)
            
            if check_result.returncode != 0:
                # Extract organism and strain from GenBank
                org_cmd = f"grep -A1 'ORGANISM' {gb_path} | tail -1 | tr -d '[:space:]'"
                org_result = run_command(org_cmd, shell=True, check=False)
                organism = org_result.stdout.strip() if org_result.returncode == 0 else accession
                
                strain_cmd = f"grep -E 'strain=\"[^\"]+\"' {gb_path} | head -1 | sed 's/.*strain=\"\\([^\"]*\\)\".*/\\1/'"
                strain_result = run_command(strain_cmd, shell=True, check=False)
                strain = strain_result.stdout.strip() if strain_result.returncode == 0 else ""
                
                # Create abbreviated name
                abbr_name = accession
                if strain:
                    # Extract first letters of each word in organism
                    organism_abbr = ''.join([word[0].upper() for word in organism.split() if word[0].isalpha()])
                    abbr_name = f"{organism_abbr}-{strain}"
                
                # Add to config
                with open(config_path, 'a') as f:
                    f.write(f"\n# {accession}\n")
                    f.write(f"{accession}.genome: {abbr_name}\n")
                    f.write(f"{accession}.chromosomes: {accession}\n")
                    f.write(f"{accession}.codonTable: Standard\n")
            
            # Build the database with more forgiving parameters for viral genomes
            cmd = f"{java_path} -jar {snpeff_jar} build -gff3 -v -noCheckProtein -noCheckCds -noLog -treatAllAsProteinCoding {accession}"
            result = run_command(cmd, shell=True, check=False)
    else:
        # Use standard GenBank approach with forgiving parameters
        logger.info(f"Using standard GenBank approach for genome {accession}")
        cmd = f"{java_path} -jar {snpeff_jar} build -genbank -v -noCheckProtein -noCheckCds -noLog -treatAllAsProteinCoding {accession}"
        result = run_command(cmd, shell=True, check=False)
    
    # Check if build was successful
    success = result.returncode == 0
    if success:
        logger.info(f"Successfully added {accession} to snpEff database")
    else:
        logger.error(f"Failed to add {accession} to snpEff database")
        logger.error(f"STDOUT: {result.stdout}")
        logger.error(f"STDERR: {result.stderr}")
    
    return success

def fix_vcf_for_snpeff(vcf_path: str) -> str:
    """
    Fix VCF file to remove problematic IUB ambiguity code lines that cause snpEff to fail.
    
    Args:
        vcf_path: Path to VCF file to fix
        
    Returns:
        Path to fixed VCF file
    """
    # First check if our helper script exists
    fix_script = find_script("fix_vcf_for_snpeff.py")
    
    if fix_script:
        # Use the helper script
        logger.info(f"Using fix_vcf_for_snpeff.py script to fix VCF: {vcf_path}")
        cmd = f"python3 {fix_script} {vcf_path}"
        run_command(cmd, shell=True)
    else:
        # Direct approach
        logger.info(f"Using grep-based filtering to remove problematic lines from VCF: {vcf_path}")
        # Create a temporary file for the fixed VCF
        fixed_vcf = f"{vcf_path}.fixed"
        
        # Extract header lines (starting with #)
        grep_header_cmd = f"grep '^#' {vcf_path} > {fixed_vcf}"
        run_command(grep_header_cmd, shell=True)
        
        # Extract content lines, filtering out problematic ones
        # - Skip lines with IUB ambiguity codes in the ALT field (tab + IUB code + tab)
        # - Skip lines with empty ALT field (tab + tab)
        grep_content_cmd = f"grep -v '^#' {vcf_path} | grep -v -E '\\t[RYSWKMBDHV]\\t' | grep -v '\\t\\t' >> {fixed_vcf}"
        run_command(grep_content_cmd, shell=True, check=False)
        
        # Count filtered lines by comparing original with fixed
        orig_count = int(run_command(f"grep -v '^#' {vcf_path} | wc -l", shell=True).stdout.strip())
        fixed_count = int(run_command(f"grep -v '^#' {fixed_vcf} | wc -l", shell=True).stdout.strip())
        filtered_count = orig_count - fixed_count
        
        logger.info(f"Removed {filtered_count} problematic variant lines from VCF file")
        
        # Replace original with fixed version
        os.rename(fixed_vcf, vcf_path)
    
    return vcf_path

def annotate_variants(vcf_file: str, accession: str, snpeff_jar: str, java_path: str = "java", min_depth: int = 200, skip_parse: bool = False) -> Dict[str, str]:
    """
    Run snpEff annotation and process results.
    
    Args:
        vcf_file: Path to filtered VCF file
        accession: Reference genome accession
        snpeff_jar: Path to snpEff.jar
        java_path: Path to Java executable
        min_depth: Minimum read depth for filtering in parsing step
        skip_parse: Skip the parsing step
        
    Returns:
        Dictionary mapping output file types to paths
    """
    # Extract sample name and directory from VCF file
    sample_name = os.path.basename(vcf_file).replace("_vars.filt.vcf", "")
    output_dir = os.path.dirname(vcf_file)
    
    # Output file paths
    ann_vcf = os.path.join(output_dir, f"{sample_name}.snpEFF.ann.vcf")
    ann_tsv = os.path.join(output_dir, f"{sample_name}.snpEFF.ann.tsv")
    summary_html = os.path.join(output_dir, f"{sample_name}_summary.html")
    summary_genes = os.path.join(output_dir, f"{sample_name}_summary.genes.txt")
    
    logger.info(f"Annotating variants for {sample_name}")
    
    # Fix VCF file to remove problematic IUB ambiguity codes
    try:
        fix_vcf_for_snpeff(vcf_file)
    except Exception as e:
        logger.warning(f"Could not fix VCF file: {str(e)}")
    
    # Make a copy of the VCF file for safety
    safe_vcf = f"{vcf_file}.safe"
    shutil.copy2(vcf_file, safe_vcf)
    
    # Run snpEff annotation
    cmd = f"{java_path} -jar -Xmx4g {snpeff_jar} -v {accession} {safe_vcf} -s {summary_html} > {ann_vcf}"
    logger.info(f"Running command: {cmd}")
    
    try:
        result = run_command(cmd, shell=True)
        logger.info(f"SnpEff annotation completed successfully")
    except subprocess.SubprocessError as e:
        logger.error(f"SnpEff annotation failed: {str(e)}")
        # Create empty output files to prevent errors downstream
        open(ann_vcf, 'w').close()
        open(ann_tsv, 'w').close()
        return {
            'ann_vcf': ann_vcf,
            'ann_tsv': ann_tsv,
            'summary_html': summary_html,
            'summary_genes': summary_genes
        }
    
    # Process the VCF file into TSV format
    logger.info(f"Converting VCF to TSV format for {sample_name}")
    
    # Create temporary files
    tmp_base = os.path.join(output_dir, f"{sample_name}.ann.base.vcf")
    tmp_info = os.path.join(output_dir, f"{sample_name}.snpEFF.ann.tmp")
    
    try:
        # Extract annotations
        extract_cmd = (
            f"grep -v '^##' {ann_vcf} | "
            f"tail -n+2 | "
            f"cut -f8 | "
            f"sed 's/|/\\t/g' | "
            f"cut -f1-16 | "
            f"sed '1i INFO\\tEFFECT\\tPUTATIVE_IMPACT\\tGENE_NAME\\tGENE_ID\\tFEATURE_TYPE\\tFEATURE_ID\\tTRANSCRIPT_TYPE\\tEXON_INTRON_RANK\\tHGVSc\\tHGVSp\\tcDNA_POSITION_AND_LENGTH\\tCDS_POSITION_AND_LENGTH\\tPROTEIN_POSITION_AND_LENGTH\\tDISTANCE_TO_FEATURE\\tERROR' > {tmp_info}"
        )
        run_command(extract_cmd, shell=True, check=False)
        
        # Extract base VCF information
        base_cmd = f"grep -v '^##' {ann_vcf} | cut -f1-7 > {tmp_base}"
        run_command(base_cmd, shell=True, check=False)
        
        # Combine into final TSV
        combine_cmd = f"paste {tmp_base} {tmp_info} > {ann_tsv}"
        run_command(combine_cmd, shell=True, check=False)
    except Exception as e:
        logger.error(f"Error converting VCF to TSV: {str(e)}")
        logger.warning("Creating empty TSV file")
        header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tEFFECT\tPUTATIVE_IMPACT\tGENE_NAME\tGENE_ID\tFEATURE_TYPE\tFEATURE_ID\tTRANSCRIPT_TYPE\tEXON_INTRON_RANK\tHGVSc\tHGVSp\tcDNA_POSITION_AND_LENGTH\tCDS_POSITION_AND_LENGTH\tPROTEIN_POSITION_AND_LENGTH\tDISTANCE_TO_FEATURE\tERROR"
        with open(ann_tsv, 'w') as f:
            f.write(header + "\n")
    finally:
        # Cleanup temporary files
        if os.path.exists(tmp_base):
            os.remove(tmp_base)
        if os.path.exists(tmp_info):
            os.remove(tmp_info)
        if os.path.exists(safe_vcf):
            os.remove(safe_vcf)
    
    # Run parse script if available and not skipped
    if not skip_parse and os.path.exists(ann_tsv) and os.path.getsize(ann_tsv) > 0:
        logger.info(f"Parsing annotations with minimum depth {min_depth}")
        
        # Find the parse script
        parse_script = find_script("parse_snpEff_annotated_vcf_for_collaborators.pl")
        
        if parse_script:
            try:
                cmd = f"perl {parse_script} {ann_tsv} {min_depth}"
                run_command(cmd, shell=True)
                logger.info(f"Successfully parsed annotations with min depth {min_depth}")
            except Exception as e:
                logger.error(f"Error parsing annotations: {str(e)}")
        else:
            logger.warning("Parse script not found, skipping parsing step")
    
    return {
        'ann_vcf': ann_vcf,
        'ann_tsv': ann_tsv,
        'summary_html': summary_html,
        'summary_genes': summary_genes
    }

def convert_gb_to_gff3(genbank_file: str, gff_output: str) -> None:
    """
    Convert GenBank file to GFF3 with advanced handling for viral genomes.
    
    Args:
        genbank_file: Path to input GenBank file
        gff_output: Path to output GFF3 file
    """
    # Find the converter script
    converter_script = find_script("convert_gb_to_gff_advanced.py")
    
    if not converter_script:
        raise FileNotFoundError("convert_gb_to_gff_advanced.py script not found")
    
    # Run the converter
    cmd = f"python3 {converter_script} {genbank_file} {gff_output} --verbose"
    run_command(cmd, shell=True)
    
    # Verify the output file
    if not os.path.exists(gff_output) or os.path.getsize(gff_output) == 0:
        raise RuntimeError(f"Failed to create GFF3 file: {gff_output}")
    
    logger.info(f"Successfully converted GenBank to GFF3: {gff_output}")

def main():
    """Main function to run the pipeline."""
    args = parse_args()
    
    if not args.command:
        logger.error("No command specified. Use --help to see available commands.")
        sys.exit(1)
    
    try:
        if args.command == "convert":
            # Convert GenBank to GFF3
            convert_gb_to_gff3(args.genbank, args.output)
            
        elif args.command == "add-genome":
            # Add genome to snpEff
            # Download genome if not provided
            if args.reference:
                fasta_path = args.reference
            else:
                output_dir = os.path.dirname(os.path.abspath(args.snpeff_jar))
                fasta_path = download_reference_genome(args.accession, output_dir)
            
            # Check if already in database
            in_database = check_snpeff_database(args.accession, args.snpeff_jar, args.java_path)
            
            if in_database:
                logger.info(f"Genome {args.accession} is already in snpEff database")
            else:
                # Add to database
                success = add_genome_to_snpeff(
                    args.accession, 
                    fasta_path, 
                    args.snpeff_jar, 
                    args.java_path,
                    args.force_gff3
                )
                
                if not success:
                    raise RuntimeError(f"Failed to add genome {args.accession} to snpEff database")
            
        elif args.command == "annotate":
            # Annotate VCF file
            # Validate input file
            if not os.path.exists(args.vcf):
                raise FileNotFoundError(f"Input VCF file not found: {args.vcf}")
            
            # Set output directory
            output_dir = args.output_dir if args.output_dir else os.path.dirname(args.vcf)
            os.makedirs(output_dir, exist_ok=True)
            
            # Check if genome is in snpEff database
            in_database = check_snpeff_database(args.accession, args.snpeff_jar, args.java_path)
            
            # Add genome to snpEff if needed
            if not in_database:
                logger.info(f"Genome {args.accession} not found in snpEff database")
                if args.add_genome:
                    # Download genome
                    fasta_path = download_reference_genome(args.accession, output_dir)
                    
                    # Add to database
                    success = add_genome_to_snpeff(
                        args.accession, 
                        fasta_path, 
                        args.snpeff_jar, 
                        args.java_path
                    )
                    
                    if not success:
                        raise RuntimeError(f"Failed to add genome {args.accession} to snpEff database")
                else:
                    raise RuntimeError(f"Genome {args.accession} not in snpEff database. Use --add-genome to add it.")
            
            # Run snpEff annotation
            annotation_files = annotate_variants(
                args.vcf, 
                args.accession, 
                args.snpeff_jar, 
                args.java_path,
                args.min_depth,
                args.skip_parse
            )
            
            logger.info("Annotation completed successfully")
            logger.info(f"Output files:")
            for file_type, file_path in annotation_files.items():
                logger.info(f"  - {file_type}: {file_path}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()