#!/usr/bin/env python3
"""
Viral Pipeline Troubleshooting Module

This module helps diagnose and resolve common issues in the viral genomics pipeline,
particularly when samples fail to produce expected variants or show poor mapping.
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def run_command(cmd, shell=False, check=True):
    """Run a command and return result."""
    logger.info(f"Running command: {cmd}")
    result = subprocess.run(
        cmd if shell else cmd.split(),
        shell=shell,
        capture_output=True,
        text=True,
        check=check
    )
    return result

def check_mapping_stats(sample_dir, sample_name=None):
    """Check mapping statistics for BAM files."""
    logger.info("=== Checking Mapping Statistics ===")
    
    # Find BAM files
    bam_pattern = "*.final.bam" if not sample_name else f"{sample_name}*.final.bam"
    bam_files = list(Path(sample_dir).rglob(bam_pattern))
    
    if not bam_files:
        logger.warning("No final BAM files found. Looking for any BAM files...")
        bam_files = list(Path(sample_dir).rglob("*.bam"))[:3]  # Limit to first 3
    
    if not bam_files:
        logger.error("No BAM files found!")
        return False
    
    for bam_file in bam_files:
        logger.info(f"\nAnalyzing: {bam_file.name}")
        logger.info(f"Size: {bam_file.stat().st_size / 1024 / 1024:.1f} MB")
        
        # Get mapping stats
        try:
            result = run_command(f"samtools flagstat {bam_file}", shell=True)
            for line in result.stdout.split('\n'):
                if 'mapped' in line or 'properly paired' in line:
                    logger.info(f"  {line.strip()}")
        except:
            logger.warning(f"  Could not get flagstat for {bam_file}")
    
    return True

def check_variant_files(sample_dir, sample_name=None):
    """Check variant files and count variants."""
    logger.info("\n=== Checking Variant Files ===")
    
    vcf_pattern = "*.vcf" if not sample_name else f"{sample_name}*.vcf"
    vcf_files = list(Path(sample_dir).rglob(vcf_pattern))
    
    if not vcf_files:
        logger.error("No VCF files found!")
        return False
    
    for vcf_file in vcf_files:
        logger.info(f"\nVCF File: {vcf_file.name}")
        try:
            # Count variants
            result = run_command(f"grep -v '^#' {vcf_file} | wc -l", shell=True)
            variant_count = int(result.stdout.strip())
            logger.info(f"  Variants: {variant_count}")
            
            if variant_count > 0:
                logger.info("  First few variants:")
                result = run_command(f"grep -v '^#' {vcf_file} | head -3 | cut -f1-5", shell=True)
                for line in result.stdout.strip().split('\n'):
                    logger.info(f"    {line}")
        except:
            logger.warning(f"  Could not analyze {vcf_file}")
    
    return True

def run_assembly_diagnosis(r1_file, r2_file, output_dir, threads=4):
    """Run MEGAHIT assembly to identify the actual organism."""
    logger.info("\n=== Running Assembly Diagnosis ===")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Run MEGAHIT
    megahit_dir = os.path.join(output_dir, "megahit_assembly")
    cmd = (f"megahit -1 {r1_file} -2 {r2_file} -o {megahit_dir} "
           f"-t {threads} --presets meta-sensitive --min-contig-len 500")
    
    try:
        # Remove existing directory if it exists
        if os.path.exists(megahit_dir):
            import shutil
            shutil.rmtree(megahit_dir)
            
        result = run_command(cmd, shell=True)
        
        # Check results
        contigs_file = os.path.join(megahit_dir, "final.contigs.fa")
        if os.path.exists(contigs_file):
            logger.info(f"Assembly completed. Contigs saved to: {contigs_file}")
            
            # Get assembly stats
            try:
                result = run_command(f"seqkit stats {contigs_file}", shell=True)
                logger.info(f"Assembly statistics:\n{result.stdout}")
            except:
                logger.warning("Could not get assembly statistics")
            
            return contigs_file
        else:
            logger.error("Assembly completed but no contigs file found")
            return None
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Assembly failed: {e}")
        return None

def blast_contigs_for_identification(contigs_file, output_dir, num_results=20):
    """BLAST contigs to identify organism."""
    logger.info("\n=== Running BLAST for Organism Identification ===")
    
    blast_output = os.path.join(output_dir, "blast_identification.tsv")
    
    # Use blastn against nt database
    cmd = (f"blastn -query {contigs_file} -db nt -remote "
           f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' "
           f"-max_target_seqs {num_results} -max_hsps 1 > {blast_output}")
    
    try:
        result = run_command(cmd, shell=True)
        logger.info(f"BLAST results saved to: {blast_output}")
        
        # Parse and display results
        logger.info("\n=== BLAST Results Summary ===")
        with open(blast_output, 'r') as f:
            lines = f.readlines()
            
        if not lines:
            logger.warning("No BLAST hits found")
            return blast_output
            
        # Group by organism
        organisms = {}
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) >= 13:
                contig = parts[0]
                accession = parts[1]
                identity = float(parts[2])
                description = parts[12]
                
                # Extract organism name
                organism = "Unknown"
                if "virus" in description.lower():
                    words = description.split()
                    for i, word in enumerate(words):
                        if "virus" in word.lower():
                            organism = " ".join(words[max(0, i-2):i+1])
                            break
                
                if organism not in organisms:
                    organisms[organism] = []
                organisms[organism].append({
                    'contig': contig,
                    'accession': accession,
                    'identity': identity,
                    'description': description
                })
        
        # Display results
        logger.info("\nTop organisms identified:")
        for organism, hits in organisms.items():
            logger.info(f"\n{organism}:")
            best_hits = sorted(hits, key=lambda x: x['identity'], reverse=True)[:3]
            for hit in best_hits:
                logger.info(f"  {hit['contig']} -> {hit['accession']} ({hit['identity']:.1f}% identity)")
                logger.info(f"    {hit['description']}")
        
        # Suggest reference genomes
        logger.info("\n=== Suggested Reference Genomes ===")
        seen_accessions = set()
        for organism, hits in organisms.items():
            best_hit = max(hits, key=lambda x: x['identity'])
            if best_hit['identity'] > 90 and best_hit['accession'] not in seen_accessions:
                logger.info(f"For {organism}: {best_hit['accession']}")
                seen_accessions.add(best_hit['accession'])
        
        return blast_output
        
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST failed: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(
        description="Troubleshoot viral genomics pipeline issues",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Check mapping stats for a sample
    python troubleshoot_pipeline.py --sample-dir /path/to/sample --check-mapping

    # Run full diagnosis including assembly and BLAST
    python troubleshoot_pipeline.py --sample-dir /path/to/sample --full-diagnosis
    
    # Assembly diagnosis only
    python troubleshoot_pipeline.py --r1 reads_R1.fq.gz --r2 reads_R2.fq.gz --assembly-diagnosis
        """
    )
    
    parser.add_argument("--sample-dir", help="Sample directory to analyze")
    parser.add_argument("--sample-name", help="Specific sample name to focus on")
    parser.add_argument("--r1", help="R1 FASTQ file for assembly")
    parser.add_argument("--r2", help="R2 FASTQ file for assembly")
    parser.add_argument("--output", default="troubleshoot_results", help="Output directory")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    
    # Analysis options
    parser.add_argument("--check-mapping", action="store_true", help="Check mapping statistics")
    parser.add_argument("--check-variants", action="store_true", help="Check variant files")
    parser.add_argument("--assembly-diagnosis", action="store_true", help="Run assembly diagnosis")
    parser.add_argument("--full-diagnosis", action="store_true", help="Run complete diagnosis")
    
    args = parser.parse_args()
    
    if not any([args.check_mapping, args.check_variants, args.assembly_diagnosis, args.full_diagnosis]):
        parser.print_help()
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    success = True
    
    if args.full_diagnosis or args.check_mapping:
        if not args.sample_dir:
            logger.error("--sample-dir required for mapping checks")
            sys.exit(1)
        success &= check_mapping_stats(args.sample_dir, args.sample_name)
    
    if args.full_diagnosis or args.check_variants:
        if not args.sample_dir:
            logger.error("--sample-dir required for variant checks")
            sys.exit(1)
        success &= check_variant_files(args.sample_dir, args.sample_name)
    
    if args.full_diagnosis or args.assembly_diagnosis:
        if not args.r1 or not args.r2:
            if args.sample_dir:
                # Try to find cleaned reads
                cleaned_dir = os.path.join(args.sample_dir, "cleaned_seqs")
                if os.path.exists(cleaned_dir):
                    r1_files = list(Path(cleaned_dir).glob("*_R1*.fastq.gz"))
                    r2_files = list(Path(cleaned_dir).glob("*_R2*.fastq.gz"))
                    if r1_files and r2_files:
                        args.r1 = str(r1_files[0])
                        args.r2 = str(r2_files[0])
                        logger.info(f"Using cleaned reads: {args.r1}, {args.r2}")
            
            if not args.r1 or not args.r2:
                logger.error("--r1 and --r2 required for assembly diagnosis")
                sys.exit(1)
        
        contigs_file = run_assembly_diagnosis(args.r1, args.r2, args.output, args.threads)
        if contigs_file:
            blast_contigs_for_identification(contigs_file, args.output)
    
    if success:
        logger.info("\nTroubleshooting completed successfully!")
        logger.info(f"Results saved to: {args.output}")
    else:
        logger.error("Some troubleshooting steps failed")
        sys.exit(1)

if __name__ == "__main__":
    main()