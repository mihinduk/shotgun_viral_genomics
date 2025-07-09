#!/usr/bin/env python3
"""
Viral Assembly Script - Performs de novo assembly of viral genomes
to help identify appropriate reference genomes
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

def run_spades_assembly(r1, r2, output_dir, threads=4):
    """Run SPAdes assembler for viral genomes."""
    logger.info("Running SPAdes assembly...")
    
    # Check if SPAdes is available
    check_spades = run_command("which spades.py", shell=True, check=False)
    if check_spades.returncode != 0:
        logger.error("SPAdes not found. Please install it or use --assembler megahit")
        return None
    
    spades_dir = os.path.join(output_dir, "spades_assembly")
    os.makedirs(spades_dir, exist_ok=True)
    
    # Use rnaviral mode for RNA viruses - best for your POWV samples
    cmd = f"spades.py --rnaviral -1 {r1} -2 {r2} -o {spades_dir} -t {threads} --careful"
    
    try:
        result = run_command(cmd, shell=True)
        
        # Find the best contigs
        contigs_file = os.path.join(spades_dir, "contigs.fasta")
        scaffolds_file = os.path.join(spades_dir, "scaffolds.fasta")
        
        # Prefer scaffolds if available
        if os.path.exists(scaffolds_file):
            logger.info(f"Assembly completed. Scaffolds saved to: {scaffolds_file}")
            return scaffolds_file
        elif os.path.exists(contigs_file):
            logger.info(f"Assembly completed. Contigs saved to: {contigs_file}")
            return contigs_file
        else:
            logger.error("SPAdes assembly completed but no contigs file found")
            return None
            
    except subprocess.CalledProcessError as e:
        logger.error(f"SPAdes assembly failed: {e}")
        return None

def run_megahit_assembly(r1, r2, output_dir, threads=4):
    """Run MEGAHIT assembler as alternative."""
    logger.info("Running MEGAHIT assembly...")
    
    # Check if MEGAHIT is available
    check_megahit = run_command("which megahit", shell=True, check=False)
    if check_megahit.returncode != 0:
        logger.error("MEGAHIT not found. Please install it")
        return None
    
    megahit_dir = os.path.join(output_dir, "megahit_assembly")
    
    # MEGAHIT is good for metagenomes and low-coverage viruses
    cmd = f"megahit -1 {r1} -2 {r2} -o {megahit_dir} -t {threads} --presets meta-sensitive --min-contig-len 500"
    
    try:
        # Remove existing directory if it exists
        if os.path.exists(megahit_dir):
            import shutil
            shutil.rmtree(megahit_dir)
            
        result = run_command(cmd, shell=True)
        
        # Find the contigs
        contigs_file = os.path.join(megahit_dir, "final.contigs.fa")
        if os.path.exists(contigs_file):
            logger.info(f"Assembly completed. Contigs saved to: {contigs_file}")
            return contigs_file
        else:
            logger.error("MEGAHIT assembly completed but no contigs file found")
            return None
            
    except subprocess.CalledProcessError as e:
        logger.error(f"MEGAHIT assembly failed: {e}")
        return None

def blast_contigs(contigs_file, output_dir, num_results=10):
    """BLAST contigs against NCBI to find best reference."""
    logger.info("Running BLAST to identify best reference genome...")
    
    blast_output = os.path.join(output_dir, "blast_results.txt")
    
    # Use blastn against nt database
    cmd = (f"blastn -query {contigs_file} -db nt -remote "
           f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' "
           f"-max_target_seqs {num_results} -max_hsps 1 > {blast_output}")
    
    try:
        result = run_command(cmd, shell=True)
        logger.info(f"BLAST results saved to: {blast_output}")
        
        # Parse and display top hits
        logger.info("\nTop BLAST hits:")
        with open(blast_output, 'r') as f:
            lines = f.readlines()[:num_results]
            for line in lines:
                parts = line.strip().split('\t')
                if len(parts) >= 13:
                    logger.info(f"  {parts[1]} - {parts[2]}% identity - {parts[12]}")
                    
        return blast_output
        
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST failed: {e}")
        return None

def check_contig_stats(contigs_file):
    """Get basic statistics about assembled contigs."""
    logger.info("Analyzing contig statistics...")
    
    cmd = f"seqkit stats {contigs_file}"
    try:
        result = run_command(cmd, shell=True)
        logger.info(f"\nContig statistics:\n{result.stdout}")
    except:
        logger.warning("Could not generate contig statistics")

def main():
    parser = argparse.ArgumentParser(description="Viral genome assembly pipeline")
    parser.add_argument("--r1", required=True, help="Forward reads (R1)")
    parser.add_argument("--r2", required=True, help="Reverse reads (R2)")
    parser.add_argument("--output", default="assembly_results", help="Output directory")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--assembler", choices=["megahit"], 
                       default="megahit", help="Which assembler to use (SPAdes removed due to stability issues)")
    parser.add_argument("--skip-blast", action="store_true", help="Skip BLAST search")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Check if reads are cleaned or need cleaning
    if "cleaned" not in args.r1:
        logger.info("Note: Using raw reads. For best results, use cleaned reads from the main pipeline.")
    
    contigs_file = None
    
    # Run assembly
    contigs_file = run_megahit_assembly(args.r1, args.r2, args.output, args.threads)
    
    if contigs_file and os.path.exists(contigs_file):
        # Get contig statistics
        check_contig_stats(contigs_file)
        
        # Run BLAST unless skipped
        if not args.skip_blast:
            blast_contigs(contigs_file, args.output)
        
        logger.info(f"\nAssembly complete! Contigs saved to: {contigs_file}")
        logger.info("You can use these contigs as a reference genome or find the best match from BLAST results.")
    else:
        logger.error("Assembly failed - no contigs produced")
        sys.exit(1)

if __name__ == "__main__":
    main()