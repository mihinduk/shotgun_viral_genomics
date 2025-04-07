#!/usr/bin/env python3
"""
Shotgun Viral Genomics Pipeline

A comprehensive pipeline for processing shotgun sequencing data from viral samples,
performing quality control, mapping to reference genomes, variant calling, and annotation.

Copyright (c) 2024
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Dict, Union, Any

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
        description="Shotgun Viral Genomics Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    
    # Input options
    parser.add_argument("--r1", required=True, help="Forward reads (R1) FASTQ files (accepts wildcards like *_R1.fastq.gz or *_R1_001.fastq.gz)")
    parser.add_argument("--r2", required=True, help="Reverse reads (R2) FASTQ files (accepts wildcards like *_R2.fastq.gz or *_R2_001.fastq.gz)")
    
    # Reference genome options
    ref_group = parser.add_argument_group("Reference Genome")
    ref_group.add_argument("--accession", help="Reference genome accession number (e.g., NC_045512.2)")
    ref_group.add_argument("--reference", help="Path to local reference genome (when not using accession)")
    ref_group.add_argument("--skip-download", action="store_true", help="Skip genome download (use local reference)")
    ref_group.add_argument("--force-download", action="store_true", help="Force download even if reference exists locally")
    
    # Pipeline control
    pipeline_group = parser.add_argument_group("Pipeline Control")
    pipeline_group.add_argument("--skip-qc", action="store_true", help="Skip quality control step")
    pipeline_group.add_argument("--skip-stats", action="store_true", help="Skip read statistics step")
    pipeline_group.add_argument("--skip-mapping", action="store_true", help="Skip read mapping step")
    pipeline_group.add_argument("--skip-variants", action="store_true", help="Skip variant calling step")
    pipeline_group.add_argument("--skip-annotation", action="store_true", help="Skip variant annotation step")
    
    # Output options
    output_group = parser.add_argument_group("Output")
    output_group.add_argument("--outdir", default=".", help="Output directory")
    output_group.add_argument("--min-depth", type=int, default=200, help="Minimum read depth for final variant reporting")
    output_group.add_argument("--keep-tmp", action="store_true", help="Keep temporary files")
    
    # Performance options
    perf_group = parser.add_argument_group("Performance")
    perf_group.add_argument("--threads", type=int, default=1, help="Number of CPU threads to use")
    
    # SnpEff options
    snpeff_group = parser.add_argument_group("SnpEff")
    snpeff_group.add_argument("--add-to-snpeff", action="store_true", help="Attempt to add genome to snpEff if not present")
    snpeff_group.add_argument("--snpeff-jar", default="snpEff.jar", help="Path to snpEff.jar")
    
    return parser.parse_args()

def download_reference_genome(accession: str, output_dir: str, force: bool = False) -> str:
    """
    Download a reference genome from NCBI using Entrez Direct utilities.
    
    Args:
        accession: The genome accession number (e.g., NC_045512.2)
        output_dir: Directory to save the downloaded genome
        force: Force download even if the file exists
        
    Returns:
        Path to the downloaded genome FASTA file
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{accession}.fasta")
    
    # Check if file already exists
    if os.path.exists(output_file) and not force:
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

def run_command(cmd: Union[str, List[str]], shell: bool = False, check: bool = True) -> subprocess.CompletedProcess:
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

def check_snpeff_database(accession: str, snpeff_jar: str) -> bool:
    """
    Check if a genome is in the snpEff database.
    
    Args:
        accession: Genome accession number
        snpeff_jar: Path to snpEff.jar
        
    Returns:
        True if the genome is in the database, False otherwise
    """
    logger.info(f"Checking if {accession} is in snpEff database")
    
    cmd = f"//usr/lib/jvm/java-11-openjdk-11.0.20.0.8-1.el7_9.x86_64/bin/java -jar {snpeff_jar} databases | grep {accession}"
    result = run_command(cmd, shell=True, check=False)
    
    return result.returncode == 0

def add_genome_to_snpeff(accession: str, fasta_path: str, snpeff_jar: str) -> bool:
    """
    Add a genome to the snpEff database.
    
    Args:
        accession: Genome accession number
        fasta_path: Path to the genome FASTA file
        snpeff_jar: Path to snpEff.jar
        
    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Adding genome {accession} to snpEff database")
    
    # Run snpEff build command
    cmd = f"//usr/lib/jvm/java-11-openjdk-11.0.20.0.8-1.el7_9.x86_64/bin/java -jar {snpeff_jar} build -genbank -v -noCheckProtein {accession}"
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

def calculate_read_stats(fastq_pattern: str, output_file: str, threads: int = 1) -> None:
    """
    Calculate read statistics using seqkit.
    
    Args:
        fastq_pattern: Pattern to match FASTQ files
        output_file: Output file for statistics
        threads: Number of threads to use
    """
    logger.info("Calculating read statistics")
    
    # Support both naming conventions (*_R1.fastq.gz and *_R1_001.fastq.gz)
    cmd = f"seqkit stats {fastq_pattern} -T -j {threads} > {output_file}"
    run_command(cmd, shell=True)
    
    logger.info(f"Read statistics saved to {output_file}")

def clean_reads(output_dir: str, threads: int = 1) -> Dict[str, Dict[str, str]]:
    """
    Clean reads using fastp.
    
    Args:
        output_dir: Output directory for cleaned reads
        threads: Number of threads to use
        
    Returns:
        Dictionary mapping sample names to input and output files
    """
    logger.info("Cleaning reads with fastp")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get R1 fastq files - support both naming conventions
    r1_files = [f for f in os.listdir('.') if re.match(r'.*_R1(?:_001)?\.fastq\.gz$', f)]
    
    if not r1_files:
        raise FileNotFoundError("No R1 FASTQ files found in the current directory")
    
    cleaned_files = {}
    
    for r1_file in r1_files:
        # Extract sample name
        match = re.match(r'(.+)_R1(?:_001)?\.fastq\.gz$', r1_file)
        if not match:
            logger.warning(f"Skipping file with unusual naming pattern: {r1_file}")
            continue
            
        sample_name = match.group(1)
        
        # Find matching R2 file - handle both naming conventions
        r2_file = r1_file.replace('_R1', '_R2')
        
        if not os.path.exists(r2_file):
            logger.warning(f"No matching R2 file found for {r1_file}")
            continue
        
        logger.info(f"Processing sample: {sample_name}")
        
        # Output file paths - standardized naming convention regardless of input
        r1_out = os.path.join(output_dir, f"{sample_name}_R1.qc.fastq.gz")
        r2_out = os.path.join(output_dir, f"{sample_name}_R2.qc.fastq.gz")
        html_report = os.path.join(output_dir, f"{sample_name}_fastp_report.html")
        json_report = os.path.join(output_dir, f"{sample_name}_fastp_report.json")
        qc_report = os.path.join(output_dir, f"{sample_name}_fastp_QC.txt")
        
        # Run fastp
        cmd = (
            f"fastp -w {threads} -q 30 "
            f"-i {r1_file} -I {r2_file} "
            f"-o {r1_out} -O {r2_out} "
            f"-h {html_report} -j {json_report}"
        )
        
        # Capture fastp output for QC report
        result = run_command(cmd, shell=True)
        
        # Extract and save QC information
        qc_info = []
        qc_info.append(f"{sample_name}")
        
        output_lines = result.stderr.split('\n')
        capture = False
        for line in output_lines:
            if "Filtering result:" in line:
                capture = True
                qc_info.append(line)
            elif capture and line.strip():
                qc_info.append(line)
            elif capture and "Insert size peak" in line:
                qc_info.append(line)
                break
        
        with open(qc_report, 'w') as f:
            f.write('\n'.join(qc_info))
        
        # Add files to result dictionary
        cleaned_files[sample_name] = {
            'r1_in': r1_file,
            'r2_in': r2_file,
            'r1_out': r1_out,
            'r2_out': r2_out,
            'qc_report': qc_report
        }
    
    logger.info(f"Cleaned {len(cleaned_files)} samples")
    return cleaned_files

def map_and_call_variants(
    reference: str, 
    output_dir: str, 
    threads: int = 1
) -> Dict[str, Dict[str, str]]:
    """
    Map reads to reference and call variants.
    
    Args:
        reference: Path to reference genome
        output_dir: Base output directory
        threads: Number of threads to use
        
    Returns:
        Dictionary mapping sample names to output files
    """
    logger.info(f"Mapping reads to reference and calling variants: {reference}")
    
    # Create output directories
    mapping_dir = os.path.join(output_dir, "mapping")
    variants_dir = os.path.join(output_dir, "variants")
    os.makedirs(mapping_dir, exist_ok=True)
    os.makedirs(variants_dir, exist_ok=True)
    
    # Index reference genome
    logger.info(f"Indexing reference genome: {reference}")
    run_command(f"bwa index {reference}", shell=True)
    
    # Get cleaned R1 files (using standardized naming)
    r1_files = [f for f in os.listdir(output_dir) if f.endswith('_R1.qc.fastq.gz')]
    
    if not r1_files:
        raise FileNotFoundError(f"No cleaned R1 files found in {output_dir}")
    
    result_files = {}
    
    for r1_file in r1_files:
        # Extract sample name
        sample_name = r1_file.replace('_R1.qc.fastq.gz', '')
        r2_file = r1_file.replace('_R1', '_R2')
        
        r1_path = os.path.join(output_dir, r1_file)
        r2_path = os.path.join(output_dir, r2_file)
        
        if not os.path.exists(r2_path):
            logger.warning(f"No matching R2 file found for {r1_path}")
            continue
        
        logger.info(f"Processing sample: {sample_name}")
        
        # Prepare output file paths
        sam_file = os.path.join(mapping_dir, f"{sample_name}.sam")
        fixmate_file = os.path.join(mapping_dir, f"{sample_name}.fixmate.bam")
        bam_file = os.path.join(mapping_dir, f"{sample_name}.bam")
        dedupe_file = os.path.join(mapping_dir, f"{sample_name}.dedupe.bam")
        realign_file = os.path.join(mapping_dir, f"{sample_name}.lofreq.realign.bam")
        indel_file = os.path.join(mapping_dir, f"{sample_name}.lofreq.indel.bam")
        final_bam = os.path.join(mapping_dir, f"{sample_name}.lofreq.final.bam")
        final_bai = os.path.join(mapping_dir, f"{sample_name}.lofreq.final.bam.bai")
        vars_file = os.path.join(variants_dir, f"{sample_name}_vars.vcf")
        
        # Alignment steps
        
        # 1. BWA MEM alignment
        logger.info(f"Aligning reads for {sample_name}")
        run_command(
            f"bwa mem -t {threads} {reference} {r1_path} {r2_path} > {sam_file}",
            shell=True
        )
        
        # 2. Fix mate information
        logger.info(f"Fixing mate information for {sample_name}")
        run_command(
            f"samtools fixmate -O bam,level=1 -m --threads {threads} {sam_file} {fixmate_file}",
            shell=True
        )
        
        # 3. Sort BAM
        logger.info(f"Sorting BAM file for {sample_name}")
        run_command(
            f"samtools sort --threads {threads} -O bam {fixmate_file} > {bam_file}",
            shell=True
        )
        
        # 4. Mark duplicates
        logger.info(f"Marking duplicates for {sample_name}")
        run_command(
            f"samtools markdup --threads {threads} -S {bam_file} {dedupe_file}",
            shell=True
        )
        
        # 5. LoFreq Viterbi realignment
        logger.info(f"LoFreq Viterbi realignment for {sample_name}")
        run_command(
            f"lofreq viterbi -f {reference} {dedupe_file} | samtools sort - --threads {threads} > {realign_file}",
            shell=True
        )
        
        # 6. LoFreq indel quality calibration
        logger.info(f"LoFreq indel quality calibration for {sample_name}")
        run_command(
            f"lofreq indelqual --dindel -f {reference} {realign_file} | samtools sort - --threads {threads} > {indel_file}",
            shell=True
        )
        
        # 7. LoFreq alignment quality calibration
        logger.info(f"LoFreq alignment quality calibration for {sample_name}")
        run_command(
            f"lofreq alnqual -b {indel_file} {reference} > {final_bam}",
            shell=True
        )
        
        # 8. Index final BAM
        logger.info(f"Indexing final BAM for {sample_name}")
        run_command(f"samtools index {final_bam}", shell=True)
        
        # 9. Call variants with LoFreq
        logger.info(f"Calling variants for {sample_name}")
        run_command(
            f"lofreq call-parallel --pp-threads {threads} --force-overwrite --no-default-filter --call-indels -f {reference} -o {vars_file} {final_bam}",
            shell=True
        )
        
        # Store output files
        result_files[sample_name] = {
            'sam': sam_file,
            'fixmate': fixmate_file,
            'bam': bam_file,
            'dedupe': dedupe_file,
            'realign': realign_file,
            'indel': indel_file,
            'final_bam': final_bam,
            'final_bai': final_bai,
            'vars': vars_file
        }
    
    logger.info(f"Mapping and variant calling completed for {len(result_files)} samples")
    return result_files

def filter_variants(variants_dir: str) -> Dict[str, str]:
    """
    Filter variant calls from LoFreq.
    
    Args:
        variants_dir: Directory containing variant files
        
    Returns:
        Dictionary mapping sample names to filtered variant files
    """
    logger.info("Filtering variant calls")
    
    # Get VCF files
    vcf_files = [f for f in os.listdir(variants_dir) if f.endswith('_vars.vcf')]
    
    if not vcf_files:
        raise FileNotFoundError(f"No variant VCF files found in {variants_dir}")
    
    filtered_files = {}
    
    for vcf_file in vcf_files:
        # Extract sample name
        sample_name = vcf_file.replace('_vars.vcf', '')
        vcf_path = os.path.join(variants_dir, vcf_file)
        
        filtered_path = os.path.join(variants_dir, f"{sample_name}_vars.filt.vcf")
        
        logger.info(f"Filtering variants for {sample_name}")
        
        # Run LoFreq filter
        run_command(
            f"lofreq filter -i {vcf_path} -o {filtered_path} -v 75",
            shell=True
        )
        
        filtered_files[sample_name] = filtered_path
    
    logger.info(f"Variant filtering completed for {len(filtered_files)} samples")
    return filtered_files

def annotate_variants(variants_dir: str, accession: str, snpeff_jar: str) -> Dict[str, Dict[str, str]]:
    """
    Annotate filtered variants using snpEff.
    
    Args:
        variants_dir: Directory containing filtered variant files
        accession: Reference genome accession
        snpeff_jar: Path to snpEff.jar
        
    Returns:
        Dictionary mapping sample names to annotation files
    """
    logger.info(f"Annotating variants with snpEff using database: {accession}")
    
    # Get filtered VCF files
    filt_files = [f for f in os.listdir(variants_dir) if f.endswith('_vars.filt.vcf')]
    
    if not filt_files:
        raise FileNotFoundError(f"No filtered variant VCF files found in {variants_dir}")
    
    annotation_files = {}
    
    for filt_file in filt_files:
        # Extract sample name
        sample_name = filt_file.replace('_vars.filt.vcf', '')
        filt_path = os.path.join(variants_dir, filt_file)
        
        # Output file paths
        ann_vcf = os.path.join(variants_dir, f"{sample_name}.snpEFF.ann.vcf")
        ann_tsv = os.path.join(variants_dir, f"{sample_name}.snpEFF.ann.tsv")
        summary_html = os.path.join(variants_dir, f"{sample_name}_summary.html")
        summary_genes = os.path.join(variants_dir, f"{sample_name}_summary.genes.txt")
        
        logger.info(f"Annotating variants for {sample_name}")
        
        # Run snpEff
        cmd = f"///usr/lib/jvm/java-11-openjdk-11.0.20.0.8-1.el7_9.x86_64/bin/java -jar -Xmx4g {snpeff_jar} {accession} {filt_path} -s {summary_html} > {ann_vcf}"
        run_command(cmd, shell=True)
        
        # Process the VCF file into TSV format
        logger.info(f"Converting VCF to TSV for {sample_name}")
        
        # Create temporary files
        tmp_base = os.path.join(variants_dir, f"{sample_name}.ann.base.vcf")
        tmp_info = os.path.join(variants_dir, f"{sample_name}.snpEFF.ann.tmp")
        
        # Extract annotations
        run_command(
            f"grep -v '^##' {ann_vcf} | "
            f"tail -n+2 | "
            f"cut -f8 | "
            f"sed 's/|/\\t/g' | "
            f"cut -f1-16 | "
            f"sed '1i INFO\\tEFFECT\\tPUTATIVE_IMPACT\\tGENE_NAME\\tGENE_ID\\tFEATURE_TYPE\\tFEATURE_ID\\tTRANSCRIPT_TYPE\\tEXON_INTRON_RANK\\tHGVSc\\tHGVSp\\tcDNA_POSITION_AND_LENGTH\\tCDS_POSITION_AND_LENGTH\\tPROTEIN_POSITION_AND_LENGTH\\tDISTANCE_TO_FEATURE\\tERROR' > {tmp_info}",
            shell=True
        )
        
        # Extract base VCF information
        run_command(
            f"grep -v '^##' {ann_vcf} | cut -f1-7 > {tmp_base}",
            shell=True
        )
        
        # Combine into final TSV
        run_command(
            f"paste {tmp_base} {tmp_info} > {ann_tsv}",
            shell=True
        )
        
        # Cleanup temporary files
        os.remove(tmp_base)
        os.remove(tmp_info)
        
        # Check for ERROR_CHROMOSOME_NOT_FOUND
        check_cmd = f"grep -q 'ERROR_CHROMOSOME_NOT_FOUND' {ann_tsv}"
        result = run_command(check_cmd, shell=True, check=False)
        
        if result.returncode == 0:
            logger.error(f"ERROR: Genome {accession} is NOT properly formatted for use with snpEff annotation.")
            raise RuntimeError(f"Genome {accession} is NOT properly formatted for use with snpEff annotation.")
        
        annotation_files[sample_name] = {
            'ann_vcf': ann_vcf,
            'ann_tsv': ann_tsv,
            'summary_html': summary_html,
            'summary_genes': summary_genes
        }
    
    logger.info(f"Variant annotation completed for {len(annotation_files)} samples")
    return annotation_files

def parse_annotations(variants_dir: str, min_depth: int) -> Dict[str, str]:
    """
    Parse annotated variants using the Perl script.
    
    Args:
        variants_dir: Directory containing annotation files
        min_depth: Minimum read depth for reporting
        
    Returns:
        Dictionary mapping sample names to parsed annotation files
    """
    logger.info(f"Parsing annotated variants with minimum depth {min_depth}")
    
    # Get annotation TSV files
    ann_files = [f for f in os.listdir(variants_dir) if f.endswith('.snpEFF.ann.tsv')]
    
    if not ann_files:
        raise FileNotFoundError(f"No annotation TSV files found in {variants_dir}")
    
    parsed_files = {}
    perl_script = "./parse_snpEff_annotated_vcf_for_collaborators.pl"
    
    # Check if Perl script exists
    if not os.path.exists(perl_script):
        raise FileNotFoundError(f"Perl script not found: {perl_script}")
    
    for ann_file in ann_files:
        # Extract sample name
        sample_name = ann_file.replace('.snpEFF.ann.tsv', '')
        ann_path = os.path.join(variants_dir, ann_file)
        
        # Parse with Perl script
        logger.info(f"Parsing annotations for {sample_name} with minimum depth {min_depth}")
        
        run_command(
            f"perl {perl_script} {ann_path} {min_depth}",
            shell=True
        )
        
        # The output file should be named by the script
        output_file = f"{ann_path.rstrip('.tsv')}_{min_depth}.tsv"
        parsed_files[sample_name] = output_file
    
    logger.info(f"Annotation parsing completed for {len(parsed_files)} samples")
    return parsed_files

def main():
    """Main function to run the pipeline."""
    args = parse_args()
    
    # Create output directory structure
    output_dir = args.outdir
    cleaned_dir = os.path.join(output_dir, "cleaned_seqs")
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(cleaned_dir, exist_ok=True)
    
    # Create temp directory
    temp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(temp_dir, exist_ok=True)
    
    try:
        # Step 1: Prepare reference genome
        reference_path = None
        if args.reference:
            reference_path = args.reference
            logger.info(f"Using provided reference genome: {reference_path}")
        elif args.accession:
            # Download reference genome if needed
            if not args.skip_download:
                reference_path = download_reference_genome(
                    args.accession, 
                    cleaned_dir,
                    args.force_download
                )
            else:
                # Use existing reference
                reference_path = os.path.join(cleaned_dir, f"{args.accession}.fasta")
                if not os.path.exists(reference_path):
                    raise FileNotFoundError(f"Reference genome not found: {reference_path}")
        else:
            raise ValueError("Either --accession or --reference must be provided")
        
        # Check if reference genome is in snpEff database
        accession = args.accession
        if not accession and reference_path:
            # Try to extract accession from reference file name
            accession = os.path.basename(reference_path).split('.')[0]
        
        if accession:
            in_database = check_snpeff_database(accession, args.snpeff_jar)
            if not in_database:
                logger.warning(f"Genome {accession} not found in snpEff database")
                if args.add_to_snpeff:
                    logger.info(f"Attempting to add {accession} to snpEff database")
                    success = add_genome_to_snpeff(accession, reference_path, args.snpeff_jar)
                    if not success:
                        logger.error(f"Failed to add {accession} to snpEff database. Annotation will likely fail.")
                else:
                    logger.warning("Use --add-to-snpeff to attempt adding the genome to snpEff")
        
        # Step 2: Calculate read statistics
        if not args.skip_stats:
            stats_file = os.path.join(output_dir, "input_stats.txt")
            calculate_read_stats(args.r1, stats_file, args.threads)
        
        # Step 3: Clean reads
        if not args.skip_qc:
            cleaned_files = clean_reads(cleaned_dir, args.threads)
        
        # Step 4: Map reads and call variants
        if not args.skip_mapping:
            variant_files = map_and_call_variants(
                reference_path, 
                cleaned_dir, 
                args.threads
            )
        
        # Step 5: Filter variants
        if not args.skip_variants:
            variants_dir = os.path.join(cleaned_dir, "variants")
            filtered_files = filter_variants(variants_dir)
        
        # Step 6: Annotate variants
        if not args.skip_annotation and accession:
            variants_dir = os.path.join(cleaned_dir, "variants")
            annotation_files = annotate_variants(variants_dir, accession, args.snpeff_jar)
            
            # Step 7: Parse annotations
            parsed_files = parse_annotations(variants_dir, args.min_depth)
        
        logger.info("Pipeline completed successfully")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}", exc_info=True)
        sys.exit(1)
    finally:
        # Cleanup
        if not args.keep_tmp and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

if __name__ == "__main__":
    main()
