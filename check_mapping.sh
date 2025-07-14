#!/bin/bash
# Script to check mapping statistics and diagnose why no variants are found

if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample_directory>"
    echo "Example: $0 /scratch/sahlab/kathie/NovaSeq_N1027_POWV_MA51240"
    exit 1
fi

SAMPLE_DIR=$1

echo "=== Checking Mapping Statistics for $SAMPLE_DIR ==="
echo

# Find BAM files
echo "Looking for BAM files..."
find "$SAMPLE_DIR" -name "*.final.bam" -o -name "*.bam" | head -10

echo
echo "=== Checking mapping rates ==="

# Check each BAM file
for BAM in $(find "$SAMPLE_DIR" -name "*.final.bam" -o -name "*.lofreq.final.bam" | head -3); do
    if [ -f "$BAM" ]; then
        echo
        echo "File: $(basename $BAM)"
        echo "Size: $(ls -lh $BAM | awk '{print $5}')"
        
        # Get basic stats
        echo "Running samtools flagstat..."
        samtools flagstat "$BAM" 2>/dev/null | grep -E "mapped|properly paired" || echo "  Error running samtools"
        
        # Check coverage
        echo "Checking coverage depth..."
        REF=$(find "$SAMPLE_DIR" -name "*.fasta" | grep -E "HM44056[012]|NC_" | head -1)
        if [ -f "$REF" ]; then
            echo "  Using reference: $(basename $REF)"
            samtools depth "$BAM" 2>/dev/null | awk '{sum+=$3; count++} END {if(count>0) print "  Average depth: " sum/count; else print "  No coverage data"}' || echo "  Error calculating depth"
        fi
    fi
done

echo
echo "=== Checking VCF files ==="

# Check VCF files
for VCF in $(find "$SAMPLE_DIR" -name "*_vars.vcf" -o -name "*_vars.filt.vcf" | head -5); do
    if [ -f "$VCF" ]; then
        echo
        echo "File: $(basename $VCF)"
        VARIANTS=$(grep -v "^#" "$VCF" 2>/dev/null | wc -l)
        echo "  Number of variants: $VARIANTS"
        
        if [ $VARIANTS -gt 0 ]; then
            echo "  First few variants:"
            grep -v "^#" "$VCF" | head -3 | cut -f1-5
        fi
    fi
done

echo
echo "=== Checking cleaned reads ==="

# Check cleaned reads
for FASTQ in $(find "$SAMPLE_DIR/cleaned_seqs" -name "*_1.fq.gz" -o -name "*_R1*.gz" | head -3); do
    if [ -f "$FASTQ" ]; then
        echo
        echo "Cleaned file: $(basename $FASTQ)"
        echo "  Size: $(ls -lh $FASTQ | awk '{print $5}')"
        
        # Quick read count
        echo "  Counting reads (first 100k)..."
        zcat "$FASTQ" 2>/dev/null | head -400000 | grep "^@" | wc -l | awk '{print "  Reads in sample: " $1 " (of first 100k lines)"}'
    fi
done

echo
echo "=== Summary ==="
echo "Check the mapping rates above. If they're very low (<1%), the reference genome is likely wrong."
echo "If mapping is good but no variants, check:"
echo "  1. Coverage depth - too low coverage means no confident variant calls"
echo "  2. Sample purity - clonal populations may have no variants"
echo "  3. Reference match - perfect match to reference = no variants"