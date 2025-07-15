#!/bin/bash
# Viral Contamination Diagnostic Module
# Runs mapping check, assembly, and viral BLAST for sample quality assessment
# Usage: ./viral_diagnostic.sh <R1> <R2> <accession> <sample_name> [threads]

# Check arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> <sample_name> [threads]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 SMS_14 4"
    echo ""
    echo "This diagnostic script will:"
    echo "  1. Check mapping statistics to reference"
    echo "  2. Perform de novo assembly with MEGAHIT"
    echo "  3. BLAST contigs against viral nt database"
    echo "  4. Generate contamination assessment report"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
SAMPLE_NAME=$4
THREADS=${5:-4}
EXTREME_MEMORY_FLAG=""

# Check for extreme memory flag
if [[ "$6" == "--extremely-large-files" ]] || [[ "$*" == *"--extremely-large-files"* ]]; then
    EXTREME_MEMORY_FLAG="--extremely-large-files"
    echo "EXTREME MEMORY MODE: Using high memory settings for large files"
fi

# Set up environment and paths
echo "========================================="
echo "VIRAL CONTAMINATION DIAGNOSTIC MODULE"
echo "========================================="
echo "Sample: $SAMPLE_NAME"
echo "R1: $R1"
echo "R2: $R2"
echo "Reference: $ACCESSION"
echo "Threads: $THREADS"
echo "Started at: $(date)"
echo "========================================="

# Create output directory (remove if exists to avoid conflicts)
DIAGNOSTIC_DIR="./diagnostic_${SAMPLE_NAME}"
if [ -d "$DIAGNOSTIC_DIR" ]; then
    echo "Removing existing diagnostic directory..."
    rm -rf "$DIAGNOSTIC_DIR"
fi
mkdir -p "$DIAGNOSTIC_DIR"
cd "$DIAGNOSTIC_DIR"

# Set up conda environment commands  
echo "Setting up environment commands..."

# Aggressively clear mamba locks and processes
echo "Cleaning up any existing mamba processes and locks..."
pkill -f "mamba" 2>/dev/null || true
rm -rf /home/mihindu/.cache/mamba/proc* 2>/dev/null || true
rm -rf /home/mihindu/.cache/mamba/locks* 2>/dev/null || true
sleep 2

# Try using conda instead of mamba to avoid lock issues entirely
source /ref/sahlab/software/anaconda3/bin/activate base
eval "$(conda shell.bash hook)"

# Check if viral_genomics tools are available via anaconda3 conda
if conda run -n viral_genomics which bwa >/dev/null 2>&1; then
    echo "Using anaconda3 conda for viral_genomics environment"
    GENOMICS_CMD="conda run -n viral_genomics"
else
    echo "Using miniforge3 mamba for viral_genomics environment (anaconda3 doesn't have tools)"
    export MAMBA_NO_BANNER=1
    GENOMICS_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"
fi

# Initialize base filenames for consistent use throughout script
R1_BASE=$(basename "$R1" .fastq.gz)
R2_BASE=$(basename "$R2" .fastq.gz)
echo "Base filenames: $R1_BASE and $R2_BASE"

# Step 1: Quick mapping check
echo ""
echo "Step 1: Mapping Statistics Check"
echo "================================"

# Download reference if needed
if [ ! -f "${ACCESSION}.fasta" ]; then
    echo "Downloading reference genome: $ACCESSION"
    $GENOMICS_CMD efetch -db nucleotide -id "$ACCESSION" -format fasta > "${ACCESSION}.fasta"
fi

# Index reference
echo "Indexing reference genome..."
$GENOMICS_CMD bwa index "${ACCESSION}.fasta"

# Quick alignment for mapping stats using cleaned reads
echo "Performing quick alignment for mapping statistics using cleaned reads..."
# Look for cleaned reads in the parent cleaned_seqs directory
CLEANED_R1="../cleaned_seqs/${R1_BASE}.qc.fastq.gz"
CLEANED_R2="../cleaned_seqs/${R2_BASE}.qc.fastq.gz"

if [ -f "$CLEANED_R1" ] && [ -f "$CLEANED_R2" ]; then
    echo "  Using cleaned reads: $CLEANED_R1 and $CLEANED_R2"
    $GENOMICS_CMD bwa mem -t "$THREADS" "${ACCESSION}.fasta" "$CLEANED_R1" "$CLEANED_R2" | \
        $GENOMICS_CMD samtools view -@ "$THREADS" -bS - | \
        $GENOMICS_CMD samtools sort -@ "$THREADS" -o "${SAMPLE_NAME}_quick.bam" -
else
    echo "  Warning: Cleaned reads not found, using original reads"
    echo "  Looking for: $CLEANED_R1 and $CLEANED_R2"
    $GENOMICS_CMD bwa mem -t "$THREADS" "${ACCESSION}.fasta" "../$R1" "../$R2" | \
        $GENOMICS_CMD samtools view -@ "$THREADS" -bS - | \
        $GENOMICS_CMD samtools sort -@ "$THREADS" -o "${SAMPLE_NAME}_quick.bam" -
fi

$GENOMICS_CMD samtools index "${SAMPLE_NAME}_quick.bam"

# Generate mapping statistics
echo "Generating mapping statistics..."
$GENOMICS_CMD samtools flagstat "${SAMPLE_NAME}_quick.bam" > "${SAMPLE_NAME}_mapping_stats.txt"
$GENOMICS_CMD samtools idxstats "${SAMPLE_NAME}_quick.bam" > "${SAMPLE_NAME}_idxstats.txt"

# Calculate mapping percentage
TOTAL_READS=$($GENOMICS_CMD samtools view -c "${SAMPLE_NAME}_quick.bam")
MAPPED_READS=$($GENOMICS_CMD samtools view -c -F 4 "${SAMPLE_NAME}_quick.bam")

# Extract duplication rate from fastp report
# First try to get it from the original main pipeline fastp report
# Remove the _R1 suffix to get the sample base name
SAMPLE_BASE=$(echo "$R1_BASE" | sed 's/_R1$//')
ORIGINAL_FASTP_JSON="../cleaned_seqs/${SAMPLE_BASE}_fastp_report.json"
echo "  Looking for fastp report: $ORIGINAL_FASTP_JSON"
if [ -f "$ORIGINAL_FASTP_JSON" ]; then
    # Extract duplication rate from main pipeline fastp JSON
    echo "  Found fastp report, extracting duplication rate..."
    # Extract duplication rate from the "duplication": {"rate": 0.67605} structure
    RAW_DUP_RATE=$(grep -A2 '"duplication"' "$ORIGINAL_FASTP_JSON" | grep -o '"rate": *[0-9.]*' | grep -o '[0-9.]*')
    echo "  Raw duplication rate extracted: '$RAW_DUP_RATE'"
    
    if [ -n "$RAW_DUP_RATE" ] && [ "$RAW_DUP_RATE" != "0" ]; then
        DUPLICATION_RATE=$(awk "BEGIN {printf \"%.2f\", $RAW_DUP_RATE * 100}")
        echo "  Using duplication rate from main pipeline: ${DUPLICATION_RATE}%"
    else
        DUPLICATION_RATE="0.00"
        echo "  Could not extract duplication rate, using ${DUPLICATION_RATE}%"
    fi
elif [ -f "${SAMPLE_NAME}_fastp_report.json" ] && [ -s "${SAMPLE_NAME}_fastp_report.json" ]; then
    # Fallback: try to extract from local fastp JSON (if it's not empty)
    DUPLICATION_RATE=$(grep -o '"dup_rate":[0-9.]*' "${SAMPLE_NAME}_fastp_report.json" 2>/dev/null | cut -d: -f2 || echo "0")
    DUPLICATION_RATE=$(awk "BEGIN {printf \"%.2f\", $DUPLICATION_RATE * 100}")
    echo "  Using duplication rate from local fastp: ${DUPLICATION_RATE}%"
else
    # Default fallback
    DUPLICATION_RATE="0.00"
    echo "  Warning: Could not extract duplication rate, using ${DUPLICATION_RATE}%"
fi

# Calculate duplication metrics based on fastp duplication rate
if [ "$TOTAL_READS" -gt 0 ]; then
    # Use awk for calculations
    MAPPING_PERCENT=$(awk "BEGIN {printf \"%.2f\", $MAPPED_READS * 100 / $TOTAL_READS}")
    
    # Calculate unique reads based on fastp duplication rate
    DUPLICATE_READS=$(awk "BEGIN {printf \"%.0f\", $TOTAL_READS * $DUPLICATION_RATE / 100}")
    UNIQUE_READS=$((TOTAL_READS - DUPLICATE_READS))
    
    if [ "$UNIQUE_READS" -gt 0 ]; then
        # Calculate unique mapped reads correctly
        # Since BWA doesn't mark duplicates, we estimate based on the duplication rate
        MAPPED_UNIQUE_READS=$(awk "BEGIN {printf \"%.0f\", $MAPPED_READS * (100 - $DUPLICATION_RATE) / 100}")
        # Calculate the actual percentage of unique reads that mapped
        DEDUPLICATED_MAPPING_PERCENT=$(awk "BEGIN {printf \"%.2f\", $MAPPED_UNIQUE_READS * 100 / $UNIQUE_READS}")
    else
        MAPPED_UNIQUE_READS=0
        DEDUPLICATED_MAPPING_PERCENT=0
    fi
else
    MAPPING_PERCENT=0
    DUPLICATE_READS=0
    UNIQUE_READS=0
    MAPPED_UNIQUE_READS=0
    DEDUPLICATED_MAPPING_PERCENT=0
fi

echo "Mapping Results:"
echo "  Total reads: $TOTAL_READS"
echo "  Duplicate reads: $DUPLICATE_READS (${DUPLICATION_RATE}%)"
echo "  Unique reads: $UNIQUE_READS"
echo "  Mapped reads (all): $MAPPED_READS (${MAPPING_PERCENT}%)"
echo "  Mapped unique reads: $MAPPED_UNIQUE_READS (${DEDUPLICATED_MAPPING_PERCENT}% of unique)"

# Step 2: De novo assembly
echo ""
echo "Step 2: De Novo Assembly"
echo "========================"

# Check for existing cleaned reads first (using pipeline naming convention)
# Use the base filenames initialized at script start
EXISTING_R1="../cleaned_seqs/${R1_BASE}.qc.fastq.gz"
EXISTING_R2="../cleaned_seqs/${R2_BASE}.qc.fastq.gz"

if [ -f "$EXISTING_R1" ] && [ -f "$EXISTING_R2" ]; then
    echo "Found existing cleaned reads from main pipeline, reusing them..."
    echo "  Using: $EXISTING_R1"
    echo "  Using: $EXISTING_R2"
    
    # Create symlinks using the same base filenames for consistency
    ln -sf "$EXISTING_R1" "${R1_BASE}.qc.fastq.gz"
    ln -sf "$EXISTING_R2" "${R2_BASE}.qc.fastq.gz"
    
    # Create a placeholder fastp report
    echo "Reused existing cleaned reads from main pipeline" > "${SAMPLE_NAME}_fastp_reused.txt"
    touch "${SAMPLE_NAME}_fastp_report.html"
    touch "${SAMPLE_NAME}_fastp_report.json"
    
else
    echo "No existing cleaned reads found, cleaning reads with fastp..."
    echo "  Looking for: $EXISTING_R1"
    echo "  Looking for: $EXISTING_R2"
    
    # Clean reads with identical settings to main pipeline
    $GENOMICS_CMD fastp -i "../$R1" -I "../$R2" \
          -o "${SAMPLE_NAME}_R1.qc.fastq.gz" \
          -O "${SAMPLE_NAME}_R2.qc.fastq.gz" \
          -h "${SAMPLE_NAME}_fastp_report.html" \
          -j "${SAMPLE_NAME}_fastp_report.json" \
          --thread "$THREADS" \
          --qualified_quality_phred 15 \
          --unqualified_percent_limit 40 \
          --length_required 36
fi

# Run MEGAHIT assembly
echo "Running MEGAHIT assembly..."

# Remove existing assembly directory if it exists
if [ -d "assembly_${SAMPLE_NAME}" ]; then
    echo "Removing existing assembly directory..."
    rm -rf "assembly_${SAMPLE_NAME}"
fi

eval "$(conda shell.bash hook)"
conda activate viral_assembly
# Build MEGAHIT command with conditional extreme memory settings
MEGAHIT_CMD="megahit -1 \"${R1_BASE}.qc.fastq.gz\" -2 \"${R2_BASE}.qc.fastq.gz\" -o \"assembly_${SAMPLE_NAME}\" --presets meta-sensitive --min-contig-len 500 -t \"$THREADS\""

if [ -n "$EXTREME_MEMORY_FLAG" ]; then
    echo "Adding extreme memory settings for MEGAHIT..."
    MEGAHIT_CMD="$MEGAHIT_CMD --memory 0.9"  # Use 90% of available memory
fi

echo "MEGAHIT command: $MEGAHIT_CMD"
eval $MEGAHIT_CMD

# Check if assembly succeeded
if [ ! -f "assembly_${SAMPLE_NAME}/final.contigs.fa" ]; then
    echo "ERROR: Assembly failed - no contigs produced"
    exit 1
fi

# Get assembly statistics
CONTIG_COUNT=$(grep -c ">" "assembly_${SAMPLE_NAME}/final.contigs.fa")
echo "Assembly completed: $CONTIG_COUNT contigs generated"

# Step 3: Viral BLAST analysis
echo ""
echo "Step 3: Viral Contamination BLAST"
echo "================================="

# Filter and sort contigs by length (>1000bp, largest first)
echo "Filtering contigs >1000bp and sorting by size for BLAST analysis..."
seqkit seq -m 1000 "assembly_${SAMPLE_NAME}/final.contigs.fa" | \
seqkit sort --by-length --reverse > "${SAMPLE_NAME}_contigs_filtered.fa"

FILTERED_COUNT=$(grep -c ">" "${SAMPLE_NAME}_contigs_filtered.fa")
echo "Filtered contigs for BLAST: $FILTERED_COUNT"

if [ "$FILTERED_COUNT" -gt 0 ]; then
    # Check for local contamination database first
    LOCAL_DB_PATH="/ref/sahlab/data/microbes_db/microbial_contaminants"
    
    # Add header to BLAST results
    echo -e "Query_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tMismatches\tGap_Opens\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score\tSubject_Title" > "${SAMPLE_NAME}_blast_all.tsv"
    
    # Use Blast environment for BLAST commands
    echo "Using Blast environment for BLAST..."
    
    eval "$(conda shell.bash hook)"
    conda activate Blast
    
    if [ -f "${LOCAL_DB_PATH}.nhr" ] || [ -f "${LOCAL_DB_PATH}.00.nhr" ]; then
        echo "Using local contamination database for fast BLAST..."
        echo "Note: Contigs are sorted by size (largest first) for priority analysis"
        
        blastn -query "${SAMPLE_NAME}_contigs_filtered.fa" \
               -db "$LOCAL_DB_PATH" \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
               -max_target_seqs 5 \
               -max_hsps 1 \
               -evalue 1e-10 \
               -num_threads "$THREADS" \
               >> "${SAMPLE_NAME}_blast_all.tsv"
    else
        echo "Local database not found, using remote BLAST (slower, may not process all contigs)..."
        echo "Note: Contigs are sorted by size (largest first) for priority analysis"
        echo "WARNING: Remote BLAST may timeout with many contigs. Consider building local database."
        
        # For remote BLAST, only process top 20 contigs to avoid timeout
        head -40 "${SAMPLE_NAME}_contigs_filtered.fa" > "${SAMPLE_NAME}_top20_contigs.fa"
        CONTIG_COUNT_TOP20=$(grep -c ">" "${SAMPLE_NAME}_top20_contigs.fa")
        echo "Processing top $CONTIG_COUNT_TOP20 contigs only (to avoid timeout)..."
        
        blastn -query "${SAMPLE_NAME}_top20_contigs.fa" \
               -remote -db nt \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
               -max_target_seqs 5 \
               -max_hsps 1 \
               -evalue 1e-10 \
               >> "${SAMPLE_NAME}_blast_all.tsv"
        
        rm -f "${SAMPLE_NAME}_top20_contigs.fa"
    fi
    
    # Filter for viral hits from the results
    echo "Filtering BLAST results for viral sequences..."
    if [ -s "${SAMPLE_NAME}_blast_all.tsv" ]; then
        # Add header to viral BLAST results
        echo -e "Query_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tMismatches\tGap_Opens\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score\tSubject_Title" > "${SAMPLE_NAME}_viral_blast.tsv"
        
        # Filter for viral hits (skip header line)
        tail -n +2 "${SAMPLE_NAME}_blast_all.tsv" | grep -i -E "(virus|viral|phage|viroid)" >> "${SAMPLE_NAME}_viral_blast.tsv" || true
        
        VIRAL_HITS=$(tail -n +2 "${SAMPLE_NAME}_viral_blast.tsv" | wc -l)
        echo "Found $VIRAL_HITS viral hits"
        
        # Create a top hits summary (best hit per contig, sorted by contig order)
        echo "Creating top hits summary..."
        echo -e "Contig_ID\tContig_Length\tSubject_ID\tPercent_Identity\tAlignment_Length\tQuery_Coverage\tE_value\tKingdom/Type\tSubject_Title" > "${SAMPLE_NAME}_top_hits.tsv"
        
        # For each contig, get the best hit and calculate query coverage
        # Parse contig lengths directly from MEGAHIT headers (e.g., ">k141_237 flag=1 multi=588.0000 len=12809")
        tail -n +2 "${SAMPLE_NAME}_blast_all.tsv" | \
        awk -F'\t' '
        BEGIN {
            # Read contig lengths from the filtered fasta file
            while ((getline line < "'${SAMPLE_NAME}_contigs_filtered.fa'") > 0) {
                if (line ~ /^>/) {
                    # Extract contig ID and length from MEGAHIT header
                    # Format: >k141_237 flag=1 multi=588.0000 len=12809
                    contig_id = substr(line, 2);
                    split(contig_id, id_parts, " ");
                    contig_id = id_parts[1];
                    
                    # Extract length from len= field
                    if (match(line, /len=([0-9]+)/, len_match)) {
                        contig_lengths[contig_id] = len_match[1];
                    }
                }
            }
            close("'${SAMPLE_NAME}_contigs_filtered.fa'");
        }
        {
            if (!seen[$1] || $11 < best_eval[$1]) {
                seen[$1] = 1;
                best_eval[$1] = $11;
                best_line[$1] = $0;
            }
        } 
        END {
            for (contig in best_line) {
                split(best_line[contig], fields, "\t");
                query_start = fields[7];
                query_end = fields[8];
                alignment_length = fields[4];
                coverage_length = query_end - query_start + 1;
                
                # Get contig length from parsed headers
                contig_length = contig_lengths[contig];
                if (contig_length == "") contig_length = "Unknown";
                
                # Calculate query coverage percentage
                if (contig_length != "Unknown" && contig_length > 0) {
                    query_coverage = sprintf("%.1f%%", coverage_length * 100 / contig_length);
                } else {
                    query_coverage = "Unknown";
                }
                
                desc = fields[13];
                kingdom = "Unknown";
                if (match(desc, "virus") || match(desc, "viral") || match(desc, "phage") || match(desc, "viroid")) kingdom = "Virus";
                else if (match(desc, "ycoplasma")) kingdom = "Mycoplasma";
                else if (match(desc, "acteria")) kingdom = "Bacteria";
                else if (match(desc, "ungi")) kingdom = "Fungi";
                else if (match(desc, "omo sapiens") || match(desc, "uman")) kingdom = "Human";
                else if (match(desc, "lant")) kingdom = "Plant";
                
                printf "%s\t%s\t%s\t%.2f%%\t%s\t%s\t%s\t%s\t%s\n", fields[1], contig_length, fields[2], fields[3], alignment_length, query_coverage, fields[11], kingdom, fields[13];
            }
        }' >> "${SAMPLE_NAME}_top_hits.tsv"
        
        echo "Top hits summary created: ${SAMPLE_NAME}_top_hits.tsv"
    else
        echo -e "Query_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tMismatches\tGap_Opens\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score\tSubject_Title" > "${SAMPLE_NAME}_viral_blast.tsv"
        echo -e "Contig_ID\tContig_Length\tSubject_ID\tPercent_Identity\tAlignment_Length\tQuery_Coverage\tE_value\tKingdom/Type\tSubject_Title" > "${SAMPLE_NAME}_top_hits.tsv"
        echo "No BLAST results obtained"
    fi
    
    # Alternative: Local viral database if available
    # blastn -query "${SAMPLE_NAME}_contigs_filtered.fa" \
    #        -db /ref/common/blastdb/viral_nt \
    #        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    #        -max_target_seqs 5 \
    #        -max_hsps 1 \
    #        -evalue 1e-5 \
    #        > "${SAMPLE_NAME}_viral_blast.tsv"
    
    echo "BLAST analysis completed"
else
    echo "No contigs >1000bp found - skipping BLAST analysis"
    touch "${SAMPLE_NAME}_viral_blast.tsv"
fi

# Step 4: Generate diagnostic report
echo ""
echo "Step 4: Generating Diagnostic Report"
echo "===================================="

cat > "${SAMPLE_NAME}_diagnostic_report.txt" << EOF
========================================
VIRAL CONTAMINATION DIAGNOSTIC REPORT
========================================
Sample: $SAMPLE_NAME
Reference: $ACCESSION
Analysis Date: $(date)
Working Directory: $(pwd)

========================================
MAPPING STATISTICS
========================================
Total Reads: $TOTAL_READS
Duplicate Reads: $DUPLICATE_READS (${DUPLICATION_RATE}%)
Unique Reads: $UNIQUE_READS

Raw Mapping: $MAPPED_READS reads (${MAPPING_PERCENT}% of total)
Deduplicated Mapping: $MAPPED_UNIQUE_READS reads (${DEDUPLICATED_MAPPING_PERCENT}% of unique)

Interpretation (based on deduplicated mapping):
- >70%: Good mapping to reference, likely correct organism
- 30-70%: Moderate mapping, possible mixed infection or variant
- <30%: Poor mapping, likely wrong reference or contamination

Note: High duplication rates (>50%) are expected in deep sequencing
of cell cultures for mutation detection and do not indicate quality issues.

========================================
ASSEMBLY STATISTICS
========================================
Total Contigs: $CONTIG_COUNT
Contigs >1000bp: $FILTERED_COUNT

========================================
BLAST RESULTS SUMMARY
========================================
EOF

# Add top hits summary to report
if [ -s "${SAMPLE_NAME}_top_hits.tsv" ] && [ $(wc -l < "${SAMPLE_NAME}_top_hits.tsv") -gt 1 ]; then
    echo "TOP HITS BY CONTIG (largest contigs first):" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    
    # Add top hits table to report
    tail -n +2 "${SAMPLE_NAME}_top_hits.tsv" | head -10 | \
    awk -F'\t' '{printf "%-15s %8s bp  %6s  %-12s  %s\n", $1, $2, $3, $6, substr($4,1,60)"..."}' >> "${SAMPLE_NAME}_diagnostic_report.txt"
    
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    
    # Summary by kingdom
    echo "CONTAMINATION SUMMARY:" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    tail -n +2 "${SAMPLE_NAME}_top_hits.tsv" | \
    awk -F'\t' '{kingdom[$6]++} END {for (k in kingdom) printf "  %s: %d contigs\n", k, kingdom[k]}' | \
    sort >> "${SAMPLE_NAME}_diagnostic_report.txt"
    
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
fi

# Add viral hits summary
if [ -s "${SAMPLE_NAME}_viral_blast.tsv" ] && [ $(tail -n +2 "${SAMPLE_NAME}_viral_blast.tsv" | wc -l) -gt 0 ]; then
    echo "VIRAL HITS DETAIL:" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    
    # Extract unique viral organisms from BLAST results
    tail -n +2 "${SAMPLE_NAME}_viral_blast.tsv" | \
    awk -F'\t' '{print $13}' | \
    sed 's/.*\[\([^]]*\)\].*/\1/' | \
    sort | uniq -c | sort -nr | head -10 >> "${SAMPLE_NAME}_diagnostic_report.txt"
    
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
fi

# Add file references
cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF

DETAILED RESULTS FILES:
- Top hits summary: ${SAMPLE_NAME}_top_hits.tsv
- All BLAST results: ${SAMPLE_NAME}_blast_all.tsv  
- Viral BLAST results: ${SAMPLE_NAME}_viral_blast.tsv
- Assembly contigs: assembly_${SAMPLE_NAME}/final.contigs.fa
EOF

if [ ! -s "${SAMPLE_NAME}_blast_all.tsv" ] || [ $(wc -l < "${SAMPLE_NAME}_blast_all.tsv") -le 1 ]; then
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "No significant BLAST hits found (E-value < 1e-10)" >> "${SAMPLE_NAME}_diagnostic_report.txt"
fi

# Add recommendations
cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF

========================================
RECOMMENDATIONS
========================================
EOF

# Generate recommendations based on deduplicated mapping percentage
if (( $(awk "BEGIN {print ($DEDUPLICATED_MAPPING_PERCENT < 30)}") )); then
    cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF
WARNING: Low deduplicated mapping percentage (${DEDUPLICATED_MAPPING_PERCENT}%)
- Check BLAST results for actual organism identity
- Consider using different reference genome
- Investigate potential sample contamination or mislabeling
- Raw mapping: ${MAPPING_PERCENT}% (includes ${DUPLICATION_RATE}% duplicates)
EOF
elif (( $(awk "BEGIN {print ($DEDUPLICATED_MAPPING_PERCENT < 70)}") )); then
    cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF
CAUTION: Moderate deduplicated mapping percentage (${DEDUPLICATED_MAPPING_PERCENT}%)
- Review BLAST results for mixed infections
- Check for strain variants or recombinants
- Consider manual review of mapping results
- Raw mapping: ${MAPPING_PERCENT}% (includes ${DUPLICATION_RATE}% duplicates)
EOF
else
    cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF
GOOD: High deduplicated mapping percentage (${DEDUPLICATED_MAPPING_PERCENT}%)
- Mapping suggests correct reference genome
- Proceed with standard analysis pipeline
- Monitor for any unusual variant patterns
- Raw mapping: ${MAPPING_PERCENT}% (includes ${DUPLICATION_RATE}% duplicates)
EOF
fi

# Step 5: Cleanup and final summary
echo ""
echo "Cleaning up temporary files..."
rm -f "${ACCESSION}.fasta."* # Remove BWA index files
rm -f "${SAMPLE_NAME}_quick.bam" "${SAMPLE_NAME}_quick.bam.bai"

echo ""
echo "========================================="
echo "GENERATING QC VISUALIZATION REPORT"
echo "========================================="

# Generate presentation-ready QC report
# Find the pipeline directory (where this script is located)
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QC_SCRIPT_PATH="${PIPELINE_DIR}/qc_with_simple_plots.py"

if [ -f "$QC_SCRIPT_PATH" ]; then
    echo "Creating QC visualization report..."
    # Run QC script from parent directory so it can find diagnostic_${SAMPLE_NAME}
    cd ..
    python "$QC_SCRIPT_PATH" "diagnostic_${SAMPLE_NAME}"
    QC_EXIT_CODE=$?
    cd "diagnostic_${SAMPLE_NAME}"
    
    if [ $QC_EXIT_CODE -eq 0 ]; then
        echo "✅ QC report generated: diagnostic_${SAMPLE_NAME}/diagnostic_${SAMPLE_NAME}_presentation_ready_report.html"
        echo "🌐 Open this file in a web browser for presentation-ready results"
    else
        echo "⚠️  QC report generation failed - check diagnostic output files"
    fi
else
    echo "⚠️  QC script not found at: $QC_SCRIPT_PATH"
    echo "   Manual QC generation: python ${PIPELINE_DIR}/qc_with_simple_plots.py diagnostic_${SAMPLE_NAME}"
fi

echo ""
echo "========================================="
echo "DIAGNOSTIC ANALYSIS COMPLETE"
echo "========================================="
echo "Report saved to: ${DIAGNOSTIC_DIR}/${SAMPLE_NAME}_diagnostic_report.txt"
echo "Key files generated:"
echo "  - ${SAMPLE_NAME}_mapping_stats.txt (mapping statistics)"
echo "  - ${SAMPLE_NAME}_fastp.html (read quality report)"
echo "  - assembly_${SAMPLE_NAME}/final.contigs.fa (assembled contigs - sorted by size)"
echo "  - ${SAMPLE_NAME}_top_hits.tsv (top hit per contig with kingdom classification)"
echo "  - ${SAMPLE_NAME}_blast_all.tsv (all BLAST results with headers)"
echo "  - ${SAMPLE_NAME}_viral_blast.tsv (viral BLAST results only)"
echo "  - ${SAMPLE_NAME}_diagnostic_report.txt (comprehensive summary report)"
echo ""
echo "Quick Summary:"
echo "  Mapping to $ACCESSION: ${MAPPING_PERCENT}%"
echo "  Contigs assembled: $CONTIG_COUNT"
echo "  Contigs BLASTed: $FILTERED_COUNT"
echo ""

# Display the report
echo "========================================="
echo "DIAGNOSTIC REPORT PREVIEW"
echo "========================================="
cat "${SAMPLE_NAME}_diagnostic_report.txt"

echo ""
echo "Diagnostic analysis completed at: $(date)"
echo "========================================="