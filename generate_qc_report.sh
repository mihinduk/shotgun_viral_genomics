#!/bin/bash
# Generate QC visualizations for viral culture analysis
# Usage: ./generate_qc_report.sh [diagnostic_directories...]

set -euo pipefail

# Check if diagnostic directories provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <diagnostic_directory1> [diagnostic_directory2] ..."
    echo ""
    echo "Examples:"
    echo "  $0 diagnostic_SMS_10"
    echo "  $0 diagnostic_SMS_* "
    echo "  $0 /path/to/diagnostic_SMS_1 /path/to/diagnostic_SMS_2"
    echo ""
    echo "This script creates:"
    echo "  - Overall quality assessment plots"
    echo "  - Contamination detection summaries"  
    echo "  - Individual sample QC reports"
    echo "  - Publication-ready figures for presentations"
    exit 1
fi

echo "========================================="
echo "VIRAL CULTURE QC REPORT GENERATOR"
echo "========================================="
echo "Started at: $(date)"
echo "Processing $(echo $# | wc -w) diagnostic directories"

# Set up Python environment
echo "Setting up Python environment..."
source /ref/sahlab/software/anaconda3/bin/activate

# Use conda python
PYTHON_CMD="python"

# Check for required Python packages
echo "Checking Python dependencies..."
MISSING_PACKAGES=""

$PYTHON_CMD -c "import pandas" 2>/dev/null || MISSING_PACKAGES="$MISSING_PACKAGES pandas"
$PYTHON_CMD -c "import matplotlib" 2>/dev/null || MISSING_PACKAGES="$MISSING_PACKAGES matplotlib"
$PYTHON_CMD -c "import seaborn" 2>/dev/null || MISSING_PACKAGES="$MISSING_PACKAGES seaborn"
$PYTHON_CMD -c "import numpy" 2>/dev/null || MISSING_PACKAGES="$MISSING_PACKAGES numpy"

if [ -n "$MISSING_PACKAGES" ]; then
    echo "Missing Python packages:$MISSING_PACKAGES"
    echo "Try: conda install$MISSING_PACKAGES"
    exit 1
fi

echo "  ✓ All dependencies found"

# Create output directory with timestamp
OUTPUT_DIR="qc_report_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

echo "Output directory: $OUTPUT_DIR"
echo ""

# Validate diagnostic directories
VALID_DIRS=""
for dir in "$@"; do
    if [ -d "$dir" ]; then
        # Check if it has diagnostic files
        SAMPLE_NAME=$(basename "$dir" | sed 's/^diagnostic_//')
        REPORT_FILE="$dir/${SAMPLE_NAME}_diagnostic_report.txt"
        
        if [ -f "$REPORT_FILE" ]; then
            echo "  ✓ Valid diagnostic directory: $dir"
            VALID_DIRS="$VALID_DIRS $dir"
        else
            echo "  ✗ Missing diagnostic report: $REPORT_FILE"
        fi
    else
        echo "  ✗ Directory not found: $dir"
    fi
done

if [ -z "$VALID_DIRS" ]; then
    echo "Error: No valid diagnostic directories found"
    exit 1
fi

echo ""
echo "Generating QC visualizations..."

# Run the Python visualization script
$PYTHON_CMD create_qc_visualizations.py $VALID_DIRS -o "$OUTPUT_DIR"

# Create summary HTML report
echo ""
echo "Creating HTML summary report..."

cat > "$OUTPUT_DIR/qc_report.html" << EOF
<!DOCTYPE html>
<html>
<head>
    <title>Viral Culture QC Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2E86AB; border-bottom: 2px solid #2E86AB; }
        h2 { color: #A23B72; margin-top: 30px; }
        .summary-box { 
            background-color: #f5f5f5; 
            padding: 15px; 
            border-left: 4px solid #2E86AB; 
            margin: 20px 0; 
        }
        .image-container { 
            text-align: center; 
            margin: 20px 0; 
            border: 1px solid #ddd; 
            padding: 10px; 
        }
        .image-container img { 
            max-width: 100%; 
            height: auto; 
        }
        .footer { 
            margin-top: 50px; 
            padding-top: 20px; 
            border-top: 1px solid #ddd; 
            color: #666; 
            font-size: 12px; 
        }
    </style>
</head>
<body>
    <h1>Viral Culture Quality Control Report</h1>
    
    <div class="summary-box">
        <strong>Report Generated:</strong> $(date)<br>
        <strong>Samples Analyzed:</strong> $(echo $VALID_DIRS | wc -w)<br>
        <strong>Analysis Type:</strong> Viral contamination and quality assessment
    </div>
    
    <h2>Overall Quality Assessment</h2>
    <p>This plot shows key quality metrics across all samples including mapping percentages, 
    duplication rates, assembly quality, and sequencing depth.</p>
    
    <div class="image-container">
        <img src="culture_quality_overview.png" alt="Culture Quality Overview">
    </div>
    
    <h2>Contamination Detection Summary</h2>
    <p>Summary of contamination types detected across samples using local BLAST database 
    against known cell culture contaminants.</p>
    
    <div class="image-container">
        <img src="contamination_summary.png" alt="Contamination Summary">
    </div>
    
    <h2>Individual Sample Reports</h2>
    <p>Detailed QC analysis for each sample:</p>
    
EOF

# Add individual sample images to HTML
for dir in $VALID_DIRS; do
    SAMPLE_NAME=$(basename "$dir" | sed 's/^diagnostic_//')
    echo "    <h3>Sample: $SAMPLE_NAME</h3>" >> "$OUTPUT_DIR/qc_report.html"
    echo "    <div class=\"image-container\">" >> "$OUTPUT_DIR/qc_report.html"
    echo "        <img src=\"${SAMPLE_NAME}_detailed_qc.png\" alt=\"${SAMPLE_NAME} QC Report\">" >> "$OUTPUT_DIR/qc_report.html"
    echo "    </div>" >> "$OUTPUT_DIR/qc_report.html"
done

cat >> "$OUTPUT_DIR/qc_report.html" << EOF
    
    <div class="footer">
        <p>Generated by Viral Genomics Pipeline QC Module</p>
        <p>For questions about interpretation, consult the pipeline documentation.</p>
    </div>
</body>
</html>
EOF

echo "  ✓ HTML report created: $OUTPUT_DIR/qc_report.html"

# Create README
cat > "$OUTPUT_DIR/README.txt" << EOF
Viral Culture QC Report
======================
Generated: $(date)

This directory contains quality control visualizations for viral culture analysis.

FILES:
------
qc_report.html                  - Complete HTML report (open in web browser)
culture_quality_overview.png    - Overall quality metrics across samples
contamination_summary.png       - Contamination detection summary
[sample]_detailed_qc.png        - Individual sample QC reports

INTERPRETATION:
--------------
Mapping Quality:
- Green (>70%): Excellent - high confidence in organism identity
- Orange (30-70%): Moderate - possible mixed infection or variant
- Red (<30%): Poor - wrong reference or heavy contamination

Duplication Rate:
- Green (<60%): Good - normal for viral sequencing
- Orange (60-80%): Moderate - acceptable for deep sequencing
- Red (>80%): High - may indicate PCR bias or low complexity

Contamination:
- Review detected contaminants for culture purity assessment
- Mycoplasma detection is critical for cell culture work
- High bacterial contamination may indicate culture problems

USAGE FOR PRESENTATIONS:
-----------------------
All images are publication-ready (300 DPI) and suitable for:
- Scientific presentations
- Reports and publications
- Quality control documentation
- Troubleshooting culture issues

For best results in presentations:
1. Use culture_quality_overview.png for overall summary
2. Use individual sample plots for detailed discussion
3. Reference contamination_summary.png for purity assessment
EOF

echo ""
echo "========================================="
echo "QC REPORT GENERATION COMPLETE"
echo "========================================="
echo "Report location: $OUTPUT_DIR/"
echo ""
echo "Key files for your presentation:"
echo "  📊 culture_quality_overview.png - Main quality summary"
echo "  🦠 contamination_summary.png - Contamination overview"
echo "  📄 qc_report.html - Complete interactive report"
echo ""
echo "To view the complete report:"
echo "  open $OUTPUT_DIR/qc_report.html"
echo ""
echo "All plots are publication-ready (300 DPI) for presentations!"
echo "========================================="