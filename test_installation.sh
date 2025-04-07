#!/bin/bash
set -e

# Test script for verifying the Shotgun Viral Genomics Pipeline installation
echo "Testing Shotgun Viral Genomics Pipeline installation..."

# Check if conda environment is activated
if [[ -z "$CONDA_PREFIX" || "$CONDA_DEFAULT_ENV" != "viral_genomics" ]]; then
    echo "Error: viral_genomics conda environment is not activated."
    echo "Please run: conda activate viral_genomics"
    exit 1
fi

# Check dependencies
echo "Checking dependencies..."
commands=(
    "seqkit --version"
    "fastp --version"
    "bwa"
    "samtools --version"
    "lofreq version"
    "perl --version"
    "python --version"
)

for cmd in "${commands[@]}"; do
    echo "  Testing: $cmd"
    if ! eval "$cmd" &> /dev/null; then
        echo "Error: Failed to run '$cmd'"
        echo "Please check your installation."
        exit 1
    fi
    echo "    ✓ Passed"
done

# Check for snpEff
echo "  Testing: snpEff"
if command -v snpEff &> /dev/null; then
    echo "    ✓ snpEff found in PATH"
elif [ -f snpEff.jar ]; then
    echo "    ✓ snpEff.jar found in current directory"
else
    echo "    ⚠ Warning: snpEff not found in PATH or current directory"
    echo "      You'll need to specify the path with --snpeff-jar when running the pipeline"
fi

# Check for viral_pipeline.py
echo "  Testing: viral_pipeline.py"
if [ -f viral_pipeline.py ] && [ -x viral_pipeline.py ]; then
    echo "    ✓ viral_pipeline.py is executable"
else
    echo "    ⚠ Warning: viral_pipeline.py not found or not executable"
    echo "      Run: chmod +x viral_pipeline.py"
fi

# Check for shell scripts
echo "  Testing: Shell scripts"
for script in isolate_cleaning.sh isolate_variant_caller.sh covid_snpotate.sh covid_snp_filter.sh; do
    if [ -f "$script" ] && [ -x "$script" ]; then
        echo "    ✓ $script is executable"
    else
        echo "    ⚠ Warning: $script not found or not executable"
        echo "      Run: chmod +x $script"
    fi
done

# Check for Perl script
echo "  Testing: Perl script"
if [ -f parse_snpEff_annotated_vcf_for_collaborators.pl ] && [ -x parse_snpEff_annotated_vcf_for_collaborators.pl ]; then
    echo "    ✓ parse_snpEff_annotated_vcf_for_collaborators.pl is executable"
else
    echo "    ⚠ Warning: parse_snpEff_annotated_vcf_for_collaborators.pl not found or not executable"
    echo "      Run: chmod +x parse_snpEff_annotated_vcf_for_collaborators.pl"
fi

echo ""
echo "Installation test complete!"