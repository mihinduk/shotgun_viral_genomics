#!/bin/bash
set -e

# Script to fix path references and make scripts compatible with the pipeline

echo "Fixing scripts for compatibility with the pipeline..."

# Update path in isolate_cleaning.sh to support both naming patterns
sed -i.bak 's/for i in \*_R1.fastq.gz; do/for i in $(find . -name "*_R1*.fastq.gz"); do/g' isolate_cleaning.sh
sed -i.bak 's/F=`basename $i _R1.fastq.gz`;/F=`basename $i | sed -E "s\/_R1(.*?)\.fastq\.gz//g"`;/g' isolate_cleaning.sh

# Update path in isolate_variant_caller.sh to support both naming patterns
sed -i.bak 's/for i in \*_R1.qc.fastq.gz; do/for i in $(find . -name "*_R1*.qc.fastq.gz"); do/g' isolate_variant_caller.sh
sed -i.bak 's/F=`basename $i _R1.qc.fastq.gz`;/F=`basename $i | sed -E "s\/_R1(.*?)\.qc\.fastq\.gz//g"`;/g' isolate_variant_caller.sh

# Update path in covid_snpotate.sh to make accession configurable
sed -i.bak 's/java -jar -Xmx4g -jar \/home\/kathiem\/snpEff\/snpEff.jar NC_045512.2/java -jar -Xmx4g $SNPEFF_JAR $ACCESSION/g' covid_snpotate.sh

# Add header to covid_snpotate.sh to accept parameters
cat > covid_snpotate.sh.new << 'EOF'
#!/bin/bash
set -e

# Get parameters
ACCESSION=${1:-NC_045512.2}
SNPEFF_JAR=${2:-snpEff.jar}

echo "Using accession: $ACCESSION"
echo "Using snpEff jar: $SNPEFF_JAR"

EOF

# Append original content
cat covid_snpotate.sh >> covid_snpotate.sh.new
mv covid_snpotate.sh.new covid_snpotate.sh
chmod +x covid_snpotate.sh

# Make all scripts executable
chmod +x *.sh *.pl

echo "Scripts fixed successfully!"