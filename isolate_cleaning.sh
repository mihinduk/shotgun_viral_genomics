#!/bin/bash
set -e
set -u
set -o pipefail

# Script to clean Illumina NextSeq runs
# Dependencies:
# fastp: https://github.com/OpenGene/fastp

# Create output directory
mkdir -p ./cleaned_seqs
OUT=./cleaned_seqs

# Clean sequences with fastp
for i in *_R1.fastq.gz; do

	F=`basename $i _R1.fastq.gz`;

	fastp -w 32 -q 30 -i "$F"_R1.fastq.gz -I "$F"_R2.fastq.gz -o $OUT/"$F"_R1.qc.fastq.gz -O $OUT/"$F"_R2.qc.fastq.gz \
		-h $OUT/fastp_report.html;

done
