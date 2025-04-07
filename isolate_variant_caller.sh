#!/bin/bash
set -e
set -u
set -o pipefail

# Script to map (along with sorting and duplicate removal) cleaned reads to a reference genome.
# Also performs variant calling with lofreq

# Dependencies:
# bwa: https://github.com/lh3/bwa
# samtools: https://github.com/samtools/samtools
# lofreq: https://csb5.github.io/lofreq/

# Help function
helpFunction()
{
        echo "Script to call variants from viral isolate sequencing."
        echo ""
        echo "Syntax: scriptTemplate [-h | -r | -t]"
        echo "options:"
	echo "-h Display help"
        echo "-r Reference Genome (FASTA format)"
	echo "-t Number of threads (default: 1)"
        echo ""
        exit 1 # Exit script after printing help
}

# Retrieve and set options/flags
while getopts "hr:t:" opt; do
        case $opt in
                h) # display Help
                helpFunction
                ;;

                r) # reference genome
                REF="$OPTARG"
                echo "The reference genome is $OPTARG"
		echo ""
                ;;

		t) # number of threads
		THREADS="$OPTARG"
		echo "The number of threads is set to $OPTARG"
		echo ""
		;;

                \?) # incorrect option
                echo "Usage: cmd [-h] [-r] [-t]"
                exit
                ;;
        esac
done

# BWA index genome
	bwa index "$REF"

mkdir -p ./mapping
mkdir -p ./variants

OUTm=./mapping
OUTv=./variants

for i in *_R1.qc.fastq.gz; do

	F=`basename $i _R1.qc.fastq.gz`;

	# Align
	bwa mem -t "$THREADS" "$REF" "$F"_R1.qc.fastq.gz "$F"_R2.qc.fastq.gz > $OUTm/"$F".sam;

	# Fix mates
        samtools fixmate -O bam,level=1 -m --threads "$THREADS" $OUTm/"$F".sam $OUTm/"$F".fixmate.bam;

	# Sort
	samtools sort --threads "$THREADS" -O bam $OUTm/"$F".fixmate.bam > $OUTm/"$F".bam;

	# Mark PCR duplicates
	samtools markdup --threads "$THREADS" -S $OUTm/"$F".bam $OUTm/"$F".dedupe.bam;

	# LoFreq Realign
	lofreq viterbi -f "$REF" $OUTm/"$F".dedupe.bam | samtools sort - --threads "$THREADS" > $OUTm/"$F".lofreq.realign.bam;

	# Insert indel qualities
	lofreq indelqual --dindel -f "$REF" $OUTm/"$F".lofreq.realign.bam | samtools sort - --threads "$THREADS" > $OUTm/"$F".lofreq.indel.bam;

	lofreq alnqual -b $OUTm/"$F".lofreq.indel.bam "$REF" > $OUTm/"$F".lofreq.final.bam;

	samtools index $OUTm/"$F".lofreq.final.bam;

	# Call variants with LoFreq
	lofreq call-parallel --pp-threads ${THREADS} --force-overwrite --no-default-filter --call-indels -f "$REF" -o $OUTv/"$F"_vars.vcf $OUTm/"$F".lofreq.final.bam;

done

echo "Run complete"
echo "The reference used was: $REF"
