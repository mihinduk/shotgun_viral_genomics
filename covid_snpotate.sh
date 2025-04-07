#!/bin/bash

# Run snpEFF using defined reference genome data

for i in *_vars.filt.vcf; do 

	F=`basename $i _vars.filt.vcf`;

	echo "Processing: $F";

	java -jar -Xmx4g -jar /home/kathiem/snpEff/snpEff.jar NC_045512.2 $i -s "$F"_summary.html > "$F".snpEFF.ann.vcf;

	grep -v "^##" "$F".snpEFF.ann.vcf | \
		tail -n+2 | \
		cut -f8 | \
		sed 's/|/\t/g' | \
		cut -f1-16 | \
		sed '1i INFO\tEFFECT\tPUTATIVE_IMPACT\tGENE_NAME\tGENE_ID\tFEATURE_TYPE\tFEATURE_ID\tTRANSCRIPT_TYPE\tEXON_INTRON_RANK\tHGVSc\tHGVSp\tcDNA_POSITION_AND_LENGTH\tCDS_POSITION_AND_LENGTH\tPROTEIN_POSITION_AND_LENGTH\tDISTANCE_TO_FEATURE\tERROR' > "$F".snpEFF.ann.tmp;

	grep -v "^##" "$F".snpEFF.ann.vcf | \
		cut -f1-7 > "$F".ann.base.vcf;

	paste "$F".ann.base.vcf "$F".snpEFF.ann.tmp > "$F".snpEFF.ann.tsv;

	rm "$F".snpEFF.ann.tmp;
	rm "$F".ann.base.vcf;

done

echo "Processing complete"
