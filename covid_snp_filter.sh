#!/bin/bash

# Filter vcf files

for i in *_vars.vcf; do

	F=`basename $i _vars.vcf`;

	echo "Filtering: $F";

	lofreq filter -i $i -o "$F"_vars.filt.vcf \
		-v 75;

	echo "Done"

done

