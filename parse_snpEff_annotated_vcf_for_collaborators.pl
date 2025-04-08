#!/usr/bin/env perl
# parse_snpEff_annotated_vcf_for_collaborators.pl
# 2023_12_01
use warnings;
use strict;

# This script takes the annotated tsv output from snpEff and parses it to generate a tab-delimited text file of hits that reach your depth requirement

# $ARGV[0] = annotated VCF from snpEff (snpEFF.ann.tsv)
# $ARGV[1] = depth requirement

my $depth =$ARGV[1];
my $doc = $ARGV[0];
$doc =~ s/.tsv//; 
my $outfile = $doc. "_$depth" .".tsv"; 
my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$EFFECT,$PUTATIVE_IMPACT,$GENE_NAME,$GENE_ID,$FEATURE_TYPE,$FEATURE_ID,$TRANSCRIPT_TYPE,
	$EXON_INTRON_RANK,$HGVSc,$HGVSp,$cDNA_POSITION_AND_LENGTH,$CDS_POSITION_AND_LENGTH,$PROTEIN_POSITION_AND_LENGTH,$DISTANCE_TO_FEATURE,$ERROR);

# Parse the VCF
open(IN, $ARGV[0]) or die "Couldn't find the annotated VCF file $!";
open (OUT, ">$outfile") or die "Couldn't create a file for your clean taxonomy $!";
print OUT "CHROM\tPOS\tID\tREF\tALT\tQUAL\tINFO\tTotal_Depth\tAllele_Frequency\tstrand_bias\tDP4\tEFFECT\tPUTATIVE_IMPACT\tGENE_NAME\tGENE_ID\tFEATURE_TYPE\tFEATURE_ID\tTRANSCRIPT_TYPE\t";
print OUT "HGVSc\tHGVSp\tcDNA_POSITION_AND_LENGTH\tCDS_POSITION_AND_LENGTH\tPROTEIN_POSITION_AND_LENGTH\tERROR\n";
	while (my $line = <IN>) {
		next if ($line =~/^#/);
		chomp ($line); 
		($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$EFFECT,$PUTATIVE_IMPACT,$GENE_NAME,$GENE_ID,$FEATURE_TYPE,$FEATURE_ID,$TRANSCRIPT_TYPE,$EXON_INTRON_RANK,$HGVSc,$HGVSp,$cDNA_POSITION_AND_LENGTH,$CDS_POSITION_AND_LENGTH,$PROTEIN_POSITION_AND_LENGTH,$DISTANCE_TO_FEATURE,$ERROR) = split(/\t/, $line, 23);
		$ERROR=~ s/,[ATGC]//; #print "\$ERROR = $ERROR\n";
		my ($DP, $AF, $SB, $DP4, $other) = split(/;/, $INFO, 5); 
		$DP=~ s/DP=//; $AF=~ s/AF=//; $SB=~ s/SB=//; $DP4=~ s/DP4=//; #print "\$DP = $DP\t\$AF = $AF\t\$SB = $SB\t\$DP4 = $DP4\t\$other = $other\n";
		next if ($DP <$depth);
		print OUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$INFO\t$DP\t$AF\t$SB\t$DP4\t$EFFECT\t$PUTATIVE_IMPACT\t$GENE_NAME\t$GENE_ID\t$FEATURE_TYPE\t$FEATURE_ID\t$TRANSCRIPT_TYPE\t";
		print OUT "$HGVSc\t$HGVSp\t'$cDNA_POSITION_AND_LENGTH'\t'$CDS_POSITION_AND_LENGTH'\t'$PROTEIN_POSITION_AND_LENGTH'\t$ERROR\n";
	}

close(IN);
close(OUT);

print "\nFons vitae caritas. Love is the fountain of life.\noutfile = $outfile\n";