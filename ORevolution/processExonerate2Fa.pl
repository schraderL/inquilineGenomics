#!/usr/bin/perl -w

use strict;
use warnings;

#my $usage = "\nThis script can be used to convert a exonerate output to a gff with the original genome coordinates and a fasta file.\n\n

#usage: perl \n\n";

#unless (scalar @ARGV == 2){die $usage}


my $exonerate = shift; #provide exonerate raw output file

open(FILE1, $exonerate) or die "Cannot open $exonerate\n";
my $switch=0;
while (<FILE1>)
	{
	if (/^###---FASTASTART---/)	
	{$switch=1;next}
	if (/^###---FASTAEND---/)	
	{$switch=0;next}
	if (/^#/){next}
	if (/-- completed exonerate analysis/)
	{last}
	if ($switch==1){
#		if (/^>/){
	print $_;
#		}
	}
	}
close(FILE1);	
	


