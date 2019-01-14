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
	
	if (/^# --- START OF GFF DUMP ---/)
	{$switch=1;
	next}
	if (/^# --- END OF GFF DUMP ---/)
	{$switch=0;
	next}
	if ($switch==1){
	if (/^#/)
		{
		#print $_;
		next
		}	
		chomp; 
		my @tab = split('\t',$_);
		my @coordinates = split /[:-]+/, $tab[0];
		$tab[3] = ($coordinates[1] + $tab[3])-1;
		$tab[4] = ($coordinates[1] + $tab[4])-1;
		print $coordinates[0],"\t";
		shift @tab;
		print join("\t",@tab),"\n";
}
	}
	
close(FILE1);	
	


