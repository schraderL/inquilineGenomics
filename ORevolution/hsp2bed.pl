#!/usr/bin/perl -w
#usage: perl hsp2bed.pl Tsep.hsp.out > Tsep.hsp.bed
use strict;
use warnings;

my $usage = "\nThis script can be used to select HSP hits from the output of identify_genes.pl and reformat it to bed format.\n\n

usage: perl hsp2bed.pl <hsp.out> <cutoff_in_aa>  > hsp.out.bed\n\n";

unless (scalar @ARGV == 2){die $usage}


my $hsp = shift; #provide hsp file
my $cutoff = shift; #provide cutoff (in aa of OR receptor)
open(FILE1, $hsp) or die "Cannot open $hsp\n";

while (<FILE1>)
	{	
	chomp; 
	my @tab = split('\s',$_);
	if (scalar(@tab) >= 7){
		if ($tab[6]-$tab[5]>= $cutoff){
			print $tab[0],"\t",$tab[2],"\t",$tab[3],"\t",$tab[4],"\t",$tab[7],"\t",$tab[1],"\t",$tab[5],"\t",$tab[6],"\t";
			print "\n";
		}
	}
}
	
close(FILE1);	
	


