#!/usr/bin/perl -w
#usage: perl bls2gff.pl <bls tab output> > <gff3>

use strict;
use warnings;

my $usage = '\nThis script can be used to convert blast output to gff.\n\n

usage: perl bls2gff.pl <bls tab file> > <gff> \n\n';

unless (scalar @ARGV == 1){die $usage}

my $bls = $ARGV[0];
my ($strand,$start,$stop);
my $counter =1;
my $previous="NONE";
open IN, "< $bls" or die "Can't open $bls!";
while (<IN>)	{

	if (/^#/){
		next;
	}
	chomp;
	my @tab = split('\t',$_); #11 elements	
	#	queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore
	#	tab[0] = AechOr158PC-RA:exon2:scaffold7:2428414-2428644
	if ($previous eq $tab[0]){
		$counter++
		}else{
		$counter=1;
		}
	my $att = "ID=".$tab[0].":hit".$counter.":".$tab[6]."-".$tab[7].";";


	$previous = $tab[0];
	if ($tab[8]>$tab[9]){
		$strand = "-";
		$start = $tab[9];
		$stop = $tab[8];
	}else{
		$strand = "+";
		$start = $tab[8];
		$stop = $tab[9];
	}
		
	
	print $tab[1],"\tblast\tHSP\t",$start,"\t",$stop,"\t",$tab[2],"\t",$strand,"\t.\t",$att,"\n";
}
close (IN);
