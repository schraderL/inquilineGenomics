#!/usr/bin/perl -w

use strict;
use warnings;

my $exonerate = shift; #provide exonerate2EVM converted file raw output file

open(FILE1, $exonerate) or die "Cannot open $exonerate\n";
while (<FILE1>)
	{
	if (/^$/){
		print $_;
		next;}
	chomp;
	my @tab = split('\t',$_);
	my @coordinates = split /[:-]+/, $tab[0];
	$tab[3] = ($coordinates[1] + $tab[3])-1;
	$tab[4] = ($coordinates[1] + $tab[4])-1;
	print $coordinates[0],"\t";
	shift @tab;
	print join("\t",@tab),"\n";
	}
	

close(FILE1);

