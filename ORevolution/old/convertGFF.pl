#!/usr/bin/perl -w

use strict;
use warnings;

#my $usage = "\nThis script can be used to convert a exonerate gff back to the original genome coordinates.\n\n

#usage: perl \n\n";

#unless (scalar @ARGV == 2){die $usage}


my $gff = shift; #provide exonerate raw gff file

open(FILE1, $gff) or die "Cannot open $gff\n";
my $temmsmsmsm;
while (<FILE1>)
	{
	if (/^#/)
		{
		print $_;
		next
		}
	chomp;
	my @tab = split('\t',$_);
	my @coordinates = split /[:-]+/, $tab[0];
	$tab[3] = $coordinates[1] + $tab[3];
	$tab[4] = $coordinates[1] + $tab[4];
	print $coordinates[0],"\t";
	shift @tab;
	print join("\t",@tab),"\n";
}

close(FILE1);
