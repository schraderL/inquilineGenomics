#!/usr/bin/perl -w
#usage: perl GOannotation2topGO.pl Acromyrmex_insinuator.function.GO > gene2GO.WE.list
use strict;
use warnings;
my $annotation = shift; #provide sequence file

open(FILE1, $annotation) or die "Cannot open $annotation\n";

while (<FILE1>)
	{	
	chomp; 
	my @tab = split("\t",$_);
	my ($seq) = $tab[0];
	my ($count) = $tab[1];

	print $seq,"\t";
		
	for (my $i=2; $i <= $count+1;  $i++){
		$tab[$i] =~ /(GO:\d+)/;
		print $1;
		unless ($i == $count+1){
			print ", ";
			}
		}	
	print "\n";
}
	
close(FILE1);	
	


