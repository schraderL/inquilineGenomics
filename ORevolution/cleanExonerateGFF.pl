#!/usr/bin/perl -w

use strict;
use warnings;

#my $usage = "\nThis script can be used to convert a exonerate output to a gff with the original genome coordinates and a fasta file.\n\n

#usage: perl \n\n";

#unless (scalar @ARGV == 2){die $usage}


my $exonerate = shift; #provide exonerate raw output file

open(FILE1, $exonerate) or die "Cannot open $exonerate\n";
my $switch=0;
my @atts;
my ($num,$ID,$att);
my $geneCount=0;
while (<FILE1>)
	
	{
	chomp;	
	my @tab = split('\t',$_);
	next if ($tab[2] eq "similarity");
	if ($tab[2] eq "gene") {
		$geneCount++;
		$att=pop @tab;
	 	@atts = split(' ; ',$att);
	 	$atts[0] =~ s/gene_id //;
		$atts[1] =~ s/sequence //;
		$ID = $atts[1]."-".$geneCount;
		print join("\t",@tab),"\tID=",$ID,"-gene;","\n";	
		$tab[2]="mRNA";
		print join("\t",@tab),"\tParent=",$ID,"-gene;ID=",$ID,";\n";	
	}
	if ($tab[2] ne "mRNA") {
		for my $i (0 .. 7) {
			print $tab[$i],"\t"	
		}
		print "Parent=",$ID,";","\n";		
	}
	}
	
close(FILE1);	
	


