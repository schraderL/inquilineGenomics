
#!/usr/bin/perl -w
#usage: perl addGeneGFF.pl <gff3 file> > wgene.gff3

use strict;
use warnings;

my $usage = '\nThis script can be used to add "gene" entries to a gff .\n\n

usage: perl addGeneGFF.pl <gff3 file> > wgene.gff3 \n\n';

unless (scalar @ARGV == 1){die $usage}

my $gff = $ARGV[0];

open IN, "< $gff" or die "Can't open $gff!";
while (<IN>)	{

	if (/^#/){
		print $_;
		next;
	}
	chomp;
	my @tab = split('\t',$_);
	
	if ($tab[2] eq "gene"){
		my @att = split(';',$tab[8]);
		$tab[8] = $att[0]."-gene";
		print join("\t", @tab),"\n";
		$tab[2] = "mRNA";
		#my $ID = $tab[8];
		$tab[8] =~ s/ID/Parent/g;
		print join("\t", @tab),";",$att[0],"\n";#,"-gene\n";
	}else { 
		print join("\t", @tab),"\n";
	}
}
close (IN);
