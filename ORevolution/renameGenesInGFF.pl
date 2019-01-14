#!/usr/bin/perl -w
use strict;
use warnings;
use File::Slurp;

my $usage = "\nUse this script to rename gff entries based on a tsv with the following format: <old_name>\t<new_name>.\n\n

#usage: perl cleanMcKenzieGFF.pl <notReannotated.McKenzie.gff>\n\n";
unless (scalar @ARGV == 2){die $usage}

my $gff = $ARGV[0]; #provide gff file
my $rename = $ARGV[1]; #provide gff file
my (%genes);
open(FILE1, $rename) or die "Cannot open $rename\n";
while (<FILE1>)
	{
		chomp;
		my @tab = split ('\t');
		$genes{$tab[0]} = $tab[1];
	}
close(FILE1);

my $fullGff = read_file($gff);
foreach my $gene ( keys %genes )
{
	if ($gene ne $genes{$gene}){
  $fullGff =~ s/($gene)/$genes{$1}-/g;

}
}

print $fullGff;
