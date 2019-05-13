#!/usr/bin/perl -w
use strict;
use warnings;

my $usage = "\nThis script can be used to check if entries in a protein fasta file contain an M at the start and a * at the end.\n\n

		usage: checkFastaGeneModel.pl fasta.fa \n\n";


use Getopt::Std;
use Getopt::Long;
use Bio::SeqIO;
use Bio::DB::SeqFeature;

my $file   = shift;
my $inseq = Bio::SeqIO->new(-file   => "$file");

while (my $seq = $inseq->next_seq) {
	my $first = substr($seq->seq,0,1);
	my $last = substr($seq->seq,-1);
	if ($first eq "M" && $last eq "\*"){
	print ">",$seq->display_id,"\n",$seq->seq,"\n";
	}
    
}
