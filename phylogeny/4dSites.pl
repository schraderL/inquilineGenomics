#!/usr/bin/perl
use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use POSIX;


my $infile = shift;
my $align = Bio::AlignIO->new(-file   => $infile,
                              -format => shift);

my $outfile = shift;
my $out = Bio::AlignIO->new(-file => ">$outfile",
                            -format => "fasta");


#print "Cobs\tPbar\tlength\ttrimmed.to\n";

my $result;
my $length;
my $good;
my %seq;
my $i = 1;
my $j = 0;
  #4d codons
  my %codons = (
		'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
	);

	while ( my $aln = $align->next_aln) {
		for (my $i = 1 ; $i < $aln->length; $i += 3) {
			my $codon = $aln->slice($i,$i+2);
			chomp;
			#check for each codon if alignment of the first two bases is the same for each species
			unless (substr(uc($codon->consensus_iupac()),0,2) =~/[^AGTC]/){
				$good = "F";
				#check if the codon is in the 4d codon table above
				foreach my $CDS ($codon->each_seq) {
					#this might be incorrect, because $good will be true if even a single sequence is a 4D codon
					if (exists($codons{$CDS->seq})){$good = "T"};
				}
				#if the codon is identical in the first 2 bases and the third base is not changing the coded aa continue
		 		if ($good eq "T"){
						$j++;
					foreach my $CDS ($codon->each_seq) {
						#extract the last base (ie the 4d site) of the codon
						$seq{$CDS->display_id} .= substr($CDS->seq,2,3);#." $i ";
					}
				}
	  	}
		}

		# SOME ALIGNMENT STATISTICS
		print $infile,"\t";
		print "Aln length:\t",$aln->length,"\t";
		print "Degenerate sites:\t",$j,"\n";
		if ($aln->length % 3 !=0) {print "ERROR $infile: Alignment length is not a multiple of 3: ",$aln->length,"\n"}
		foreach my $seqSTATS ($aln->each_seq) {
			my $aa = $seqSTATS->translate->seq;
			if (3*length($aa) != $aln->length){
			print STDERR "ERROR  $infile: The length of sequence ",$seqSTATS->display_id," does not match the length of the alignment!",length($aa),"-",$aln->length,"/3 != 0\n" ;
			}

	  	}
	}
#my $degenAln = Bio::SimpleAlign->new(-seqs => %seq);
my $degenSite = Bio::SimpleAlign->new(-format => "fasta");
foreach my $key ( keys %seq )
{
my $seqs = Bio::LocatableSeq->new(-seq => $seq{$key},
                    -id  => $key."/4D sites",
                    -start => 1,
                    -end   => length($seq{$key}));
$degenSite->addSeq($seqs);
}
$out->write_aln($degenSite);
