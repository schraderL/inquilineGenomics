
#!/usr/bin/perl -w
#usage: perl translateCDS2AA.pl <CDS.fa> > <AA.fa>
use strict;
use warnings;
use Bio::SeqIO;
my $file = $ARGV[0];
my $inseq = Bio::SeqIO->new(-file   => $file);
my $counter=0;
while (my $seq = $inseq->next_seq) {
	$counter++;
	print ">",$seq->display_id,"-",$counter," ",$seq->desc,"\n";
	print $seq->translate->seq,"\n";
}