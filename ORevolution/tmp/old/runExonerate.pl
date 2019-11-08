
#!/usr/bin/perl -w
#usage: perl ../scripts/runExonerate.pl config_$species"2"
use strict;
use warnings;
use Bio::DB::Fasta;
use File::Basename;
my $usage = "\nThis script can be used to annotate target loci with exonerate.\n\n

usage: perl ../scripts/runExonerate.pl config_species2 \n\n";

unless (scalar @ARGV == 1){die $usage}


# extract the genomic locus at which a target gene is suspected.

my $config = $ARGV[0];
my %config;
open IN, "< $config" or die "Can't open $config!";
while (<IN>)	{
    next if (/^#/);
    next if (/^$/);
	chomp;
	my @config = split /=/;
	$config{$config[0]}=$config[1];
}
close (IN);

################################
# example config file:
#BED File containing potential loci
#BED=Tsep.NCBI.bed
#target genome
#GENOME=/usr/local/home/lschrader/data/OR/genomes/NCBI/Tsep/GCF_001594115.1_Tsep1.0_genomic.fna
#target folder for target loci fasta
#OUTFOLDER=out
#query fasta file
#QUERY=AechOR.fa
#out folder for query fasta
#OUTQUERY=AechFa
################################


my $bed = $config{BED}; #provide hsp file
my $genome = $config{GENOME}; #provide genome file
my $outfolder = $config{OUTFOLDER}; #provide outfolder
my $queryFa = $config{QUERY};
my $outQuery = $config{OUTQUERY};
`mkdir -p $outQuery`;
my $range = $config{RANGE};;
my $db = Bio::DB::Fasta->new( $queryFa );
my ($start,$stop);
my $counter = 1;
#my $cutoff = shift; #provide cutoff (in aa of OR receptor)
open(FILE1, $bed) or die "Cannot open $bed\n";
while (<FILE1>)
	{
	chomp;
	my @tab = split('\t',$_);
	unless ($outfolder =~ /\/$/){$outfolder.="/"}
	#print $outfolder;
	`mkdir -p $outfolder`;
	if ($tab[1] < $tab[2]){
		$start = $tab[1]-$range;
		$stop  = $tab[2]+$range;
	}else{
		$start = $tab[2]-$range;
		$stop  = $tab[1]+$range;
	}
	if ($start < 0){$start = 0}
	if ($stop < 0){$stop = 0}
	`samtools faidx  $genome $tab[0]:$start-$stop > $outfolder$tab[3].$counter.target.fa`;
	getFasta($tab[3]);
  $counter++;
	}

close(FILE1);



#get $tab[3] as fasta entry
sub getFasta {
	my $query = shift;
	open(my $fh, ">$outQuery/" . basename($query).".fa") or die "Could not open file '$query'.fa $!";
    my $sequence = $db->seq($query);
    if  (!defined( $sequence )) {
            die "Sequence $query not found. \n"
    }
    print $fh ">$query\n", "$sequence\n";
}
