#!/usr/local/bin/perl -w
#Lukas Schrader
 use strict;
 use warnings;
# use  List::MoreUtils;
 use syntax 'junction';
 use Data::Dumper qw(Dumper);

my $usage = "\nquality control for MCL gene family clustering  \n\n
			usage: perl filterMCLcluster.pl IPR.annotations.tsv final.seq.mci.I15";

my $famID=1;
my $isAnn=0;
my $infile1 = $ARGV[0]; # annotation file
my $infile2 = $ARGV[1]; # mcl.output
my (%HoA,%HoA2,%seen);
my ($count,$bestIPR,@uniqueIPR,$unique,@members,@ann,@IPR);

open(FILE1, $infile1) or die "Cannot open $infile1\n";
open(FILE2, $infile2) or die "Cannot open $infile1\n";

#Store IPR annotations in HoA, where each element of the hash is referenced by the geneID
while(<FILE1>){
	chomp;
	@ann = split("\t");
	@IPR = @ann[ 2 .. $#ann ];
	$HoA2{$ann[0]} = [@IPR];
	}
close(FILE1);


# Store all the members of a Cluster in a hash of arrays, with each array containing every gene in the cluster.
while(<FILE2>){
	chomp;
	@members = split("\t");
	$HoA{$famID} = [@members];
	$famID++;

}
close(FILE2);

print "ClusterID\tFamSize\tAnn\tNo.Top\tNo.Unique\tBest\tUnique\n";

for my $family ( sort {$a<=>$b} keys %HoA ) {
	my %f;
    foreach my $gene (@{ $HoA{$family} }){
		if (exists $HoA2{$gene}){
				#Count how many genes have annotations
				$isAnn++;
				#Count the occurence of every IPR term in the annotation of all members of a Cluster
			    for (@{$HoA2{$gene}}) {
	 		    	$f{$_}++;
					}
			    #foreach my $string (@{$HoA2{$gene}}) {
	 			#	next unless $seen{$string}++;
				#    $count++;
				#	}
			    }
		}
	#Calculate total Cluster size
	my $FamSize = scalar keys ($HoA{$family});
	print "Cluster$family\t$FamSize\t$isAnn\t";
	$isAnn=0;
	$unique=0;
	my $max=0;
	#calculate the highest frequency of IPR annotations
	foreach (sort {$b<=>$a} values %f) {
		$max=$_;
		last;}
	#Get all the unique annotations per Cluster
	foreach (sort keys %f) {
		if ($f{$_}==1){
			$unique++;
			push @uniqueIPR, $_;
		}
		if ($f{$_}==$max){$bestIPR = $_;}
	   }

	print "$max\t$unique\t$bestIPR\t";
	$bestIPR="NA";
	if (@uniqueIPR){
          print "$uniqueIPR[0]";
		}else{print "NA"};
	@uniqueIPR=();


#########################################
#	foreach (sort keys %f) {
	## To print IPR frequencies
		#print "$f{$_};";
#	   }

##To Print IPR terms
#print "\t";
#	foreach (sort keys %f) {
#  		print "$f{$_} ";
#  		print substr($_,0,9),"; ";
#
#  		}
#########################################

   	print "\n";
  		}
