#!/usr/bin/perl -w
use strict;
use warnings;

my $usage = "\nThis script can be used to rename the IDs generated from webapollo.\n\n

#usage: perl cleanWebapolloGFF.pl <manual.curation.gff3>\n\n";
unless (scalar @ARGV == 1){die $usage}


my $gff = shift; #provide hsp file
my (%id2par,%exons);
my $exonID;
#my $counter=1;
open(FILE1, $gff) or die "Cannot open gff\n";
while (<FILE1>)
	{
	chomp;

	if (/\#/){print $_,"\n"; next;};

	#print "AAAAAAAAAAAA";
	my @tab = split('\s',$_);
	if ($tab[2] =~ /CDS/){next;};
	my ($name) = $tab[8] =~ /Name=(.*?)(-|$)/;

	my ($id) = $tab[8] =~ /ID=(.*?)(;|$)/;
	$id2par{$id} = $name;
	my ($date) = $tab[8] =~ /date_creation=(.*?)(;|$)/;
	my ($Parent) = $tab[8] =~ /Parent=(.*?)(;|$)/;
	my ($owner) = $tab[8] =~ /owner=(.*?)(;|$)/;
	my ($dlm) = $tab[8] =~ /date_last_modified=(.*?)(;|$)/;

	#print $_;
	pop @tab;
	print join("\t",@tab);
	print "\t";

	 if ($tab[2] eq "gene"){
		 if (defined $name){print "ID=",$name,";";}
		 if (defined $name){print "Name=",$name;}
		 		 #if (defined $owner){print $owner,"\t";}
		 #if (defined $dlm){print $dlm;}
	 }
	 if (defined $Parent){
		 if (exists $exons{$id2par{$Parent}."-".$tab[2]}){$exons{$id2par{$Parent}."-".$tab[2]}++;}
		 else{$exons{$id2par{$Parent}."-".$tab[2]}=1;}
		 $name=$id2par{$Parent}."-".$tab[2]."-".$exons{$id2par{$Parent}."-".$tab[2]};
		 print "ID=",$name,";";
		 print "Name=",$name,";";
		 $id2par{$id} = $name;
		 print "Parent=",$id2par{$Parent};
		 if ($tab[2] =~ /exon/){
			 my $CDSline = "\n".join("\t",@tab)."\t"."ID=".$name.";"."Name=".$name.";"."Parent=".$id2par{$Parent};
			 my $cl2 = $CDSline =~ s/exon/CDS/g;
			 print $CDSline;

		 	};
	 	}
	 #if (defined $date){ print $date,"\t";}
	 #if (defined $Parent){}
	# if (defined $id){print $id,"\t";}
	# if (defined $owner){print $owner,"\t";}
	# if (defined $dlm){print $dlm;}
	 print "\n";
	}
	#else{print $_;}


close(FILE1);
