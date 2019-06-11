use strict;
use warnings;

die "Usage: $0 hash_file gff_file\n" if @ARGV < 2;

my ($hashfile, $gfffile) = @ARGV;

open(my $hfile, '<', $hashfile) or die "Cannot open $hashfile: $!";
my %names;
while (my $line = <$hfile>) {
		 $line =~ /^#/ and next;
		chomp($line);
		my @tab = split /\t/, $line;
		my $oldname = $tab[6];
		my $newname = $tab[0];
    $names{$oldname} = $newname;
}
close $hfile;

my $regex = join '|', sort { length $b <=> length $a } keys %names;
$regex = qr/$regex/;

open(my $gfile, '<', $gfffile) or die "Cannot open $gfffile: $!";


while (my $line = <$gfile>) {
    chomp($line);
    $line =~ s/($regex)/$names{$1}/g;
    print $line, "\n";
}

close $gfile;
