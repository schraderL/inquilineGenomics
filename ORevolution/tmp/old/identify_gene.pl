#!/usr/bin/perl -w
use strict;

my $config = $ARGV[0];
my $blast_out = $ARGV[1];

my %config;
open IN, "< $config" or die "Can't open $config!";
while (<IN>)	{
    next if (/^#/);
	chomp;
	my @config = split /=/;
	$config{$config[0]}=$config[1];
}
close (IN);

my %genome_length = &genome_length($config{DATABASE});

my %refseq = &refseq($config{INPUT_GENEFAM});

unless (-e "$blast_out")	{
	system "tblastn -db $config{DATABASE} -query $config{INPUT_GENEFAM} -out $blast_out -num_threads $config{NUMPROC} -evalue $config{EXPECT_VAL} -outfmt 6";
}

my %hsps; my %contigs; my $hsp = 0;
open IN, "< $blast_out" or die "Can't open $blast_out!";
while (<IN>)	{
	my @hsp = split /\s+/;
	$hsps{$hsp}->{REF_GENE} = $hsp[0];
	if ($hsp[8] < $hsp[9])	{
		$hsps{$hsp}->{CONTIG} = $hsp[1]." +";
	}	else	{
		$hsps{$hsp}->{CONTIG} = $hsp[1]." -";
	}
	$hsps{$hsp}->{Q_START} = $hsp[6];
	$hsps{$hsp}->{Q_END} = $hsp[7];
	$hsps{$hsp}->{S_START} = $hsp[8];
	$hsps{$hsp}->{S_END} = $hsp[9];
	$hsps{$hsp}->{E_VALUE} = $hsp[10];
	$hsps{$hsp}->{BIT_SCORE} = $hsp[11];
	$contigs{$hsps{$hsp}->{CONTIG}}->{$hsp} = 1;
	$hsp++;
}
close (IN);

&filter_hsp(\%hsps, $config{PERC_CUTOFF}, $config{TOP_PERC_IN}, \%contigs);

my %gene_builds;
foreach my $contig (keys %contigs)	{
	&merge_hsp($contig, $contigs{$contig}, \%hsps, \%gene_builds, $config{MIN_HSP_DIST}, $config{MAX_OVERLAP_HSP});
	&filter_gene_build(\%gene_builds, $contig, $config{MIN_OVERLAP_GENE});
	&check_gene_build(\%gene_builds, $contig, $config{MIN_GENE_DIST}, $config{MAX_OVERLAP_HSP});
	(my $contig_name = $contig) =~ s/[\s\+\-]//g;
	&report_gene_build(\%gene_builds, $contig, $config{DATABASE}, $config{PGENE_EXT}, $genome_length{$contig_name}, \%refseq);
}

exit;

sub genome_length	{
	my $database = $_[0];

	my %genome_length; my $contig;
	open IN, "< $database" or die "Can't open $database!";
	while (<IN>)	{
		if (/>/)	{
			($contig = $_) =~ s/[\s\>]//g;
		}	else	{
			s/\s//g;
			$genome_length{$contig} .= $_;
		}
	}
	close (IN);
	
	foreach my $contig (keys %genome_length)	{
		$genome_length{$contig} = length($genome_length{$contig});
	}
	
	return %genome_length;
}

sub refseq	{
	my $input_genefam = $_[0];

	my %refseq; my $refseq;
	open IN, "< $input_genefam" or die "Can't open $input_genefam";
	while (<IN>)	{
		if (/>/)	{
			($refseq = $_) =~ s/[\s\>]//g;
		}
		$refseq{$refseq} .= $_;
	}
	close (IN);
	
	return %refseq;
}

sub filter_hsp	{
	(my $hsps, my $perc_cutoff, my $top_perc_in, my $contigs) = @_;
	
	my @bit_score;
	foreach my $hsp (keys %{$hsps}) {
		push @bit_score, $hsps->{$hsp}->{BIT_SCORE};
	}
	@bit_score = sort {$b <=> $a} @bit_score;
	my $bit_score;
	my $in = int((scalar @bit_score)*$top_perc_in*0.01);
	foreach my $hsp (0..$in)	{
		$bit_score += $bit_score[$hsp];
	}
	$bit_score *= $perc_cutoff*0.01/$in;
	foreach my $hsp (keys %{$hsps}) {
		if ($hsps->{$hsp}->{BIT_SCORE} < $bit_score) {
			delete $contigs->{$hsps->{$hsp}->{CONTIG}}->{$hsp};
			delete $hsps->{$hsp};
		}
	}

    return;
}

sub merge_hsp	{
	(my $contig, my $contigs, my $hsps, my $gene_builds, my $min_hsp_dist, my $max_overlap_hsp) = @_;
	
	my %ref_gene;
	foreach my $hsp (keys %{$contigs})	{
		$ref_gene{$hsps->{$hsp}->{REF_GENE}}->{$hsp} = 1;
	}
	foreach my $ref_gene (keys %ref_gene)	{
		my @ordered_ref_gene;
		if ($contig =~ /\+$/)	{
			@ordered_ref_gene = sort {$hsps->{$a}->{S_START} <=> $hsps->{$b}->{S_START}} keys %{$ref_gene{$ref_gene}};
		}	else	{
			@ordered_ref_gene = sort {$hsps->{$b}->{S_START} <=> $hsps->{$a}->{S_START}} keys %{$ref_gene{$ref_gene}};
		}
		my %current_gene;
		foreach my $hsp (@ordered_ref_gene)	{
			if (defined($current_gene{BIT_SCORE}))	{
				my $s_flag = &s_flag($contig, $hsps->{$hsp}->{S_START}, $current_gene{S_END}, 0);
				if ($s_flag <= $min_hsp_dist)	{
					if ($current_gene{Q_END} <= $hsps->{$hsp}->{Q_START} + $max_overlap_hsp) {
						$current_gene{Q_END} = $hsps->{$hsp}->{Q_END};
						$current_gene{S_END} = $hsps->{$hsp}->{S_END};
						$current_gene{BIT_SCORE} += $hsps->{$hsp}->{BIT_SCORE};
						next;
					}
				}
				my $gene_build_name = $current_gene{S_START}."\t".$current_gene{S_END}."\t".$current_gene{REF_GENE}."\t".$current_gene{Q_START}."\t".$current_gene{Q_END}."\t".$current_gene{BIT_SCORE};
				%{$gene_builds->{$contig}->{$gene_build_name}} = %current_gene;
				$gene_builds->{COUNT}++;
			}	
			%current_gene = %{$hsps->{$hsp}};
		}
		my $gene_build_name = $current_gene{S_START}."\t".$current_gene{S_END}."\t".$current_gene{REF_GENE}."\t".$current_gene{Q_START}."\t".$current_gene{Q_END}."\t".$current_gene{BIT_SCORE};
		%{$gene_builds->{$contig}->{$gene_build_name}} = %current_gene;
		$gene_builds->{COUNT}++;
	}

    return;
}

sub s_flag	{
	(my $contig, my $gene1, my $gene2, my $adjustment) = @_;

	my $s_flag;
	if ($contig =~ /\+$/)	{
		$s_flag = $gene1 - $gene2 + $adjustment;
	}	else	{
		$s_flag = $gene2 - $gene1 + $adjustment;
	}

	return $s_flag;
} 

sub filter_gene_build	{
	(my $gene_builds, my $contig, my $min_overlap_gene) = @_;
	
	my @gene_build;
	if ($contig =~ /\+$/)	{
		@gene_build = sort {$gene_builds->{$contig}->{$a}->{S_START} <=> $gene_builds->{$contig}->{$b}->{S_START}} keys %{$gene_builds->{$contig}};
	}	else	{
		@gene_build = sort {$gene_builds->{$contig}->{$b}->{S_START} <=> $gene_builds->{$contig}->{$a}->{S_START}} keys %{$gene_builds->{$contig}};
	}
	my $current_gene;
	foreach my $gene_build (@gene_build)	{
		if (defined($current_gene))	{
			my $current_gene_q_length = abs($gene_builds->{$contig}->{$current_gene}->{Q_END} - $gene_builds->{$contig}->{$current_gene}->{Q_START}) + 1;
			my $current_gene_s_length = abs($gene_builds->{$contig}->{$current_gene}->{S_END} - $gene_builds->{$contig}->{$current_gene}->{S_START}) + 1;
			my $gene_build_q_length = abs($gene_builds->{$contig}->{$gene_build}->{Q_END} - $gene_builds->{$contig}->{$gene_build}->{Q_START}) + 1;
			my $gene_build_s_length = abs($gene_builds->{$contig}->{$gene_build}->{S_END} - $gene_builds->{$contig}->{$gene_build}->{S_START}) + 1;
			my $s_flag = &s_flag($contig, $gene_builds->{$contig}->{$current_gene}->{S_END}, $gene_builds->{$contig}->{$gene_build}->{S_START}, 1);
			if ($s_flag >= ($current_gene_s_length + $gene_build_s_length)*$min_overlap_gene/200)	{
				if ($gene_builds->{$contig}->{$gene_build}->{BIT_SCORE} > $gene_builds->{$contig}->{$current_gene}->{BIT_SCORE})	{
					delete $gene_builds->{$contig}->{$current_gene};
					$current_gene = $gene_build;
				}	elsif ($gene_builds->{$contig}->{$gene_build}->{BIT_SCORE} < $gene_builds->{$contig}->{$current_gene}->{BIT_SCORE})	{
					delete $gene_builds->{$contig}->{$gene_build};
				}	else	{
					if ($current_gene_q_length >= $gene_build_q_length)	{
						delete $gene_builds->{$contig}->{$gene_build};
					}	else	{
						delete $gene_builds->{$contig}->{$current_gene};
						$current_gene = $gene_build;
					}
				}
				$gene_builds->{COUNT}--;
				next;
			}
		}
		$current_gene = $gene_build;
	}

    return;
}

sub check_gene_build	{
	(my $gene_builds, my $contig, my $min_gene_dist, my $max_overlap_hsp) = @_;

	my @gene_build;
	if ($contig =~ /\+$/)	{
		@gene_build = sort {$gene_builds->{$contig}->{$a}->{S_START} <=> $gene_builds->{$contig}->{$b}->{S_START}} keys %{$gene_builds->{$contig}};
	}	else	{
		@gene_build = sort {$gene_builds->{$contig}->{$b}->{S_START} <=> $gene_builds->{$contig}->{$a}->{S_START}} keys %{$gene_builds->{$contig}};
	}
	my $current_gene;
	foreach my $gene_build (@gene_build)	{
		if (defined($current_gene))	{
			my $s_flag = &s_flag($contig, $gene_builds->{$contig}->{$gene_build}->{S_START}, $gene_builds->{$contig}->{$current_gene}->{S_END}, 0);
			if ($s_flag <= $min_gene_dist)	{
				if ($gene_builds->{$contig}->{$current_gene}->{Q_END} <= $gene_builds->{$contig}->{$gene_build}->{Q_START} + $max_overlap_hsp) {
					$gene_builds->{$contig}->{$current_gene}->{Q_END} = $gene_builds->{$contig}->{$gene_build}->{Q_END};
					$gene_builds->{$contig}->{$current_gene}->{S_END} = $gene_builds->{$contig}->{$gene_build}->{S_END};
					if ($gene_builds->{$contig}->{$current_gene}->{BIT_SCORE} < $gene_builds->{$contig}->{$gene_build}->{BIT_SCORE})	{
						$gene_builds->{$contig}->{$current_gene}->{REF_GENE} = $gene_builds->{$contig}->{$gene_build}->{REF_GENE}; 
					}
					$gene_builds->{$contig}->{$current_gene}->{BIT_SCORE} += $gene_builds->{$contig}->{$gene_build}->{BIT_SCORE};
					delete $gene_builds->{$contig}->{$gene_build};
					$gene_builds->{COUNT}--;
					next;
				}
			}
			my $gene_build_name = $gene_builds->{$contig}->{$current_gene}->{S_START}."\t".$gene_builds->{$contig}->{$current_gene}->{S_END}."\t".$gene_builds->{$contig}->{$current_gene}->{REF_GENE}."\t".$gene_builds->{$contig}->{$current_gene}->{Q_START}."\t".$gene_builds->{$contig}->{$current_gene}->{Q_END}."\t".$gene_builds->{$contig}->{$current_gene}->{BIT_SCORE};
			unless ($gene_build_name eq $current_gene)	{
				%{$gene_builds->{$contig}->{$gene_build_name}} = %{$gene_builds->{$contig}->{$current_gene}};
				delete $gene_builds->{$contig}->{$current_gene};
			}
		}
		$current_gene = $gene_build;
	}
	my $gene_build_name = $gene_builds->{$contig}->{$current_gene}->{S_START}."\t".$gene_builds->{$contig}->{$current_gene}->{S_END}."\t".$gene_builds->{$contig}->{$current_gene}->{REF_GENE}."\t".$gene_builds->{$contig}->{$current_gene}->{Q_START}."\t".$gene_builds->{$contig}->{$current_gene}->{Q_END}."\t".$gene_builds->{$contig}->{$current_gene}->{BIT_SCORE};
	unless ($gene_build_name eq $current_gene)      {
		%{$gene_builds->{$contig}->{$gene_build_name}} = %{$gene_builds->{$contig}->{$current_gene}};
		delete $gene_builds->{$contig}->{$current_gene};
	}

    return;
}

sub report_gene_build	{
	(my $gene_builds, my $contig, my $database, my $pgene_ext, my $contig_length, my $refseqs) = @_;
	
	my @gene_build;
	if ($contig =~ /\+$/)   {
		@gene_build = sort {$gene_builds->{$contig}->{$a}->{S_START} <=> $gene_builds->{$contig}->{$b}->{S_START}} keys %{$gene_builds->{$contig}};
	}       else    {
		@gene_build = sort {$gene_builds->{$contig}->{$b}->{S_START} <=> $gene_builds->{$contig}->{$a}->{S_START}} keys %{$gene_builds->{$contig}};
	}
	my $strand;
	if ($contig =~ /\+$/)	{
		$strand = "plus";
	}	else	{
		$strand = "minus";
	}
	my @range;
	my $left; my $right;
	my $current_gene;
	foreach my $gene_build (@gene_build)	{
		print "$gene_build\n";
		if (defined($current_gene))	{
			if (defined($range[1]))	{
				$range[0] = $range[1];
			}	else	{
				$range[0] = $pgene_ext;
			}
			my $s_flag = &s_flag($contig, $gene_builds->{$contig}->{$gene_build}->{S_START}, $gene_builds->{$contig}->{$current_gene}->{S_END}, 0);
			if ($s_flag <= $pgene_ext*2)	{
				$range[1] = int($s_flag/2);
			}	else	{
				$range[1] = $pgene_ext;
			}
			if ($contig =~ /\+$/)	{
				$left = $gene_builds->{$contig}->{$current_gene}->{S_START} - $range[0];
				if ($left < 1)	{
					$left = 1;
				}
				$right = $gene_builds->{$contig}->{$current_gene}->{S_END} + $range[1];
				if ($right > $contig_length)	{
					$right = $contig_length;
				}
			}	else	{
				$left = $gene_builds->{$contig}->{$current_gene}->{S_START} + $range[0];
				if ($left > $contig_length)	{
					$left = $contig_length;
				}
				$right = $gene_builds->{$contig}->{$current_gene}->{S_END} - $range[1];
				if ($right < 1)	{
					$right = 1;
				}
			}
            print "$contig\t$current_gene\n";
		}
		$current_gene = $gene_build;
	}
	if (defined($range[1]))	{
		$range[0] = $range[1];
	}	else	{
		$range[0] = $pgene_ext;
	}
	$range[1] = $pgene_ext;
	if ($contig =~ /\+$/)	{
		$left = $gene_builds->{$contig}->{$current_gene}->{S_START} - $range[0];
		if ($left < 1)	{
			$left = 1;
		}
		$right = $gene_builds->{$contig}->{$current_gene}->{S_END} + $pgene_ext;
		if ($right > $contig_length)	{
			$right = $contig_length;
		}
	}	else	{
		$left = $gene_builds->{$contig}->{$current_gene}->{S_START} + $range[0];
		if ($left > $contig_length)	{
			$left = $contig_length;
		}
		$right = $gene_builds->{$contig}->{$current_gene}->{S_END} - $pgene_ext;
		if ($right < 1)	{
			$right = 1;
		}
	}
    print "$contig\t$current_gene\n";

    return;
}
