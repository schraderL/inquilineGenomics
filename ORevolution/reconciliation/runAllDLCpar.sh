#!/usr/bin/bash
export base=/usr/local/home/lschrader/data/inqGen18/ORphylogeny/reconciledClades
for clade in $(ls $base/cladeData/*.tsv|perl -pe 's/.*\/(.*)\.tsv/$1/g')
do
  echo $clade

  cd $base
  mkdir $base/$clade
  cd $base/$clade
  cut -f 1 $base/cladeData/$clade.tsv|seqkit grep -rf - $base/data/ORs.cds.fa |perl -pe 's/-mRNA.*//g'> $clade.cds.fa
  #grep ">" $clade.cds.fa -c

  ## Check if codons are in order
  crooked=$(seqkit fx2tab $clade.cds.fa |awk -F $'\t' '{OFS = FS} {if (length($2)%3!=0) print $1,$2}'|wc -l)

  if [ "$crooked" != 0 ];
  then
    seqkit fx2tab $clade.cds.fa |awk -F $'\t' '{print $0,length($2)}'|awk -F $'\t' '{OFS = FS} {if ($3%3!=0) print $1,$2}'|seqkit tab2fx > $clade.tmp.fa
    TransDecoder.LongOrfs -t $clade.tmp.fa
    cat $clade.tmp.fa.transdecoder_dir/longest_orfs.cds |perl -pe 's/\>.* (.*)/\>$1/g'> $clade.trimmed.fa
    seqkit fx2tab $clade.cds.fa |awk -F $'\t' '{print $0,length($2)}'|awk -F $'\t' '{OFS = FS} {if ($3%3==0) print $1,$2}'|seqkit tab2fx > $clade.good.fa
    cat $clade.good.fa $clade.trimmed.fa > $clade.cds.trimmed.fa
    cut -f 1 -d " " $clade.cds.trimmed.fa|cut -f 1 -d ":" > $clade.cds.edit2.fa
  else
  echo "all good"
    cut -f 1 -d " " $clade.cds.fa|cut -f 1 -d ":" > $clade.cds.edit2.fa
  fi

  prank -t=$base/cladeData/$clade.tre -d=$clade.cds.edit2.fa -o=$clade.cds.out -codon -iterate=3

  ## Compute gene tree phylogeny with FastTree
  #bash
  FastTree -nt $clade.cds.out.best.fas > $clade.cds.out.FastTree.fas

  #midpoint reroot tree
  # outgroup.lst
  cat $clade.cds.out.FastTree.fas|gotree reroot midpoint -i -  > $clade.cds.out.FastTree.rooted.fas

  ## Create #smap# file for dlcpar
  cat $base/cladeData/$clade.tsv |sed 1d|cut -f 1,6 |perl -pe 's/-.*\t/\*\t/g' |sort|uniq|awk -F"\t" '{OFS=FS}{print $1,toupper($2)}' > $clade.smap

  ## Run DLCpar
  ## reformat trees so dlcpar search does accept them

  mkdir $base/$clade/DLC
  cd $base/$clade/DLC

  ln -s ../$clade.cds.out.FastTree.rooted.fas .
  cat $clade.cds.out.FastTree.rooted.fas> tmp2.tre
  python $base/scripts/nodelabel.py tmp2.tre|perl -pe 's/\;/root\;/g' > $clade.FT.tre
  nice dlcpar search -S ../$clade.smap -s $base/data/species.tree $clade.FT.tre
  nice dlcpar events --3t -S ../$clade.smap -s $base/data/species.tree $clade.FT.tre.dlcsearch.coal.tree > $clade.recon.tsv
  # mowgli
  mkdir $base/$clade/mowgli
  cd $base/$clade/mowgli
  cat ../DLC/$clade.FT.tre|perl -pe 's/Or-(...).*?:/_$1:/g'|tr [a-z] [A-Z] > tmp.tre

  ~/software/Mowgli-v2/Mowgli_linux_i386_64_static -s $base/data/species.tree -g tmp.tre  -t=-1 -b
  cd $base

done
