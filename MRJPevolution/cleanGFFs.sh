#! /usr/bin/bash
species=$1
curatedGff=$2
gFa=$3


grep -shoP  "ID\=.{36}" $curatedGff |cut -f 2 -d "="|sort|uniq |awk '{printf($0"\tel%03d\n", NR)}'|tr " " "\t" > ids.tsv
awk -F'\t' 'FNR==NR{a[$1]=$2;next} {for (i in a)gsub(i, a[i]);print}' ids.tsv $curatedGff  > tmp.gff
awk '{if($3=="gene" || $3=="pseudogene" ) print $0}' tmp.gff |perl -pe 's/.*ID=(.*?);.*Name=(.*?);.*/$1\t$2/g' > names1.tsv
awk -F'\t' 'FNR==NR{a[$1]=$2;next} {for (i in a)gsub(i, a[i]);print}' names1.tsv tmp.gff  > tmp2.gff
awk '{if($3=="mRNA"||$3=="transcript") print $0}' tmp.gff |perl -pe 's/.*ID=(.*?);.*Name=(.*?);.*/$1\t$2/g' > names2.tsv
awk -F'\t' 'FNR==NR{a[$1]=$2;next} {for (i in a)gsub(i, a[i]);print}' names2.tsv tmp2.gff  > tmp3.gff
cat tmp3.gff|perl -pe 's/(.*)date_last_modified=.*?;(.*)/$1$2/g'|perl -pe 's/(.*)owner=.*?;(.*)/$1$2/g'|perl -pe 's/(.*)date_creation=.*?/$1/g' > tmp4.gff

# clean gff
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp4.gff

# split pseudogenes and regular genes
awk '{if($3=="pseudogene" ) print $0}' tmp4.gff|perl -pe 's/.*ID=(.*?);.*/$1/g' |grep -f - tmp4.gff > tmp.pseudogenes.gff
awk '{if($3=="pseudogene" ) print $0}' tmp4.gff|perl -pe 's/.*ID=(.*?);.*/$1/g' |grep -f - -v tmp4.gff > tmp.genes.gff

# create clean gene file
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.genes.gff
/usr/local/home/lschrader/software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -addintrons tmp.genes_clean.gff > $species.MRJPs.gff3

# create clean pseudogene file
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.pseudogenes.gff
/usr/local/home/lschrader/software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -addintrons tmp.pseudogenes_clean.gff > $species.pseudogenes.MRJPs.gff3
rm tmp.gff tmp2.gff tmp3.gff names1.tsv names2.tsv ids.tsv tmp4.gff tmp4_clean.gff tmp.genes_clean.gff tmp.pseudogenes_clean.gff tmp.pseudogenes.gff tmp.genes.gff
/usr/local/home/lschrader/software/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl $species.MRJPs.gff3 $gFa |sed "s/ \../\t/g"> $species.MRJPs.pep
gffread $species.MRJPs.gff3 -x $species.MRJPs.cds.fa -g $gFa

awk '{if($3=="pseudogene" ) print $0}' $species.pseudogenes.MRJPs.gff3|perl -pe 's/.*ID=(.*?);.*/$1/g' > pseudogenes.lst
while read gene; do
  grep $gene $species.pseudogenes.MRJPs.gff3|grep exon > tmp.$gene.gff
  bedtools getfasta -fi $gFa -bed tmp.$gene.gff -fo - | grep -v ">" |perl -pe 'chomp'|perl -lpe "s/^/>${gene}-pseudogene\n/g"
done <pseudogenes.lst > $species.pseudogenes.MRJPs.exon.fa
rm tmp.*.gff pseudogenes.lst
