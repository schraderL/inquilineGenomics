#
base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/p450/tree

##################################################################
# TAG: 1. align proper genes (peptides)
##################################################################

cd $base/

for i in $(ls ../*.p450/*.p450.GeMoMa.fa)
do
  species=$( echo $i|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|perl -pe 's/(.).*\_(...).*/$1$2/g')
  cat $i |sed "s/>/>$species-/g"
done > p450.pep
##################################################################
# TAG: 2. prepare fastas
##################################################################

for i in $(readlink -f ../*.p450/)
do
  awk '{print $0,FILENAME}' $i/*.p450.GeMoMa.details.tsv
done | sed "s|/corefac/cse/lukas/inqGen18/geneFamilies/specificFamilies/p450/||g"|cut -f 1 -d "/" > tmp.details.tsv

# create clean summary file with names p450
cat tmp.details.tsv |perl -pe 's/(.*)\t\s+(.).*\_(...).*/$2$3-$1/g' > p450.details.tsv

# subset to only full-length p450
awk '{if($2=="p450" && $3=="M*") print $1}' p450.details.tsv > p450.complete.lst
awk '{if($2=="p450") print $1}' p450.details.tsv > p450.incomplete.lst
faSomeRecords p450.pep p450.complete.lst p450.complete.pep
faSomeRecords p450.pep p450.incomplete.lst p450.incomplete.pep


##################################################################
# TAG: 3. align all p450s containing a p450 domain (peptides)
##################################################################
unset MAFFT_BINARIES
prank -o=p450.incomplete -d=p450.incomplete.pep &
cat p450.incomplete.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln
nice FastTreeMP tmp.aln > p450.incomplete.fast.tre
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0027731


##################################################################
# TAG: 3. align all COMPLETE p450s containing a p450 domain (peptides)
##################################################################
unset MAFFT_BINARIES
prank -o=p450.complete -d=p450.complete.pep &
cat p450.complete.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.cmp.aln
nice FastTreeMP tmp.cmp.aln > p450.complete.fast.tre
