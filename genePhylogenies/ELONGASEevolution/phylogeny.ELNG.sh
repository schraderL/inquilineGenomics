#
base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/ELG/tree

##################################################################
# TAG: 1. align proper genes (peptides)
##################################################################

cd $base/


for i in $(ls ../*.ELG/*.ELG.GeMoMa.fa)
do
  species=$( echo $i|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|perl -pe 's/(.).*\_(...).*/$1$2/g')
  cat $i |sed "s/>/>$species-/g"
done > ELGs.pep
unset MAFFT_BINARIES
prank -o=ELG.pep.aln -d=ELGs.pep
cat ELG.pep.aln.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln

antAln=$(readlink -f ELG.pep.aln.best.fas)

FastTree tmp.aln > ELG.fast.tre

##################################################################
# TAG: 2. Complete genes only
##################################################################


for i in $(readlink -f ../*.ELG/)
do
  awk '{print $0,FILENAME}' $i/*.ELG.GeMoMa.details.tsv
done | sed "s|/corefac/cse/lukas/inqGen18/geneFamilies/specificFamilies/ELG/||g"|cut -f 1 -d "/" > tmp.details.tsv

# create clean summary file with names ELGs
cat tmp.details.tsv |perl -pe 's/(.*)\t\s+(.).*\_(...).*/$2$3-$1/g' > ELG.details.tsv

# subset to only full-length ELGs
awk '{if($2=="ELO" && $3=="M*") print $1}' ELG.details.tsv > ELG.complete.lst
faSomeRecords ELGs.pep ELG.complete.lst ELG.complete.pep

unset MAFFT_BINARIES
prank -o=ELG.complete -d=ELG.complete.pep
cat ELG.complete.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln

antAln=$(readlink -f ELG.complete.best.fas)

FastTree tmp.aln > ELG.complete.fast.tre

##################################################################
# TAG: 2. Genes with ELO domain only
##################################################################

# subset to only full-length ELGs
awk '{if($2=="ELO") print $1}' ELG.details.tsv > ELG.ELO.lst
faSomeRecords ELGs.pep ELG.ELO.lst ELG.ELO.pep

unset MAFFT_BINARIES
prank -o=ELG.ELO -d=ELG.ELO.pep
cat ELG.complete.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln

antAln=$(readlink -f ELG.complete.best.fas)

FastTree tmp.aln > ELG.complete.fast.tre
