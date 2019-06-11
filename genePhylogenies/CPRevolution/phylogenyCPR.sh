#
base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/CPR/tree

##################################################################
# TAG: 1. align proper genes (peptides)
##################################################################

cd $base/

for i in $(ls ../*.CPR/*.CPR.GeMoMa.fa)
do
  species=$( echo $i|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|perl -pe 's/(.).*\_(...).*/$1$2/g')
  cat $i |sed "s/>/>$species-/g"
done > CPR.pep
unset MAFFT_BINARIES
prank -o=CPR.pep.aln -d=CPR.pep

cat CPR.pep.aln.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln

antAln=$(readlink -f CPR.pep.aln.best.fas)

FastTree tmp.aln > CPR.best.tre


##################################################################
# TAG: 2. prepare fastas
##################################################################

cd $base/

for i in $(readlink -f ../*.CPR/)
do
  awk '{print $0,FILENAME}' $i/*.CPR.GeMoMa.details.tsv
done | sed "s|/corefac/cse/lukas/inqGen18/geneFamilies/specificFamilies/CPR/||g"|cut -f 1 -d "/" > tmp.details.tsv

# create clean summary file with names CPR
cat tmp.details.tsv |perl -pe 's/(.*)\t\s+(.).*\_(...).*/$2$3-$1/g' > CPR.details.tsv

# subset to only full-length CPR
awk '{if($2=="Chitin_bind_4" && $3=="M*") print $1}' CPR.details.tsv > CPR.complete.lst
awk '{if($2=="Chitin_bind_4") print $1}' CPR.details.tsv > CPR.incomplete.lst
faSomeRecords CPR.pep CPR.complete.lst CPR.complete.pep
faSomeRecords CPR.pep CPR.incomplete.lst CPR.incomplete.pep


##################################################################
# TAG: 3. align all CPRs containing a Chitin_bind_4 domain (peptides)
##################################################################
unset MAFFT_BINARIES
prank -o=CPR.incomplete -d=CPR.incomplete.pep
cat CPR.incomplete.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln
FastTree tmp.aln > CPR.incomplete.fast.tre
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0027731


##################################################################
# TAG: 3. align all COMPLETE CPRs containing a Chitin_bind_4 domain (peptides)
##################################################################
unset MAFFT_BINARIES
prank -o=CPR.complete -d=CPR.complete.pep
cat CPR.complete.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.cmp.aln
nice FastTreeMP tmp.cmp.aln > CPR.complete.fast.tre
