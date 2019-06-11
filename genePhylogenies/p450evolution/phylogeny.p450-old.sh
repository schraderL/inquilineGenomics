#
base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/p450/tree

##################################################################
# TAG: 1. align proper genes (peptides)
##################################################################

cd $base/

for i in $(ls ../*.p450/*.p450.GeMoMa.fa);
do
     species=$( echo $i|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|perl -pe 's/(.).*\_(...).*/$1$2/g');
     cat $i |sed "s/>/>$species-/g";
done > p450.pep


unset MAFFT_BINARIES
prank -o=GRs.pep.aln -d=GRs.pep

cat GRs.pep.aln.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln

antAln=$(readlink -f GRs.pep.aln.best.fas)

FastTree tmp.aln > GRs.best.tre

##################################################################
# TAG: 2. raxML phylogeny
##################################################################

nice raxmlHPC-PTHREADS -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE -s $antAln -n r1 -T 20
#/corefac/cse/lukas/inqGen18/geneFamilies/specificFamilies/MRJPs/aln/mergeAln/RAxML_bipartitionsBranchLabels.r1


##################################################################
# TAG: 3. align only complete GRs (peptides)
##################################################################

cd $base/

for i in $(readlink -f ../*.p450/);
do
  awk '{print $0,FILENAME}' $i/*.p450.GeMoMa.details.tsv;
done | sed "s|/corefac/cse/lukas/inqGen18/geneFamilies/specificFamilies/p450/||g"|cut -f 1 -d "/" > tmp.details.tsv

 # create clean summary file with names p450
cat tmp.details.tsv |perl -pe 's/(.*)\t\s+(.).*\_(...).*/$2$3-$1/g' > GRs.details.tsv

# subset to only full-length GRs
awk '{if($2=="7_TM7" && $3=="M*") print $1}' GRs.details.tsv > GRs.complete.lst
faSomeRecords GRs.pep GRs.complete.lst GRs.complete.pep

unset MAFFT_BINARIES
prank -o=GRs.complete -d=GRs.complete.pep
cat GRs.complete.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln

antAln=$(readlink -f GRs.complete.best.fas)

FastTree tmp.aln > GRs.complete.fast.tre
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0027731
##################################################################
# TAG: 3. raxML phylogeny
##################################################################

nice raxmlHPC-PTHREADS -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE -s $antAln -n r1 -T 20
#/corefac/cse/lukas/inqGen18/geneFamilies/specificFamilies/MRJPs/aln/mergeAln/RAxML_bipartitionsBranchLabels.r1
