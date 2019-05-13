############################################
# Prepare environment
############################################

base=/corefac/cse/lukas/inqGen18/ORphylogeny/
assignment=$base/data/SubfamAssignment.tab
cd $base
cat $assignment|sed 1d|perl -pe 's/(.*)\t(.*)\t/>$2_$1\n/g' > $assignment.fa
makeblastdb -dbtype prot -in $assignment.fa

#
# results/
# ├── aln
# │   ├── subfamilyAlignments
# │   ├── subfamilyFastas
# │   └── subfamilyTree
# └── species
#     ├── Acha
#     │   ├── fasta
#     │   └── lists
#     ├── Aech
#     │   ├── fasta
#     │   └── lists
#     ├── Ahey
#     │   ├── fasta
#     │   └── lists
#     ├── Ains
#     │   ├── fasta
#     │   └── lists
#     └── Parg
#         ├── fasta
#         └── lists

############################################
# retrieve OR annotations
############################################

cd $base/data/ORsets
for i in $(readlink -f /corefac/cse/lukas/inqGen18/ORs/*/final/finalSet/*OR.fa|grep -v "NCBI")
do
ln -s $i .
done

############################################
# Blast each OR annotation against assignment Table (pers. comm. McKenzie)
############################################

for i in $(readlink -f $base/data/ORsets/*.fa)
do
  species=$(echo $i|rev|cut -f 1 -d "/"|rev)
  spAbr=$(echo $species|cut -f 1 -d ".")
  nice blastp -query $i -num_threads 20 -db $assignment.fa -outfmt 6 -evalue 1e-50|sort -gk11,11 |sort -buk 1,1 > $base/$spAbr.OR.fa.bls
done


############################################
# Create fasta for each subfamily
############################################

for i in $(readlink -f $base/data/ORsets/*.fa)
do
  species=$(echo $i|rev|cut -f 1 -d "/"|rev)
  spAbr=$(echo $species|cut -f 1 -d ".")
  mkdir -p $base/results/species/$spAbr/
  mkdir -p $base/results/species/$spAbr/lists/
  mkdir -p $base/results/species/$spAbr/fasta/
  cd $base/results/species/$spAbr/lists/
  cut -f 1 -d "_" $base/$spAbr.OR.fa.bls|awk -F"\t" '{print>$2}'
  for subfam in 9E A B C D E F G H I J K L M N NV O Orco P Q R S T U V W X XA Z ZA
  do
    faSomeRecords $base/data/ORsets/$spAbr.OR.fa $subfam ../fasta/$spAbr.$subfam.fa
  done
done



############################################
# run alignments for each subfamily
############################################
mkdir $base/results/aln/
mkdir $base/results/aln/subfamilyFastas/
mkdir $base/results/aln/subfamilyAlignments/
mkdir $base/results/aln/subfamilyTree/

unset MAFFT_BINARIES
cd $base/results/aln/subfamilyAlignments/
#for subfam in 9E A B C D E F G H I J K L M N NV O Orco P Q R S T U V W X XA Z ZA
for subfam in 9E A B C D E F G H I J K L M N NV O Orco P Q R S T U V W X XA Z ZA
do
 cat $base/results/species/*/fasta/*.$subfam.fa > $base/results/aln/subfamilyFastas/$subfam.fa
 #time prank -o=$subfam.aln -d=$base/results/aln/subfamilyFastas/$subfam.fa
 time linsi $base/results/aln/subfamilyFastas/$subfam.fa > $subfam.mafft.aln
 #nice FastTreeMP $base/results/aln/subfamilyAlignments/$subfam.mafft.aln > $subfam.tre

mkdir alnRaw
mv *.aln alnRaw/

cd $base/results/aln/subfamilyAlignments/
# add subfam as name in mafft
for subfam in 9E A B C D E F G H I J K L M N NV O Orco P Q R S T U V W X XA Z ZA
do
 cat alnRaw/$subfam.mafft.aln|sed -re "s/>(.*) />\1 $subfam /g" > $subfam.mafft.aln
done

############################################
# merge alignments for all subfamilies
############################################
# https://mafft.cbrc.jp/alignment/software/merge.html
# could check out http://mergealign.appspot.com/
#https://mafft.cbrc.jp/alignment/software/makemergetable.rb
cd $base/results/aln/subfamilyTree/
alignments=$(readlink -f $base/results/aln/subfamilyAlignments/*.mafft.aln |tr "\n" " ")
cat $alignments > cat.ORsubfam.aln
ruby $base/scripts/makemergetable.rb $alignments > merge.table
mafft --localpair --maxiterate 100 --merge merge.table cat.ORsubfam.aln > merge.ORsubfam.aln

############################################
# Run tree inference with fastTree
############################################

nice FastTreeMP merge.ORsubfam.aln > ant.OR.tre

############################################
# Run tree inference with Raxml
############################################
cd $base/results/aln/subfamilyTree/raxml

nice raxmlHPC-PTHREADS -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE -s ../merge.ORsubfam.aln -n r1 -T 20
#sequence duplicates and undetermined columns removed is printed to file ../merge.ORsubfam.aln.reduced
