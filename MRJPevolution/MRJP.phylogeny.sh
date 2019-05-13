##################################################################
# TAG: 1. prepare environment
##################################################################
scripts=~/data/inqGen18/geneFamilies/specificFamilies/scripts
base=~/data/inqGen18/geneFamilies/specificFamilies/MRJPs

cd $base/webapolloGFFs

##################################################################
# TAG: 2. clean GFFs
##################################################################

# make folders for each species
cd $base/cleanGFFs/Acep
bash $scripts/cleanGFFs.sh Acep $base/webapolloGFFs/cleaned_Acep_MRJP_apollo_Annotations.gff3 /corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_cephalotes/genome/Atta_cephalotes.v2.0.fa
cd $base/cleanGFFs/Acol
bash $scripts/cleanGFFs.sh Acol $base/webapolloGFFs/cleaned_Acol_MRJP_apollo_Annotations.gff3 /corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_colombica/genome/Atta_colombica.v2.0.fa
cd $base/cleanGFFs/Parg
bash $scripts/cleanGFFs.sh Parg $base/webapolloGFFs/cleaned_Parg_MRJP_apollo_Annotations.gff3 /corefac/cse/lukas/inqGen18/inquilines_v2.1/Pseudoatta_argentina.2.1/genome/Pseudoatta_argentina.v2.1.fa
cd $base/cleanGFFs/Acha
bash $scripts/cleanGFFs.sh Acha $base/webapolloGFFs/cleaned_Acha_MRJP_apollo_Annotations.gff3 /corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_charruanus.2.1/genome/Acromyrmex_charruanus.v2.1.fa
cd $base/cleanGFFs/Ahey
bash $scripts/cleanGFFs.sh Ahey $base/webapolloGFFs/cleaned_Ahey_MRJP_apollo_Annotations.gff3 /corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_heyeri.2.1/genome/Acromyrmex_heyeri.v2.1.fa
cd $base/cleanGFFs/Ains
bash $scripts/cleanGFFs.sh Ains $base/webapolloGFFs/cleaned_Ains_MRJP_apollo_Annotations.gff3 /corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_insinuator.2.1/genome/Acromyrmex_insinuator.v2.1.fa
cd $base/cleanGFFs/Aech
bash $scripts/cleanGFFs.sh Aech $base/webapolloGFFs/cleaned_Aech_MRJP_apollo_Annotations.gff3 /corefac/cse/lukas/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/genome/Acromyrmex_echinatior.v2.0.fa

##################################################################
# TAG: 3. align proper genes (peptides)
##################################################################

cd $base/aln
#ls  ../cleanGFFs/*/*pep
ln -s ../cleanGFFs/*/*pep .
for i in $(ls *.MRJPs.pep)
do
  species=$(echo $i|cut -f 1 -d ".")
  cat $i |sed "s/>/>$species-/g"
done > MRJPs.pep
unset MAFFT_BINARIES
prank -o=MRJPs.pep.aln -d=MRJPs.pep
#cat MRJPs.pep.aln.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln

antAln=$(readlink -f MRJPs.pep.aln.best.fas)

##################################################################
# TAG: 4. Create A. mellifera MRJP/yellow phylogeny
##################################################################
cd $base/Amel.reference
mkdir gb
esearch -db protein -query "(yellow OR MRJP) AND (Apis mellifera[organism]) AND srcdb+refseq[prop]" < /dev/null|elink -target gene|efetch -format tabular > Amel.y.tsv

for i in $(cut -f 3 Amel.y.tsv|sed 1d); do
  esearch -db protein -query "($i) AND (Apis mellifera[organism]) AND srcdb+refseq[prop]" < /dev/null|efetch -format gb > $base/gb/Amel.$i.gb

for i in $(cut -f 3 Amel.y.tsv|sed 1d); do
# convert gb entry to fasta with gene="<>" tag and all splice forms
cat gb/Amel.$i.gb|\
awk '/^ACCESSION   / {printf(">%s\t",$2);next;} \
     /[[:space:]]*\/gene\=/ {{match($0,"gene=.*",a)}{print a[0];next;}} \
     /^ORIGIN/ {inseq=1;next;} /^\/\// {inseq=0;} {if(inseq==0) next;\
      gsub(/[0-9 ]/,"",$0); printf("%s\n",toupper($0));}' > fa/Amel.$i.fa

# select longest isoform from all splice forms
cat fa/Amel.$i.fa|\
awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' |\
awk -F '\t'  '{printf("%s\t%d\n",$0,length($3));}' |\
sort -k2,2 -k4,4nr |\
sort -k2,2 -u -s |\
awk '{print $1";symbol="$2";length="$4"\n"$3}'> lfa/Amel.$i.longestIsoform.fa
done
cd $base/aln

cat ../lfa/* > Amel.fa
unset MAFFT_BINARIES
prank -o=Amel.mrjp.aln -d=Amel.fa

AmelAln=$(readlink -f Amel.mrjp.aln.best.fas)

##################################################################
# TAG: 5. Combine Ant and Amel alignments
##################################################################

cd $base/aln/mergeAln
prank -d1=$antAln -d2=$AmelAln -o=merged.aln
cat merged.aln.fas |perl -pe 's/(>.*?-.*?-.*?)-.*/$1/g'|perl -pe 's/.*gene=\"(.*?)\";.*/>$1/g' > output.aln.fas
FastTree output.aln.fas > output.aln.tre

# RaxML on aa aln
nice raxmlHPC-PTHREADS -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE -s output.aln.fas -n r1 -T 20


#FastTree tmp.aln > MRJP.best.tre

##################################################################
# TAG: 3. raxML phylogeny
##################################################################

nice raxmlHPC-PTHREADS -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE -soutput.aln.fas -n r1 -T 20
#/corefac/cse/lukas/inqGen18/geneFamilies/specificFamilies/MRJPs/aln/mergeAln/RAxML_bipartitionsBranchLabels.r1
