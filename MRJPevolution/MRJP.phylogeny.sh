##################################################################
# TAG: 1. prepare environment
##################################################################
scripts=~/data/inqGen18/geneFamilies/specificFamilies/scripts
base=~/data/inqGen18/geneFamilies/specificFamilies/

cd $base/cleanGFFs

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

cd $base/fasta
ln -s ../cleanGFFs/*/*pep .
for i in $(ls *.MRJPs.pep)
do
  species=$(echo $i|cut -f 1 -d ".")
  cat $i |sed "s/>/>$species-/g"
done > MRJPs.pep
unset MAFFT_BINARIES
prank MRJPs.pep

cat output.best.fas |perl -pe 's/(>.*?)_.*/$1/g' > tmp.aln
FastTree tmp.aln > MRJP.best.tre


##################################################################
# TAG: 3. raxML phylogeny
##################################################################

nice raxmlHPC-PTHREADS -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE -soutput.aln.fas -n r1 -T 20
#/corefac/cse/lukas/inqGen18/geneFamilies/specificFamilies/MRJPs/aln/mergeAln/RAxML_bipartitionsBranchLabels.r1
