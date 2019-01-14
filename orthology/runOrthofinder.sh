########################################
# 1. Infer orthology with ORTHOFINDER
########################################

# Run orthofinder to infer orthology between species. Using all 7 reannotated attines and the 4 new genomes.
base=/usr/local/home/lschrader/data/inqGen18/orthologs/orthofinder/
bd=/usr/local/home/lschrader/data/inqGen18/


########################################
# get soft links to original peptide files
for file in $(find $bd/inquilines_v2.1/  -name "*.pep.fa")
do
short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
echo $short with $file
ln -s $file data/$short.pep.fa
done
########################################

########################################
# clean up peptide files to make their structure identical between inquilines and attines
mkdir $base/cleanData
cd $base/cleanData
for file in $(find $base/data/ -name "*.pep.fa")
do
  echo "cleaning $file"
  name=$(echo $file|rev|cut -f 1 -d "/"|rev)
  perl -pe 's/(>.*?)[\s|-].*/$1$2/g' $file > $name
done
########################################

########################################
# run Orthoginder on cleaned pep files
cd $base
# Use files with "TE proteins"
# consider removing all those proteins that have been identified as being TE-derived.

nice orthofinder -t 32 -a 10 -f cleanData/
 mv cleanData/Results_Sep03/ $
