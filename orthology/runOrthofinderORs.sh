########################################
# 1. Infer orthology with ORTHOFINDER
########################################

# Run orthofinder to infer orthology between ORs in attines
base=/usr/local/home/lschrader/data/inqGen18/ORphylogeny/orthofinder/
bd=/usr/local/home/lschrader/data/inqGen18/ORs/


########################################
# get soft links to original peptide files
for file in $(find $bd/*/final/finalSet/  -name "*.fa"|grep "NCBI" -v)
do
ln -s $file $base/data/.
done
########################################

################################################################################
# RUN ON ALL ATTINES
################################################################################

########################################
# run Orthoginder on OR pep files
cd $base
nice orthofinder -t 32 -a 10 -f data/
mv data/Results_Apr19/ .
