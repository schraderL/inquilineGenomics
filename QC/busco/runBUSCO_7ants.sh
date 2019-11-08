###########################################
# QC of genome assemblies and annotations
###########################################
# We ran BUSCO on the genome assembly and the set of predicted proteins for each of the attine genomes.

#######################
# BUSCO
#######################

# prepare BUSCO env
export AUGUSTUS_CONFIG_PATH="/usr/local/home/lschrader/software/augustus-3.2.3/config/"

# prepare environment
export base="/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7antsHYM"
export baseDB="/corefac/cse/lukas/inqGen18/QC/BUSCO_inquilinesv2.1/"
#######################
# run for each species
#######################

cd $base

export genomeBase=~/data/genomes/reannotations_of_7_ants/

# genome mode
export BUSCO_CONFIG_FILE=$base/genome.config.ini

for genome in $(readlink -f  $genomeBase/*/genome/*.v2.0.fa)
do
  species=$(echo $genome|perl -pe 's/.*\/(.).*\_(...).*/$1$2/g')
  out=$(echo $genome|rev|cut -f 1 -d "/" |rev)
  nice python ~/software/buscoV3/scripts/run_BUSCO.py --in $genome --out $out --lineage $baseDB/datasets/hymenoptera_odb9 --mode genome --cpu 1 -f
done

# protein mode
export BUSCO_CONFIG_FILE=$base/protein.config.ini

for protein in $(readlink -f  $genomeBase/*/annotation/gene_annotation/*.pep.fa)
do
  species=$(echo $protein|perl -pe 's/.*\/(.).*\_(...).*/$1$2/g')
  out=$(echo $protein|rev|cut -f 1 -d "/" |rev)
  nice  python ~/software/buscoV3/scripts/run_BUSCO.py --in $protein --out $out.prot --lineage $baseDB/datasets/hymenoptera_odb9 --mode prot --cpu 20
done
