###########################################
# QC of genome assemblies and annotations
###########################################
# We ran BUSCO on the genome assembly and the set of predicted proteins for each of the 7 reannotated genomes.

#######################
# BUSCO
#######################

# prepare BUSCO env
export AUGUSTUS_CONFIG_PATH="/usr/local/home/lschrader/software/augustus-3.2.3/config/"

# prepare environment
export base="/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants"

#######################
# run for each species
#######################

for species in Acep Acol Aech Ccos Tcor Tsep Tzet
do
  cd $base
  mkdir $species
  mkdir $base/$species/genome
  mkdir $base/$species/protein
  ## Run BUSCO
  # genome
  cd $base/$species/genome
  export BUSCO_CONFIG_FILE=$base/configs/$species.genome.config.ini
  #nice python ~/software/buscoV3/scripts/run_BUSCO.py -f > $species.genome.busco.out
  # protein
  cd $base/$species/protein
  export BUSCO_CONFIG_FILE=$base/configs/$species.protein.config.ini
  nice python ~/software/buscoV3/scripts/run_BUSCO.py -f > $species.protein.busco.out
done


cd $base
mkdir $base/output/

for species in Acep Acol Aech Ccos Tcor Tsep Tzet
do
 cp $base/$species/protein/run_"$species"_protein/short_summary_"$species"_protein.txt ./output/
 cp $base/$species/genome/run_"$species"_genome/short_summary_"$species"_genome.txt ./output/
done

$base/$species/protein/run_"$species"_genome/short_summary_"$species"_genome.txt
