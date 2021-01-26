
# QC of genome assemblies and annotations

 We ran BUSCO on the genome assembly and the set of predicted proteins for each of the 4 genomes.


## BUSCO

```bash

# prepare BUSCO env
export AUGUSTUS_CONFIG_PATH="/corefac/cse/lukas/software/augustus-3.2.3/config/"

# prepare environment
export base="/corefac/cse/lukas/inqGen18/reannotation/"


#######################
# run for each species
#######################

cd $base
cp /usr/local/home/lschrader/data/inqGen18/QC/BUSCO_inquilinesv2.1/protein.config.ini .

# protein mode
export BUSCO_CONFIG_FILE=$base/protein.config.ini
for species in /corefac/cse/lukas/inqGen18/reannotation/results/*/GAF/*.GeMoMa.longestIsoform.pep.fa
do
  out=$(echo $species|rev|cut -f 1 -d "/" |rev)
  echo "nice  python ~/software/buscoV3/scripts/run_BUSCO.py --in $species --out $out.BUSCO.prot --lineage /usr/local/home/lschrader/data/inqGen18/QC/BUSCO_inquilinesv2.1/datasets/hymenoptera_odb9 --mode prot --cpu 10"
done
#nice  python ~/software/buscoV3/scripts/run_BUSCO.py --in /corefac/cse/lukas/inqGen18/reannotation/results/Acromyrmex_echinatior.v2.0/GAF/Aech.GeMoMa.longestIsoform.pep.fa --out Aech.GeMoMa.longestIsoform.pep.fa.BUSCO.prot --lineage /usr/local/home/lschrader/data/inqGen18/QC/BUSCO_inquilinesv2.1/datasets/hymenoptera_odb9 --mode prot --cpu 10

```
