###########################################
# QC of genome assemblies and annotations
###########################################
# We ran BUSCO on the genome assembly and the set of predicted proteins for each of the 4 genomes.

#######################
# BUSCO
#######################

# prepare BUSCO env
export AUGUSTUS_CONFIG_PATH="/usr/local/home/lschrader/software/augustus-3.2.3/config/"

# prepare environment
export base="/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_inquilinesv2.1"

#######################
# run for each species
#######################

cd $base
#gunzip ~/data/genomes/inquilines_2.0/*/genome/*.fa.gz

genomeBase=~/data/inqGen18/inquilines_v2.1/
ACHAG=$genomeBase/Acromyrmex_charruanus.2.1/genome/Acromyrmex_charruanus.v2.1.fa
AHEYG=$genomeBase/Acromyrmex_heyeri.2.1/genome/Acromyrmex_heyeri.v2.1.fa
AINSG=$genomeBase/Acromyrmex_insinuator.2.1/genome/Acromyrmex_insinuator.v2.1.fa
PARGG=$genomeBase/Pseudoatta_argentina.2.1/genome/Pseudoatta_argentina.v2.1.fa

ACHA=$genomeBase/Acromyrmex_charruanus.2.1/annotation/gene_annotation/Acromyrmex_charruanus.v2.1.pep.fa
AHEY=$genomeBase/Acromyrmex_heyeri.2.1/annotation/gene_annotation/Acromyrmex_heyeri.v2.1.pep.fa
AINS=$genomeBase/Acromyrmex_insinuator.2.1/annotation/gene_annotation/Acromyrmex_insinuator.v2.1.pep.fa
PARG=$genomeBase/Pseudoatta_argentina.2.1/annotation/gene_annotation/Pseudoatta_argentina.v2.1.pep.fa

for species in $AINSG $AHEYG $PARGG $ACHAG $AINS $AHEY $PARG $ACHA
do
  ls -lh $species
done

# genome mode
export BUSCO_CONFIG_FILE=$base/genome.config.ini
$ACHAG
for species in $AINSG $AHEYG $PARGG
do
  out=$(echo $species|rev|cut -f 1 -d "/" |rev)
  nice  python ~/software/buscoV3/scripts/run_BUSCO.py --in $species --out $out --lineage $base/datasets/hymenoptera_odb9 --mode genome --cpu 1 -f
done

# protein mode
export BUSCO_CONFIG_FILE=$base/protein.config.ini
for species in $AINS $AHEY $PARG $ACHA
do
  out=$(echo $species|rev|cut -f 1 -d "/" |rev)
  nice  python ~/software/buscoV3/scripts/run_BUSCO.py --in $species --out $out.prot --lineage $base/datasets/hymenoptera_odb9 --mode prot --cpu 20
done
