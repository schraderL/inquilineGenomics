
base=/usr/local/home/lschrader/data/inqGen18/WGA/progressiveCactus
bd=/usr/local/home/lschrader/data/inqGen18/
cd $bd
mkdir /usr/local/home/lschrader/data/inqGen18/WGA/
mkdir $base
cd $base
mkdir data/

########################################
# get soft links to original genome files for inquiline genomes
rm $base/data/genomeFile.txt
for file in $(find $bd/inquilines_v2.1/*/genome/  -name "*v2.1.fa"|egrep -v "fungal|mitome|bacterial")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  SHORT=$(echo $short|tr [a-z] [A-Z])
  echo $SHORT $file >> $base/data/genomeFile.txt
  ln -s $file $base/data/$short.genome.fa
done
########################################

########################################
# get soft links to original genome files for 7 ants
for file in $(find /usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/*/genome/ -name "*.fa")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  SHORT=$(echo $short|tr [a-z] [A-Z])
  echo $SHORT $file >> $base/data/genomeFile.txt
  ln -s $file $base/data/$short.genome.fa
done

########################################

treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")



echo "
# Sequence data for progressive alignment of 11 genomes
# Aech is flagged as good assemblies.
# <this file> /usr/local/home/lschrader/data/WGA/progressiveCactus/attines.txt
# <workDir> = /usr/local/home/lschrader/data/inqGen18/WGA/progressiveCactus/attines
# <out HAL> = /usr/local/home/lschrader/data/inqGen18/WGA/progressiveCactus/attines/attines.hal
# tree from $treeFile
$treeOriginal" > $base/tmp.txt
cat $base/tmp.txt $base/data/genomeFile.txt |sed 's/^AECH/*AECH/g'> $base/attines.txt
rm tmp.txt


source /usr/local/home/lschrader/software/cactus/progressiveCactus/environment
cd /usr/local/home/lschrader/software/cactus/progressiveCactus/
nice bin/runProgressiveCactus.sh $base/attines.txt \
$base/attines \
$base/attines/attines.hal \
--database kyoto_tycoon --maxThreads 20 --overwrite

cd $base/attines
source /usr/local/home/lschrader/software/cactus/progressiveCactus/environment
halValidate attines.hal
halSummarizeMutations attines.hal > mutations.attines.tsv
halStats  --allCoverage attines.hal > stats.attines.tsv

#hal2maf Ahey.Ains.Aech.Acol.hal Ahey.Ains.Aech.Acol.maf

for node in AHEY ACHA AECH AINS Anc09
do
  nice halBranchMutations attines.hal $node --refFile $node.ins.bed --parentFile $node.del.bed
done

for node in ACOL ACEP Anc06 Anc07 Anc05 Anc08 Anc06 Anc04
do
  nice halBranchMutations attines.hal $node --refFile $node.ins.bed --parentFile $node.del.bed
done



#nice halPhyloPTrain.py attines.hal AECH neutralRegions.bed neutralModel.mod --numProc 12
