# see /Users/lukas/CSE/inquilineGenomics/Phylogeny4/workflow.phylogeny4.sh

#######################
# Prepare environment
#######################

#orthologBase=/usr/local/home/lschrader/data/inqGen18/orthologs/orthofinder/Results_Jul09/
orthologBase=/usr/local/home/lschrader/data/inqGen18/orthologs/orthofinder/Results_Sep03/
base=~/data/inqGen18/phylogeny/SCO
data=~/data/inqGen18/phylogeny/data
scripts=~/data/inqGen18/phylogeny/scripts/

cd $data

##################################
# prepare 4 inquiline genome files
#######################

# peptide
bd=/usr/local/home/lschrader/data/inqGen18/
for file in $(find $bd/inquilines_v2.1/  -name "*.pep.fa")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file $short.pep.fa
done

# cds
for file in $(find $bd/inquilines_v2.1/  -name "*.cds.fa")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file $short.cds.fa
done

##################################
# prepare 7 attine genome files
#######################

# peptide
r7=~/data/genomes/reannotations_of_7_ants
for file in $(find $r7/*/  -name "*.pep.fa")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file $short.pep.fa
done

# cds
for file in $(find $r7/*/  -name "*.cds.fa")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file $short.cds.fa
done
##################################


#######################
# 1 Select SCOs
#######################

cd $base
# cat all fastas but rename to standard ids (ie. ACEP12345)
mkdir tmp
cat $data/*pep.fa|perl -pe 's/(>[A-Z]{4}[0-9]{5}).*/$1/g' > tmp/pep.all.fa
cat $data/*cds.fa|perl -pe 's/(>[A-Z]{4}[0-9]{5}).*/$1/g' > tmp/cds.all.fa
mkdir pep.OG
mkdir cds.OG
mkdir OGlst
SCOG=$(readlink -f $orthologBase/SingleCopyOrthogroups.txt)
grep -f $SCOG $orthologBase/Orthogroups.txt > tmp/SCO.tsv

# create OG group files for fasta retrieval. Rename genes to standard convention
while read p; do
  OG=$(echo $p|cut -f 1 -d ":")
  echo $p|cut -f 2 -d ":"|tr " " "\n"|sed '/^\s*$/d' > OGlst/$OG.lst
done < tmp/SCO.tsv

# get all PEP and CDS per OG
ls OGlst/*.lst|parallel --nice 10 faSomeRecords tmp/cds.all.fa {} cds.OG/{/.}.cds.fa
ls OGlst/*.lst|parallel --nice 10 faSomeRecords tmp/pep.all.fa {} pep.OG/{/.}.pep.fa
#######################


#######################
# 2 Align with prank
#######################
software=/usr/local/home/lschrader/software
scripts=~/data/inqGen18/phylogeny/scripts/
cd $base
mkdir prank
mkdir prank/pep
mkdir prank/cds
mkdir prank/guidetrees


unset MAFFT_BINARIES

# prepare species trees
echo "(CCOS,(TZET,(TCOR,(TSEP,((ACOL,ACEP),((AHEY,(ACHA,PARG)),(AECH,AINS)))))));" > tmp/unscaledSpeciesTree.tre
ls OGlst/*.lst|parallel --nice 10 perl $scripts/renameTree.2.1.pl {} tmp/unscaledSpeciesTree.tre  prank/guidetrees/{/.}.pep.tre

# align sequences with prank
ls pep.OG/*.fa|parallel --nice 10 prank -t=prank/guidetrees/{/.}.tre -d={} -o=prank/pep/{/.}.out -quiet

# back-translate with pal2nal
ls prank/pep/*.fas|sed s/.pep.out.best.fas//g|parallel --nice 10 "perl $software/pal2nal.v14/pal2nal.pl  {}.pep.out.best.fas cds.OG/{/}.cds.fa -nogap -nomismatch -output fasta > prank/cds/{/}.cds.aln.fa"

# remove files with sequence length of 0
for file in prank/cds/*
do
two=$(head -n2 $file |tail -n 1)
 if [ -z $two ]
 then
 echo "deleting $file"
# rm $file
 fi
done
#######################

################################################################################
# 2.3 Test for recombination in alignments
################################################################################
mkdir $base/recombination/
cd $base/recombination/
mkdir $base/prank/cds.recombining/
ls $base/prank/cds/*.cds.aln.fa| parallel --nice 10 'Phi -f {} |grep -Po "OG.*?\.|PHI.*"'|awk 'ORS=NR%2?" ":"\n"' > Phi.test.tsv
perl -pe 's/\..*\://g' Phi.test.tsv |awk '{if ($2< 0.01) print $1}'|grep "\-\-" -v > Phi.significant.lst


while read p;
do
    mv $base/prank/cds/$p* $base/prank/cds.recombining/
done < $base/recombination/Phi.significant.lst


###########################################################################
# 2.4 Extract 4d sites
###########################################################################
cd $base
# Extract 4-fold degenerate sites from untrimmed CDS alignment
mkdir 4d/
mkdir 4d/OGs
ls $base/prank/cds/*.cds.aln.fa|parallel --nice 10  "perl $scripts/4dSites.pl {} fasta 4d/OGs/{/.}.4d" > 4d/CDSaln4d.statistics.tsv

# remove empty files
find $base/4d/OGs/* -size 0 -delete

################################################################################
# 2.5 Identify saturated 4d alignments
################################################################################


Rscript $scripts/analyseAln.R $base/4d/OGs/ "*.4d$" QC4d
Rscript $scripts/postProcess4dQC.R  QC4ddists.titvs.aln.Rfile QC4d.raw.tsv QC4d.
# creates 4 files
#QC4ddists.titvs.aln.Rfile  QC4d.pdf  QC4d.raw.tsv  QC4.saturated.tsv

# QC4.saturated.tsv contains all OGs that have a correlation coefficient of the uncorrected genetic distance to the TN93 genetic distance (< 0.6241, 1.5 IQR for all correlations)
mkdir $base/4d/saturated.4d
sed 1d $base/4d/QC4d/QC4.saturated.tsv|cut -f 3|xargs -I -- mv $base/4d/OGs/-- $base/4d/saturated.4d/

################################################################################
# 2.6 concatenate Alignments
################################################################################
mkdir $base/4d/4dtmp
cd $base
ls $base/4d/OGs/*|parallel --nice 10 'perl -pe "s/(>....).*/\$1/g" {} > ./4d/4dtmp/{/}'
python3 $software/amas-0.98/amas/AMAS.py concat -f fasta -d dna -i 4d/4dtmp/*.4d
mv concatenated.out 4d/4dAll.fa
mv partitions.txt 4d/4dAll.partitions.txt
#prank -convert -d=4d/4dAll.fa -o=4d/4dAll -f=phylipi
rm -rf 4d/4dtmp/


################################################################################
# 2.7 summary statistics
################################################################################
cd $base/4d/QC4d
python /usr/local/home/lschrader/software/amas-0.98/amas/AMAS.py summary -i $base/4d/4dAll.fa -f fasta -d dna
# creates $base/4d/QC4d/summary.txt
