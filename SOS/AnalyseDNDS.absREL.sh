base=/usr/local/home/lschrader/data/inqGen18/SOS

# HYPHY 2.3.14.20190217beta(MP) for Linux on x86_64


cd $base/data
ln -s /corefac/cse/lukas/inqGen18/orthologs/orthofinder/cleanData.Acromyrmex/Results_Feb15/SingleCopyOrthogroups.txt .
ln -s /corefac/cse/lukas/inqGen18/orthologs/orthofinder/cleanData.Acromyrmex/Results_Feb15/Orthogroups.txt .
ln -s ~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre .

SCO=$base/data/SingleCopyOrthogroups.txt
OGs=$base/data/Orthogroups.txt

tree="(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);"


########################################
# TAG Retrieve fasta files
########################################

########################################
# Acromyrmex
########################################
bd=/usr/local/home/lschrader/data/inqGen18/
cd $base/data/pep
for file in $(find $bd/inquilines_v2.1/  -name "*.pep.fa")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file $short.pep.fa
done

cd $base/data/cds
for file in $(find $bd/inquilines_v2.1/  -name "*.cds.fa")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file $short.cds.fa
done

########################################
# Acol & Aech
########################################

r7=~/data/genomes/reannotations_of_7_ants
cd $base/data/pep
file=/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Atta_colombica/annotation/gene_annotation/Atta_colombica.v2.0.pep.fa
echo Acol with $file
ln -s $file Acol.pep.fa

cd $base/data/cds
file=/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Atta_colombica/annotation/gene_annotation/Atta_colombica.v2.0.cds.fa
echo Acol with $file
ln -s $file Acol.cds.fa


r7=~/data/genomes/reannotations_of_7_ants
cd $base/data/pep
file=/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/annotation/gene_annotation/Acromyrmex_echinatior.v2.0.pep.fa
echo Aech with $file
ln -s $file Aech.pep.fa

cd $base/data/cds
file=/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/annotation/gene_annotation/Acromyrmex_echinatior.v2.0.cds.fa
echo Aech with $file
ln -s $file Aech.cds.fa


#######################
# TAG 1 Select SCOs
#######################

cd $base
# cat all fastas but rename to standard ids (ie. ACEP12345)
mkdir tmp
cat $base/data/pep/*pep.fa|perl -pe 's/(>[A-Z]{4}[0-9]{5}).*/$1/g' > tmp/pep.all.fa
cat $base/data/cds/*cds.fa|perl -pe 's/(>[A-Z]{4}[0-9]{5}).*/$1/g' > tmp/cds.all.fa
mkdir pep.OG
mkdir cds.OG
mkdir OGlst
grep -f $SCO $OGs > tmp/SCO.tsv

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
# TAG Align with prank
#######################
software=/usr/local/home/lschrader/software
scripts=$base/scripts/

cd $base
mkdir prank
mkdir prank/pep
mkdir prank/cds
mkdir prank/guidetrees

unset MAFFT_BINARIES

# prepare species trees with! branch lenghts.
echo $tree|perl -pe 's/\:[0-9]\.[0-9]{6}//g' > tmp/unscaledSpeciesTree.tre
echo $tree > tmp/unscaledSpeciesTree.branchLengths.tre
# label nodes
# (((AHEY,(ACHA,PARG)P1)A1,(AECH,AINS)A2)AO,ACOL)LC;
# P1= origin of parasitism
# A1= origin of Acromyrmex heyeri/A. charruanus/P.argentina
# A2= origin of A. echinatior/A. insunuator
# AO= origin of Acromyrmex (letter "o")
# LC = origin of leaf cutters

ls OGlst/*.lst|parallel --nice 10 perl $scripts/renameTree.2.1.pl {} tmp/unscaledSpeciesTree.branchLengths.tre  prank/guidetrees/{/.}.pep.tre

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
#deleting prank/cds/OG0003876.cds.aln.fa
#deleting prank/cds/OG0007381.cds.aln.fa
#######################


################################################################################
# 2.3 Test for recombination in alignments
################################################################################
mkdir $base/recombination/
cd $base/recombination/
mkdir $base/prank/cds.recombining/
ls $base/prank/cds/*.cds.aln.fa| parallel --nice 10 'Phi -f {} |grep -Po "OG.*?\.|PHI.*"'|awk 'ORS=NR%2?" ":"\n"' > Phi.test.tsv &
perl -pe 's/\..*\://g' Phi.test.tsv |awk '{if ($2< 0.01) print $1}'|grep "\-\-" -v > Phi.significant.lst

while read p;
do
    mv $base/prank/cds/$p* $base/prank/cds.recombining/
done < $base/recombination/Phi.significant.lst

################################################################################
# TAG Run ABSrel
################################################################################

#http://hyphy.org/resources/json-fields.pdf
#https://github.com/veg/hyphy/issues/752
mkdir $base/absREL
cd $base/absREL

HMdir=/corefac/cse/lukas/software/hyphy-2.3.14
bf=$HMdir/res/TemplateBatchFiles/SelectionAnalyses/aBSREL.bf

ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "mkdir $base/absREL/{/.}/"
ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "ln -s $base/prank/cds/{}.cds.aln.fa $base/absREL/{/}/{}.aln"
# started April 12, 12:16
ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "HYPHYMP LIBPATH=$HMdir/res/ $bf 'Universal' '$base/absREL/{/}/{}.aln' '$base/prank/guidetrees/{/}.pep.tre' 'All' '$base/absREL/{/}/{}.absREL.out'"

# Generates HYPHY commands:
#HYPHYMP
# LIBPATH=/corefac/cse/lukas/software/hyphy-2.3.14/res/
# /corefac/cse/lukas/software/hyphy-2.3.14/res/TemplateBatchFiles/SelectionAnalyses/aBSREL.bf
# 'Universal'
# '/usr/local/home/lschrader/data/inqGen18/SOS/absREL/OG0001201/OG0001201.aln'
# '/usr/local/home/lschrader/data/inqGen18/SOS/prank/guidetrees/OG0001201.pep.tre'
# 'All'
# '/usr/local/home/lschrader/data/inqGen18/SOS/absREL/OG0001201/OG0001201.absREL.out'
HYPHYMP LIBPATH=/corefac/cse/lukas/software/hyphy-2.3.14/res/ /corefac/cse/lukas/software/hyphy-2.3.14/res/TemplateBatchFiles/SelectionAnalyses/aBSREL.bf 'Universal' '/usr/local/home/lschrader/data/inqGen18/SOS/absREL/OG0001202/OG0001202.aln' '/usr/local/home/lschrader/data/inqGen18/SOS/prank/guidetrees/OG0001202.pep.tre' 'All' '/usr/local/home/lschrader/data/inqGen18/SOS/absREL/OG0001202/OG0001202.absREL.out'


################################################################################
# TAG Load all json files
################################################################################

scp -r -o ProxyCommand="ssh -W %h:%p zcm795@ssh-bio-labnet.science.ku.dk" lschrader@pallas.bio.ku.dk:/usr/local/home/lschrader/data/inqGen18/SOS/absREL/*/*.json /Users/lukas/sciebo/inquilineGenomics18/SOS/results/

################################################################################
# TAG Run RELAX: Test=all Parasites; Reference=all other
################################################################################
mkdir $base/prank/gtRELAX
ls prank/guidetrees/*|parallel --nice 10 "cat {}| sed 's/\(AINS[0-9][0-9][0-9][0-9][0-9]:[0-9]\.[0-9]*\)/\1\{T\}/g' |sed 's/\(PARG[0-9][0-9][0-9][0-9][0-9]:[0-9]\.[0-9]*\)/\1\{T\}/g' |sed 's/\(ACHA[0-9][0-9][0-9][0-9][0-9]:[0-9]\.[0-9]*\)/\1\{T\}/g' | sed 's/\()P1\:[0-9]\.[0-9]*\)/\1\{T\}/g ' > $base/prank/gtRELAX/{/.}.tre"
HMdir=/corefac/cse/lukas/software/hyphy-2.3.14
bf=$HMdir/res/TemplateBatchFiles/SelectionAnalyses/RELAX.bf

ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "mkdir $base/RELAX/{/.}/"
ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "ln -s $base/prank/cds/{}.cds.aln.fa $base/RELAX/{/}/{}.aln"
ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "HYPHYMP LIBPATH=$HMdir/res/ $bf 'Universal' '$base/RELAX/{/}/{}.aln' '$base/prank/gtRELAX/{/}.pep.tre' 'T' 'Minimal'"


################################################################################
# TAG Run RELAX: Test=all Parasites; Reference=AHEY,AECH,A1,A2
################################################################################
mkdir $base/prank/gtRELAX2
mkdir $base/RELAX2
# Test=all Parasites; Reference=AHEY,AECH,A1,A2
ls prank/guidetrees/*|parallel --nice 10 "cat {}| sed 's/\()A1\:[0-9]\.[0-9]*\)/\1\{R\}/g'| sed 's/\()A2\:[0-9]\.[0-9]*\)/\1\{R\}/g'  | sed 's/\(AECH[0-9][0-9][0-9][0-9][0-9]:[0-9]\.[0-9]*\)/\1\{R\}/g' | sed 's/\(AHEY[0-9][0-9][0-9][0-9][0-9]:[0-9]\.[0-9]*\)/\1\{R\}/g' | sed 's/\(AINS[0-9][0-9][0-9][0-9][0-9]:[0-9]\.[0-9]*\)/\1\{T\}/g' |sed 's/\(PARG[0-9][0-9][0-9][0-9][0-9]:[0-9]\.[0-9]*\)/\1\{T\}/g' |sed 's/\(ACHA[0-9][0-9][0-9][0-9][0-9]:[0-9]\.[0-9]*\)/\1\{T\}/g' | sed 's/\()P1\:[0-9]\.[0-9]*\)/\1\{T\}/g '> $base/prank/gtRELAX2/{/.}.tre"

HMdir=/corefac/cse/lukas/software/hyphy-2.3.14
bf=$HMdir/res/TemplateBatchFiles/SelectionAnalyses/RELAX.bf

ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "mkdir $base/RELAX2/{/.}/"
ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "ln -s $base/prank/cds/{}.cds.aln.fa $base/RELAX2/{/}/{}.aln"
ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "HYPHYMP LIBPATH=$HMdir/res/ $bf 'Universal' '$base/RELAX2/{/}/{}.aln' '$base/prank/gtRELAX2/{/}.pep.tre' 'T' 'R' 'Minimal'"

################################################################################
# TAG Run RELAX: Test=all Parasites; Reference=AHEY,AECH,A1,A2, Full model
################################################################################
mkdir $base/RELAX3
# Test=all Parasites; Reference=AHEY,AECH,A1,A2; Model=Full
#$base/prank/gtRELAX2/*.tre

HMdir=/corefac/cse/lukas/software/hyphy-2.3.14
bf=$HMdir/res/TemplateBatchFiles/SelectionAnalyses/RELAX.bf

ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "mkdir $base/RELAX3/{/.}/"
ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "ln -s $base/prank/cds/{}.cds.aln.fa $base/RELAX3/{/}/{}.aln"
ls $base/prank/cds/*.aln.fa|rev|cut -f 1 -d "/"|rev|cut -f 1 -d "."|parallel --nice 10 "HYPHYMP LIBPATH=$HMdir/res/ $bf 'Universal' '$base/RELAX3/{/}/{}.aln' '$base/prank/gtRELAX2/{/}.pep.tre' 'T' 'R' 'All'"
