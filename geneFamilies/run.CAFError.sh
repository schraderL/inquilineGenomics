# ************************************************************************
# SUMMARY: Runs 1 and 2 yielded ~same results
# ************************************************************************
:'
# =======================================================================
# Final error estimates by species:
# PARG                             0.129265136719
# ACOL                             0.0
# AINS                             0.0430883789063
# AHEY                             0.0478759765625
# ACHA                             0.0287255859375
# AECH                             0.0526635742187
# Score with individual errors:    27970.041073
# Lambda with individual errors:   0.00048795590293
# =======================================================================
# ************************************
# Score with no errormodel:        30193.176929
# Lambda with no errormodel:       0.00095486428317
# ************************************
# Global Error Estimation:         0.0478759765625
# Score with global errormodel:    28615.716574
# Lambda with global errormodel:   0.00052122244019
# ************************************
'

###############################
#TAG: Prepare template CAFError script
###############################
mkdir ~/data/inqGen18/geneFamilies/CAFE/
base=~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction
bd=~/data/inqGen18/geneFamilies
software=~/data/software
mkdir $base
mkdir $base/ctl
cd $base/

echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

load -i <CLUSTERS> -t 20 -l <OUTPUTSUMMARY> -p 0.05 -filter

tree <TREE>

#Estimate lambda for all species using a single rate
# Error model estimate restricted to lambda. Mu is not included.
lambda -s -t <LAMBDATREE>
report <OUTPUTREPORT>" \
> $base/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe

###############################
#TAG: Prepare subsetted MCL cluster file
###############################

# input files
originalClusters=$bd/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

###############################
#TAG: Filter genes with >100 copies in one species from Cluster file
###############################
cd $base
rawclusters=$bd/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
cat $rawclusters |awk  -F $'\t' 'BEGIN {OFS=FS} {if (($3 <100 && $4 <100 && $5 <100 && $6 <100 && $7 <100 && $8 <100 && $9 <100 && $10 <100) || $1=="Description") print $0}' \
> $bd/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv

############################################################################################################################

###############################
#RUN: 1: Run errorTest.smallSet.filter100
###############################
cd $base
runName=RUN1
############################################
#TAG: INPUT files
############################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree:
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run4.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error prediction cafe script
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe
'

# input files
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
ls $treeFile

###############################
# TAG: Define runName
###############################
runName=errorTest.smallSet.filter100

  runbase=$base/$runName
  cscript=$runbase/$runName.cafe
  outsummary=$runbase/$runName.out
  outreport=$runbase/report/$runName.report.out

# prepare environment
#mkdir $base/CAFE/CAFE.errorprediction/run1 # First test run. Ignore.
mkdir $runbase
mkdir $runbase/tmp
mkdir $runbase/err
cd $runbase

## Run the error prediction script, predicting error rates for each species

###############################
# TAG: prepare tree from MCMCtree output
###############################
treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# select only Acro and Acol

#(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526
tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
echo $tree

###############################
# TAG: Prepare cafe script from template.
###############################
cat $base/ctl/errorTest.template.cafe | \
  sed "s|<CLUSTERS>|$clusters|g" | \
  sed "s|<TREE>|$tree|g" | \
  sed "s|<LAMBDATREE>|$LAMBDATREE|g" | \
  sed "s|<OUTPUTSUMMARY>|$outsummary|g" | \
  sed "s|<OUTPUTREPORT>|$outreport|g" \
> $cscript
less $cscript

###############################
#TAG: Run CAFError
###############################
nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1 > caferror.run.out

###############################
#TAG:  Gather best fitting error models
###############################

cd $runbase
grep "Final error estimates by species:" $runbase/tmp/caferrorLog.txt -A 200 > $runbase/FinalErrorEstimates.txt

# get final estimates for each species
for i in AECH AHEY AINS ACHA ACOL PARG
do
  #grep $i tmp|perl -pe 's/# ([A-Z]{4}).?+([0-9].+$)/$1\t$2/g'
  grep $i $runbase/FinalErrorEstimates.txt |perl -pe 's/# ([A-Z]{4}).+?([0-9].+$)/$1\tcafe_errormodel_$2.txt/g'
done > $runbase/FinalErrorEstimates.perSpecies.txt

# Retrieve all scores for each species
mkdir $runbase/err/
cd $runbase/tmp/
for i in AECH AHEY AINS ACHA ACOL PARG
do
 grep "Search Result"  *$i* -A 1|grep -E "Score:.*"|sed s/cafe_//g|sed s/'_.*Score: '/' '/g > ../err/$i.err
done

cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt

# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s $runbase/FinalErrorEstimates.txt .
ln -s $runbase/FinalErrorEstimates.perSpecies.txt .

############################################################################################################################
############################################################################################################################

###############################
#RUN: 2: Run errorTest.smallSet.filter100
###############################
cd $base
runName=RUN2
############################################
#TAG: INPUT files
############################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree:
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run1.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error prediction cafe script
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe
'

# input files
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
ls $treeFile

###############################
# TAG: Define runName
###############################
runName=errorTest.smallSet.filter100.RUN2

  runbase=$base/$runName
  cscript=$runbase/$runName.cafe
  outsummary=$runbase/$runName.out
  outreport=$runbase/report/$runName.report.out

# prepare environment
#mkdir $base/CAFE/CAFE.errorprediction/run1 # First test run. Ignore.
mkdir $runbase
mkdir $runbase/tmp
mkdir $runbase/err
cd $runbase

## Run the error prediction script, predicting error rates for each species

###############################
# TAG: prepare tree from MCMCtree output
###############################
treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# select only Acro and Acol
#((((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665):0.020045);
#(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526
tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
echo $tree

###############################
# TAG: Prepare cafe script from template.
###############################
cat $base/ctl/errorTest.template.cafe | \
  sed "s|<CLUSTERS>|$clusters|g" | \
  sed "s|<TREE>|$tree|g" | \
  sed "s|<LAMBDATREE>|$LAMBDATREE|g" | \
  sed "s|<OUTPUTSUMMARY>|$outsummary|g" | \
  sed "s|<OUTPUTREPORT>|$outreport|g" \
> $cscript
less $cscript

###############################
#TAG: Run CAFError
###############################
nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1 > caferror.run.out

###############################
#TAG:  Gather best fitting error models
###############################

cd $runbase
grep "Final error estimates by species:" $runbase/tmp/caferrorLog.txt -A 200 > $runbase/FinalErrorEstimates.txt

# get final estimates for each species
for i in AECH AHEY AINS ACHA ACOL PARG
do
  #grep $i tmp|perl -pe 's/# ([A-Z]{4}).?+([0-9].+$)/$1\t$2/g'
  grep $i $runbase/FinalErrorEstimates.txt |perl -pe 's/# ([A-Z]{4}).+?([0-9].+$)/$1\tcafe_errormodel_$2.txt/g'
done > $runbase/FinalErrorEstimates.perSpecies.txt

# Retrieve all scores for each species
mkdir $runbase/err/
cd $runbase/tmp/
for i in AECH AHEY AINS ACHA ACOL PARG
do
 grep "Search Result"  *$i* -A 1|grep -E "Score:.*"|sed s/cafe_//g|sed s/'_.*Score: '/' '/g > ../err/$i.err
done

cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt

# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s $runbase/FinalErrorEstimates.txt .
ln -s $runbase/FinalErrorEstimates.perSpecies.txt .

############################################################################################################################
# =======================================================================
# Final error estimates by species:
# PARG                             0.129265136719
# ACOL                             0.0
# AINS                             0.0430883789063
# AHEY                             0.0478759765625
# ACHA                             0.0287255859375
# AECH                             0.0526635742187
# Score with individual errors:    27970.041436
# Lambda with individual errors:   0.00048794472171
# =======================================================================
# ************************************
# Score with no errormodel:        30193.176198
# Lambda with no errormodel:       0.00095488221842
# ************************************
# Global Error Estimation:         0.0478759765625
# Score with global errormodel:    28615.714074
# Lambda with global errormodel:   0.00052123236134
# ************************************
# =======================================================================
# Caferror finished at:            11.24.2018 | 05:37:19
# Runtime:                         41.0628216346 minutes


############################################################################################################################

###############################
#RUN: 3: Run errorTest.smallSet.filter100 : simple repetition of first two runs
###############################
cd $base
runName=RUN3
############################################
#TAG: INPUT files
############################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree:
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run1.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error prediction cafe script
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe
'

# input files
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
ls $treeFile

###############################
# TAG: Define runName
###############################
runName=errorTest.smallSet.filter100.RUN3

  runbase=$base/$runName
  cscript=$runbase/$runName.cafe
  outsummary=$runbase/$runName.out
  outreport=$runbase/report/$runName.report.out

# prepare environment
#mkdir $base/CAFE/CAFE.errorprediction/run1 # First test run. Ignore.
mkdir $runbase
mkdir $runbase/tmp
mkdir $runbase/err
cd $runbase

## Run the error prediction script, predicting error rates for each species

###############################
# TAG: prepare tree from MCMCtree output
###############################
treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# select only Acro and Acol
#((((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665):0.020045);
#(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526
tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
echo $tree

###############################
# TAG: Prepare cafe script from template.
###############################
cat $base/ctl/errorTest.template.cafe | \
  sed "s|<CLUSTERS>|$clusters|g" | \
  sed "s|<TREE>|$tree|g" | \
  sed "s|<LAMBDATREE>|$LAMBDATREE|g" | \
  sed "s|<OUTPUTSUMMARY>|$outsummary|g" | \
  sed "s|<OUTPUTREPORT>|$outreport|g" \
> $cscript
head $cscript

###############################
#TAG: Run CAFError
#
###############################
nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1 > caferror.run.out

###############################
#TAG:  Gather best fitting error models
###############################

cd $runbase
grep "Final error estimates by species:" $runbase/tmp/caferrorLog.txt -A 200 > $runbase/FinalErrorEstimates.txt

# get final estimates for each species
for i in AECH AHEY AINS ACHA ACOL PARG
do
  #grep $i tmp|perl -pe 's/# ([A-Z]{4}).?+([0-9].+$)/$1\t$2/g'
  grep $i $runbase/FinalErrorEstimates.txt |perl -pe 's/# ([A-Z]{4}).+?([0-9].+$)/$1\tcafe_errormodel_$2.txt/g'
done > $runbase/FinalErrorEstimates.perSpecies.txt

# Retrieve all scores for each species
mkdir $runbase/err/
cd $runbase/tmp/
for i in AECH AHEY AINS ACHA ACOL PARG
do
 grep "Search Result"  *$i* -A 1|grep -E "Score:.*"|sed s/cafe_//g|sed s/'_.*Score: '/' '/g > ../err/$i.err
done

cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt

# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s $runbase/FinalErrorEstimates.txt .
ln -s $runbase/FinalErrorEstimates.perSpecies.txt .

# =======================================================================
# Final error estimates by species:
# PARG                             0.129265136719
# ACOL                             0.0
# AINS                             0.0430883789063
# AHEY                             0.0478759765625
# ACHA                             0.0287255859375
# AECH                             0.0526635742187
# Score with individual errors:    27970.041248
# Lambda with individual errors:   0.00048795778972
# =======================================================================
# ************************************
# Score with no errormodel:        30193.178423
# Lambda with no errormodel:       0.00095486984295
# ************************************
# Global Error Estimation:         0.0478759765625
# Score with global errormodel:    28615.713342
# Lambda with global errormodel:   0.00052122001394
# ************************************
# =======================================================================
############################################################################################################################
