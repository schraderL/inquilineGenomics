# ************************************************************************
# SUMMARY: Runs 3 and 4 yielded ~same results
# ************************************************************************
:'
The two runs (named Run 3 and Run 4) on the subsetted dataset (Tsep+leafcutters) yielded very similar results.
After rerunning mcmctree with longer chains, I decided to rerun CAFE errorprediction as well, using as input the tree produced by MCMCtree run4.


#

 paste /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest3.smallSet/tmp/caferrorLog.txt /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest4.smallSet/tmp/caferrorLog.txt /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest5.smallSet/tmp/caferrorLog.txt|tail -n 40

# Final error estimates by species:
# ************************************************************************
# species                          errorTest3    errorTest4           error Test5
# ************************************************************************
# PARG                             0.1290625     0.130859375          0.129265136719
# ACOL                             0.023046875   0.026171875          0.03830078125
# AINS                             0.041484375   0.041875             0.0430883789063
# TSEP                             0.0           5.20417042793e-18    0.0
# AHEY                             0.04609375    0.047109375          0.0478759765625
# ACEP                             0.059921875   0.0628125            0.0622387695313
# ACHA                             0.032265625   0.03140625           0.0335131835937
# AECH                             0.04609375    0.047109375          0.0478759765625
# Score with individual errors:    32397.068767  32396.043792         32412.317419
# Lambda with individual errors:   0.0003796688  0.0003768867         0.00036820716549
# ************************************************************************
# Score with no errormodel:        35262.219287 35262.217702          35294.582562
# Lambda with no errormodel:       0.0007641781 0.0007641845          0.00075852947559
# ************************************************************************
# Global Error Estimation:         0.04609375   0.05234375            0.0478759765625
# Score with global errormodel:    33051.31241  33063.771311          33055.56009
# Lambda with global errormodel:   0.0003928821 0.00037832964         0.00038526594786
'

############################################
##CAFE gene family size evolution: Error rate estimation
############################################
base=/usr/local/home/lschrader/data/inqGen18/geneFamilies
bd=/usr/local/home/lschrader/data/inqGen18
software=~/data/software

mkdir $base/CAFE
mkdir $base/CAFE/CAFE.errorprediction
mkdir $base/CAFE/CAFE.errorprediction/ctl

###############################
#TAG: Prepare template CAFError script
###############################

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
> /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe

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
clusters=$base/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run1.tre

############################################
# Prepare Environment
############################################
base=/usr/local/home/lschrader/data/inqGen18/geneFamilies
bd=/usr/local/home/lschrader/data/inqGen18
software=~/data/software

###############################
#TAG: Prepare subsetted MCL cluster file
###############################

# input files
originalClusters=$base/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run1.tre

# Prepare Cluster file containing only Tsep and the leafcutters
# In R:
:'
a<-read.csv("/usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv",T,sep="\t")
b<-a[,c(1:8,10,12)]
b$Sum<-rowSums(b[,3:10])
set<-subset(b,Sum!=0)
set<-set[,-11]
dim(set)
write.table(set,"/usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv",quote=F,row.names=F,sep="\t",)
'

###############################
#TAG: Filter genes with >100 copies in one species from Cluster file
###############################
cd $base
rawclusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
cat $rawclusters |awk  -F $'\t' 'BEGIN {OFS=FS} {if ($3 <100 || $4 <100|| $5 <100|| $6 <100|| $7 <100|| $8 <100|| $9 <100|| $10 <100) print $0}'

############################################################################################################################
###############################
#RUN: 1: Run errorTest3.smallSet
###############################

runName=errorTest3.smallSet
###############################
  runbase=$base/CAFE/CAFE.errorprediction/$runName
  cscript=$runbase/$runName.cafe
  outsummary=$runbase/$runName.out
  outreport=$runbase/report/$runName.report.out
  clusters=$base/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv

# prepare environment
mkdir $runbase
mkdir $runbase/tmp
cd $runbase

## Run the error prediction script, predicting error rates for each species
# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# subset tree to include only Tsep and leafcutters
  subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

# Prepare cafe script from template.
cat $base/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe | \
  sed "s|<CLUSTERS>|$clusters|g" | \
  sed "s|<TREE>|$tree|g" | \
  sed "s|<LAMBDATREE>|$LAMBDATREE|g" | \
  sed "s|<OUTPUTSUMMARY>|$outsummary|g" | \
  sed "s|<OUTPUTREPORT>|$outreport|g" \
> $cscript


# Run CAFerror script to estimate error rates.
nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1

# Retrieve error model likelihoods for each branch
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
 grep "Search Result"  *$i* -A 1|grep -E "Score:.*"|sed s/cafe_//g|sed s/'_.*Score: '/' '/g > ../err/$i.err
done


# Copy best error model files to folder bestFitModels
cd $runbase/err/
mkdir $runbase/bestFitModels
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  fit=$(cat $i.err |sort -k 2,2|cut -f 1 -d " "|head -n1)
  #ls "../tmp/cafe_errormodel_"$fit".txt"
  cp "../tmp/cafe_errormodel_"$fit".txt"  $runbase"/bestFitModels/cafe_errormodel_"$fit"_"$i".txt"
done

cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt

# =======================================================================
# Final error estimates by species:
# PARG                             0.1290625
# ACOL                             0.023046875
# AINS                             0.041484375
# TSEP                             0.0
# AHEY                             0.04609375
# ACEP                             0.059921875
# ACHA                             0.032265625
# AECH                             0.04609375
# Score with individual errors:    32397.068767
# Lambda with individual errors:   0.00037966884999
# =======================================================================
# ************************************
# Score with no errormodel:        35262.219287
# Lambda with no errormodel:       0.00076417812417
# ************************************
# Global Error Estimation:         0.04609375
# Score with global errormodel:    33051.31241
# Lambda with global errormodel:   0.00039288212261

# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s ../bestFitModels/ .


############################################################################################################################

###############################
#RUN: 2: Run errorTest4.smallSet
###############################

# input files
#clusters=$base/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run1.tre


# Prepare Cluster file containing only Tsep and the leafcutters
# In R:
:'
a<-read.csv("/usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv",T,sep="\t")
b<-a[,c(1:8,10,12)]
b$Sum<-rowSums(b[,3:10])
set<-subset(b,Sum!=0)
set<-set[,-11]
dim(set)
write.table(set,"/usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv",quote=F,row.names=F,sep="\t",)
'

clusters=$base/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv

runName=errorTest4.smallSet
###############################
  runbase=$base/CAFE/CAFE.errorprediction/$runName
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
# retrieve tree from MCMCtree output
treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526
tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

# Prepare cafe script from template.
cat $base/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe | \
  sed "s|<CLUSTERS>|$clusters|g" | \
  sed "s|<TREE>|$tree|g" | \
  sed "s|<LAMBDATREE>|$LAMBDATREE|g" | \
  sed "s|<OUTPUTSUMMARY>|$outsummary|g" | \
  sed "s|<OUTPUTREPORT>|$outreport|g" \
> $cscript

# Run CAFError
nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1 > caferror.run.out

# Retrieve all scores for each species
cd $runbase/tmp/
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
 grep "Search Result"  *$i* -A 1|grep -E "Score:.*"|sed s/cafe_//g|sed s/'_.*Score: '/' '/g > ../err/$i.err
done


# Best error models
cd $runbase/err/
mkdir $runbase/bestFitModels
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  fit=$(cat $i.err |sort -k 2,2|cut -f 1 -d " "|head -n1)
#  ls "../tmp/cafe_errormodel_"$fit".txt"
  cp "../tmp/cafe_errormodel_"$fit".txt"  $runbase"/bestFitModels/cafe_errormodel_"$fit"_"$i".txt"
done

cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt

# Final error estimates by species:
# PARG                             0.130859375
# ACOL                             0.026171875
# AINS                             0.041875
# TSEP                             5.20417042793e-18
# AHEY                             0.047109375
# ACEP                             0.0628125
# ACHA                             0.03140625
# AECH                             0.047109375
# Score with individual errors:    32396.043792
# Lambda with individual errors:   0.00037688672535
# =======================================================================
# ************************************
# Score with no errormodel:        35262.217702
# Lambda with no errormodel:       0.00076418450812
# ************************************
# Global Error Estimation:         0.05234375
# Score with global errormodel:    33063.771311
# Lambda with global errormodel:   0.00037832964566
# ************************************

# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s ../bestFitModels/ .


############################################################################################################################

############################################################################################################################

###############################
#RUN: 3: Run errorTest5.smallSet
###############################
cd $base
# input files
#clusters=$base/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
ls $treeFile

# Prepare Cluster file containing only Tsep and the leafcutters
# In R:
:'
a<-read.csv("/usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv",T,sep="\t")
b<-a[,c(1:8,10,12)]
b$Sum<-rowSums(b[,3:10])
set<-subset(b,Sum!=0)
set<-set[,-11]
dim(set)
write.table(set,"/usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv",quote=F,row.names=F,sep="\t",)
'

clusters=$base/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv

runName=errorTest5.smallSet
###############################
  runbase=$base/CAFE/CAFE.errorprediction/$runName
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
# retrieve tree from MCMCtree output
treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526
tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
echo $tree
# Prepare cafe script from template.
cat $base/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe | \
  sed "s|<CLUSTERS>|$clusters|g" | \
  sed "s|<TREE>|$tree|g" | \
  sed "s|<LAMBDATREE>|$LAMBDATREE|g" | \
  sed "s|<OUTPUTSUMMARY>|$outsummary|g" | \
  sed "s|<OUTPUTREPORT>|$outreport|g" \
> $cscript
less $cscript
# Run CAFError
nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1 > caferror.run.out

# Retrieve all scores for each species
cd $runbase/tmp/
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
 grep "Search Result"  *$i* -A 1|grep -E "Score:.*"|sed s/cafe_//g|sed s/'_.*Score: '/' '/g > ../err/$i.err
done


# Best error models
cd $runbase/err/
mkdir $runbase/bestFitModels
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  fit=$(cat $i.err |sort -k 2,2|cut -f 1 -d " "|head -n1)
#  ls "../tmp/cafe_errormodel_"$fit".txt"
  cp "../tmp/cafe_errormodel_"$fit".txt"  $runbase"/bestFitModels/cafe_errormodel_"$fit"_"$i".txt"
done

cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt

# Final error estimates by species:
# PARG                             0.130859375
# ACOL                             0.026171875
# AINS                             0.041875
# TSEP                             5.20417042793e-18
# AHEY                             0.047109375
# ACEP                             0.0628125
# ACHA                             0.03140625
# AECH                             0.047109375
# Score with individual errors:    32396.043792
# Lambda with individual errors:   0.00037688672535
# =======================================================================
# ************************************
# Score with no errormodel:        35262.217702
# Lambda with no errormodel:       0.00076418450812
# ************************************
# Global Error Estimation:         0.05234375
# Score with global errormodel:    33063.771311
# Lambda with global errormodel:   0.00037832964566
# ************************************

# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s ../bestFitModels/ .


############################################################################################################################

############################################################################################################################
############################################
#TAG: OUTPUT files
############################################

tree /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest*.smallSet/results/

:'
/usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest3.smallSet/results/
├── bestFitModels -> ../bestFitModels/
└── caferrorLog.txt -> ../tmp/caferrorLog.txt
/usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest4.smallSet/results/
├── bestFitModels -> ../bestFitModels/
└── caferrorLog.txt -> ../tmp/caferrorLog.txt
/usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest4.smallSet/results/
'
############################################################################################################################
