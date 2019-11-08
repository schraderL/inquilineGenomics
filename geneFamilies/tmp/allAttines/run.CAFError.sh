# ************************************************************************
# SUMMARY: Runs 3 and 4 yielded ~same results
# ************************************************************************
:'
The two runs (named Run 3 and Run 4) on the subsetted dataset (Tsep+leafcutters) yielded very similar results.
After rerunning mcmctree with longer chains, I decided to rerun CAFE errorprediction as well, using as input the tree produced by MCMCtree run4.


#
paste /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100.replicate/tmp/caferrorLog.txt \
/usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp/caferrorLog.txt |tail -n 40

# Final error estimates by species:     # Final error estimates by species:
# PARG                             0.128737792969       # PARG                             0.128737792969
# ACOL                             0.0238403320312      # ACOL                             0.0238403320312
# AINS                             0.0429125976562      # AINS                             0.0429125976562
# TSEP                             0.0  # TSEP                             0.0
# AHEY                             0.0476806640625      # AHEY                             0.0476806640625
# ACEP                             0.0619848632812      # ACEP                             0.0619848632812
# ACHA                             0.0333764648437      # ACHA                             0.0333764648437
# AECH                             0.0476806640625      # AECH                             0.0476806640625
# Score with individual errors:    31909.320767 # Score with individual errors:    31909.319042
# Lambda with individual errors:   0.00035527061921     # Lambda with individual errors:   0.00035527214121
# =======================================================================       # =======================================================================
# ************************************  # ************************************
# Score with no errormodel:        34888.448602 # Score with no errormodel:        34888.448932
# Lambda with no errormodel:       0.00074289397613     # Lambda with no errormodel:       0.00074290996097
# ************************************  # ************************************
# Global Error Estimation:         0.0476806640625      # Global Error Estimation:         0.0476806640625
# Score with global errormodel:    32567.729975 # Score with global errormodel:    32567.727621
# Lambda with global errormodel:   0.00036497273575     # Lambda with global errormodel:   0.00036497083717
'



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

###############################
#TAG: Prepare subsetted MCL cluster file
###############################

# input files
originalClusters=$base/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

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
cat $rawclusters |awk  -F $'\t' 'BEGIN {OFS=FS} {if (($3 <100 && $4 <100 && $5 <100 && $6 <100 && $7 <100 && $8 <100 && $9 <100 && $10 <100) || $1=="Description") print $0}' \
> ~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv

############################################################################################################################

###############################
#RUN: 1: Run errorTest.smallSet.filter100
###############################
cd $base

############################################
#TAG: Prepare Environment
############################################
base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction
bd=/usr/local/home/lschrader/data/inqGen18
software=~/data/software

mkdir $base/
mkdir $base/ctl


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
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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
subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
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
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  #grep $i tmp|perl -pe 's/# ([A-Z]{4}).?+([0-9].+$)/$1\t$2/g'
  grep $i $runbase/FinalErrorEstimates.txt |perl -pe 's/# ([A-Z]{4}).+?([0-9].+$)/$1\tcafe_errormodel_$2.txt/g'
done > $runbase/FinalErrorEstimates.perSpecies.txt

# Retrieve all scores for each species
mkdir $runbase/err/
cd $runbase/tmp/
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
 grep "Search Result"  *$i* -A 1|grep -E "Score:.*"|sed s/cafe_//g|sed s/'_.*Score: '/' '/g > ../err/$i.err
done

cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt

# Final error estimates by species:
# PARG                             0.128737792969
# ACOL                             0.0238403320312
# AINS                             0.0429125976562
# TSEP                             0.0
# AHEY                             0.0476806640625
# ACEP                             0.0619848632812
# ACHA                             0.0333764648437
# AECH                             0.0476806640625
# Score with individual errors:    31909.319042
# Lambda with individual errors:   0.00035527214121
# =======================================================================
# ************************************
# Score with no errormodel:        34888.448932
# Lambda with no errormodel:       0.00074290996097
# ************************************
# Global Error Estimation:         0.0476806640625
# Score with global errormodel:    32567.727621
# Lambda with global errormodel:   0.00036497083717
# ************************************

# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s $runbase/FinalErrorEstimates.txt .
ln -s $runbase/FinalErrorEstimates.perSpecies.txt .


############################################################################################################################

###############################
#RUN: 2: Run errorTest.smallSet.filter100 (repetition of RUN1)
###############################
############################################
#TAG: Prepare Environment
############################################
base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction
bd=/usr/local/home/lschrader/data/inqGen18
software=~/data/software

cd $base
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
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
ls $treeFile

###############################
# TAG: Define runName
###############################
runName=errorTest.smallSet.filter100.replicate

  runbase=$base/$runName
  cscript=$runbase/$runName.cafe
  outsummary=$runbase/$runName.out
  outreport=$runbase/report/$runName.report.out

# prepare environment
mkdir $runbase
mkdir $runbase/tmp
mkdir $runbase/err
cd $runbase

## Run the error prediction script, predicting error rates for each species

###############################
# TAG: prepare tree from MCMCtree output
###############################
treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
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
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  #grep $i tmp|perl -pe 's/# ([A-Z]{4}).?+([0-9].+$)/$1\t$2/g'
  grep $i $runbase/FinalErrorEstimates.txt |perl -pe 's/# ([A-Z]{4}).+?([0-9].+$)/$1\tcafe_errormodel_$2.txt/g'
done > $runbase/FinalErrorEstimates.perSpecies.txt

# Retrieve all scores for each species
mkdir $runbase/err/
cd $runbase/tmp/
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
 grep "Search Result"  *$i* -A 1|grep -E "Score:.*"|sed s/cafe_//g|sed s/'_.*Score: '/' '/g > ../err/$i.err
done

cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt

# Final error estimates by species:
# PARG                             0.128737792969
# ACOL                             0.0238403320312
# AINS                             0.0429125976562
# TSEP                             0.0
# AHEY                             0.0476806640625
# ACEP                             0.0619848632812
# ACHA                             0.0333764648437
# AECH                             0.0476806640625
# Score with individual errors:    31909.320767
# Lambda with individual errors:   0.00035527061921
# =======================================================================
# ************************************
# Score with no errormodel:        34888.448602
# Lambda with no errormodel:       0.00074289397613
# ************************************
# Global Error Estimation:         0.0476806640625
# Score with global errormodel:    32567.729975
# Lambda with global errormodel:   0.00036497273575

# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s $runbase/FinalErrorEstimates.txt .
ln -s $runbase/FinalErrorEstimates.perSpecies.txt .


############################################################################################################################


############################################################################################################################
############################################
#TAG: OUTPUT files
############################################

tree ~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100.replicate/results/

:'
'
############################################################################################################################
