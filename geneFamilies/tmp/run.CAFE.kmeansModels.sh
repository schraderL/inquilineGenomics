##########################################################################################################
# SUMMARY : Results of RUN1
##########################################################################################################

:'
Running CAFE as follows:

- MCMCtree ultrametric tree with long chain
- Filtered MCL output (families < 100)
- Errormodels computed for this setup.

'

##########################################################################################################
#TAG: CAFE gene family size evolution: kmeans models
##########################################################################################################

############################################################################################################################
# RUN: 1: Run models with 7 parameters
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN1
mkdir ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans
base=~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/$run
bd=~/data/inqGen18/geneFamilies
mkdir $base
mkdir $base/ctl
cd $base/

##############################################################
#TAG: INPUT files
##############################################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
errorbase=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
##############################################################
#TAG: Setup for subsetted data set (Tsep+leafcutters)
##############################################################

# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# subset tree to include only Tsep and leafcutters
  subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"

##############################################################
#TAG: Prepare CAFE clustering LAMBDA trees
##############################################################

# Folder with kmeans clustered data
cd ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/data
cut -f 1,7,8,16 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.7.txt

cd $base/ctl
VARIABLETREE="(<M01>,(((<M02>,(<M03>,<M04>)<M05>)<M06>,(<M07>,<M08>)<M09>)<M10>,(<M11>,<M12>)<M13>)<M14>);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d ,)
  CLUSTER=$(echo $KMEANS | cut -f2 -d ,)
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s/<$MODEL>/$CLUSTER/g")
done <$base/ctl/kmeans.7.txt

##############################################################
#TAG: Prepare CAFE script RUN1
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript

# section to run different lambdamu models
echo "
#calculate lambdamu
#(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)));
lambdamu -s -t $VARIABLETREE
report MC07
" >> $cafescript
echo $cafescript

##############################################################
# TAG: Run cafe RUN1
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################


###########################################################################################################################
# RUN: 2: Run models with 7 parameters (repetition of RUN1)
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN2
mkdir ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans
base=~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/$run
bd=~/data/inqGen18/geneFamilies
mkdir $base
mkdir $base/ctl
cd $base/

##############################################################
#TAG: INPUT files
##############################################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
errorbase=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
##############################################################
#TAG: Setup for subsetted data set (Tsep+leafcutters)
##############################################################

# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# subset tree to include only Tsep and leafcutters
  subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"

##############################################################
#TAG: Prepare CAFE clustering LAMBDA trees
##############################################################

# Folder with kmeans clustered data
cd ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/data
cut -f 1,7,8,16 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.7.txt

cd $base/ctl
VARIABLETREE="(<M01>,(((<M02>,(<M03>,<M04>)<M05>)<M06>,(<M07>,<M08>)<M09>)<M10>,(<M11>,<M12>)<M13>)<M14>);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d ,)
  CLUSTER=$(echo $KMEANS | cut -f2 -d ,)
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s/<$MODEL>/$CLUSTER/g")
done <$base/ctl/kmeans.7.txt

##############################################################
#TAG: Prepare CAFE script RUN1
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript

# section to run different lambdamu models
echo "
#calculate lambdamu
#(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)));
lambdamu -s -t $VARIABLETREE
report MC07
" >> $cafescript
echo $cafescript

##############################################################
# TAG: Run cafe RUN2
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################

###########################################################################################################################
# RUN: 3: Run models with 8 parameters
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN3
mkdir ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans
base=~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/$run
bd=~/data/inqGen18/geneFamilies
mkdir $base
mkdir $base/ctl
cd $base/

##############################################################
#TAG: INPUT files
##############################################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
errorbase=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
##############################################################
#TAG: Setup for subsetted data set (Tsep+leafcutters)
##############################################################

# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# subset tree to include only Tsep and leafcutters
  subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"

##############################################################
#TAG: Prepare CAFE clustering LAMBDA trees
##############################################################

# Folder with kmeans clustered data
cd ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/data
#cut -f 1,7,8,16 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.7.txt
cut -f 1,7,8,17 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.8.txt

cd $base/ctl
VARIABLETREE="(<M01>,(((<M02>,(<M03>,<M04>)<M05>)<M06>,(<M07>,<M08>)<M09>)<M10>,(<M11>,<M12>)<M13>)<M14>);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d ,)
  CLUSTER=$(echo $KMEANS | cut -f2 -d ,)
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s/<$MODEL>/$CLUSTER/g")
done <$base/ctl/kmeans.8.txt

##############################################################
#TAG: Prepare CAFE script RUN1
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript

# section to run different lambdamu models
echo "
#calculate lambdamu
#(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)));
lambdamu -s -t $VARIABLETREE
report MC08
" >> $cafescript
echo $cafescript

##############################################################
# TAG: Run cafe RUN3
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################

###########################################################################################################################
# RUN: 4: Run models with 8 parameters (repeat RUN3)
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN4
mkdir ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans
base=~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/$run
bd=~/data/inqGen18/geneFamilies
mkdir $base
mkdir $base/ctl
cd $base/

##############################################################
#TAG: INPUT files
##############################################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
errorbase=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
##############################################################
#TAG: Setup for subsetted data set (Tsep+leafcutters)
##############################################################

# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# subset tree to include only Tsep and leafcutters
  subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"

##############################################################
#TAG: Prepare CAFE clustering LAMBDA trees
##############################################################

# Folder with kmeans clustered data
cd ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/data
#cut -f 1,7,8,16 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.7.txt
cut -f 1,7,8,17 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.8.txt

cd $base/ctl
VARIABLETREE="(<M01>,(((<M02>,(<M03>,<M04>)<M05>)<M06>,(<M07>,<M08>)<M09>)<M10>,(<M11>,<M12>)<M13>)<M14>);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d ,)
  CLUSTER=$(echo $KMEANS | cut -f2 -d ,)
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s/<$MODEL>/$CLUSTER/g")
done <$base/ctl/kmeans.8.txt

##############################################################
#TAG: Prepare CAFE script RUN1
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript

# section to run different lambdamu models
echo "
#calculate lambdamu
#(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)));
lambdamu -s -t $VARIABLETREE
report MC08
" >> $cafescript
echo $cafescript

##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################

###########################################################################################################################
# RUN: 5: Run models with 10 parameters
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN5
mkdir ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans
base=~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/$run
bd=~/data/inqGen18/geneFamilies
mkdir $base
mkdir $base/ctl
cd $base/

##############################################################
#TAG: INPUT files
##############################################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
errorbase=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
##############################################################
#TAG: Setup for subsetted data set (Tsep+leafcutters)
##############################################################

# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# subset tree to include only Tsep and leafcutters
  subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"

##############################################################
#TAG: Prepare CAFE clustering LAMBDA trees
##############################################################

# Folder with kmeans clustered data
cd ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/data
#cut -f 1,7,8,16 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.7.txt
cut -f 1,7,8,19 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.8.txt

cd $base/ctl
VARIABLETREE="(<M01>,(((<M02>,(<M03>,<M04>)<M05>)<M06>,(<M07>,<M08>)<M09>)<M10>,(<M11>,<M12>)<M13>)<M14>);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d ,)
  CLUSTER=$(echo $KMEANS | cut -f2 -d ,)
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s/<$MODEL>/$CLUSTER/g")
done <$base/ctl/kmeans.8.txt

##############################################################
#TAG: Prepare CAFE script RUN1
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript

# section to run different lambdamu models
echo "
#calculate lambdamu
#(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)));
lambdamu -s -t $VARIABLETREE
report MC08
" >> $cafescript
echo $cafescript

##############################################################
# TAG: Run cafe RUN5
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################

###########################################################################################################################
# RUN: 6: Run models with 10 parameters (REPEAT RUN 5)
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN6
mkdir ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans
base=~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/$run
bd=~/data/inqGen18/geneFamilies
mkdir $base
mkdir $base/ctl
cd $base/

##############################################################
#TAG: INPUT files
##############################################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
errorbase=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
##############################################################
#TAG: Setup for subsetted data set (Tsep+leafcutters)
##############################################################

# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# subset tree to include only Tsep and leafcutters
  subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"

##############################################################
#TAG: Prepare CAFE clustering LAMBDA trees
##############################################################

# Folder with kmeans clustered data
cd ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/data
#cut -f 1,7,8,16 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.7.txt
cut -f 1,7,8,19 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.10.txt


cd $base/ctl
VARIABLETREE="(<M01>,(((<M02>,(<M03>,<M04>)<M05>)<M06>,(<M07>,<M08>)<M09>)<M10>,(<M11>,<M12>)<M13>)<M14>);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d ,)
  CLUSTER=$(echo $KMEANS | cut -f2 -d ,)
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s/<$MODEL>/$CLUSTER/g")
done <$base/ctl/kmeans.10.txt

##############################################################
#TAG: Prepare CAFE script RUN6
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript

# section to run different lambdamu models
echo "
#calculate lambdamu
#(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)));
lambdamu -s -t $VARIABLETREE
report MC08
" >> $cafescript
echo $cafescript

##############################################################
# TAG: Run cafe RUN3
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################

###########################################################################################################################
# RUN: 7: Run models with 10 parameters (REPEAT RUN 5 and 6)
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN7
mkdir ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans
base=~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/$run
bd=~/data/inqGen18/geneFamilies
mkdir $base
mkdir $base/ctl
cd $base/

##############################################################
#TAG: INPUT files
##############################################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
errorbase=$bd/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
##############################################################
#TAG: Setup for subsetted data set (Tsep+leafcutters)
##############################################################

# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# subset tree to include only Tsep and leafcutters
  subtree=$(echo $treeOriginal |perl -pe 's/\(.*?\(.*?\(.*?(\(.+\)).*?\).*?\).*?\)/$1/g')
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"

##############################################################
#TAG: Prepare CAFE clustering LAMBDA trees
##############################################################

# Folder with kmeans clustered data
cd ~/data/inqGen18/geneFamilies/CAFE/CAFE.kmeans/data
#cut -f 1,7,8,16 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.7.txt
cut -f 1,7,8,19 clustered.rates.CAFE.R1.tsv |sed 1,2d|cut -f 1,4|tr "\t" ","|tr "\n" "\t"|tr "\t" "\n"> $base/ctl/kmeans.10.txt


cd $base/ctl
VARIABLETREE="(<M01>,(((<M02>,(<M03>,<M04>)<M05>)<M06>,(<M07>,<M08>)<M09>)<M10>,(<M11>,<M12>)<M13>)<M14>);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d ,)
  CLUSTER=$(echo $KMEANS | cut -f2 -d ,)
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s/<$MODEL>/$CLUSTER/g")
done <$base/ctl/kmeans.10.txt

##############################################################
#TAG: Prepare CAFE script RUN6
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript

# section to run different lambdamu models
echo "
#calculate lambdamu
#(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)));
lambdamu -s -t $VARIABLETREE
report MC08
" >> $cafescript
echo $cafescript

##############################################################
# TAG: Run cafe RUN3
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################
