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
#TAG: CAFE gene family size evolution: 2-rate models
##########################################################################################################

############################################################################################################################
# RUN: 1: Run null model all two parameter models, updated MCMCtree phylogeny and updated error estimates
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN1
base=~/data/inqGen18/geneFamilies/CAFE/CAFE/$run
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
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(1,1)1)1);
#(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)));
report null
lambdamu -s -t (2,(((1,(1,1)1)1,(1,1)1)1,(1,1)1)1);
report r01
lambdamu -s -t (1,(((2,(1,1)1)1,(1,1)1)1,(1,1)1)1);
report r02
lambdamu -s -t (1,(((1,(2,1)1)1,(1,1)1)1,(1,1)1)1);
report r03
lambdamu -s -t (1,(((1,(1,2)1)1,(1,1)1)1,(1,1)1)1);
report r04
lambdamu -s -t (1,(((1,(1,1)2)1,(1,1)1)1,(1,1)1)1);
report r05
lambdamu -s -t (1,(((1,(1,1)1)2,(1,1)1)1,(1,1)1)1);
report r06
lambdamu -s -t (1,(((1,(1,1)1)1,(2,1)1)1,(1,1)1)1);
report r07
lambdamu -s -t (1,(((1,(1,1)1)1,(1,2)1)1,(1,1)1)1);
report r08
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)2)1,(1,1)1)1);
report r09
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)2,(1,1)1)1);
report r10
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(2,1)1)1);
report r11
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(1,2)1)1);
report r12
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(1,1)2)1);
report r13
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(1,1)1)2);
report r14
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

############################################################################################################################
# RUN: 2: Run null model all two parameter models, updated MCMCtree phylogeny and updated error estimates (replicate RUN1)
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN2
base=~/data/inqGen18/geneFamilies/CAFE/CAFE/$run
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
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(1,1)1)1);
#(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)));
report null
lambdamu -s -t (2,(((1,(1,1)1)1,(1,1)1)1,(1,1)1)1);
report r01
lambdamu -s -t (1,(((2,(1,1)1)1,(1,1)1)1,(1,1)1)1);
report r02
lambdamu -s -t (1,(((1,(2,1)1)1,(1,1)1)1,(1,1)1)1);
report r03
lambdamu -s -t (1,(((1,(1,2)1)1,(1,1)1)1,(1,1)1)1);
report r04
lambdamu -s -t (1,(((1,(1,1)2)1,(1,1)1)1,(1,1)1)1);
report r05
lambdamu -s -t (1,(((1,(1,1)1)2,(1,1)1)1,(1,1)1)1);
report r06
lambdamu -s -t (1,(((1,(1,1)1)1,(2,1)1)1,(1,1)1)1);
report r07
lambdamu -s -t (1,(((1,(1,1)1)1,(1,2)1)1,(1,1)1)1);
report r08
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)2)1,(1,1)1)1);
report r09
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)2,(1,1)1)1);
report r10
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(2,1)1)1);
report r11
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(1,2)1)1);
report r12
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(1,1)2)1);
report r13
lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(1,1)1)2);
report r14
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
