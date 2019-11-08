##########################################################################################################
# SUMMARY
##########################################################################################################

:'
Running CAFE as follows:

- MCMCtree ultrametric tree with long chain
- Filtered MCL output (families < 100)
- Errormodels computed for this setup.
- Clustered lambda and mu based on clustering of RUN1 /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/RUN1/treeRates.RUN1.tsv
'

##########################################################################################################
#TAG: CAFE gene family size evolution: kmeans models
##########################################################################################################

# upload clustered rates
scp -o ProxyCommand="ssh -W %h:%p zcm795@ssh-bio-labnet.science.ku.dk" /Users/lukas/sciebo/inquilineGenomics18/geneFamilies/CAFE/results/Acro.Atta/clustered.rates.CAFE.R1.tsv  lschrader@pallas.bio.ku.dk:~/data/inqGen18/geneFamilies/CAFE/
scp -o ProxyCommand="ssh -W %h:%p zcm795@ssh-bio-labnet.science.ku.dk" /Users/lukas/sciebo/inquilineGenomics18/geneFamilies/CAFE/results/Acro.Atta/clustered.rates.CAFE.RUN5.bestFits.tsv  lschrader@pallas.bio.ku.dk:~/data/inqGen18/geneFamilies/CAFE/

############################################################################################################################
# RUN: 1: Run models with 5 parameters
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN1C
base=~/data/inqGen18/geneFamilies/CAFE/$run
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
~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL PARG
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
cut -f 8,13 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.R1.tsv |sed 1d > $base/ctl/kmeans.5.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done <$base/ctl/kmeans.5.txt

##############################################################
#TAG: Prepare CAFE script RUN1
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads
# run with 5 clusters
tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript
# section to run different lambdamu models
echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $VARIABLETREE
report m5a

lambdamu -s -t $VARIABLETREE
report m5b


" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run


cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

##############################################################
# TAG: Gather results
##############################################################
#(((AHEY1,(ACHA2,PARG3)4)5,(AECH6,AINS7)8)9,ACOL10);

grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run

############################################################################################################################


############################################################################################################################
# RUN: 2: Run models with 5 parameters, multiple repetitions
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN2C
base=~/data/inqGen18/geneFamilies/CAFE/$run
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
~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL PARG
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
cut -f 8,13 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.R1.tsv |sed 1d > $base/ctl/kmeans.5.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done <$base/ctl/kmeans.5.txt

##############################################################
#TAG: Prepare CAFE script for run
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads
# run with 5 clusters
tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript
# section to run different lambdamu models
echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $VARIABLETREE
report m5a

lambdamu -s -t $VARIABLETREE
report m5b

lambdamu -s -t $VARIABLETREE
report m5c

lambdamu -s -t $VARIABLETREE
report m5d

lambdamu -s -t $VARIABLETREE
report m5e

lambdamu -s -t $VARIABLETREE
report m5f

" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run


cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

##############################################################
# TAG: Gather results
##############################################################
#(((AHEY1,(ACHA2,PARG3)4)5,(AECH6,AINS7)8)9,ACOL10);

grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run

############################################################################################################################


############################################################################################################################
# RUN: 3: Run models with 5 parameters, multiple repetitions
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN3C
base=~/data/inqGen18/geneFamilies/CAFE/$run
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
~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL PARG
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
cut -f 8,13 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.R1.tsv |sed 1d > $base/ctl/kmeans.5.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done <$base/ctl/kmeans.5.txt

##############################################################
#TAG: Prepare CAFE script for run
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads
# run with 5 clusters
tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript
# section to run different lambdamu models
echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $VARIABLETREE
report m5a

lambdamu -s -t $VARIABLETREE
report m5b

lambdamu -s -t $VARIABLETREE
report m5c

lambdamu -s -t $VARIABLETREE
report m5d

lambdamu -s -t $VARIABLETREE
report m5e

lambdamu -s -t $VARIABLETREE
report m5f

lambdamu -s -t $VARIABLETREE
report m5g

lambdamu -s -t $VARIABLETREE
report m5h

lambdamu -s -t $VARIABLETREE
report m5i


" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run


cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

##############################################################
# TAG: Gather results
##############################################################
#(((AHEY1,(ACHA2,PARG3)4)5,(AECH6,AINS7)8)9,ACOL10);

grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run

############################################################################################################################

############################################################################################################################
# RUN: 4C: Run models with 6 parameters, multiple repetitions
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN4C
base=~/data/inqGen18/geneFamilies/CAFE/$run
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
~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL PARG
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
# -f 14 = 6 clusters
cut -f 8,14 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.R1.tsv |sed 1d > $base/ctl/kmeans.6.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done <$base/ctl/kmeans.6.txt

##############################################################
#TAG: Prepare CAFE script for run
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads
# run with 6 clusters
tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript
# section to run different lambdamu models
echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $VARIABLETREE
report m6a

lambdamu -s -t $VARIABLETREE
report m6b

lambdamu -s -t $VARIABLETREE
report m6c

lambdamu -s -t $VARIABLETREE
report m6d

lambdamu -s -t $VARIABLETREE
report m6e

lambdamu -s -t $VARIABLETREE
report m6f

lambdamu -s -t $VARIABLETREE
report m6g

lambdamu -s -t $VARIABLETREE
report m6h

lambdamu -s -t $VARIABLETREE
report m6i


" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run


cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

##############################################################
# TAG: Gather results
##############################################################
#(((AHEY1,(ACHA2,PARG3)4)5,(AECH6,AINS7)8)9,ACOL10);
#(((2,(3,6)4)5,(1,4)5)5,5);

grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run

############################################################################################################################

############################################################################################################################
# RUN: 5C: Run models with 6 parameters, multiple repetitions, use clustered.rates.CAFE.RUN5.bestFits.tsv as input clusters
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN5C
base=~/data/inqGen18/geneFamilies/CAFE/$run
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
~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL PARG
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
# -f 1,18 = labels & 6 clusters
cut -f 1,18 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.RUN5.bestFits.tsv |sed 1d > $base/ctl/kmeans.6.tmp.txt
# sort by length of node name, so that replacement works properly
awk -vOFS='\t' '{ print length($1), $0 }' $base/ctl/kmeans.6.tmp.txt | sort -r -k1,1n -k2,2 | cut -f2-|tac > $base/ctl/kmeans.6.txt
rm $base/ctl/kmeans.6.tmp.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done < $base/ctl/kmeans.6.txt

##############################################################
#TAG: Prepare CAFE script for run
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads
# run with 6 clusters
tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript
# section to run different lambdamu models
echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $VARIABLETREE
report m6a

lambdamu -s -t $VARIABLETREE
report m6b

lambdamu -s -t $VARIABLETREE
report m6c

lambdamu -s -t $VARIABLETREE
report m6d

lambdamu -s -t $VARIABLETREE
report m6e

lambdamu -s -t $VARIABLETREE
report m6f

lambdamu -s -t $VARIABLETREE
report m6g

lambdamu -s -t $VARIABLETREE
report m6h

lambdamu -s -t $VARIABLETREE
report m6i


" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run


cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

##############################################################
# TAG: Gather results
##############################################################
#(((AHEY1,(ACHA2,PARG3)4)5,(AECH6,AINS7)8)9,ACOL10);
#(((2,(3,6)4)5,(1,4)5)5,5);

grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run

############################################################################################################################

############################################################################################################################
# RUN: 6C: Run models with 4 parameters, multiple repetitions, use clustered.rates.CAFE.RUN5.bestFits.tsv as input clusters
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN6C
base=~/data/inqGen18/geneFamilies/CAFE/$run
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
~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL PARG
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
# -f 1,18 = labels & 4 clusters
cut -f 1,16 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.RUN5.bestFits.tsv |sed 1d > $base/ctl/kmeans.tmp.txt
# sort by length of node name, so that replacement works properly
awk -vOFS='\t' '{ print length($1), $0 }' $base/ctl/kmeans.tmp.txt | sort -r -k1,1n -k2,2 | cut -f2-|tac > $base/ctl/kmeans.4.txt
rm $base/ctl/kmeans.tmp.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done < $base/ctl/kmeans.4.txt

##############################################################
#TAG: Prepare CAFE script for run
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads
# run with 6 clusters
tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript
# section to run different lambdamu models
echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $VARIABLETREE
report m4a

lambdamu -s -t $VARIABLETREE
report m4b

lambdamu -s -t $VARIABLETREE
report m4c

lambdamu -s -t $VARIABLETREE
report m4d

lambdamu -s -t $VARIABLETREE
report m4e

lambdamu -s -t $VARIABLETREE
report m4f

lambdamu -s -t $VARIABLETREE
report m4g

lambdamu -s -t $VARIABLETREE
report m4h

lambdamu -s -t $VARIABLETREE
report m4i


" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run


cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

##############################################################
# TAG: Gather results
##############################################################
#(((AHEY1,(ACHA2,PARG3)4)5,(AECH6,AINS7)8)9,ACOL10);
#(((2,(3,6)4)5,(1,4)5)5,5);

grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run

############################################################################################################################

############################################################################################################################
# RUN: 7C: Run models with 7 parameters, multiple repetitions, use clustered.rates.CAFE.RUN5.bestFits.tsv as input clusters
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN7C
base=~/data/inqGen18/geneFamilies/CAFE/$run
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
~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL PARG
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
# -f 1,18 = labels & 7 clusters
cut -f 1,19 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.RUN5.bestFits.tsv |sed 1d > $base/ctl/kmeans.tmp.txt
# sort by length of node name, so that replacement works properly
awk -vOFS='\t' '{ print length($1), $0 }' $base/ctl/kmeans.tmp.txt | sort -r -k1,1n -k2,2 | cut -f2-|tac > $base/ctl/kmeans.7.txt
rm $base/ctl/kmeans.tmp.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done < $base/ctl/kmeans.7.txt

##############################################################
#TAG: Prepare CAFE script for run
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads
# run with 6 clusters
tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript
# section to run different lambdamu models
echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $VARIABLETREE
report m7a

lambdamu -s -t $VARIABLETREE
report m7b

lambdamu -s -t $VARIABLETREE
report m7c

lambdamu -s -t $VARIABLETREE
report m7d

lambdamu -s -t $VARIABLETREE
report m7e

lambdamu -s -t $VARIABLETREE
report m7f

lambdamu -s -t $VARIABLETREE
report m7g

lambdamu -s -t $VARIABLETREE
report m7h

lambdamu -s -t $VARIABLETREE
report m7i


" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run


cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

##############################################################
# TAG: Gather results
##############################################################
#(((AHEY1,(ACHA2,PARG3)4)5,(AECH6,AINS7)8)9,ACOL10);
#(((2,(3,6)4)5,(1,4)5)5,5);

grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run

############################################################################################################################

############################################################################################################################
# RUN: 8C: Run manually curated models with 6 parameters, multiple repetitions
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN8C
base=~/data/inqGen18/geneFamilies/CAFE/$run
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
~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL PARG
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
# -f 14 = 6 clusters
cut -f 8,14 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.R1.tsv |sed 1d > $base/ctl/kmeans.6.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done <$base/ctl/kmeans.6.txt

# create tree with
## cluster AECH [1]
## cluster ACOL [2]
## cluster AHEY & ACHA [3]
## cluster AINS & ACHA/PARG [4]
## cluster AHEY/PARG/ACHA,ACRO,AECH/AINS [5]
## cluster PARG [6]
VARIABLETREE="(((3,(3,6)4)5,(1,4)5)5,2);"

##############################################################
#TAG: Prepare CAFE script for run
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads
# run with 6 clusters
tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript
# section to run different lambdamu models
echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $VARIABLETREE
report m6a

lambdamu -s -t $VARIABLETREE
report m6b

lambdamu -s -t $VARIABLETREE
report m6c

lambdamu -s -t $VARIABLETREE
report m6d

lambdamu -s -t $VARIABLETREE
report m6e

lambdamu -s -t $VARIABLETREE
report m6f

lambdamu -s -t $VARIABLETREE
report m6g

lambdamu -s -t $VARIABLETREE
report m6h

lambdamu -s -t $VARIABLETREE
report m6i


" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run


cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

##############################################################
# TAG: Gather results
##############################################################

grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run

############################################################################################################################


############################################################################################################################
# RUN: 9C: Run manually curated models with 6 parameters, multiple repetitions
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN9C
base=~/data/inqGen18/geneFamilies/CAFE/$run
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
~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
## Run 5, run on MCMCtree run4
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/results/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=~/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv
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

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')

##############################################################
#TAG: Retrieve error models
##############################################################

for i in AECH AHEY AINS ACHA ACOL PARG
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
# -f 14 = 6 clusters
cut -f 8,14 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.R1.tsv |sed 1d > $base/ctl/kmeans.6.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done <$base/ctl/kmeans.6.txt

# create tree with
## cluster AECH [1]
## cluster AINS[2]
## cluster AHEY & ACHA [3]
## cluster ACHA/PARG [4]
## cluster AHEY/PARG/ACHA,ACRO,AECH/AINS,ACOL [5]
## cluster PARG [6]
VARIABLETREE="(((3,(3,6)4)5,(1,2)5)5,5);"

##############################################################
#TAG: Prepare CAFE script for run
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads
# run with 6 clusters
tree $tree
load -i $clusters -t 20 -l $out -p 0.05 -filter

#assign error models
" >  $cafescript

# section to load error models
cat $base/ctl/errorModels.txt >> $cafescript
# section to run different lambdamu models
echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $VARIABLETREE
report m6a

lambdamu -s -t $VARIABLETREE
report m6b

lambdamu -s -t $VARIABLETREE
report m6c

lambdamu -s -t $VARIABLETREE
report m6d

lambdamu -s -t $VARIABLETREE
report m6e

lambdamu -s -t $VARIABLETREE
report m6f

lambdamu -s -t $VARIABLETREE
report m6g

lambdamu -s -t $VARIABLETREE
report m6h

lambdamu -s -t $VARIABLETREE
report m6i


" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run


cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

##############################################################
# TAG: Gather results
##############################################################

grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run

############################################################################################################################
