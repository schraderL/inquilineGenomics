############################################################################################################################
# RUN: 1: Run models with 5 parameters, different groupings
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN1DC
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
# "classic" clusters
lambdamu -s -t (((3,(3,1)5)2,(4,5)2)2,2);
#report m5a

# ACHA/PARG and AINS in separate clusters
# ACHA, AHEY, AECH in one
lambdamu -s -t (((3,(3,1)4)2,(3,5)2)2,2);
lambdamu -s -t (((3,(3,1)4)2,(3,5)2)2,2);
lambdamu -s -t (((3,(3,1)4)2,(3,5)2)2,2);
lambdamu -s -t (((3,(3,1)4)2,(3,5)2)2,2);
report m5b

# ACHA/PARG and AINS in one
# ACHA single
# AHEY and AECH in one
lambdamu -s -t (((3,(4,1)5)2,(3,5)2)2,2);
lambdamu -s -t (((3,(4,1)5)2,(3,5)2)2,2);
lambdamu -s -t (((3,(4,1)5)2,(3,5)2)2,2);
report m5c

# ACHA/PARG and AINS in one
# ACHA, AHEY and AECH in one
# ACRO split from main
lambdamu -s -t (((3,(3,1)5)2,(3,5)2)4,2);
lambdamu -s -t (((3,(3,1)5)2,(3,5)2)4,2);
lambdamu -s -t (((3,(3,1)5)2,(3,5)2)4,2);
report m5d




" >> $cafescript


echo $cafescript

##############################################################
# TAG: Run cafe
##############################################################
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
