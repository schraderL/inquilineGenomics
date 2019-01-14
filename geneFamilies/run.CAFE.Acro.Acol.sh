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

mkdir ~/data/inqGen18/geneFamilies/CAFE/
run=RUN1
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
#TAG: Setup for subsetted data set (Acro + Acol)
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
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $LAMBDATREE
report null
" >> $cafescript

#echo $LAMBDATREE|grep -o 1|wc -l
for i in {1..10}
do
  printf '%s\n' "report r$i"
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
done  >> $cafescript

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

############################################################################################################################


############################################################################################################################
# RUN: 2: Run null model, full model, 4-parameter model
# updated MCMCtree phylogeny and updated error estimates
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################

mkdir ~/data/inqGen18/geneFamilies/CAFE/
run=RUN2
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
#TAG: Setup for subsetted data set (Acro + Acol)
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
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $LAMBDATREE
report null
" >> $cafescript


FULLTREE=$LAMBDATREE
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {1..10}
do
  FULLTREE=$(echo $FULLTREE|sed -e "s/1/$i/2")
done

echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for FULLTREE
lambdamu -s -t $FULLTREE
report full
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

############################################################################################################################
# RUN: 3: Run null model, all two parameter models, full model
# updated MCMCtree phylogeny and updated error estimates
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################

mkdir ~/data/inqGen18/geneFamilies/CAFE/
run=RUN3
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
#TAG: Setup for subsetted data set (Acro + Acol)
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
#TAG: Prepare CAFE script RUN3
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
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $LAMBDATREE
report null
" >> $cafescript


FULLTREE=$LAMBDATREE
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {1..10}
do
  FULLTREE=$(echo $FULLTREE|sed -e "s/1/$i/2")
done

echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for FULLTREE
lambdamu -s -t $FULLTREE
report full
" >> $cafescript

#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
CLUSTERTREE="(((4,(4,2)3)1,(4,3)1)1,1);"

echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for clusters
lambdamu -s -t $CLUSTERTREE
report cluster01
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
mv $base/full.cafe $base/$run
mv $base/cluster01.cafe $base/$run

: '
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
(((4,(4,2)3)1,(4,3)1)1,1);
0.00036173817362  1 internal branches
0.00029090061278  2 PARG
0.00034671847156  3 PARG/ACHA,AINS
0.00131463284504  4 AECH,AHEY,ACHA
& Score: 26377.1
Mu :
0.00015588206004  1 internal branches
0.00707610535453  2 PARG
0.00251104450323  3 PARG/ACHA,AINS
0.00112650941324  4 AECH,AHEY,ACHA
 & Score: 26377.1
'
############################################################################################################################

############################################################################################################################
# RUN: 4: repeat run3 but compute full model 4 times
# updated MCMCtree phylogeny and updated error estimates
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################

mkdir ~/data/inqGen18/geneFamilies/CAFE/
run=RUN4
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
#TAG: Setup for subsetted data set (Acro + Acol)
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
#TAG: Prepare CAFE script RUN4
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
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu
lambdamu -s -t $LAMBDATREE
report null
" >> $cafescript


FULLTREE=$LAMBDATREE
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {1..10}
do
  FULLTREE=$(echo $FULLTREE|sed -e "s/1/$i/2")
done

echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for FULLTREE
lambdamu -s -t $FULLTREE
report full1

#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for FULLTREE
lambdamu -s -t $FULLTREE
report full2

#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for FULLTREE
lambdamu -s -t $FULLTREE
report full3

#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for FULLTREE
lambdamu -s -t $FULLTREE
report full4

#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for FULLTREE
lambdamu -s -t $FULLTREE
report full5

" >> $cafescript

#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
CLUSTERTREE="(((4,(4,2)3)1,(4,3)1)1,1);"

echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for clusters
lambdamu -s -t $CLUSTERTREE
report cluster01
" >> $cafescript

echo $cafescript

##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
#(((AHEY1,(ACHA2,PARG3)4)5,(AECH6,AINS7)8)9,ACOL10);
#grep "Search" $run.txt -A 3|grep "Lambda :" -A 2 > rates.$run.txt
grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/full.cafe $base/$run
mv $base/cluster01.cafe $base/$run
