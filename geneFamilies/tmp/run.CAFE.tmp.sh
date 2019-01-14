##########################################################################################################
# SUMMARY : Results of run2
##########################################################################################################

:'
(TSEP1,(((AHEY2,(ACHA3,PARG4)5)6,(AECH7,AINS8)9)10,(ACEP11,ACOL12)13)14);
00 Lambda : 0.00033392813177                  & Score: 32385.5Mu : 0.00041921778892                  & Score: 32385.5
01 Lambda : 0.00034573011849,0.00031085704148 & Score: 31808.4Mu : 0.00064750896817,0.00000000000003 & Score: 31808.4
02 Lambda : 0.00028944125786,0.00141000486012 & Score: 32222.5Mu : 0.00039928828085,0.00080485366097 & Score: 32222.5
03 Lambda : 0.00031921477799,0.00104400226672 & Score: 32344.4Mu : 0.00040786128739,0.00062330194260 & Score: 32344.4
04 Lambda : 0.00032347373162,0.00030481467410 & Score: 31701.1Mu : 0.00034821272570,0.00644397960341 & Score: 31701.1
05 Lambda : 0.00033875133930,0.00001870567946 & Score: 32303.7Mu : 0.00039571120838,0.00220549195794 & Score: 32303.7
06 Lambda : 0.00090562411531,0.00014131516418 & Score: 35520.4Mu : 0.00033130632816,0.00034024106462 & Score: 35520.4
07 Lambda : 0.00031957993246,0.00170316109102 & Score: 32313.4Mu : 0.00040849239965,0.00113899295692 & Score: 32313.4
08 Lambda : 0.00032448694014,0.00067807148701 & Score: 32302.1Mu : 0.00040583327714,0.00252570310286 & Score: 32302.1
09 Lambda : 0.00034842051880,0.00017132614847 & Score: 32370.9Mu : 0.00043144098315,0.00034294923085 & Score: 32370.9
10 Lambda : 0.00034796915619,0.00023718954603 & Score: 32361.2Mu : 0.00045161191926,0.00024576648782 & Score: 32361.2
11 Lambda : 0.00033619653538,0.00021382706210 & Score: 32259.2Mu : 0.00039432628370,0.00214580665661 & Score: 32259.2
12 Lambda : 0.00032805034643,0.00049534147227 & Score: 32353.4Mu : 0.00040099724110,0.00104288681873 & Score: 32353.4
13 Lambda : 0.00039937697586,0.00008854958156 & Score: 32286.9Mu : 0.00042722247769,0.00036417417934 & Score: 32286.9
14 Lambda : 0.00034272683597,0.00024113519350 & Score: 32341.5Mu : 0.00036641681285,0.00109125867368 & Score: 32341.5
'

##########################################################################################################
#TAG: CAFE gene family size evolution: Error rate estimation
##########################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
base=~/data/inqGen18/geneFamilies/CAFE/CAFE
bd=~/data/inqGen18/geneFamilies
mkdir $base
mkdir $base/ctl
cd $base/

##############################################################
#TAG: INPUT files
##############################################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree:
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run1.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest3.smallSet/results/bestFitModels
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=$bd/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run1.tre


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
errorbase=$bd/CAFE/CAFE.errorprediction/errorTest3.smallSet/results/bestFitModels

# retrieve error models
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
file=$(ls $errorbase/*$i*)
  echo "errormodel -model $file -sp $i"
done > $base/ctl/errorModels.txt

############################################################################################################################
# RUN: 1: Run null model all two parameter models "Run2"
############################################################################################################################

cafescript=$base/ctl/run.subsettedMCLclusters.errorcorrected.cafe
out=$base/run2.txt
##############################################################
#TAG: Prepare CAFE script RUN 1
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

tree $tree
load -i $clusters -t 20 -l run2.txt -p 0.05 -filter

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

##############################################################
# TAG: Run cafe RUN 1
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > run2.out

# Gather results
grep "Search" run2.txt -A 1|grep "Lambda :" > rates.run1.txt

############################################################################################################################


############################################################################################################################
# RUN: 2: Run null model all two parameter models (repeat RUN1) "run3"
############################################################################################################################

cafescript=$base/ctl/run2.subsettedMCLclusters.errorcorrected.cafe
run=run3
out=$base"/"$run".txt"
##############################################################
#TAG: Prepare CAFE script RUN 2
##############################################################

# section to load tree and cluster file
echo "
#! /corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

tree $tree
load -i $clusters -t 20 -l run3.txt -p 0.05 -filter

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
# TAG: Run cafe RUN 2
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################


############################################################################################################################
# RUN: 3: Run null model all two parameter models (repeat RUN1) "run4"
############################################################################################################################

cafescript=$base/ctl/run4.subsettedMCLclusters.errorcorrected.cafe
run=run4
out=$base"/"$run".txt"
##############################################################
#TAG: Prepare CAFE script RUN 2
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
# TAG: Run cafe RUN 2
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
# RUN: 5: Run null model all two parameter models, updated MCMCtree phylogeny and updated error estimates
############################################################################################################################
base=~/data/inqGen18/geneFamilies/CAFE/CAFE
bd=~/data/inqGen18/geneFamilies
cd $base
run=run5
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe

out=$base"/"$run".txt"
##############################################################
#TAG: INPUT files
##############################################################
:'
# Ultrametric, MCMCtree dated, phylogenetic tree:
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run1.tre
# MCL clustered genes:
## inflation parameter I = 1.5
~/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Error models from CAFError
~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest3.smallSet/results/bestFitModels
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=$bd/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre


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
errorbase=$bd/CAFE/CAFE.errorprediction/errorTest5.smallSet/bestFitModels

# retrieve error models
for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
file=$(ls $errorbase/*$i*)
  echo "errormodel -model $file -sp $i"
done > $base/ctl/errorModelsR5.txt


##############################################################
#TAG: Prepare CAFE script RUN 2
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
cat $base/ctl/errorModelsR5.txt >> $cafescript

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
# TAG: Run cafe RUN 2
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

# Gather results
grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################













########################################################################################################################
#TEMP First test runs
########################################################################################################################
#  Run CAFE
/corefac/cse/lukas/software/CAFE/release/cafe
echo "
# $base
#				3		 4     5         6         7	    8    9    10    11   12       13       14
tree (TSEP:0136.666,(((AHEY:0024.759,(ACHA:0016.195,PARG:0016.195):0008.563):0027.448,(AECH:0009.565,AINS:0009.565):0042.642):0064.592,(ACEP:0020.253,ACOL:0020.253):0096.546):0019.867);
#load -i $clusters -t 7 -l $out -p 0.05 -filter
load -i /usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/Tsep.Atta.Acro.largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv -t 10 -l run1.txt -p 0.05 -filter

#assign error models
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies//CAFE/CAFE.errorprediction/errorTest3.smallSet/bestFitModels/cafe_errormodel_0.050703125_AECH.txt -sp AECH
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies//CAFE/CAFE.errorprediction/errorTest3.smallSet/bestFitModels/cafe_errormodel_0.050703125_AHEY.txt -sp AHEY
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies//CAFE/CAFE.errorprediction/errorTest3.smallSet/bestFitModels/cafe_errormodel_0.041484375_AINS.txt -sp AINS
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies//CAFE/CAFE.errorprediction/errorTest3.smallSet/bestFitModels/cafe_errormodel_0.032265625_ACHA.txt -sp ACHA
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies//CAFE/CAFE.errorprediction/errorTest3.smallSet/bestFitModels/cafe_errormodel_0.023046875_ACOL.txt -sp ACOL
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies//CAFE/CAFE.errorprediction/errorTest3.smallSet/bestFitModels/cafe_errormodel_0.059921875_ACEP.txt -sp ACEP
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies//CAFE/CAFE.errorprediction/errorTest3.smallSet/bestFitModels/cafe_errormodel_0.1290625_PARG.txt -sp PARG
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies//CAFE/CAFE.errorprediction/errorTest3.smallSet/bestFitModels/cafe_errormodel_0.0_TSEP.txt -sp TSEP

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

lambdamu -s -t (1,(((1,(1,1)1)1,(1,1)1)1,(1,1)1)1);
report r14

#calculate lambdamu for manually selected clusters
lambdamu -s -t (1,(((2,(3,4)1)1,(1,3)1)1,(1,1)1)1);
report manualClusters

#Two clusters based on kmeans
# 2  2  1  2  2  2  1  1  1  1  1  1
#lambdamu -s -t ((((2,2)1,(2,(2,2)1)1)1,1)1,1)
#report M2K

#Three clusters based on kmeans
# 2  2  1  3  2  3  1  1  1  1  1  1
#lambdamu -s -t ((((2,2)1,(3,(2,3)1)1)1,1)1,1)
#report M3K

#Four clusters based on kmeans
# 4  4  3  1  4  1  2  2  3  3  3  3
#lambdamu -s -t ((((4,4)3,(1,(4,1)2)2)3,3)3,3)
#report M4Krepeat
#report M4Krepeat

#Four "good" clusters based on kmeans
# 2  1  3  2  4  1  3  3  3  3  3  3
#lambdamu -s -t ((((2,1)3,(2,(1,4)3)3)3,3)3,3)
#report M4aK


#Five clusters based on kmeans
#2  2  4  3  2  1  5  5  4  4  4  4
#lambdamu -s -t ((((2,2)4,(3,(2,1)5)5)4,4)4,4)
#report M5K

#Five "good" clusters based on kmeans #actually not good because Aech and Ains are confused
#1  2  2  3  4  5  2  2  2  2  2  2
#lambdamu -s -t ((((1,2)2,(3,(4,5)2)2)2,2)2,2)
#report M5aK


#############################################
#BEST
#############################################
#Five actually "good" clusters based on kmeans
#2  1  5  3  1  4  5  5  5  5  5  5
#lambdamu -s -t ((((2,1)5,(3,(1,4)5)5)5,5)5,5)
#report M5bK
#############################################


#Five alternative clusters joining Aech & Parasites but splitting large Cluster
#2  1  5  3  1  4  5  5  5  5  5  5
lambdamu -s -t ((((1,1)5,(3,(1,4)2)2)5,5)5,5)
report M5cK
#

#Six clusters based on kmeans
# 2  5  4  1  5  1  6  3  4  3  4  4
#lambdamu -s -t ((((2,5)4,(1,(5,1)6)3)4,3)4,4)
#report M6K

#Seven clusters based on kmeans
# 1  1  7  4  1  4  2  6  7  6  3  5
#lambdamu -s -t ((((2,5)4,(1,(5,1)6)3)4,3)4,4)
#report M6bK

#Six "good" clusters based on kmeans
#3  6  2  5  6  1  4  4  2  2  2  2
#lambdamu -s -t ((((3,6)2,(5,(6,1)4)4)2,2)2,2)
#report M6aK

#Seven "good" clusters based on kmeans
# 4  7  3  1  7  2  6  5  3  5  3  3
lambdamu -s -t ((((4,7)3,(1,(7,2)6)5)3,5)3,3)
M7K
\
> /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe
