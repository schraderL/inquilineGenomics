
# CAFE runs

<!-- TOC -->

- [CAFE runs](#cafe-runs)
  - [RUN1: null model and 2-parameter models](#run1-null-model-and-2-parameter-models)
    - [Preapare Environment](#preapare-environment)
    - [INPUT files](#input-files)
    - [Setup data](#setup-data)
    - [Check trees](#check-trees)
    - [Retrieve error models](#retrieve-error-models)
    - [Prepare CAFE script RUN1](#prepare-cafe-script-run1)
    - [Check CAFE script](#check-cafe-script)
    - [Run CAFE](#run-cafe)
    - [Gather results](#gather-results)
  - [RUN2: Run all 2-parameter models](#run2-run-all-2-parameter-models)
    - [Preapare Environment](#preapare-environment-1)
    - [INPUT files](#input-files-1)
    - [Setup data](#setup-data-1)
    - [Check trees](#check-trees-1)
    - [Retrieve error models](#retrieve-error-models-1)
    - [Prepare CAFE script](#prepare-cafe-script)
    - [Check CAFE script](#check-cafe-script-1)
    - [Run CAFE](#run-cafe-1)
    - [Gather results](#gather-results-1)
  - [RUN3: 6-parameter models and 2-parameter models](#run3-6-parameter-models-and-2-parameter-models)
    - [Preapare Environment](#preapare-environment-2)
    - [INPUT files](#input-files-2)
    - [Setup data](#setup-data-2)
    - [Check trees](#check-trees-2)
    - [Retrieve error models](#retrieve-error-models-2)
    - [Prepare Cluster Tree](#prepare-cluster-tree)
    - [Prepare CAFE script RUN3](#prepare-cafe-script-run3)
    - [Check CAFE script](#check-cafe-script-2)
    - [Run CAFE](#run-cafe-2)
    - [Gather results](#gather-results-2)
  - [RUN5: Repetition 6-parameter models and 2-parameter models](#run5-repetition-6-parameter-models-and-2-parameter-models)
    - [Preapare Environment](#preapare-environment-3)
    - [INPUT files](#input-files-3)
    - [Setup data](#setup-data-3)
    - [Check trees](#check-trees-3)
    - [Retrieve error models](#retrieve-error-models-3)
    - [Prepare Cluster Tree](#prepare-cluster-tree-1)
    - [Prepare CAFE script](#prepare-cafe-script-1)
    - [Check CAFE script](#check-cafe-script-3)
    - [Run CAFE](#run-cafe-3)
    - [Gather results](#gather-results-3)
  - [RUN6: Second repetition 6-parameter models and 2-parameter models](#run6-second-repetition-6-parameter-models-and-2-parameter-models)
    - [Preapare Environment](#preapare-environment-4)
    - [INPUT files](#input-files-4)
    - [Setup data](#setup-data-4)
    - [Check trees](#check-trees-4)
    - [Retrieve error models](#retrieve-error-models-4)
    - [Prepare Cluster Tree](#prepare-cluster-tree-2)
    - [Prepare CAFE script](#prepare-cafe-script-2)
    - [Check CAFE script](#check-cafe-script-4)
    - [Run CAFE](#run-cafe-4)
    - [Gather results](#gather-results-4)
- [Combine all runs](#combine-all-runs)
- [Supplementary](#supplementary)
  - [RUN: 4-cluster runs](#run-4-cluster-runs)
  - [Run cafe RUN4](#run-cafe-run4)

<!-- /TOC -->


# NOTE:
We removed OG1 with over 100 copies in each species (except Parg).

rep| OG|  ACHA|ACOL|AECH|AHEY|AINS|PARG
|-|-|-|-|-|-|-|-|
ACOLG00000014530.1|  OG00001 |145| 179| 154| 157| 146| 56

This is a orthogroup of odorant receptors!
<font size=1>
> \>**ACOLG00000014530.1** gene=ACOLG00000014530 seq_id=scaffold305 type=cds
MDFFDTRYFRINKIFLSFIGLWPYQTSFIKFLTQSFAILGVAIMCAPQIAYMFKHVDDLD
NMFELMPILAGTVICIAKIISLTCNSEMFKGLLQHMQDDWNNLLTSEETQILTHYAEKSR
TLMLAYSISVIGFVFCYALLPLTEPVFDIILPMNETRPRKLPHLADFVILDQEKYYYTLL
LILYVSYVVCVSIAVAADILYIFLVEHICGMYGVLCHRLRNLATHDNLRWINGNYIHKEI
GRYVQRCIQLHERIRLFIEMMESTISLFLFFDIGLGFLLHTSSCIMIIVRMGSSEIMRYV
ALMLMQSCRLFFNSWAGQEVTDHSVEVSIAAYDGIWYNASVKVQKLLLFLIARSQKVSQI
TIAKLYVINLEGFSKFMRTSVYYITIMVSLNSDM*


rep|OG|ACHA|ACOL|AECH|AHEY|AINS|PARG
|-|-|-|-|-|-|-|-|
ACOLG00000014530.1|OG00001|145|179|154|157|146|56
ACHAG00000000065.1|OG00002|66|73|80|73|76|54
AHEYG00000008852.1|OG00003|70|73|69|71|70|68
AINSG00000010943.1|OG00004|67|68|83|70|79|53
ACOLG00000010388.1|OG00005|73|76|72|70|71|51
ACOLG00000004415.1|OG00006|68|43|62|68|62|69
AECHG00000001331.1|OG00007|52|51|52|53|51|51
AHEYG00000005716.2|OG00008|49|55|54|47|54|33
ACOLG00000012226.1|OG00009|45|48|47|46|44|45
AINSG00000007534.1|OG00010|46|47|50|51|46|18
PARGG00000005389.1|OG00011|41|44|42|41|41|40
ACOLG00000008449.1|OG00012|41|44|43|44|43|25
AHEYG00000002132.1|OG00013|40|43|42|40|44|30
ACOLG00000011799.2|OG00014|36|35|35|35|35|35
AECHG00000009397.1|OG00015|34|35|35|35|33|35
AECHG00000004244.1|OG00016|35|32|33|40|33|32
AHEYG00000015543.1|OG00017|33|31|32|31|31|31
AECHG00000011290.1|OG00018|31|33|31|31|31|31
AECHG00000015667.1|OG00019|30|29|33|32|29|23
ACHAG00000016535.1|OG00020|30|22|32|33|31|27
AECHG00000011108.1|OG00021|29|36|29|25|32|24
ACHAG00000010571.1|OG00022|29|28|31|26|35|25
PARGG00000005421.1|OG00023|24|28|25|30|31|25
AHEYG00000004737.1|OG00024|27|17|33|33|30|21
ACHAG00000012955.1|OG00025|26|25|27|26|27|26
ACHAG00000010690.1|OG00026|25|26|24|26|23|26
ACOLG00000006158.2|OG00027|33|15|29|25|25|15
AHEYG00000002292.1|OG00028|23|23|23|22|23|23
AECHG00000004413.1|OG00029|21|20|22|22|21|20
AINSG00000015957.1|OG00030|28|15|19|27|20|15
AHEYG00000002164.1|OG00031|21|20|25|24|22|10
PARGG00000011028.1|OG00032|22|19|22|20|21|17
AECHG00000009109.1|OG00033|20|19|19|20|19|20
AECHG00000008853.1|OG00034|20|18|19|21|19|20
AINSG00000001936.1|OG00035|20|17|21|20|18|18
AECHG00000000343.1|OG00036|18|18|18|18|18|18
ACOLG00000012015.2|OG00037|17|17|18|19|18|18
ACOLG00000008481.1|OG00038|17|18|18|19|17|15
ACOLG00000008321.1|OG00039|17|16|18|16|17|17
AECHG00000009141.1|OG00040|17|17|17|18|17|14
AHEYG00000008820.1|OG00041|19|15|13|18|14|17
ACHAG00000008002.1|OG00042|14|1|22|13|32|10
ACHAG00000004528.2|OG00043|14|20|15|16|15|11
PARGG00000003299.1|OG00044|14|17|16|16|14|13
ACHAG00000005795.1|OG00045|14|14|14|14|15|15
AECHG00000007394.1|OG00046|16|14|14|16|13|13
PARGG00000004118.1|OG00047|14|15|14|14|14|14
AHEYG00000003365.1|OG00048|14|14|14|14|14|15
AINSG00000002076.2|OG00049|14|14|14|14|14|14
AHEYG00000014557.1|OG00050|14|14|14|14|14|14
AINSG00000010284.1|OG00051|14|14|14|14|14|14
ACOLG00000001248.1|OG00052|15|11|14|14|14|15
AECHG00000011733.1|OG00053|16|10|14|14|15|11
AINSG00000011866.1|OG00054|12|13|15|12|14|13
AINSG00000007374.2|OG00055|13|13|13|13|13|13
ACOLG00000014685.1|OG00056|13|13|13|13|13|13
AINSG00000016270.1|OG00057|13|13|13|13|13|12
ACOLG00000012111.2|OG00058|13|12|14|13|12|12
ACHAG00000014801.1|OG00059|11|8|15|16|14|11
AECHG00000004072.1|OG00060|13|11|13|13|13|12
AECHG00000004317.1|OG00061|12|13|14|12|13|10
AHEYG00000010512.1|OG00062|14|12|12|12|12|12
AHEYG00000010932.1|OG00063|12|12|14|12|13|10
ACOLG00000008513.1|OG00064|12|12|12|12|12|12
ACHAG00000012196.1|OG00065|12|12|12|12|12|12
AECHG00000015195.1|OG00066|12|12|12|12|12|12
ACHAG00000000001.2|OG00067|12|11|12|13|12|11
ACHAG00000007706.1|OG00068|16|12|13|12|12|6

</font>



## RUN1: null model and 2-parameter models

--------------------
**RUN NAME = RUN1**

see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN1/`
see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN1/ctl/RUN1.subsettedMCLclusters.errorcorrected.cafe`

--------------------

Run the null model (1 lambda and mu rate) and all 2-parameter models.

### Preapare Environment
```bash
############################################################################################################################
#RUN: 1: Run null model all two parameter models, updated MCMCtree phylogeny and updated error estimates
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################

#mkdir /corefac/cse/lukas/inqGen18/reannotation/CAFE/
run=RUN1
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/
mkdir $base
mkdir $base/ctl
cd $base/
```
###  INPUT files

>**Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):**
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

>**MCL clustered genes:**
*inflation parameter I = 1.5*
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv

>**Error models from CAFError**
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt

```bash
# cluster file
clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt)
errorbase=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/tmp)

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
```
###  Setup data
```bash
#retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
#restricted tree to ACOL+Acromyrmex
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
#multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
#prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
```

### Check trees
`echo $tree`
`echo $LAMBDATREE`

### Retrieve error models
```bash
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


```

### Prepare CAFE script RUN1
```bash
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


```
### Check CAFE script
`less $cafescript`

### Run CAFE

```bash
##############################################################
# TAG: Run cafe
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out
```

### Gather results
```bash

grep "Search" $run.txt -A 1|grep "Lambda :" > rates.$run.txt
mkdir $base/$run
mv $base/null.cafe $base/$run
mv $base/r*.cafe $base/$run

############################################################################################################################
```

## RUN2: Run all 2-parameter models

--------------------
**RUN NAME = RUN2**

see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN2`
see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN2/ctl/RUN2.subsettedMCLclusters.errorcorrected.cafe`

--------------------

Repeat two-parameter models and run final 6-parameter models. Each model is repeated 5 times.

### Preapare Environment
```bash
run=RUN2
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/

mkdir $base
mkdir $base/ctl
cd $base/
```


###  INPUT files

>**Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):**
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

>**MCL clustered genes:**
*inflation parameter I = 1.5*
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv

>**Error models from CAFError**
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt


```bash

# cluster file (subsetted to Tsep+leafcutters)
clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt)
errorbase=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/tmp)

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
```
###  Setup data
```bash
# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
```

### Check trees
`echo $tree`
`echo $LAMBDATREE`

### Retrieve error models

```bash
for i in AECH AHEY AINS ACHA ACOL PARG
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"
```

### Prepare CAFE script
```bash
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
# two parameter models
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {10..1}
do
  #printf '%s\n' "report r$i"
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.a"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.b"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.c"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.d"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.e"$i
done  >> $cafescript


```
### Check CAFE script
`less $cafescript`

### Run CAFE

```bash
##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

```

### Gather results
```bash
##############################################################
# TAG: Gather Results
##############################################################
grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run
```

--------------------------

## RUN3: 6-parameter models and 2-parameter models

--------------------
**RUN NAME = RUN3**

see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN3`
see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN3/ctl/RUN3.subsettedMCLclusters.errorcorrected.cafe`

--------------------

Repeat two-parameter models and run final 6-parameter models. Each model is repeated 5 times.

### Preapare Environment
```bash
run=RUN3
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/

mkdir $base
mkdir $base/ctl
cd $base/
```


###  INPUT files

>**Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):**
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

>**MCL clustered genes:**
*inflation parameter I = 1.5*
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv

>**Error models from CAFError**
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt


```bash

# cluster file (subsetted to Tsep+leafcutters)
clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt)
errorbase=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/tmp)

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
```
###  Setup data
```bash
# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
```

### Check trees
`echo $tree`
`echo $LAMBDATREE`

### Retrieve error models

```bash
for i in AECH AHEY AINS ACHA ACOL PARG
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"
```


### Prepare Cluster Tree

Here we set a tree with 6 clusters to run in CAFE.
```bash
# Folder with kmeans clustered data
## Retrieve clusters from original runs
#cut -f 8,14 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.R1.tsv |sed 1d > $base/ctl/kmeans.6.txt
cp ~/data/inqGen18/geneFamilies/CAFE/RUN4C/ctl/kmeans.6.txt $base/ctl/kmeans.6.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done <$base/ctl/kmeans.6.txt

```

### Prepare CAFE script RUN3
```bash
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
## Run final 5-parameter model
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

" >> $cafescript

# two parameter models
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {1..5}
do
  #printf '%s\n' "report r$i"
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.a"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.b"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.c"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.d"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.e"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.f"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.g"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.h"$i
done  >> $cafescript


```
### Check CAFE script
`less $cafescript`

### Run CAFE

```bash
##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

```

### Gather results
```bash
##############################################################
# TAG: Gather Results
##############################################################
grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run
```
-------------------

--------------------------

## RUN5: Repetition 6-parameter models and 2-parameter models

--------------------
**RUN NAME = RUN5**

see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN5`
see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN5/ctl/RUN5.subsettedMCLclusters.errorcorrected.cafe`

--------------------

Repeat two-parameter models and run final 6-parameter models. Each model is repeated 5 times.

### Preapare Environment
```bash
run=RUN5
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/

mkdir $base
mkdir $base/ctl
cd $base/
```


###  INPUT files

>**Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):**
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

>**MCL clustered genes:**
*inflation parameter I = 1.5*
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv

>**Error models from CAFError**
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt


```bash

# cluster file (subsetted to Tsep+leafcutters)
clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt)
errorbase=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/tmp)

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
```
###  Setup data
```bash
# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
```

### Check trees
`echo $tree`
`echo $LAMBDATREE`

### Retrieve error models

```bash
for i in AECH AHEY AINS ACHA ACOL PARG
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"
```


### Prepare Cluster Tree

Here we set a tree with 6 clusters to run in CAFE.
```bash
# Folder with kmeans clustered data
## Retrieve clusters from original runs
#cut -f 8,14 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.R1.tsv |sed 1d > $base/ctl/kmeans.6.txt
cp ~/data/inqGen18/geneFamilies/CAFE/RUN4C/ctl/kmeans.6.txt $base/ctl/kmeans.6.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done <$base/ctl/kmeans.6.txt

```

### Prepare CAFE script
```bash
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
## Run final 5-parameter model
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

" >> $cafescript

# two parameter models
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {1..5}
do
  #printf '%s\n' "report r$i"
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.a"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.b"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.c"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.d"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.e"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.f"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.g"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.h"$i
done  >> $cafescript


```
### Check CAFE script
`less $cafescript`

### Run CAFE

```bash
##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

```

### Gather results
```bash
##############################################################
# TAG: Gather Results
##############################################################
grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run
```
-------------------
-------------------


## RUN6: Second repetition 6-parameter models and 2-parameter models

--------------------
**RUN NAME = RUN6**

see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN6`
see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN6/ctl/RUN6.subsettedMCLclusters.errorcorrected.cafe`

--------------------

Repeat two-parameter models and run final 6-parameter models. Each model is repeated 5 times.

### Preapare Environment
```bash
run=RUN6
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/

mkdir $base
mkdir $base/ctl
cd $base/
```


###  INPUT files

>**Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):**
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

>**MCL clustered genes:**
*inflation parameter I = 1.5*
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv

>**Error models from CAFError**
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt


```bash

# cluster file (subsetted to Tsep+leafcutters)
clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt)
errorbase=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/tmp)

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
```
###  Setup data
```bash
# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
```

### Check trees
`echo $tree`
`echo $LAMBDATREE`

### Retrieve error models

```bash
for i in AECH AHEY AINS ACHA ACOL PARG
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"
```


### Prepare Cluster Tree

Here we set a tree with 6 clusters to run in CAFE.
```bash
# Folder with kmeans clustered data
## Retrieve clusters from original runs
#cut -f 8,14 ~/data/inqGen18/geneFamilies/CAFE/clustered.rates.CAFE.R1.tsv |sed 1d > $base/ctl/kmeans.6.txt
cp ~/data/inqGen18/geneFamilies/CAFE/RUN4C/ctl/kmeans.6.txt $base/ctl/kmeans.6.txt

cd $base/ctl
VARIABLETREE="(((AHEY,(ACHA,PARG)ACHA/PARG)AHEY/ACHA/PARG,(AECH,AINS)AECH/AINS)ACRO,ACOL);"
while read KMEANS; do
  MODEL=$(echo $KMEANS | cut -f1 -d " ")
  CLUSTER=$(echo $KMEANS | cut -f2 -d " ")
  VARIABLETREE=$(echo $VARIABLETREE|sed  "s|$MODEL|$CLUSTER|g")
done <$base/ctl/kmeans.6.txt

```

### Prepare CAFE script
```bash
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
## Run final 5-parameter model
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

lambdamu -s -t $VARIABLETREE
report m6j


" >> $cafescript

# two parameter models
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {1..10}
do
  #printf '%s\n' "report r$i"
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.a"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.b"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.c"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.d"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.e"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.f"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.g"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.h"$i
done  >> $cafescript


```
### Check CAFE script
`less $cafescript`

### Run CAFE

```bash
##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

```

### Gather results
```bash
##############################################################
# TAG: Gather Results
##############################################################

#run=RUN7
#base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/
#cd $base
grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.txt |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run
```

## RUN7: Run 2-parameter models for Aech and Ains

--------------------
**RUN NAME = RUN7**

see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN7`
see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN7/ctl/RUN7.subsettedMCLclusters.errorcorrected.cafe`

--------------------

Repeat two-parameter models and run final 6-parameter models. Each model is repeated 5 times.

### Preapare Environment
```bash
run=RUN7
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/

mkdir $base
mkdir $base/ctl
cd $base/
```


###  INPUT files

>**Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):**
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

>**MCL clustered genes:**
*inflation parameter I = 1.5*
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv

>**Error models from CAFError**
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt


```bash

# cluster file (subsetted to Tsep+leafcutters)
clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt)
errorbase=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/tmp)

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
```
###  Setup data
```bash
# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
```

### Check trees
`echo $tree`
`echo $LAMBDATREE`

### Retrieve error models

```bash
for i in AECH AHEY AINS ACHA ACOL PARG
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"
```

### Prepare CAFE script
```bash
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
# two parameter models
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {6..7}
do
  #printf '%s\n' "report r$i"
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.a"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.b"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.c"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.d"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.e"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.f"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.g"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.h"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.i"$i
done  >> $cafescript


```
### Check CAFE script
`less $cafescript`

### Run CAFE

```bash
##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

```

### Gather results
```bash
##############################################################
# TAG: Gather Results
##############################################################
grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run
```

--------------------------

## RUN7: Run 2-parameter models for Aech and Ains

--------------------
**RUN NAME = RUN8**

see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN8`
see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN8/ctl/RUN8.subsettedMCLclusters.errorcorrected.cafe`

--------------------

Repeat two-parameter models and run final 6-parameter models. Each model is repeated 5 times.

### Preapare Environment
```bash
run=RUN8
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/

mkdir $base
mkdir $base/ctl
cd $base/
```


###  INPUT files

>**Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):**
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

>**MCL clustered genes:**
*inflation parameter I = 1.5*
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv

>**Error models from CAFError**
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt


```bash

# cluster file (subsetted to Tsep+leafcutters)
clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt)
errorbase=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/tmp)

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
```
###  Setup data
```bash
# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
```

### Check trees
`echo $tree`
`echo $LAMBDATREE`

### Retrieve error models

```bash
for i in AECH AHEY AINS ACHA ACOL PARG
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"
```

### Prepare CAFE script
```bash
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
# two parameter models
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {7..7}
do
  #printf '%s\n' "report r$i"
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.a"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.b"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.c"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.d"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.e"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.f"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.g"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.h"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.i"$i
done  >> $cafescript

for i in {6..6}
do
  #printf '%s\n' "report r$i"
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.a"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.b"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.c"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.d"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.e"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.f"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.g"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.h"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.i"$i
done  >> $cafescript

```
### Check CAFE script
`less $cafescript`

### Run CAFE

```bash
##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

```

### Gather results
```bash
##############################################################
# TAG: Gather Results
##############################################################
grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run
```

--------------------------



## RUN7: Run 2-parameter models for Aech and Ains on unfiltered gene families

--------------------
**RUN NAME = RUN8**

see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN9`
see `/corefac/cse/lukas/inqGen18/reannotation/CAFE/RUN9/ctl/RUN9.subsettedMCLclusters.errorcorrected.cafe`

--------------------

Repeat two-parameter models and run final 6-parameter models. Each model is repeated 5 times.

### Preapare Environment
```bash
run=RUN9
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/

mkdir $base
mkdir $base/ctl
cd $base/

#/corefac/cse/lukas/inqGen18/reannotation/CAFE/geneFamilies.TEfree.edit.tsv
```


###  INPUT files

>**Ultrametric, MCMCtree dated, phylogenetic tree from run 4 (sf 50000):**
~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

>**MCL clustered genes:**
*inflation parameter I = 1.5*
/corefac/cse/lukas/inqGen18/reannotation/CAFE/geneFamilies.TEfree.edit.tsv

>**Error models from CAFError**
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt


```bash

# cluster file (subsetted to Tsep+leafcutters)
clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/geneFamilies.TEfree.edit.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt)
errorbase=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/tmp)

# cafe script
cafescript=$base/ctl/$run.subsettedMCLclusters.errorcorrected.cafe
```
###  Setup data
```bash
# retrieve tree from MCMCtree output
  treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")

  #(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);
  subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
  tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
  LAMBDATREE=$(echo $subtree|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
```

### Check trees
`echo $tree`
`echo $LAMBDATREE`

### Retrieve error models

```bash
for i in AECH AHEY AINS ACHA ACOL PARG
do
  file=$(grep $i $errorfile|cut -f 2)
  echo "errormodel -model $errorbase/$file -sp $i"
done > $base/ctl/errorModels.txt

cd $base
out=$base"/"$run".txt"
```

### Prepare CAFE script
```bash
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
# two parameter models
#echo $LAMBDATREE|grep -o 1|wc -l
for i in {6..7}
do
  #printf '%s\n' "report r$i"
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.a"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.b"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.c"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.d"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.e"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.f"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.g"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.h"$i
  echo -ne "lambdamu -s -t "
  echo $LAMBDATREE|sed -e "s/1/2/$i"
  echo "report p2.i"$i
done  >> $cafescript

```
### Check CAFE script
`less $cafescript`

### Run CAFE

```bash
##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out

```

### Gather results
```bash
##############################################################
# TAG: Gather Results
##############################################################
grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run
```

--------------------------

-------------------
-------------------


# Combine all runs
```bash
bd=/corefac/cse/lukas/inqGen18/reannotation/CAFE
mkdir $bd/results
cd $bd/results
grep -r -n ".*" $bd/RUN*/tree*|perl -pe 's/(.*RUN([0-9]+).tsv)\:([0-9]+)\:(.*?)\t/R$2m$3\t$1\t$4\t/g' > allProperRuns.tsv
awk '$3 !~ /[3-6]/' allProperRuns.tsv |awk '$3 ~ /2/' > 2pModels.tsv
```
`/corefac/cse/lukas/inqGen18/reannotation/CAFE/results/2pModels.tsv`
`/corefac/cse/lukas/inqGen18/reannotation/CAFE/results/allProperRuns.tsv`

# Supplementary

## RUN: 4-cluster runs
<!--
```bash
############################################################################################################################
# RUN: 4: repeat run3 but compute full model 4 times
# updated MCMCtree phylogeny and updated error estimates
############################################################################################################################

##############################################################
#TAG: Preapare Environment
##############################################################
run=RUN4
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/$run/

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
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv
# Error models from CAFError
## Run 1, run on MCMCtree run4
/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt
'

# cluster file (subsetted to Tsep+leafcutters)
clusters=clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre

# error models
errorfile=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.perSpecies.txt)
errorbase=$(readlink -f $base/../CAFE.errorprediction/errorTest.smallSet.filter100/tmp)

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

echo $LAMBDATREE
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

#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
CLUSTERTREE="(((4,(4,2)3)1,(4,3)1)1,1);"

echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for clusters
lambdamu -s -t $CLUSTERTREE
report cluster01
" >> $cafescript

echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for clusters
lambdamu -s -t $CLUSTERTREE
report cluster02
" >> $cafescript

echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for clusters
lambdamu -s -t $CLUSTERTREE
report cluster03
" >> $cafescript

echo "
#(((AHEY,(ACHA,PARG)),(AECH,AINS)),ACOL);
#calculate lambdamu for clusters
lambdamu -s -t $CLUSTERTREE
report cluster04
" >> $cafescript

```


`less $cafescript`


## Run cafe RUN4
```bash
##############################################################
# TAG: Run cafe RUN4
##############################################################
cd $base/
nice /corefac/cse/lukas/software/CAFE/release/cafe $cafescript > $run.out


##############################################################
# TAG: Gather Results
##############################################################
grep "Search" $run.txt -A 3|grep "Lambda :" -A 2|perl -00pe 's/DONE: Lamda,Mu Search or setting, for command:\n/ /g' |perl -00pe 's/ lambdamu .*?\n//g'|grep "Lambda :"> rates.$run.txt
cat rates.$run.txt |tr "," "\t"|sed 's/Lambda : //g'|sed 's/ \& Score://g'|sed 's/Mu ://g'  |tr " " "\t" > rates.$run.tsv
grep "Lambda Tree:" $run.out |cut -f 3 -d " "|paste - rates.$run.tsv > treeRates.$run.tsv
mkdir $base/$run
mv $base/*.cafe $base/$run


```
