# CAFE error estimates

We repeated the error estimation three times in independent runs. All runs returned the same estimates for each species.

-----------------------------
<!-- TOC -->

- [CAFE error estimates](#cafe-error-estimates)
- [Summary](#summary)
  - [Prepare folders and input files](#prepare-folders-and-input-files)
    - [Prepare template CAFError script](#prepare-template-caferror-script)
    - [Prepare subsetted MCL cluster file](#prepare-subsetted-mcl-cluster-file)
    - [Filter genes with >100 copies in one species from Cluster file](#filter-genes-with-100-copies-in-one-species-from-cluster-file)
  - [RUN1](#run1)
    - [Prepare Run](#prepare-run)
    - [Setup Run environment](#setup-run-environment)
    - [Prepare the error prediction script](#prepare-the-error-prediction-script)
    - [Run the error prediction script, predicting error rates for each species](#run-the-error-prediction-script-predicting-error-rates-for-each-species)
    - [Gather best fitting error models](#gather-best-fitting-error-models)
    - [Store results](#store-results)
  - [RUN2](#run2)
    - [Prepare Run](#prepare-run-1)
    - [Setup Run environment](#setup-run-environment-1)
    - [Prepare the error prediction script](#prepare-the-error-prediction-script-1)
    - [Run the error prediction script, predicting error rates for each species](#run-the-error-prediction-script-predicting-error-rates-for-each-species-1)
    - [Gather best fitting error models](#gather-best-fitting-error-models-1)
    - [Store results](#store-results-1)
  - [RUN3](#run3)
    - [Prepare Run](#prepare-run-2)
    - [Setup Run environment](#setup-run-environment-2)
    - [Prepare the error prediction script](#prepare-the-error-prediction-script-2)
    - [Run the error prediction script, predicting error rates for each species](#run-the-error-prediction-script-predicting-error-rates-for-each-species-2)
    - [Gather best fitting error models](#gather-best-fitting-error-models-2)
    - [Store results](#store-results-2)

<!-- /TOC -->

-----------------------------

# Summary
`cat /corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/FinalErrorEstimates.txt`
**see also**
>************************************************************************
>SUMMARY: Runs 1 and 2 yielded ~same results
>************************************************************************
>Final error estimates by species:
**PARG                             0.0228515625**
**ACOL                             0.0**
**AINS                             0.0087890625**
**AHEY                             0.010546875**
**ACHA                             0.00703125**
**AECH                             0.0140625**
**Score with individual errors:    22274.478868**
Lambda with individual errors:   0.00042094533780
>
>************************************
>Score with no errormodel:        22532.956177
Lambda with no errormodel:       0.00052466465961
>************************************
>Global Error Estimation:         0.0087890625
Score with global errormodel:    22356.238701
Lambda with global errormodel:   0.00044082832371
>************************************
>
>Caferror finished at:            09.27.2020 | 11:59:01
Runtime:                         42.6983719667 minutes
>************************************

```bash
cat /corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100.run2/FinalErrorEstimates.txt
cat /corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100.run3/FinalErrorEstimates.txt
```

## Prepare folders and input files

### Prepare template CAFError script
```bash
###############################
#TAG: Prepare template CAFError script
###############################
base=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/
old=~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction
bd=~/data/inqGen18/geneFamilies
software=~/data/software
mkdir /corefac/cse/lukas/inqGen18/reannotation/CAFE/
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
> $base/ctl/errorTest.template.cafe
```
### Prepare subsetted MCL cluster file
```bash
###############################
#TAG: Prepare subsetted MCL cluster file
###############################

# input files
originalClusters=/corefac/cse/lukas/inqGen18/reannotation/MCL/geneFamilies.TEfree.tsv
# treeFile
treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
```
### Filter genes with >100 copies in one species from Cluster file
```bash
###############################
#TAG: Filter genes with >100 copies in one species from Cluster file
###############################
cd $base
rawclusters=/corefac/cse/lukas/inqGen18/reannotation/MCL/geneFamilies.TEfree.tsv

head -n 1 $rawclusters > $base/filtered.geneFamilies.TEfree.tsv
cat $rawclusters |awk  -F $'\t' 'BEGIN {OFS=FS} {if (($3 <100 && $4 <100 && $5 <100 && $6 <100 && $7 <100 && $8 <100 && $9 <100 && $10 <100) || $1=="Description") print $0}' \
>> $base/filtered.geneFamilies.TEfree.tsv

cd $base/..
ln -s $base/filtered.geneFamilies.TEfree.tsv .

```

---------------------------------------------------------

## RUN1
 **Ultrametric, MCMCtree dated, phylogenetic tree:**
`~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf500.sr17.5cal.sigma2.5.5.run1/FigTree.run4.tre`
 **MCL clustered genes:**
 *inflation parameter I = 1.5*
`/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv`
 **Error prediction cafe script**
`~/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe`

### Prepare Run
```bash
############################################
#TAG: INPUT files
############################################
# input files
clusters=/corefac/cse/lukas/inqGen18/reannotation/CAFE/CAFE.errorprediction/filtered.geneFamilies.TEfree.tsv

###############################
# TAG: Define runName
###############################
runName=errorTest.smallSet.filter100

  runbase=$base/$runName
  cscript=$runbase/$runName.cafe
  outsummary=$runbase/$runName.out
  outreport=$runbase/report/$runName.report.out
```

### Setup Run environment
```bash

# prepare environment
#mkdir $base/CAFE/CAFE.errorprediction/run1 # First test run. Ignore.
mkdir $runbase
mkdir $runbase/tmp
mkdir $runbase/err
cd $runbase
```

### Prepare the error prediction script
```bash
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
echo $LAMBDATREE

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

```

Check the script:
`less $cscript`

### Run the error prediction script, predicting error rates for each species
```bash
###############################
#TAG: Run CAFError
###############################
nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1 > caferror.run.out

```

### Gather best fitting error models
```bash
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

#cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt
```
### Store results
```bash
# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s $runbase/FinalErrorEstimates.txt .
ln -s $runbase/FinalErrorEstimates.perSpecies.txt .

```

-------------------------------

## RUN2
### Prepare Run
 Repeat same as above but change runName to `runName=RUN2`
```bash
###############################
#RUN: 2: Run errorTest.smallSet.filter100
###############################
cd $base
###############################
# TAG: Define runName
###############################

runName=errorTest.smallSet.filter100.run2
  runbase=$base/$runName
  cscript=$runbase/$runName.cafe
  outsummary=$runbase/$runName.out
  outreport=$runbase/report/$runName.report.out
```

### Setup Run environment
```bash

# prepare environment
#mkdir $base/CAFE/CAFE.errorprediction/run1 # First test run. Ignore.
mkdir $runbase
mkdir $runbase/tmp
mkdir $runbase/err
cd $runbase
```

### Prepare the error prediction script
```bash
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
echo $LAMBDATREE

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

```

Check the script:
`less $cscript`

### Run the error prediction script, predicting error rates for each species
```bash
###############################
#TAG: Run CAFError
###############################
nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1 > caferror.run.out

```

### Gather best fitting error models
```bash
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

#cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt
```
### Store results
```bash
# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s $runbase/FinalErrorEstimates.txt .
ln -s $runbase/FinalErrorEstimates.perSpecies.txt .

```


## RUN3
### Prepare Run
Repeat same as above but change runName to `runName=RUN3`
```bash
###############################
#RUN: 2: Run errorTest.smallSet.filter100
###############################
runName=errorTest.smallSet.filter100.run3
  runbase=$base/$runName
  cscript=$runbase/$runName.cafe
  outsummary=$runbase/$runName.out
  outreport=$runbase/report/$runName.report.out
```

### Setup Run environment
```bash

# prepare environment
#mkdir $base/CAFE/CAFE.errorprediction/run1 # First test run. Ignore.
mkdir $runbase
mkdir $runbase/tmp
mkdir $runbase/err
cd $runbase
```

### Prepare the error prediction script
```bash
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
echo $LAMBDATREE

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

```

Check the script:
`less $cscript`

### Run the error prediction script, predicting error rates for each species
```bash
###############################
#TAG: Run CAFError
###############################
nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1 > caferror.run.out

```

### Gather best fitting error models
```bash
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

#cd $runbase/bestFitModels
cat $runbase/tmp/caferrorLog.txt
```
### Store results
```bash
# Store results
cd $runbase
mkdir $runbase/results
cd $runbase/results
ln -s ../tmp/caferrorLog.txt .
ln -s $runbase/FinalErrorEstimates.txt .
ln -s $runbase/FinalErrorEstimates.perSpecies.txt .

```
