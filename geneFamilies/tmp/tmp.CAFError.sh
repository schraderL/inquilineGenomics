

###############################
#RUN: 3: Setup Run errorTest2
###############################

runName=errorTest2
###############################
  runbase=$base/CAFE/CAFE.errorprediction/$runName
  cscript=$runbase/$runName.cafe
  outsummary=$runbase/$runName.out
  outreport=$runbase/report/$runName.report.out

# prepare environment
#mkdir $base/CAFE/CAFE.errorprediction/run1 # First test run. Ignore.
mkdir $runbase
mkdir $runbase/report/
mkdir $runbase/tmp
cd $runbase

## Run the error prediction script, predicting error rates for each species
# retrieve tree from MCMCtree output
treeOriginal=$(cat $treeFile|grep UTREE|perl -pe 's/\[.*?\]//g'|cut -f 2 -d "="|sed "s/ //g")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526
tree=$(echo $treeOriginal|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
LAMBDATREE=$(echo $treeOriginal|perl -pe 's/([A-Z]{4})\://g'|perl -pe 's/[0-9]\.[0-9]+/1/g'|perl -pe 's/://g')
# Prepare cafe script from template.
cat $base/CAFE/CAFE.errorprediction/ctl/errorTest.template.cafe | \
  sed "s|<CLUSTERS>|$clusters|g" | \
  sed "s|<TREE>|$tree|g" | \
  sed "s|<LAMBDATREE>|$LAMBDATREE|g" | \
  sed "s|<OUTPUTSUMMARY>|$outsummary|g" | \
  sed "s|<OUTPUTREPORT>|$outreport|g" \
> $cscript

nice python $software/CAFE/cafe/caferror.py -i $cscript -d $runbase/tmp -s 1

#python caferror.py -i cafe_run2.sh -d errorout -f 1

for i in AECH AHEY AINS ACHA ACOL ACEP PARG TSEP
do
 grep "Search Result"  *$i* -A 1|grep -E "Score:.*"|sed s/cafe_//g|sed s/'_.*Score: '/' '/g > ../err/$i.err
done
