# InterproScan
```bash
# get Data
base=/home/s/schradel/data/inqGen/InterProScan
mkdir $base
mkdir $base/data
cd $base/data

scp -o ProxyCommand="ssh -W %h:%p zcm795@ssh-bio-labnet.science.ku.dk" lschrader@pallas.bio.ku.dk:/corefac/cse/lukas/inqGen18/reannotation/results/*/GAF/*.GeMoMa.longestIsoform.pep.fa .

```


```bash
#@jgpogo
mkdir $base/interpro

for file in $base/data/*
do
  pep=$(readlink -f $file)
  cat $pep |perl -pe 's/\*//g' > $pep.edit
  cd /global/homes/jg/schradel/dbs/interpro/interproscan-5.44-79.0

  ./interproscan.sh \
  --input $pep.edit \
  --disable-precalc \
  --output-dir $base/interpro \
  --formats TSV,XML,GFF3 \
  --goterms
done
cd /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro

```
