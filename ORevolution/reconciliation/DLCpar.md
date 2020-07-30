# Gene tree reconciliation for Clade285 of the 9Exon family

We used ```dlcpar 2.0.1``` for inferring a reconciled gene tree. See https://github.com/wutron/dlcpar/blob/master/MANUAL.md and https://github.com/wutron/dlcpar/blob/master/EXAMPLES.md.

Awesome tool right here: https://github.com/evolbioinfo/gotree

## Prepare environment
```bash
base=/usr/local/home/lschrader/data/inqGen18/ORphylogeny/reconciledClades
cd $base
mkdir data
mkdir cladeData

```

## retrieve CDS entries of ORs
```bash
cd $base/data
ln -s ~/data/inqGen18/ORs/*/final/finalSet/*.OR.gff3 .

for i in ~/data/genomes/reannotations_of_7_ants/*/genome/*.v2.0.fa
do
  echo $i
  short=$(echo $i|rev|cut -f 1 -d "/"|rev|perl -pe 's/(.).*_(...).*/$1$2/g')
  ln -s $i $short.genome.fa
done

for i in $(ls ~/data/inqGen18/inquilines_v2.1/*/genome/*.v2.1.fa|egrep "mitome|fungal|bacterial" -v)
do
  echo $i
  short=$(echo $i|rev|cut -f 1 -d "/"|rev|perl -pe 's/(.).*_(...).*/$1$2/g')
  ln -s $i $short.genome.fa
done

ls *.gff3|cut -f 1 -d "."|parallel "gt extractfeat -type CDS -join yes -seqid yes -target yes -retainids yes -matchdescstart -seqfile {}.genome.fa {}.OR.gff3 > {}.CDS.fa"

cat *.CDS.fa > ORs.cds.fa
```

# Run for each clade
```bash
clade=9E-Clade-285
mkdir $base/$clade
cd $base/$clade
```

## Retrieve cds sequence for particular clade

```bash
cut -f 1 $base/cladeData/$clade.tsv|seqkit grep -rf - $base/data/ORs.cds.fa |perl -pe 's/-mRNA.*//g'> $clade.cds.fa
grep ">" $clade.cds.fa -c

#seqkit fx2tab $clade.cds.fa |awk '{print $0,length($2)}'|awk -F $'\t' '{OFS = FS} {if ($3%3!=0) print $1,substr($2,1,$3-$3%3);else print $1,$2}' |seqkit tab2fx > $clade.cds.trimmed.fa
```

## Check if codons are in order
```bash
#edit
cp $clade.cds.fa $clade.cds.edit.fa

#manually fix sequences not multiple of 3
seqkit fx2tab $clade.cds.fa |awk -F $'\t' '{print $0,length($2)}'|awk -F $'\t' '{OFS = FS} {if ($3%3!=0) print $1,$3,$3%3}'

AcepOr-393-NTEfd         532    1 start -1
AechOr-434-NTEfd         532    1
AechOr-377-CTE   1129   1 end -1
AechOr-421-NTEfd         743    2˘

nano $clade.cds.edit.fa # manually edit CDS
seqkit fx2tab $clade.cds.edit.fa |awk -F $'\t' '{print $0,length($2)}'|awk -F $'\t' '{OFS = FS} {if ($3%3!=0) print $1,$3,$3%3}'

```
## Align CDS with ```prank``` in ```codon```mode

```bash
cut -f 1 -d " " $clade.cds.edit.fa > $clade.cds.edit2.fa
prank -t=$base/cladeData/$clade.tre -d=$clade.cds.edit2.fa -o=$clade.repeat.cds.out -codon -iterate=3
# died but still produced an alignment
```
## Compute rooted gene tree (rooted at Ccos OR) phylogeny with RaxML (bootstrap and ML)
```bash
raxmlHPC-PTHREADS -f a -m GTRGAMMA -p 12345 -x 12345 -# 500 -s  $clade.cds.out.best.fas -n rooted -T 4 -o CcosOr-363
```

## Compute gene tree phylogeny with FastTree
```bash
FastTree -nt $clade.cds.out.best.fas > $clade.cds.out.FastTree.fas
#FastTree -nosupport -nt $clade.cds.out.best.fas > $clade.cds.out.FastTree.fas

#manually reroot tree
# outgroup.lst
echo "CcosOr-363" > outgroup.lst
cat $clade.cds.out.FastTree.fas|gotree reroot outgroup -i - -l outgroup.lst > $clade.cds.out.FastTree.rooted.fas

```
## Create ```smap``` file for dlcpar
```bash
cat $base/cladeData/$clade.tsv |sed 1d|cut -f 1,6 |perl -pe 's/-.*\t/\*\t/g' |sort|uniq|awk -F"\t" '{OFS=FS}{print $1,toupper($2)}' > $clade.smap
```
## Run DLCpar
From the manual:

>The ilp command finds a most parsimonious gene tree-species tree reconciliation through integer linear programming. It is a useful alternative to dlcpar ilp when (1) the gene tree is too large or is highly incongruent to the species tree, as these cases may be too complex for the dynamic programming approach, or (2) you wish to limit the maximum runtime or memory when inferring an MPR.

## reformat trees so dlcpar search does accept them
```bash
mkdir DLC
cd DLC

ln -s ../RAxML_bipartitionsBranchLabels.rooted .
cat RAxML_bipartitionsBranchLabels.rooted |perl -pe 's/\[.*?\]//g'> tmp1.tre
python /corefac/cse/lukas/inqGen18/reconClade285/scripts/nodelabel.py tmp1.tre|perl -pe 's/\;/root\;/g' > $clade.raxML.tre


nice dlcpar search -S ../$clade.smap -s $base/data/species.tree $clade.raxML.tre
nice dlcpar events --3t -S ../$clade.smap -s $base/data/species.tree $clade.raxML.tre.dlcsearch.coal.tree > $clade.raxML.recon.tsv


ln -s ../$clade.cds.out.FastTree.rooted.fas .
cat $clade.cds.out.FastTree.rooted.fas> tmp2.tre
python $base/scripts/nodelabel.py tmp2.tre|perl -pe 's/\;/root\;/g' > $clade.FT.tre

nice dlcpar dp -S ../$clade.smap -s $base/data/species.tree $clade.FT.tre --output_format 3t

nice dlcpar search -S ../$clade.smap -s $base/data/species.tree $clade.FT.tre
nice dlcpar events --3t -S ../$clade.smap -s $base/data/species.tree $clade.FT.tre.dlcsearch.coal.tree > $clade.recon.tsv

nice dlcpar dp -S ../$clade.smap -s $base/data/species.tree $clade.FT.tre
```

# Get tree-relations file

```bash
PYTHONPATH=$PYTHONPATH/:/usr/local/home/lschrader/software/dlcoal-1.0/

for clade in $(ls $base/cladeData/*.tsv|perl -pe 's/.*\/(.*)\.tsv/$1/g')
do
  echo $clade
  cd $base/$clade/DLC
  ls $clade.FT.tre.dlcsearch.locus.tree| ~/software/dlcoal-1.0/bin/tree-relations     -d -s $base/data/species.tree -S ../$clade.smap -T .tree -R .recon > $clade.rel.tsv
done

grep "loss" 9E-Clade-285.rel.tsv|cut -f 2-4 > 9E-Clade-285.rel.losses.tsv

```


<!--
```bash
# dp mode is the full mode
ln -s ../RAxML_rootedTree.test .
nice dlcpar dp -S ../$clade.smap -s ../../species.tree tmp.tre  --output_format 3t
nice dlcpar search -S ../$clade.smap -s ../../species.tree tmp.tre  --output_format 3t
nice dlcpar events --3t -S ../$clade.smap -s ../../species.tree RAxML_rootedTree.test.dlcdp.coal.tree
nice dlcpar events --3t -S ../$clade.smap -s ../../species.tree tmp.tre.dlcsearch.coal.tree

nice dlcpar dp -S ../$clade.smap -s ../../species.tree ../$clade.cds.out.FastTree.rooted.fas  --output_format 3t
# ilp dies due to pulp error
#nice dlcpar ilp -S ../$clade.smap -s ../../species.tree ../$clade.cds.out.FastTree.rooted.fas  --output_format 3t
ln -s ../$clade.cds.out.FastTree.rooted.fas .
nice dlcpar search -S ../$clade.smap -s ../../species.tree V-Clade-258.cds.out.FastTree.rooted.fas
nice dlcpar dp -S ../$clade.smap -s ../../species.tree V-Clade-258.cds.out.FastTree.rooted.fas
nice dlcpar events --3t -S ../$clade.smap -s ../../species.tree V-Clade-258.cds.out.FastTree.rooted.fas.dlcsearch.coal.tree
```
#visualize results
```bash
dlcpar landscape \
    -S $clade.smap \
    -s ../species.tree RAxML_rootedTree.test.dlcdp.coal.tree \
    RAxML_rootedTree.test.dlcdp.coal.tree
```


## Compute gene tree phylogeny with RaxML (bootstrap and ML), root with ```-f I```
```bash
#raxmlHPC-PTHREADS -f I -m GTRGAMMA -p 12345 -x 12345 -# 500 -s $clade.cds.out.best.fas -n run1 -T 4
raxmlHPC -f I -m GTRGAMMA -t /corefac/cse/lukas/inqGen18/reconClade285/V-Clade-258/RAxML_bipartitionsBranchLabels.unrooted -n test -s $clade.cds.out.best.fas

```
