

# MCL gene family clustering

```bash
# 5 acromyrmex + Acol genomes

base=/corefac/cse/lukas/inqGen18/reannotation/MCL/
bd=/usr/local/home/lschrader/data/inqGen18/
mkdir $base/
mkdir $base/data
mkdir $base/blsDB

#cp -r tmp/geneFamilies/scripts/ geneFamilies/
```
# 1.0 Prepare peptide files for all genomes
```bash
# get soft links to original peptide files
cd $base
for file in $(find $bd/inquilines_v2.1/  -name "*.pep.fa")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file data/$short.pep.fa
done

for file in $(find ~/data/genomes/reannotations_of_7_ants/  -name "*Atta_colombica*.pep.fa" -o -name "*Acromyrmex_echinatior*.pep.fa" )
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file data/$short.pep.fa
done
```

# 2.0 BLAST ALL-VS-ALL with mmseqs
```bash

# clean up peptide files to make their structure identical between inquilines and attines
#cat $base/data/* |perl -pe 's/(>.*?)[\s|-].*/$1$2/g' > blsDB/allProteins.pep.fa
cat /corefac/cse/lukas/inqGen18/reannotation/results/*/GAF/*.GeMoMa.longestIsoform.pep.fa  > blsDB/allProteins.pep.fa

# Blast all 11 attine  proteomes against all ant proteomes

# Make blastdb
cd $base/blsDB/
mmseqs easy-search --start-sens 1  --sens-steps 2  -s 7 allProteins.pep.fa allProteins.pep.fa allvsall.out tmp
awk '($11 + 0) < 1E-5' allvsall.out > allvsall.cutoff.out

```

# 3.0 MCL Clustering

```bash

#http://micans.org/mcl/
#http://micans.org/mcl/man/clmdist.html
#e-value 1e-5
#alignment length cutoff of 50 %
#length cutoff of 25 aa
#  In the end, I didn't use an alignment cutoff. I only used an e-value cutoff of 1e-05.

cd $base/

mcxload -abc $base/blsDB/allvsall.cutoff.out --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o seq.mci -write-tab seq.tab

#tested different inflation parameters
# inflation of 1.5 works best
# -te <n> for multi threading
nice mcl seq.mci -I 1.5 -use-tab seq.tab -o seq.mci.I01.5 -te 10

awk -F'|' 'BEGIN{print "I01.5", "lineNum"}{print 1+gsub(/[A-Z]{4}[0-9]{5}/,"") "\t" NR}' seq.mci.I01.5


for i in ACOL AECH AINS AHEY ACHA PARG
do
  echo -e "rep\tOG\t$i" > $i.ogs.tsv
  awk -v OFS='\t' '{print $1,"OG"sprintf("%05d", NR),gsub(/'"$i"'/,"")}' seq.mci.I01.5 >> $i.ogs.tsv
done
paste *.ogs.tsv|cut -f 1-3,6,9,12,15,18 > geneFamilies.tsv
wc -l geneFamilies.tsv
awk '{if ($3==1 && $4==1 && $5==1 && $6==1 && $7==1 && $8==1) print $0}' geneFamilies.tsv |wc -l

# for i in ACOL AECH AINS AHEY ACHA PARG
# do
#   echo -e "rep\tOG\t$i" > $i.ogs.tsv
#   awk -v OFS='\t' '{print $1,"OG"sprintf("%05d", NR),gsub(/'"$i"'/,"")}' /usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl/seq.af50lf25.mci.I01.5 >> $i.ogs.tsv
# done
# paste *.ogs*|cut -f 1-3,6,9,12,15,18 > geneFamilies.tsv
# awk '{if ($3==1 && $4==1 && $5==1 && $6==1 && $7==1 && $8==1) print $0}' geneFamilies.tsv |wc -l


```

# 5.0 Incorporate TransposonPsi annotation in MCL cluster QC
```bash

##############################################################
# 5.1 Clean TE proteins from clusters
#     and at the same time compute some QC statistics for each MCL cluster
##############################################################

# cat all transposonPSI annotations
TEpsiOut=$base/../Tpsi/
cat $TEpsiOut/Acol*.allHits $TEpsiOut/Aech*.allHits $TEpsiOut/Ahey*.allHits $TEpsiOut/Ains*.allHits $TEpsiOut/Acha*.allHits $TEpsiOut/Parg*.allHits \
> $base/all.pep.TEpsi.allHits

cut -f 6 all.pep.TEpsi.allHits > all.pep.TEpsi.allHits.genes.lst
grep -v -f all.pep.TEpsi.allHits.genes.lst seq.mci.I01.5 > seq.mci.I01.5.TEfree

for i in ACOL AECH AINS AHEY ACHA PARG
do
  echo -e "rep\tOG\t$i" > $i.ogs.TEfree.tsv
  awk -v OFS='\t' '{print $1,"OG"sprintf("%05d", NR),gsub(/'"$i"'/,"")}' seq.mci.I01.5.TEfree >> $i.ogs.TEfree.tsv
done
paste *.ogs.TEfree.tsv|cut -f 1-3,6,9,12,15,18 > geneFamilies.TEfree.tsv
awk '{if ($3==1 && $4==1 && $5==1 && $6==1 && $7==1 && $8==1) print $0}' geneFamilies.TEfree.tsv |wc -l
awk '{if ($3==1 && $4==1 && $5==1 && $6==1 && $7==1 && $8==1) print $0}' geneFamilies.tsv |wc -l
wc -l geneFamilies.*


# Clean MCL output from TE genes and create QC file
#remove all TE genes from MCL clusters and produce QC files to later produce CAFE-ready gene family cluster files

#perl $base/scripts/TEpsiFilter.pl all.pep.TEpsi.allHits $base/results/mcl/seq.af50lf25.mci.I01.5 $base/mclQC/TEannotations/I01.5.qc.tsv  > $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.5

#grep ";EMPTY" $base/mclQC/TEannotations/I01.5.qc.tsv -v > $base/mclQC/TEannotations/I01.5.qc.emptiesRemoved.tsv

##############################################################
# Repeat QC for TEfiltered clusters
##############################################################

# prepare Ipr QC for filtered files
#perl $base/scripts/filterMCLcluster.pl $base/mclQC/funAnn/IPR.annotations.tsv $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.0 > $base/results/mcl.TEfiltered/QC.TEfiltered.I01.0.MCL.tsv
#perl $base/scripts/filterMCLcluster.pl $base/mclQC/funAnn/IPR.annotations.tsv $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.5 > $base/results/mcl.TEfiltered/QC.TEfiltered.I01.5.MCL.tsv
#perl $base/scripts/filterMCLcluster.pl $base/mclQC/funAnn/IPR.annotations.tsv $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I03.0 > $base/results/mcl.TEfiltered/QC.TEfiltered.I03.0.MCL.tsv

# prepare phyletic profile QC for filtered files
#perl $base/scripts/phyleticProfile.MCL.pl $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.0 $base/mclQC/funAnn/species.tsv > $base/results/mcl.TEfiltered/phyProf.TEfiltered.MCL.I01.0.tsv
#perl $base/scripts/phyleticProfile.MCL.pl $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.5 $base/mclQC/funAnn/species.tsv > $base/results/mcl.TEfiltered/phyProf.TEfiltered.MCL.I01.5.tsv
#perl $base/scripts/phyleticProfile.MCL.pl $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I03.0 $base/mclQC/funAnn/species.tsv > $base/results/mcl.TEfiltered/phyProf.TEfiltered.MCL.I03.0.tsv


##############################################################
# 5.2 Paste QC for each inflation parameter output of MCL
##############################################################

# paste phyletic profiles, QC and TE annotations for each MCL cluster
#cd $base/results/mcl.TEfiltered/
#paste phyProf.*I01.0.tsv QC.*I01.0.MCL.tsv $base/mclQC/TEannotations/I01.0.qc.emptiesRemoved.tsv > largeQC.TEpso.nofilter.MCL.I01.0.tsv
#paste phyProf.*I01.5.tsv QC.*I01.5.MCL.tsv $base/mclQC/TEannotations/I01.5.qc.emptiesRemoved.tsv > largeQC.TEpso.nofilter.MCL.I01.5.tsv
#paste phyProf.*I03.0.tsv QC.*I03.0.MCL.tsv $base/mclQC/TEannotations/I03.0.qc.emptiesRemoved.tsv > largeQC.TEpso.nofilter.MCL.I03.0.tsv


##############################################################
# 5.3. Run Rscript to prepare CAFE-ready output files
#       All MCL clusters with one or more TE proteins has been removed entirely
##############################################################
#$base/scripts/Ranalysis.TEfilteredMCL.R

```

```bash

###################################################
# MCL folder structure
###################################################


# CAFE ready file
/corefac/cse/lukas/inqGen18/reannotation/MCL/geneFamilies.TEfree.tsv
# Copy file to pallas
```
