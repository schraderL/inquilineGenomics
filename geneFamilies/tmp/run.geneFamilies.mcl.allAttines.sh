###################################################
# MCL gene family clustering
###################################################

# MOVED FOLDER TO /corefac/cse/lukas/inqGen18/tmp/geneFamilies

# 11 attine genomes

base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/
bd=/usr/local/home/lschrader/data/inqGen18/
mkdir $base/data
mkdir $base/blsDB

###################################################
# 1.0 Prepare peptide files for all 11 attines
###################################################
# get soft links to original peptide files
for file in $(find $bd/inquilines_v2.1/  -name "*.pep.fa")
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file data/$short.pep.fa
done

for file in $(find ~/data/genomes/reannotations_of_7_ants/  -name "*.pep.fa" )
do
  short=$(echo $file|rev|cut -f 1 -d "/" |rev|perl -pe 's/(.).*?_(...).*/$1$2/g')
  echo $short with $file
  ln -s $file data/$short.pep.fa
done


###################################################
# 2.0 BLAST ALL-VS-ALL
###################################################
# clean up peptide files to make their structure identical between inquilines and attines
cat $base/data/* |perl -pe 's/(>.*?)[\s|-].*/$1$2/g' > blsDB/11attines.pep.fa

# Blast all 11 attine  proteomes against all ant proteomes

# Make blastdb
cd $base/blsDB/
makeblastdb -in 11attines.pep.fa -dbtype prot

#Run Blastp in parallel
cat blsDB/11attines.pep.fa| nice parallel -j20 -k --block 5k --recstart '>' --pipe 'blastp -db blsDB/11attines.pep.fa -query - -outfmt 6 -seg yes -evalue 1e-5' > results/11attines.pep.bls

# -k = Keep sequence of output same as the order of input. Normally the output of a job will be printed as soon as the job completes.
# The block size is determined by --block.
# The strings --recstart and --recend tell GNU parallel how a record starts and/or ends.
# The block read will have the final partial record removed before the block is passed on to the job. The partial record will be prepended to next block.

###################################################
# 2.1 Calculate protein length for alignment length cut-off
###################################################

# This information can be used to apply an alignment cutoff to filter the blast results

bioawk -c fastx '{ print $name, length($seq) }' $base/blsDB/11attines.pep.fa > $base/blsDB/11attines.pep.length

###################################################
# 2.2 Filter Blast results
###################################################

#Filter blast results based on alignment length

perl $base/scripts/addGeneLength2BLAST.pl $base/blsDB/11attines.pep.length $base/results/11attines.pep.bls > $base/results/11attines.pep.extended.bls

#unfiltered
awk -v OFS='\t' '{print $1,$2,$11}' $base/results/11attines.pep.extended.bls > $base/results/11attines.pep.unfiltered.bls

# 50 % alignment filter
awk -v OFS='\t' '{if ($15>.50 && $16>.50) print $1,$2,$11}' $base/results/11attines.pep.extended.bls > $base/results/11attines.pep.alnFilter50.bls

# 50 % alignment filter & 25 length filter
awk -v OFS='\t' '{if ($13>25 && $14>25 &&$15>.50 && $16>.50)  print $1,$2,$11}' $base/results/11attines.pep.extended.bls > $base/results/11attines.pep.alnFilter50lengthFilter25.bls


###################################################
# 3.0 MCL Clustering
###################################################

#http://micans.org/mcl/
#http://micans.org/mcl/man/clmdist.html
#e-value 1e-5
#alignment length cutoff of 50 %
#length cutoff of 25 aa
#  In the end, I didn't use an alignment cutoff. I only used an e-value cutoff of 1e-05.

mcxload -abc $base/results/11attines.pep.alnFilter50lengthFilter25.bls --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o seq.af50lf25.mci -write-tab seq.af50lf25.tab

#tested different inflation parameters
# inflation of 1.5 works best
# -te <n> for multi threading
mkdir results/mcl/
cd $base/results
nice mcl seq.af50lf25.mci -I 1.0 -use-tab seq.af50lf25.tab -o $base/results/mcl/seq.af50lf25.mci.I01.0 -te 10
nice mcl seq.af50lf25.mci -I 1.5 -use-tab seq.af50lf25.tab -o $base/results/mcl/seq.af50lf25.mci.I01.5 -te 10
nice mcl seq.af50lf25.mci -I 3.0 -use-tab seq.af50lf25.tab -o $base/results/mcl/seq.af50lf25.mci.I03.0 -te 10


###################################################
# 4.0 QC of different inflation parameters
###################################################

# The QC of the protein family clustering is based on interpro annotations, transposonPSI annotations and phyletic profiles

# Calculate
cd $base/results/mcl
awk -F'|' 'BEGIN{print "I01.0", "lineNum"}{print 1+gsub(/[A-Z]{4}[0-9]{5}/,"") "\t" NR}' $base/results/mcl/seq.af50lf25.mci.I01.0 > $base/mclQC/FamSizeI1.0.tsv
awk -F'|' 'BEGIN{print "I01.5", "lineNum"}{print 1+gsub(/[A-Z]{4}[0-9]{5}/,"") "\t" NR}' mcl/seq.af50lf25.mci.I01.5 > $base/mclQC/FamSizeI1.5.tsv
awk -F'|' 'BEGIN{print "I03.0", "lineNum"}{print 1+gsub(/[A-Z]{4}[0-9]{5}/,"") "\t" NR}' mcl/seq.af50lf25.mci.I03.0 > $base/mclQC/FamSizeI3.0.tsv

paste $base/mclQC/FamSize* | head -n 17075|awk -v OFS='\t' '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19}' > $base/mclQC/qc.diffI.tsv

mkdir $base/mclQC/funAnn
cd $base/mclQC/funAnn
inqFunann=$bd/inquilines_v2.1/*/*/function_annotation
attineFunann=~/data/genomes/reannotations_of_7_ants/*/*/function_annotation

#cat all IPR annotations
cat $inqFunann/*.gene.ipr $attineFunann/*.gene.ipr > $base/mclQC/funAnn/IPR.annotations.tsv

# Run filterMCLcluster.pl to get an idea of the IPR annotations in each cluster
cat $base/results/mcl/seq.af50lf25.mci.I01.0|perl $base/scripts/filterMCLcluster.pl $base/mclQC/funAnn/IPR.annotations.tsv - > $base/mclQC/funAnn/QC.01.0.MCL.tsv
cat $base/results/mcl/seq.af50lf25.mci.I01.5|perl $base/scripts/filterMCLcluster.pl $base/mclQC/funAnn/IPR.annotations.tsv - > $base/mclQC/funAnn/QC.01.5.MCL.tsv
cat $base/results/mcl/seq.af50lf25.mci.I03.0|perl $base/scripts/filterMCLcluster.pl $base/mclQC/funAnn/IPR.annotations.tsv - > $base/mclQC/funAnn/QC.03.0.MCL.tsv

# Run phyleticProfile.MCL.pl to get phyletic profiles for each family.
cut -c 1-4 IPR.annotations.tsv |uniq|tr "\n" "\t" > species.tsv
perl $base/scripts/phyleticProfile.MCL.pl $base/results/mcl/seq.af50lf25.mci.I01.0 $base/mclQC/funAnn/species.tsv > $base/mclQC/funAnn/phyProf.MCL.I01.0.tsv
perl $base/scripts/phyleticProfile.MCL.pl $base/results/mcl/seq.af50lf25.mci.I01.5 $base/mclQC/funAnn/species.tsv > $base/mclQC/funAnn/phyProf.MCL.I01.5.tsv
perl $base/scripts/phyleticProfile.MCL.pl $base/results/mcl/seq.af50lf25.mci.I03.0 $base/mclQC/funAnn/species.tsv > $base/mclQC/funAnn/phyProf.MCL.I03.0.tsv

#paste IPR annotations overviews and phyletic profiles
cd $base/mclQC/funAnn/

paste $base/mclQC/funAnn/phyProf.MCL.I01.0.tsv $base/mclQC/funAnn/QC.01.0.MCL.tsv >  $base/mclQC/funAnn/largeQC.I01.0.tsv
paste $base/mclQC/funAnn/phyProf.MCL.I01.5.tsv $base/mclQC/funAnn/QC.01.5.MCL.tsv >  $base/mclQC/funAnn/largeQC.I01.5.tsv
paste $base/mclQC/funAnn/phyProf.MCL.I03.0.tsv $base/mclQC/funAnn/QC.03.0.MCL.tsv >  $base/mclQC/funAnn/largeQC.I03.0.tsv

# some potentially useful one-liners
head -n1 $base/results/mcl/seq.af50lf25.mci.I01.0 |tail|cut -f 40 -d ";"|sed s/","/" "/g |awk 'BEGIN { OFS = "\n" } { $1=$1; print }'|xargs -I --  grep -- $base/mclQC/funAnn/IPR.annotations.tsv



##############################################################
# 5.0 Incorporate TransposonPsi annotation in MCL cluster QC
##############################################################

base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/
bd=/usr/local/home/lschrader/data/inqGen18/
mkdir $base/mclQC/TEannotations
cd $base/mclQC/TEannotations
mkdir $base/results/mcl.TEfiltered
mkdir $base/results/mcl.TEfiltered

# see $base/mcl.overview.txt for details of all files produced here.

##############################################################
# 5.1 Clean TE proteins from clusters
#     and at the same time compute some QC statistics for each MCL cluster
##############################################################

# cat all transposonPSI annotations
cat $bd/QC/TransposonPSI/output/*.allHits > $base/mclQC/TEannotations/all.pep.TEpsi.allHits

# Clean MCL output from TE genes and create QC file
#remove all TE genes from MCL clusters and produce QC files to later produce CAFE-ready gene family cluster files
perl ../../scripts/transposonPSI/TEpsiFilter.pl ../all.pep.TEpsi.allHits ../../mcl/unfiltered/allVSall.unfiltered.mci.-- TEfiltered.allVSall.unfiltered.mci.--.qc  > TEfiltered.allVSall.unfiltered.mci.--
perl $base/scripts/TEpsiFilter.pl all.pep.TEpsi.allHits $base/results/mcl/seq.af50lf25.mci.I01.0 $base/mclQC/TEannotations/I01.0.qc.tsv  > $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.0
perl $base/scripts/TEpsiFilter.pl all.pep.TEpsi.allHits $base/results/mcl/seq.af50lf25.mci.I01.5 $base/mclQC/TEannotations/I01.5.qc.tsv  > $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.5
perl $base/scripts/TEpsiFilter.pl all.pep.TEpsi.allHits $base/results/mcl/seq.af50lf25.mci.I03.0 $base/mclQC/TEannotations/I03.0.qc.tsv  > $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I03.0

grep ";EMPTY" $base/mclQC/TEannotations/I01.0.qc.tsv -v > $base/mclQC/TEannotations/I01.0.qc.emptiesRemoved.tsv
grep ";EMPTY" $base/mclQC/TEannotations/I01.5.qc.tsv -v > $base/mclQC/TEannotations/I01.5.qc.emptiesRemoved.tsv
grep ";EMPTY" $base/mclQC/TEannotations/I03.0.qc.tsv -v > $base/mclQC/TEannotations/I03.0.qc.emptiesRemoved.tsv

##############################################################
# Repeat QC for TEfiltered clusters
##############################################################

# prepare Ipr QC for filtered files
perl $base/scripts/filterMCLcluster.pl $base/mclQC/funAnn/IPR.annotations.tsv $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.0 > $base/results/mcl.TEfiltered/QC.TEfiltered.I01.0.MCL.tsv
perl $base/scripts/filterMCLcluster.pl $base/mclQC/funAnn/IPR.annotations.tsv $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.5 > $base/results/mcl.TEfiltered/QC.TEfiltered.I01.5.MCL.tsv
perl $base/scripts/filterMCLcluster.pl $base/mclQC/funAnn/IPR.annotations.tsv $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I03.0 > $base/results/mcl.TEfiltered/QC.TEfiltered.I03.0.MCL.tsv

# prepare phyletic profile QC for filtered files
perl $base/scripts/phyleticProfile.MCL.pl $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.0 $base/mclQC/funAnn/species.tsv > $base/results/mcl.TEfiltered/phyProf.TEfiltered.MCL.I01.0.tsv
perl $base/scripts/phyleticProfile.MCL.pl $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I01.5 $base/mclQC/funAnn/species.tsv > $base/results/mcl.TEfiltered/phyProf.TEfiltered.MCL.I01.5.tsv
perl $base/scripts/phyleticProfile.MCL.pl $base/results/mcl.TEfiltered/TEfiltered.seq.af50lf25.mci.I03.0 $base/mclQC/funAnn/species.tsv > $base/results/mcl.TEfiltered/phyProf.TEfiltered.MCL.I03.0.tsv


##############################################################
# 5.2 Paste QC for each inflation parameter output of MCL
##############################################################

# paste phyletic profiles, QC and TE annotations for each MCL cluster
cd $base/results/mcl.TEfiltered/
paste phyProf.*I01.0.tsv QC.*I01.0.MCL.tsv $base/mclQC/TEannotations/I01.0.qc.emptiesRemoved.tsv > largeQC.TEpso.nofilter.MCL.I01.0.tsv
paste phyProf.*I01.5.tsv QC.*I01.5.MCL.tsv $base/mclQC/TEannotations/I01.5.qc.emptiesRemoved.tsv > largeQC.TEpso.nofilter.MCL.I01.5.tsv
paste phyProf.*I03.0.tsv QC.*I03.0.MCL.tsv $base/mclQC/TEannotations/I03.0.qc.emptiesRemoved.tsv > largeQC.TEpso.nofilter.MCL.I03.0.tsv


##############################################################
# 5.3. Run Rscript to prepare CAFE-ready output files
#       All MCL clusters with one or more TE proteins has been removed entirely
##############################################################
$base/scripts/Ranalysis.TEfilteredMCL.R



###################################################
# MCL folder structure
###################################################

~/data/inqGen18/geneFamilies/results
:'
├── mclQC                                             # QC files for MCL clusters (TEs not filtered)
│   ├── FamSizeI1.5.tsv                               # MCL family sizes for I=1.0
│   ├── funAnn
│   │   ├── IPR.annotations.tsv                       # all 11 IPR annotations
│   │   ├── largeQC.I01.5.tsv                         # full QC for unfiltered MCL clusters
│   │   ├── phyProf.MCL.I01.5.tsv                     # phyletic profile for "
│   │   ├── QC.01.0.MCL.tsv                           # IPR QC for           "
│   │   └── species.tsv                               # list of all analysed species
│   ├── qc.diffI.tsv                                  # comparison of fam sizes for I1.0, I1.5 and I3.0
│   └── TEannotations                                 # TransposonPSI results and QC
│       ├── all.pep.TEpsi.allHits                     # TransposonPSI annotations
│       ├── I01.5.qc.emptiesRemoved.tsv               # table of all clusters and their respective TE contaminations
│       └── I01.5.qc.tsv                              # tmp of "
├── results
│   ├── 11attines.pep.alnFilter50lengthFilter25.bls   # Filtered Blastp results
│   ├── 11attines.pep.bls                             # Uniltered Blastp results
│   ├── 11attines.pep.extended.bls                    # extended Blastp results (aln length)
│   ├── mcl
│   │   ├── seq.af50lf25.mci.I01.5                    # MCL Cluster results
│   ├── mcl.TEfiltered                                # TE filtered MCL Clusters (all TE proteins removed, but TE contaminated families kept) and QC
│   │   ├── largeQC.TEpso.nofilter.MCL.I01.5.tsv      # full QC of all clusters
│   │   ├── phyProf.TEfiltered.MCL.I01.5.tsv          # phyletic profile of TE filtered clusters
│   │   ├── QC.TEfiltered.I01.5.MCL.tsv               # IPR QC of TE filtered clusters
│   │   └── TEfiltered.seq.af50lf25.mci.I01.5         # MCL clusters with TE proteins removed (but TE contaminated families kept)
    ├── mcl.clean                                     # Separating TE-infested and TE-free gene families
    │   ├── BasicStats.pdf                            # Basic stats for different Inflation parameters
    │   ├── finalQC.I01.5.tsv                         # Final and complete QC overview for all clusters for preferred inflation parameter
    │   ├── hierachicalClustering.TE-free.pdf         # hierachical Clustering of Cluster sizes. TE-free clusters only
    │   ├── hierachicalClustering.TE-infested.pdf     # hierachical Clustering of Cluster sizes. TE-infested clusters only
    │   ├── largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv # cafe-ready input file
    │   ├── mds.TE-free.pdf                           # MDS analysis of Cluster sizes. TE-free clusters only
    │   └── mds.TE-infested.pdf                       # MDS analysis of Cluster sizes. TE-free clusters only
'

# CAFE ready file
$base/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.tsv
# Copy file to pallas
