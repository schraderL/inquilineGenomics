
##############################################################
# Incorporate TransposonPsi annotation in MCL cluster QC
##############################################################

base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/
bd=/usr/local/home/lschrader/data/inqGen18/
mkdir $base/mclQC/TEannotations
cd $base/mclQC/TEannotations
mkdir $base/results/mcl.TEfiltered
mkdir $base/results/mcl.TEfiltered

# see $base/mcl.overview.txt for details of all files produced here.

##############################################################
# Clean TE proteins from clusters
#    and at the same time compute some QC statistics for each MCL cluster
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
# Paste QC for each inflation parameter output of MCL
##############################################################

# paste phyletic profiles, QC and TE annotations for each MCL cluster
cd $base/results/mcl.TEfiltered/
paste phyProf.*I01.0.tsv QC.*I01.0.MCL.tsv $base/mclQC/TEannotations/I01.0.qc.emptiesRemoved.tsv > largeQC.TEpso.nofilter.MCL.I01.0.tsv
paste phyProf.*I01.5.tsv QC.*I01.5.MCL.tsv $base/mclQC/TEannotations/I01.5.qc.emptiesRemoved.tsv > largeQC.TEpso.nofilter.MCL.I01.5.tsv
paste phyProf.*I03.0.tsv QC.*I03.0.MCL.tsv $base/mclQC/TEannotations/I03.0.qc.emptiesRemoved.tsv > largeQC.TEpso.nofilter.MCL.I03.0.tsv


##############################################################
# Run Rscript to prepare CAFE-ready output files
# All MCL clusters with one or more TE proteins has been removed entirely
##############################################################
$base/scripts/Ranalysis.TEfilteredMCL.R
