This repository holds all the code that was used in the Acromyrmex inquiline comparative genomics project.

# QC
These files hold the workflows to analyse the basic quality of the genomes used in the comparison.

#####################
# inquilineGenomics18
#####################

# Proper documentation of the attine inquiline comparative genomics project.
# mac
macbase=/Users/lukas/sciebo/inquilineGenomics18/
# pallas
base=/corefac/cse/lukas/inqGen18/

#####################
# Files
#####################

# The newly generated assemblies and annotations are located on Pallas. The raw files are located at /usr/local/home/lschrader/data/genomes/inquilines_2.0/

#####################
### On Pallas
#####################

##########################################
# assembly v2.0/annotation v2.0
##########################################

Acromyrmex charruanus:  ./Acha.files.txt
Acromyrmex insinuator:  ./Ains.files.txt
Acromyrmex heyeri:      ./Ahey.files.txt
Pseudoatta argentina:   ./Parg.files.txt

##########################################
# assembly v2.1/annotation v2.1
##########################################

# clean genomes and annotations of fungal and mitochondrial scaffolds.
# mac
./QC/filterContaminants/inquilines2.0-to-2.1.sh

# pallas
$base/inquilines_v2.1/
$base/inquilines_v2.1/Acromyrmex_charruanus.2.1/Acromyrmex_charruanus.v2.1.files.tree
$base/inquilines_v2.1/Acromyrmex_heyeri.2.1/Acromyrmex_heyeri.v2.1.files.tree
$base/inquilines_v2.1/Acromyrmex_insinuator.2.1/Acromyrmex_insinuator.v2.1.files.tree
$base/inquilines_v2.1/Pseudoatta_argentina.2.1/Pseudoatta_argentina.v2.1.files.tree

#####################
# QC
#####################

#####################
# BUSCO
#####################

# run BUSCO on v2.1 assembly and annotation of the new genomes
./QC/busco/runBUSCO_inq.v2.1.sh
$base/QC/BUSCO_inquilinesv2.1/runBUSCO_inq.v2.1.sh

# run BUSCO on each assembly and annotation of the 7 reannotated ants
./QC/busco/runBUSCO_7ants.sh
$base/QC/BUSCO_7ants
#####################

#####################
# TransposonPSI
#####################
./QC/transposonPSI/runTransposonPsi.sh
ls $base/QC/TransposonPSI/inquilines_2.1/
#####################


#####################
# Orthologs
#####################

I used orthofinder to identify orthologs and orthogroups
# orthofinder
./code/orthology/runOrthofinder.sh

#####################
# OR annotations
#####################

# all scripts for the OR annotation pipeline are found at
./ORevolution

$base/ORs/

# final annotations appropriate for v2.1 are found here
ls $base/ORs/Acha/final/finalSet
ls $base/ORs/Ahey/final/finalSet
ls $base/ORs/Ains/final/finalSet
ls $base/ORs/Parg/final/finalSet


#####################
# OR phylogeny
#####################

./ORphylogeny/ORphylogeny.trimal.sh

#####################
# MRJP phylogeny
#####################

./MRJPevolution/specificGeneFamilies.Analysis.sh

#####################
# Find and align single copy orthologs
#####################

./singleCopyOrthologs/runSCOalignment.sh

#####################
# Species Phylogeny
#####################

./phylogeny/runPhylogeny.sh

#####################
# Tree divergence dating
#####################

./TreeDating/treeDating.workflow.sh

#####################
# Gene Family Evolution with CAFE
#####################

./geneFamilies/README.txt

#####################
# whole genome alignment
#####################

./WGA/run.progressiveCactus.sh
