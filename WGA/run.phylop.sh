########################################################################
# Compute neutral model for full alignment
########################################################################

base=~/data/inqGen18/WGA/PHAST
cd $base
source /usr/local/home/lschrader/software/cactus/progressiveCactus/environment

cat Acromyrmex_echinatior.v2.0.gff|awk -F'\t' -vOFS='\t' '{if ($3=="CDS") print $0;gsub("CDS", "exon", $3);print $0}' > tmp.gff
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff
gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -addintrons tmp_clean.gff > Acromyrmex_echinatior.v2.0.clean.gff
gff3ToGenePred Acromyrmex_echinatior.v2.0.clean.gff tmp.genePred
genePredToBed tmp.genePred Acromyrmex_echinatior.v2.0.clean.bed

nice halPhyloPTrain.py attines.hal AECH Acromyrmex_echinatior.v2.0.clean.bed neutralMod.progCac.mod --numProc 12
# NEUTRAL MODEL AT
## /corefac/cse/lukas/inqGen18/WGA/PHAST/neutralMod.progCac.mod

nice halTreePhyloP.py attines.hal test.mod ./outdir/ --bigWig --numProc 12

########################################################################
# Convert hal to maf, remove parasites, split by scaffold (over 500kb)
########################################################################
base=/usr/local/home/lschrader/data/inqGen18/WGA/hal2maf
ln -s ~/data/inqGen18/WGA/progressiveCactus/attines/attines.hal .

# Loop over AECH scaffolds
cut -f 1,2 Acromyrmex_echinatior.v2.0.fa.fai |sort -k2,2 -nr|awk '{if ($2>500000) print $1}' > Aech.over500kb.lst
cd $base/rawMaf/

# get all species
# TCOR,TSEP,AHEY,ACHA,PARG,AECH,AINS,ACEP,ACOL,TZET,CCOS,

cat ../Aech.over500kb.lst|parallel --nice 10 "hal2maf ../attines.hal Aech.{}.attines.maf --refGenome AECH --noAncestors --refSequence {}"

########################################################################
# Filter MAF
########################################################################
base=/usr/local/home/lschrader/data/inqGen18/WGA/hal2maf
aln=/corefac/cse/lukas/inqGen18/WGA/hal2maf/AECH
source /usr/local/home/lschrader/software/cactus/progressiveCactus/environment

cd $base/
# 1. remove duplicates and all regions not aligned across all Species
# 2. remove repeats (according to soft-masking of genomes)
# 3. save file to ./noDups/Aech.*.ND.maf
# 4. remove parasites
# 5. save file to ./noDupsNoP/Aech.*.ND.NP.maf
cat $base/Aech.over500kb.lst|grep scaffold585|parallel --nice 10 '\
maffilter input.file=/corefac/cse/lukas/inqGen18/WGA/hal2maf/rawMaf/Aech.{}.attines.maf \
input.file.compression=none input.format=Maf output.log=./maffilterLogs/{}.log \
maf.filter="\
    Subset(\
         species=(TCOR,TSEP,AHEY,ACHA,PARG,AECH,AINS,ACEP,ACOL,TZET,CCOS),\
         strict=yes,                        \
         keep=no,                           \
         remove_duplicates=yes),            \
    MaskFilter(                             \
         species=(TCOR,TSEP,AHEY,ACHA,PARG,AECH,AINS,ACEP,ACOL,TZET,CCOS),\
         window.size=10,                     \
         window.step=1,                      \
         max.masked=2),                      \
   MinBlockLength(                          \
   min_length=1),\
    Output(                                 \
         file=./noDups/Aech.{}.ND.maf,           \
         compression=none,              \
         mask=no),                     \
    Subset(\
         species=(TCOR,TSEP,AHEY,AECH,ACEP,ACOL,TZET,CCOS),\
         strict=yes,                        \
         keep=no,                           \
         remove_duplicates=no),            \
    Output(                                 \
        file=./noDupsNoP/Aech.{}.ND.NP.maf,                \
        compression=none,                   \
        mask=no)"'

########################################################################
# Run phast on each alignment using R
########################################################################


1. remove parasites from hal alignment
2. calculate global neutral model using halPhyloPTrain.py
3. run rPHAST on parasite-free alignments for each scaffold (see phyloP.nonParasites.R)
3a. read alignment, gff and tree
3b. read global neutral Model
3c. run phastCons()
3d. run phyloP()
4. run PAR analysis (see PAR.R)
4a. Identify elements rapidly evolving in inquiline
4b. read full alignment, gff and tree, consElements.bed,
4c. Calculate neutral model again?
4d. Filter well-aligned regions
4e. Run phyloP() on inquiline branches
4f. Simulate alignments
