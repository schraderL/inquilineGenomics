# Acromyrmex inquiline social parasite comparative genomics
This repository holds all the code that was used in the _Acromyrmex_ inquiline comparative genomics project.

* * *

These files hold the workflows to analyse the basic quality of the genomes used in the comparison.

##### Local and remote copies of files and folders



```bash
# mac
macbase=/Users/lukas/sciebo/inquilineGenomics18/
# pallas
base=/corefac/cse/lukas/inqGen18/
```

#### Raw assemblies & annotations (v.2.0)

The latest assemblies and annotations delivered from BGI are located on Pallas. The raw files are located at `/usr/local/home/lschrader/data/genomes/inquilines_2.0/`

```bash
Acromyrmex charruanus:  ./code/Acha.files.txt
Acromyrmex insinuator:  ./code/Ains.files.txt
Acromyrmex heyeri:      ./code/Ahey.files.txt
Pseudoatta argentina:   ./code/Parg.files.txt
```

#### assembly v2.1/annotation v2.1

Clean genomes and annotations of fungal and mitochondrial scaffolds.

```bash
./code/QC/filterContaminants/inquilines2.0-to-2.1.sh

# location on pallas
$base/inquilines_v2.1/
$base/inquilines_v2.1/Acromyrmex_charruanus.2.1/Acromyrmex_charruanus.v2.1.files.tree
$base/inquilines_v2.1/Acromyrmex_heyeri.2.1/Acromyrmex_heyeri.v2.1.files.tree
$base/inquilines_v2.1/Acromyrmex_insinuator.2.1/Acromyrmex_insinuator.v2.1.files.tree
$base/inquilines_v2.1/Pseudoatta_argentina.2.1/Pseudoatta_argentina.v2.1.files.tree
```

## QC

Code, programs and scripts used for quality control of assemblies and annotations are provided below.

#### BUSCO

BUSCO (v3) was run on v2.1 assembly and annotation of the new genomes and on the reannotated previously published attine genomes

```bash
# newly sequenced genomes
./code/QC/busco/runBUSCO_inq.v2.1.sh

# previously published attines
./code/QC/busco/runBUSCO_7ants.sh
```

#### Quast

Quast (v5.0.2) was run to assess N50 contig and scaffold sizes and other assembly parameters

```bash
./code/QC/N50/runQuast.sh
```

#### TransposonPSI

TransposonPSI (v.08222010) was run to identify putative TE-encoded proteins.

```bash
./code/QC/transposonPSI/runTransposonPsi.sh
```

### GAG

GAG v2.0.1 was run on all newly sequenced genomes (see `./code/QC/filterContaminants/inquilines2.0-to-2.1.sh`) and also on reannoated attines.

```bash
./code/QC/gag/runGAG_7attines.sh
```

## Ortholog inference

I used orthofinder (v.2.2.6) to identify orthologs and orthogroups between different species.

```bash
./code/orthology/runOrthofinder.sh
```

## Species Phylogeny

#### Orthologous 4-fold degenerate site alignment

Alignment of single copy orthologs, testing for recombination and extraction of 4d sites.

```bash
./code/singleCopyOrthologs/runSCOalignment.sh
```

#### Phylogenetic inference

Phylogenetic inference with PAUP, PhyML, MrBayes, and RaxML

```bash
./code/phylogeny/runPhylogeny.sh
```

### Tree divergence dating

```bash
./code/TreeDating/treeDating.workflow.sh
```

## Gene Family Evolution

The complete workflow for analysing gene family evolution with CAFE is described here:

```bash
./code/geneFamilies/README.md

# This readme refers to the following scripts:
  run.geneFamilies.mcl.sh
  run.CAFError.sh
  run.CAFE.sh
  analyseCAFEmodels.Rmd
  CAFE.lhtest.R
  analyseCAFEoutput.Rmd
```

## Evolutionary rate analyses
Evolutionary rate analyses and signatures of selection (SOS) on single copy orthologs with HyPhy (v2.3.14)

```
./code/SOS/AnalyseDNDS.hyphy.sh
./code/SOS/RELAX.summary.Rmd
./code/SOS/absREL.summary.Rmd
./code/SOS/analyseSOS.Rmd
```
## Gene Ontology (GO) enrichment
GO enrichment analyses were run on genes under relaxed, positive, and identified selection and on shrinking gene families.
```bash
./code/GOenrichment/RELAX.GOenrichment.Rmd
./code/GOenrichment/absREL.GOenrichment.Rmd
./code/GOenrichment/CAFE.GOenrichment.Rmd
# scripts for revigo plots
./code/GOenrichment/revigoPlots/
```

## Whole genome alignment
Whole genome alignment with progressiveCactus v0.1
```bash
./code/WGA/run.progressiveCactus.sh
```

## Synteny Analysis
Synteny was analyzed with i-adhore v3.0.1
```bash
./code/synteny/run.iadhore.sh
./code/synteny/analyse.synteny.Rmd

# analysis of annotated OR genes only
./code/synteny/run.iadhore.ORs.sh
./code/synteny/plotOR.synteny.Rmd

```

## Analysis of the attine OR gene family

#### OR annotation

The full annotation pipeline is contained in

```bash
./code/ORevolution/clean/ORannotation.190408.sh
```

The following scipts are called in the above bash script:

```bash
ORannotation.190408.sh
filterForDomain.sh
runGeMoMa.pipe.sh
runExonerate.sh
runEVM.sh
prepareEVM.sh
addIncompleteGenes.sh
rename.R
runExonBlast.sh
runGeMoMa.sh
runExonerate.pl
processExonerate2Fa.pl
cleanExonerateGFF.pl
translateCDS2AA.pl
processExonerate2GFF.pl
addGeneGFF.pl
identify_geneLS.pl
hsp2bed.pl
```

```bash
# final annotations appropriate for v2.1 and published attine genomes
./ORevolution/results/annotations/
```

#### OR phylogeny

Scripts for OR alignments, alignment QC, tree inference, tree plotting, and analyses.

```bash
./code/ORphylogeny/ORphylogeny.sh
./code/ORphylogeny/ORplotTreeClean.Rmd
./code/ORphylogeny/analyseTree.Rmd
```

#### Evolutionary rate analyses in ORs with HyPhy (v2.3.14)

```bash
./code/ORevolution/OR.SOS/hyphy.ORs.sh
./code/ORevolution/OR.SOS/summary.RELAX.OR.Rmd
./code/ORevolution/OR.SOS/summary.absREL.OR.Rmd

# Detailed analysis of clade 285
./code/ORevolution/OR.SOS/analyseClade285.Rmd
```

## Annnotation & analysis of other gene families

```bash
# MRJPs
./code/genePhylogenies/MRJPevolution/MRJP.annotation.sh
./code/genePhylogenies/MRJPevolution/MRJP.phylogeny.sh
./code/genePhylogenies/MRJPevolution/plotMRJP.tree.Rmd
./code/genePhylogenies/MRJPevolution/mrjp.analyseTree.Rmd

# Gustatory Receptors (GR)
./code/genePhylogenies/GRevolution/annotateGRs.sh
./code/genePhylogenies/GRevolution/phylogenyGRs.sh
./code/genePhylogenies/GRevolution/plotGRtree.Rmd

# Elongases
./code/genePhylogenies/ELONGASEevolution/annotate.ELNG.sh
./code/genePhylogenies/ELONGASEevolution/phylogeny.ELNG.sh
./code/genePhylogenies/ELONGASEevolution/plotELGtree.Rmd

# Cuticular Proteins (CPRs)
./code/genePhylogenies/CPRevolution/annotateCPR.sh
./code/genePhylogenies/CPRevolution/phylogenyCPR.sh
./code/genePhylogenies/CPRevolution/plotCPRtree.Rmd
```

#### Script to plot barplots for all analyzed gene families

```bash
./code/genePhylogenies/barplot.Rmd
```
