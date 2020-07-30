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
Acromyrmex charruanus:  ./Acha.files.txt
Acromyrmex insinuator:  ./Ains.files.txt
Acromyrmex heyeri:      ./Ahey.files.txt
Pseudoatta argentina:   ./Parg.files.txt
```

#### assembly v2.1/annotation v2.1

[Clean genomes and annotations of fungal and mitochondrial scaffolds.](./QC/filterContaminants/inquilines2.0-to-2.1.sh)

```bash
./QC/filterContaminants/inquilines2.0-to-2.1.sh

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

BUSCO (v3) was run on v2.1 assembly and annotation of [the new genomes](./QC/busco/runBUSCO_inq.v2.1.sh) and on the [reannotated previously published attine genomes](./QC/busco/runBUSCO_7ants.sh).

```bash
# newly sequenced genomes
./QC/busco/runBUSCO_inq.v2.1.sh

# previously published attines
./QC/busco/runBUSCO_7ants.sh
```

#### Quast

[Quast (v5.0.2) was run](./QC/N50/runQuast.sh) to assess N50 contig and scaffold sizes and other assembly parameters.

```bash
./QC/N50/runQuast.sh
```

#### TransposonPSI

[TransposonPSI (v.08222010) was run](./QC/transposonPSI/runTransposonPsi.sh) to identify putative TE-encoded proteins.

```bash
./QC/transposonPSI/runTransposonPsi.sh
```

### GAG

[GAG v2.0.1 was run](./QC/gag/runGAG_7attines.sh) on all newly sequenced genomes (see `./QC/filterContaminants/inquilines2.0-to-2.1.sh`) and also on reannoated attines.

```bash
./QC/gag/runGAG_7attines.sh
```

## Ortholog inference

I used orthofinder (v.2.2.6) to [identify orthologs and orthogroups between different species](./orthology/runOrthofinder.sh).

```bash
./orthology/runOrthofinder.sh
```

## Species Phylogeny

#### Orthologous 4-fold degenerate site alignment

[Alignment of single copy orthologs, testing for recombination and extraction of 4d sites](./singleCopyOrthologs/runSCOalignment.sh).

```bash
./singleCopyOrthologs/runSCOalignment.sh
```

#### Phylogenetic inference

[Phylogenetic inference](./phylogeny/runPhylogeny.sh) with PAUP, PhyML, MrBayes, and RaxML

```bash
./phylogeny/runPhylogeny.sh
```

### Tree divergence dating
[Divergence estimates with mcmcTree](./TreeDating/treeDating.workflow.sh)
```bash
./TreeDating/treeDating.workflow.sh
```

## Gene Family Evolution

The [complete workflow](./geneFamilies/README.md) for analysing gene family evolution with CAFE is described here:

```bash
./geneFamilies/README.md

# This readme refers to the following scripts:
  run.geneFamilies.mcl.sh
  run.CAFError.sh
  run.CAFE.sh
  analyseCAFEmodels.Rmd
  CAFE.lhtest.R
  analyseCAFEoutput.Rmd
```

## Demographic inference with MSMC2
Demographic history for A. insinuator and A. echinatior with MSMC2:

./Ne/MSMC2.refAech.md
./Ne/plotMSMC2.Rmd

## Evolutionary rate analyses
[Evolutionary rate analyses](./SOS/AnalyseDNDS.hyphy.sh) and signatures of selection (SOS) on single copy orthologs with HyPhy (v2.3.14)

```
./SOS/AnalyseDNDS.hyphy.sh
./SOS/RELAX.summary.Rmd
./SOS/absREL.summary.Rmd
./SOS/analyseSOS.Rmd
```
## Gene Ontology (GO) enrichment
GO enrichment analyses were run on genes under [relaxed, intensified](./GOenrichment/RELAX.GOenrichment.Rmd), and [positive](./GOenrichment/absREL.GOenrichment.Rmd) selection and on shrinking gene families.
```bash
./GOenrichment/RELAX.GOenrichment.Rmd
./GOenrichment/absREL.GOenrichment.Rmd
./GOenrichment/CAFE.GOenrichment.Rmd
# scripts for revigo plots
./GOenrichment/revigoPlots/
```

## Whole genome alignment
Whole genome alignment with progressiveCactus v0.1
```bash
./WGA/run.progressiveCactus.sh
```

## Synteny Analysis
Synteny was analyzed with i-adhore v3.0.1
```bash
./synteny/run.iadhore.sh
./synteny/analyse.synteny.Rmd

# analysis of annotated OR genes only
./synteny/run.iadhore.ORs.sh
./synteny/plotOR.synteny.Rmd

```

## Analysis of the attine OR gene family

#### OR annotation

The full annotation pipeline is contained in

```bash
./ORevolution/clean/ORannotation.190408.sh
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
./ORphylogeny/ORphylogeny.sh
./ORphylogeny/ORplotTreeClean.Rmd
./ORphylogeny/analyseTree.Rmd
```

#### OR gene tree reconciliation

Scripts for clade-wide OR alignments, alignment QC, tree inference, tree plotting, and analyses.

```bash
./ORevolution/reconciliation/allClades.Rmd
./ORevolution/reconciliation/DLCpar.md
./ORevolution/reconciliation/plotReconTreeDLCpar.Rmd
./ORevolution/reconciliation/runAllDLCpar.sh
```


#### Evolutionary rate analyses in ORs with HyPhy (v2.3.14)

```bash
./ORevolution/OR.SOS/hyphy.ORs.sh
./ORevolution/OR.SOS/summary.RELAX.OR.Rmd
./ORevolution/OR.SOS/summary.absREL.OR.Rmd

# Detailed analysis of clade 285
./ORevolution/OR.SOS/analyseClade285.Rmd
```

## Annnotation & analysis of other gene families

```bash
# MRJPs
./genePhylogenies/MRJPevolution/MRJP.annotation.sh
./genePhylogenies/MRJPevolution/MRJP.phylogeny.sh
./genePhylogenies/MRJPevolution/plotMRJP.tree.Rmd
./genePhylogenies/MRJPevolution/mrjp.analyseTree.Rmd

# Gustatory Receptors (GR)
./genePhylogenies/GRevolution/annotateGRs.sh
./genePhylogenies/GRevolution/phylogenyGRs.sh
./genePhylogenies/GRevolution/plotGRtree.Rmd

# Elongases
./genePhylogenies/ELONGASEevolution/annotate.ELNG.sh
./genePhylogenies/ELONGASEevolution/phylogeny.ELNG.sh
./genePhylogenies/ELONGASEevolution/plotELGtree.Rmd

# Cuticular Proteins (CPRs)
./genePhylogenies/CPRevolution/annotateCPR.sh
./genePhylogenies/CPRevolution/phylogenyCPR.sh
./genePhylogenies/CPRevolution/plotCPRtree.Rmd
```

#### Script to plot barplots for all analyzed gene families

```bash
./genePhylogenies/barplot.Rmd
```
