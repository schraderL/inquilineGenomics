
# Gene family size Evolution analysis of inquiline genomes
 We will use the following workflow to repeat the analysis:

 1. GeMoMa
 2. AGAT
 3. Orthofinder
 4. TransposonPSI?
 5. CAFE5
-----------------------------------

 <!-- TOC -->
 - [Software required](#software-required)
 - [Prepare data](#prepare-data)
 - [Homology-based gene prediction with `GeMoMa`](#homology-based-gene-prediction-with-gemoma)
   - [Process GeMoMa predictions](#process-gemoma-predictions)
     - [Combine predictions](#combine-predictions)
     - [Clean GFFs and retrieve fasta sequences](#clean-gffs-and-retrieve-fasta-sequences)
 - [QC of genome assemblies and annotations](#qc-of-genome-assemblies-and-annotations)
   - [BUSCO](#busco)
 - [TransposonPsi](#transposonpsi)
 <!-- /TOC -->


# Software required

* gemoma: www.jstacs.de/download.php?which=GeMoMa
* mmseqs: mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
* orthofinder: github.com/davidemms/OrthoFinder/releases/download/2.4.0/OrthoFinder.tar.gz
? Maybe rather use mmseqs easy-cluster

* CAFE5: github.com/hahnlab/CAFExp/releases/download/v5.0b2/CAFE5-5.0.0.tar.gz
* AGAT: github.com/NBISweden/AGAT

These software are installe here:
```
/corefac/cse/lukas/software/GeMoMa-1.7.1/GeMoMa-1.7.1.jar
/corefac/cse/lukas/software/CAFE5/bin/cafexp
/corefac/cse/lukas/software/OrthoFinder-2.4.0/orthofinder
/corefac/cse/lukas/software/mmseqs-sse41/bin/mmseqs
```

# Prepare data
```bash
base=/corefac/cse/lukas/inqGen18/reannotation
mkdir $base/data


# download all RefSeq myrmmicinae from NCBI
unzip /corefac/cse/lukas/inqGen18/reannotation/data/Myrmicinae.zip
## delete rna fasta as it is not required
rm $base/data/ncbi_dataset/data/*/rna.fna
mkdir $base/data/ncbi_dataset/data/unused
#GCF_001594045.1 Atta colombica
#GCF_000204515.1 Acromyrmex echinatior
#GCF_000143395.1 Atta cephalotes
mv GCF_001594045.1 unused
mv GCF_000204515.1 unused
mv GCF_000143395.1 unused

# retrieve genomes to annotate
cd $base/data
ln -s /corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_insinuator.2.1/genome/Acromyrmex_insinuator.v2.1.fa .
ln -s /corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_charruanus.2.1/genome/Acromyrmex_charruanus.v2.1.fa .
ln -s /corefac/cse/lukas/inqGen18/inquilines_v2.1/Pseudoatta_argentina.2.1/genome/Pseudoatta_argentina.v2.1.fa .
ln -s /corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_heyeri.2.1/genome/Acromyrmex_heyeri.v2.1.fa .
ln -s /corefac/cse/lukas/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/genome/Acromyrmex_echinatior.v2.0.fa .
ln -s /corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_colombica/genome/Atta_colombica.v2.0.fa .


# cufflinks assembled transcriptomes
cp /Users/lukas/CSE/inquiline_genomes/P.argentina/Pseudoatta_argentina.transcripts.gtf .
cp /Users/lukas/CSE/inquiline_genomes/A.charruanus/Acromyrmex_charruana.transcripts.gtf .
cp /Users/lukas/CSE/inquiline_genomes/A.heyeri/Acromyrmex_heyeri.transcripts.gtf .
cp /Users/lukas/CSE/inquiline_genomes/A.insinuator/Acromyrmex_insinuator.transcripts.gtf .

```

# Homology-based gene prediction with `GeMoMa`
We used `GeMoMa-v1.7.1` to reannotate the six genomes used in the gene family size evolution analysis. We retrieved genome annotations from all Myrmicine ants in RefSeq (except leaf-cutting ant species).


Used accessions for annotating the genomes:
>GCF_000187915.1 Pogonomyrmex barbatus
GCF_000188075.2 Solenopsis invicta
GCF_000949405.1 Vollenhovia emeryi
GCF_000956235.1 Wasmannia auropunctata
GCF_001594055.1 Trachymyrmex zeteki
GCF_001594065.1 Cyphomyrmex costatus
GCF_001594075.1 Trachymyrmex cornetzi
GCF_001594115.1 Trachymyrmex septentrionalis
GCF_003070985.1 Temnothorax curvispinosus
GCF_003260585.2 Monomorium pharaonis
```bash
gemoma=/corefac/cse/lukas/software/GeMoMa-1.7.1/GeMoMa-1.7.1.jar

# loop over all genomes to annotate
cd $base
for genome in $base/data/*.fa
do
  target=$(echo $genome|perl -pe 's/.*\/(.*)\.fa/$1/g')
  # loop over all reference genomes to annotate the genomes with
  for accession in $base/data/ncbi_dataset/data/G*/
  do
  ref=$(grep "species_name:" $accession/data_report.yaml |perl -pe 's/.*species_name: (.).* (...).*/$1$2/g')
  reffa=$(readlink -f $accession/*.fna)
  refgff=$(readlink -f $accession/genomic.gff)
  #run GeMoMa
  echo "nice java -Xmx40G -jar $gemoma CLI GeMoMaPipeline threads=25 outdir=$base/results/$target/$ref/ GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$genome i=$ref a=$refgff g=$reffa GeMoMa.pa=false"
  done
done > runGemoma.sh

# run all commands
bash runGemoma.sh > runGemoma.out 2> runGemoma.err

```

## Process GeMoMa predictions
### Combine predictions
```bash
gemoma=/corefac/cse/lukas/software/GeMoMa-1.7.1/GeMoMa-1.7.1.jar
cd $base

for in in /corefac/cse/lukas/inqGen18/reannotation/results/*
do
  allrefs=$(ls $in/*/final_annotation.gff |xargs -I -- echo "g="--" ")
  echo $allrefs
  # combine all predictions
  java -Xmx25G -jar $gemoma CLI GAF outdir=$in/GAF/ $allrefs
done
```

### Clean GFFs and retrieve fasta sequences
We retrieved the longest isoform for each gene.

```bash
# convert to proper GFF

cd $base
for in in /corefac/cse/lukas/inqGen18/reannotation/results/*
do
  genomeFa=$(echo $in|perl -pe 's/.*\/(.*)/$1.fa/g')
  short=$(echo $in|perl -pe 's/.*\/(.).*_(...).*/$1$2/g')
  SHORT=$(echo $short|tr a-z A-Z)
  echo -e $genomeFa"\t"$short"\t"$in
  echo "cd $in/GAF/"
  # clean up GFF
  echo "agat_sp_manage_IDs.pl -gff filtered_predictions.gff -o filtered_predictions_renamed.gff --tair --type_dependent --ensembl --prefix $SHORT"
  echo "agat_convert_sp_gxf2gxf.pl -g filtered_predictions_renamed.gff  -o $short.GeMoMa.gff3"

  # retrieve Fasta sequences of CDS (aa and nt)
  #agat_sp_extract_sequences.pl -g $short.GeMoMa.gff3 -f $base/data/$genomeFa  --aa -t CDS -o $short.GeMoMa.pep.fa
  #agat_sp_extract_sequences.pl -g $short.GeMoMa.gff3 -f $base/data/$genomeFa  -t CDS -o $short.GeMoMa.cds.fa

  # retrieve longest isoform gff
  echo "agat_sp_keep_longest_isoform.pl -gff $short.GeMoMa.gff3 -o $short.GeMoMa.longestIsoform.gff3"
  # retrieve longest isoform Fasta sequences of CDS (aa and nt)
  echo "seqkit seq $base/data/$genomeFa -w 100 > $base/data/$genomeFa.wrap"
  echo "agat_sp_extract_sequences.pl -g $short.GeMoMa.longestIsoform.gff3 -f $base/data/$genomeFa.wrap  --aa -t CDS -o $short.GeMoMa.longestIsoform.pep.fa"
  echo "agat_sp_extract_sequences.pl -g $short.GeMoMa.longestIsoform.gff3 -f $base/data/$genomeFa.wrap  -t CDS -o $short.GeMoMa.longestIsoform.cds.fa"
done > agatRun.sh

# check
grep ">" *.GeMoMa.longestIsoform.pep.fa -c
#14550
```


# QC of genome assemblies and annotations

 We ran BUSCO on the genome assembly and the set of predicted proteins for each of the 4 genomes.


## BUSCO

```bash

# prepare BUSCO env
export AUGUSTUS_CONFIG_PATH="/corefac/cse/lukas/software/augustus-3.2.3/config/"

# prepare environment
export base="/corefac/cse/lukas/inqGen18/reannotation/"
export BUSCO_CONFIG_FILE=$base/protein.config.ini

cd $base
cp /usr/local/home/lschrader/data/inqGen18/QC/BUSCO_inquilinesv2.1/protein.config.ini .
#######################
# run for each species
#######################

# protein mode
for species in /corefac/cse/lukas/inqGen18/reannotation/results/*/GAF/*.GeMoMa.longestIsoform.pep.fa
do
  out=$(echo $species|rev|cut -f 1 -d "/" |rev)
  echo "nice  python ~/software/buscoV3/scripts/run_BUSCO.py --in $species --out $out.BUSCO.prot --lineage /usr/local/home/lschrader/data/inqGen18/QC/BUSCO_inquilinesv2.1/datasets/hymenoptera_odb9 --mode prot --cpu 20 -f"
done > BUSCO.run.sh

# outfolder is
/usr/local/home/lschrader/data/inqGen18/QC/BUSCO/out.protein/
mkdir BUSCO
mv /usr/local/home/lschrader/data/inqGen18/QC/BUSCO/out.protein/*GeMoMa.longestIsoform.pep.fa.BUSCO.prot ./BUSCO/

```

>**ACHA**: C:97.8%[S:97.3%,D:0.5%],F:1.0%,M:1.2%,n:4415
>**AECH**: C:97.5%[S:97.0%,D:0.5%],F:1.2%,M:1.3%,n:4415
>**AHEY**: C:95.8%[S:95.3%,D:0.5%],F:2.4%,M:1.8%,n:4415
>**PARG**: C:97.2%[S:96.8%,D:0.4%],F:1.4%,M:1.4%,n:4415
>**AINS**: C:97.8%[S:97.3%,D:0.5%],F:1.0%,M:1.2%,n:4415
>**ACOL**: C:97.5%[S:97.3%,D:0.2%],F:1.2%,M:1.3%,n:4415


# TransposonPsi
```bash
base=/corefac/cse/lukas/inqGen18/reannotation

mkdir $base/Tpsi
cd $base/Tpsi
ls /corefac/cse/lukas/inqGen18/reannotation/results/*/GAF/*.GeMoMa.longestIsoform.pep.fa|parallel --nice 10 "perl ~/software/TransposonPSI_08222010/transposonPSI.pl {} prot 1>{}.PSI.out 2> {}.PSI.err"

```
<!--

```bash
https://github.com/hahnlab/CAFExp
/corefac/cse/lukas/software/GeMoMa-1.7.1/GeMoMa-1.7.1.jar
/corefac/cse/lukas/software/CAFE5/bin/cafexp
/corefac/cse/lukas/software/OrthoFinder-2.4.0/orthofinder


```
# TransposonPSI?

## Orthofinder
We ran orthofinder (v.2.4.0) to retrieve orthogroups across all available myrmicines.
```bash

mkdir $base/orthofinder/input
cd $base/orthofinder/input
ln -s $base/results/Acromyrmex_charruanus.v2.1/GAF/Acha.GeMoMa.longestIsoform.pep.fa .
ln -s $base/results/Acromyrmex_charruanus.v2.1/GAF/Acha.GeMoMa.longestIsoform.pep.fa test.fa

cd $base/orthofinder/
nice /corefac/cse/lukas/software/OrthoFinder-2.4.0/orthofinder -f $base/orthofinder/input/ -t 20

mv $base/input/OrthoFinder/ $base
```

## orthoMCL
### mmseqs
```bash
ln -s $base/results/Acromyrmex_charruanus.v2.1/GAF/Acha.GeMoMa.longestIsoform.pep.fa .
ln -s $base/results/Acromyrmex_charruanus.v2.1/GAF/Acha.GeMoMa.longestIsoform.pep.fa test.fa
/usr/local/home/lschrader/data/inqGen18/
```

# CAFEv5
## Prepare orthogroups
```bash
cd $base/tmp/
sed '1 s/\.longestIsoform//g' Orthogroups.GeneCount.tsv|awk -v OFS='\t'  'NF{NF-=1};1'|cut -f 1-7|perl -pe 's/^/\t/g'|sed '1 s/^\t/Desc\t/g' > Orthogroups.CAFE.tsv

# rename for test run
nano Orthogroups.CAFE.tsv

head -n 1 Orthogroups.CAFE.tsv > tmp.tsv
sed -n '100,8000p' Orthogroups.CAFE.tsv >> tmp.tsv
```

## Prepare tree input
```bash
#treeFile=~/data/inqGen18/phylogeny/treeDating/MCMCtree/mcmctree.sf50000.sr17.5cal.sigma2.5.5.run4/FigTree.run4.tre
subtree=$(echo "(((AHEY:0.024983,(ACHA:0.016335,PARG:0.016335):0.008648):0.027567,(AECH:0.009605,AINS:0.009605):0.042945):0.065115,ACOL:0.117665);")
# multiply branch lengths by 100 so that 0.157526 becomes 0157.526 to avoid CAFE issues (-Inf probabilites)
tree=$(echo $subtree|perl -pe 's/([0-9])\.([0-9]{3})/$1$2./g')
# prepare lambda tree with single lambda rate for all branches
LAMBDATREE=$(echo $subtree|perl -pe 's/[0-9]\.[0-9]+/1/g')
echo $LAMBDATREE > lambda.tre
echo $tree > ultrametric.tre
cat lambda.tre |perl -pe 's/:[0-9]+\.*[0-9]*/:1/g'|perl -pe 's/PARG:1/PARG:2/g'|perl -pe 's/AECH:1/AECH:3/g' > lambda3.tre
```

```bash
# run cafe with error estimation
/corefac/cse/lukas/software/CAFE5/bin/cafexp -i tmp.tsv -t ultrametric.tre -y lambda.tre -e
# run cafe with one lambda
/corefac/cse/lukas/software/CAFE5/bin/cafexp -i tmp.tsv -t ultrametric.tre -y lambda.tre
# run cafe with three lambda
/corefac/cse/lukas/software/CAFE5/bin/cafexp -i tmp.tsv -t ultrametric.tre -y lambda3.tre


results/base_clade_results.txt
```

```bash
python /corefac/cse/lukas/software/CAFE5/python_scripts/cafetutorial_draw_tree.py -i reports/summary_run1_node.txt -t '((((cat:68.7105,horse:68.7105):4.56678,cow:73.2773):20.7227,(((((chimp:4.44417,human:4.44417):6.68268,orang:11.1268):2.28586,gibbon:13.4127):7.21153,(macaque:4.56724,baboon:4.56724):16.057):16.0607,marmoset:36.6849):57.3151):38.738,(rat:36.3024,mouse:36.3024):96.4356)' -d '((((cat<0>,horse<2>)<1>,cow<4>)<3>,(((((chimp<6>,human<8>)<7>,orang<10>)<9>,gibbon<12>)<11>,(macaque<14>,baboon<16>)<15>)<13>,marmoset<18>)<17>)<5>,(rat<20>,mouse<22>)<21>)<19>' -o reports/summary_run1_tree_rapid.png -y Rapid

```


#
# A fatal error has been detected by the Java Runtime Environment:
#
#  SIGSEGV (0xb) at pc=0x00007fa6ebe997c9, pid=60940, tid=0x00007f9bc16bf700
#
# JRE version: Java(TM) SE Runtime Environment (8.0_101-b13) (build 1.8.0_101-b13)
# Java VM: Java HotSpot(TM) 64-Bit Server VM (25.101-b13 mixed mode linux-amd64 )
# Problematic frame:
# V  [libjvm.so+0x8917c9][thread 140305214883584 also had an error][thread 140305202251520 also had an error]
