# MSMC2 analysis to infer effective population size of A. echinatior and A. insinuator

The Pairwise Sequentially Markovian Coalescent (PSMC) model uses information in the complete diploid sequence of a single individual to infer the history of population size changes. PSMC uses the coalescent approach to estimate changes in population size. Each diploid genome is a collection of hundreds of thousands independent loci. Estimating TMRCA of the two alleles at each locus is used to create a TMRCA distribution across the genome. And since the rate of coalescent events is inversely proportional to effective population size (Ne), PSMC identifies periods of Ne change. For example, when many loci coalesce at the same time, it is a sign of small Ne at that time.

Based on the comparisons, they set the following (rather strict) cut off values for running a reliable PSMC analysis:

- A mean genome-wide coverage of at least 18X
- No more than 25% of missing data
- A per-site filter of ≥10 reads

The **raw read files** can be found at:
```/global/homes/jg/schradel/data/AcromyrmexPopGen/00.raw_data_pop/```

The **reference genomes** can be found at
```/global/homes/jg/schradel/data/genomes/reannotations_of_7_ants/```


# samples
| | | ||
|-|-|-|-|
|GAi-073|wHAXPI040443-85|Ae503|**new**
|GAi−020|wHAXPI040635-77|Ae374|
|GAe-067|wHAXPI040607-13|Ae459|
|GAe-077|wHAXPI040608-14|Ae077|
see   ```/global/homes/jg/schradel/data/inqGen/MSMC/samples.tsv```

## Prepare environment

```bash
base=/global/homes/jg/schradel/data/inqGen/MSMC
data=/global/homes/jg/schradel/data/AcromyrmexPopGen/00.raw_data_pop/
software=/global/homes/jg/schradel/software
trimmomatic=$software/Trimmomatic-0.39/trimmomatic-0.39.jar
bwa=$software/bwa-mem2-2.0pre1_x64-linux/

cd $base
```

## Read trimming
We use Trimmomatic 0.39 for read trimming.

```bash
mkdir $base/01.trimming
cd $base/01.trimming

# prepare output folder
trimOut=$base/01.trimming/trimmed.filtered
mkdir $trimOut

# select appropriate adapter file
## "TruSeq3 (as used by HiSeq and MiSeq machines)"
adapterFile="/global/homes/jg/schradel/software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"
# define trimming
trimSettings="LEADING:10 TRAILING:10 SLIDINGWINDOW:4:10 MINLEN:30 ILLUMINACLIP:$adapterFile:2:30:10"

# (dry) run in parallel across all cutadapt filtered files.
# not tested whether the  "> $base/trimLogs.log" works
cut -f 2 $base/samples.tsv > $base/samples.lst
readlink -f $data/*/*/*/*_1.fq.gz|grep -f $base/samples.lst|grep "wHAXPI040443-85"|perl -pe "s/_1.fq.gz//g" |nice parallel --dryrun -j 1 java -jar $trimmomatic PE {}_1.fq.gz {}_2.fq.gz $trimOut/trimmed.{/.}.1.fq.gz $trimOut/trimmed.{/.}.1.unpaired.fq.gz $trimOut/trimmed.{/.}.2.fq.gz $trimOut/trimmed.{/.}.2.unpaired.fq.gz -threads 5 $trimSettings
readlink -f $data/*/*/*/*_1.fq.gz|grep -f $base/samples.lst|grep "wHAXPI040443-85"|perl -pe "s/_1.fq.gz//g" |nice parallel -j 1 java -jar $trimmomatic PE {}_1.fq.gz {}_2.fq.gz $trimOut/trimmed.{/.}.1.fq.gz $trimOut/trimmed.{/.}.1.unpaired.fq.gz $trimOut/trimmed.{/.}.2.fq.gz $trimOut/trimmed.{/.}.2.unpaired.fq.gz -threads 5 $trimSettings
```

## Read alignment

```bash
mkdir $base/02.mapping
mkdir $base/02.mapping/mapping
cd $base/02.mapping/

# get reference 1
ln -s /global/homes/jg/schradel/data/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/genome/Acromyrmex_echinatior.v2.0.fa $base/Aech.v2.0.fa

# index reference
#$bwa/bwa-mem2 index $base/Aech.v2.0.fa
bwa index $base/Aech.v2.0.fa

for file in $(readlink -f $base/01.trimming/trimmed.filtered/*.1.fq.gz)
do
 fw=$file
 rev=$(echo $file|perl -pe "s/1.fq.gz/2.fq.gz/g")
 sample=$(echo $file|perl -pe 's/.*\_(.+?)\.[0-9].fq.gz/$1/g')
 library=$(echo $file|perl -pe 's/.*\/trimmed\.(.*)_.*.fq.gz/$1/g')
 sampleGroup=$(grep $sample $base/samples.tsv |cut -c 1-3)
 readgroup=$(echo "@RG\\tID:$sample\\tLB:$library\\tPL:ILLUMINA\\tCN:BGI-Shenzhen\\tKS:NNNNNNNN\\tSM:$sample")
 reference=$base/Aech.v2.0.fa
echo "echo \"working on $sample:\""
echo "#reads are taken from $file."
echo "#using the following settings:"
echo -e "#library=$library\tsampleGroup=$sampleGroup\tsample=$sample"
echo "#bwa: align,samtools view: extract mapped reads, samtools sort: sort"
echo "nice bwa mem -t 20 -R '$readgroup' $reference $fw $rev |nice samtools view -@ 14 -b -F 4 - |nice samtools sort -@ 14 -o $base/02.mapping/mapping/$sample.bam"

echo ""
done > $base/02.mapping/bwa.alignment.sh

# check commands
head $base/02.mapping/bwa.alignment.sh
# run commands
bash $base/02.mapping/bwa.alignment.sh
```

```bash
#average read depth
ls $base/02.mapping/mapping/*.bam| nice parallel --dryrun -j 1 "samtools depth -a  {} |  awk '{sum+=\$3} END { print \"{}\",sum/NR}'"
ls $base/02.mapping/mapping/*.bam| nice parallel -j 1 "samtools depth -a  {} |  awk '{sum+=\$3} END { print \"{}\",sum/NR}'" > coverage.tsv &

# samtools view -H mapping/wHAXPI040635-77.bam |sed "s/KS:/KS:NNNNNNNN/g"| samtools reheader - mapping/wHAXPI040635-77.bam > mapping/wHAXPI040635-77.reheader.bam
# set -d to a third and -D to twice
# i.e. -d = 9 -D = 54 if coverage = 27
# i.e. -d = 8 -D = 46 if coverage = 23
```

# Remove duplicates
```bash
mkdir $base/03.duplicates
cd $base/03.duplicates
ls $base/02.mapping/mapping/*.bam|egrep "wHAXPI040607-13" | nice parallel -j 4 java -jar /global/homes/jg/schradel/software/picard.jar MarkDuplicates \
      I={} \
      O={/.}.marked_duplicates.bam \
      M={/.}.marked_dup_metrics.txt \
      REMOVE_DUPLICATES=T \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

ls $base/03.duplicates/*.bam|parallel -j 4 samtools index {}

```

# Batch script to run bamCaller on each scaffold
```bash
#!/bin/bash
#$ -S /bin/bash # tell the SGE to use bash
#$ -cwd         #run in current working dir
#$ -pe smp 1    # how many CPUs?
#$ -l h_vmem=4G # how much RAM?
#$ -w e         # reject jobs with error
#$ -V           # export all environment variables
#$ -N bamCaller     # job name
#$ -e ./bamCaller.err
#$ -o ./bamCaller.out

# Usage
#qsub -v file="Aech.v2.0.over10000.lst" -v reference="Ains.v2.1.fa" -v sample="sample.bam" -t 1-$(wc -l < listOfScaffolds.lst) batch.bamCaller.qsub.sh

processing=$(sed "${SGE_TASK_ID}q;d" ${file})
sampleName=$(echo ${sample}|perl -pe 's/.*\/(.*?)\..*/$1/g')

samtools mpileup -r ${processing} -q 20 -Q 20 -C 50 -u -f ${reference} ${sample} | \
bcftools call -c -V indels | \
/global/homes/jg/schradel/software/msmc2/msmc-tools-master/bamCaller.py 11 ${sampleName}.${processing}_mask.bed.gz| gzip -c > ${sampleName}.${processing}.vcf.gz
```

# Batch submission for each sample
```bash
mkdir $base/04.MSMC
cd $base
samtools faidx Aech.v2.0.fa
cat Aech.v2.0.fa.fai |sort -k2,2 -nr|awk -v FS='\t'  '{if ($2 > 100000) print $0}'|cut -f 1 > Aech.v2.0.over100000.lst

cd $base/04.MSMC

#|GAe-077|wHAXPI040608-14|Ae077|
qsub -N GAe077.bamCaller -v file="../Aech.v2.0.over100000.lst" -v reference="/global/homes/jg/schradel/data/inqGen/MSMC/Aech.v2.0.fa" -v sample="../03.duplicates/wHAXPI040608-14.marked_duplicates.bam" -t 1-$(wc -l < ../Aech.v2.0.over100000.lst) ../batch.bamCaller.qsub.sh

#|GAe-067|wHAXPI040607-13|Ae459|
qsub -N GAe067.bamCaller -v file="../Aech.v2.0.over100000.lst" -v reference="/global/homes/jg/schradel/data/inqGen/PSMC/Aech.v2.0.fa" -v sample="../03.duplicates/wHAXPI040607-13.marked_duplicates.bam" -t 1-$(wc -l < ../Aech.v2.0.over100000.lst) ../batch.bamCaller.qsub.sh

#|GAi-073|wHAXPI040443-85|Ae503|
qsub -N GAi073.bamCaller -v file="../Aech.v2.0.over100000.lst" -v reference="/global/homes/jg/schradel/data/inqGen/MSMC/Aech.v2.0.fa" -v sample="../03.duplicates/wHAXPI040443-85.marked_duplicates.bam" -t 1-$(wc -l < ../Aech.v2.0.over100000.lst) ../batch.bamCaller.qsub.sh

#|GAi−020|wHAXPI040635-77|Ae374|
qsub -N GAi020.bamCaller -v file="../Aech.v2.0.over100000.lst" -v reference="/global/homes/jg/schradel/data/inqGen/MSMC/Aech.v2.0.fa" -v sample="../03.duplicates/wHAXPI040635-77.marked_duplicates.bam" -t 1-$(wc -l < ../Aech.v2.0.over100000.lst) ../batch.bamCaller.qsub.sh

mkdir vcf
mkdir maskBed

mv *.vcf.gz vcf
mv *.bed.gz maskBed
```

#Combine samples for each scaffold
```bash
mkdir $base/04.MSMC/multihetstep
cd $base/04.MSMC
# ORDER: AECH1 AECH2 AINS1 AINS2
cat ../Aech.v2.0.over100000.lst|parallel -j 16 "generate_multihetsep.py --chr {} --mask maskBed/wHAXPI040607-13.{}_mask.bed.gz --mask maskBed/wHAXPI040608-14.{}_mask.bed.gz --mask maskBed/wHAXPI040635-77.{}_mask.bed.gz --mask maskBed/wHAXPI040443-85.{}_mask.bed.gz vcf/wHAXPI040607-13.{}.vcf.gz vcf/wHAXPI040608-14.{}.vcf.gz vcf/wHAXPI040635-77.{}.vcf.gz vcf/wHAXPI040443-85.{}.vcf.gz> multihetstep/multihetsep.{}.txt"

# TEST
generate_multihetsep.py --chr scaffold7 --mask maskBed/wHAXPI040607-13.scaffold7_mask.bed.gz --mask maskBed/wHAXPI040608-14.scaffold7_mask.bed.gz --mask maskBed/wHAXPI040635-77.scaffold7_mask.bed.gz --mask maskBed/wHAXPI040443-85.scaffold7_mask.bed.gz vcf/wHAXPI040607-13.scaffold7.vcf.gz vcf/wHAXPI040608-14.scaffold7.vcf.gz vcf/wHAXPI040635-77.scaffold7.vcf.gz vcf/wHAXPI040443-85.scaffold7.vcf.gz> multihetstep/multihetsep.scaffold7.txt

```

# Run MSMC2 on all scaffolds
```
# CAREFUL WITH THIS OPTION!
-I, --pairIndices:                  this can be given in two flavors. First, you can enter a
                                    single comma-separated list like this "-I 0,1,4,5".
                                    In this case, the program will
                                    run over all pairs of haplotypes within this set of
                                    indices. This is useful for running on multiple phased
                                    diploid genomes sampled from one population.
                                    In the second flavor, you can give a list of pairs, like
                                    this: "-I 0-1,2-3,4-5". In this case, the
                                    program will run only those specified pairs. This can be
                                    used to run on a number of unphased genomes, to avoid pairs
                                    of haplotypes from different individuals. This should also
                                    be used to indicate a cross-population run, where you want
                                    to run the program only over pairs of haplotypes across
                                    population boundaries. So with two phased genomes, one from
                                    each population you'd run "-I 0-2,0-3,1-2,1-3", and the
                                    program would run only those four pairs of haplotypes.
                                    Note that if you do not use this parameter altogether,
                                    MSMC2 will run on all pairs of input haplotypes.
```

## Run MSMC2

see http://biorxiv.org/cgi/reprint/382036v1
MSMC is robust for scaffolds > 1MB, here I restrict the analysis to scaffolds over 2MB.

```bash
cd $base/04.MSMC
# delete multihetstep files without any SNPs
find $base/04.MSMC/multihetstep/multihetsep.*.txt -empty -delete

# process scaffolds in chunks of 10
hets=$(ls $base/04.MSMC/multihetstep/multihetsep.*.txt|tr "\n" " ")

# restrict to largest Scaffolds
cat $base/Aech.v2.0.fa.fai |sort -k2,2 -nr|less
hets=$(cat ../Aech.v2.0.over100000.lst |sed 's|^|/global/homes/jg/schradel/data/inqGen/MSMC/04.MSMC/multihetstep/multihetsep.|g'|sed 's|$|.txt|g' | sed -n '1,30 p'|tr "\n" " ")

# 1-4 k SNPs

# FOR AECH
mkdir $base/04.MSMC/results
# Fix mutation rate across runs? tested. Not giving sensible results.

# Default P 1*2+15*1+1*2
msmc2 -t 20 -o $base/04.MSMC/results/Aech.msmc2 -p 1*2+15*1+1*2 -I 0-1,2-3 $hets
# FOR AINS
msmc2 -t 20 -o $base/04.MSMC/results/Ains.msmc2 -p 1*2+15*1+1*2 -I 4-5,6-7 $hets

# FOR BOTH
# results for cross-coalesce likely not working because genomes are not phased
msmc2 -t 20 -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -s -p 1*2+15*1+1*2 -o $base/04.MSMC/results/Ains.Aech.msmc2 $hets

# COMBINE
combineCrossCoal.py $base/04.MSMC/results/Ains.Aech.msmc2.final.txt $base/04.MSMC/results/Ains.msmc2.final.txt \
    $base/04.MSMC/results/Aech.msmc2.final.txt > $base/04.MSMC/results/combined.msmc2.final.txt
```
# Bootstrap
```bash
mkdir $base/04.MSMC/bootstrap
mkdir $base/04.MSMC/bootstrap/results
cd $base/04.MSMC/bootstrap

count=$(echo $hets|tr " " "\n"|wc -l)

# options set to resemble used real-data scaffolds
multihetsep_bootstrap.py -s 100000 --chunks_per_chromosome 20 -n 500 --nr_chromosomes $count ./ $hets
cd $base/04.MSMC/bootstrap/
qsub -t 1-500 batch.boots.qsub.sh

#Two bootstraps did not converge: 203, 455
find $base/04.MSMC/bootstrap/results/ -empty -delete


SGE_TASK_ID=302
  boots=$(ls ./_${SGE_TASK_ID}/*multihetsep.*.txt|tr "\n" " ")
  msmc2 -t 2 -o ./results/Aech.boot${SGE_TASK_ID}.msmc2 -p 1*2+15*1+1*2 -I 0-1,2-3 ${boots}
  msmc2 -t 2 -o ./results/Ains.boot${SGE_TASK_ID}.msmc2 -p 1*2+15*1+1*2 -I 4-5,6-7 ${boots}
  msmc2 -t 2 -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -s -p 1*2+15*1+1*2 -o ./results/Ains.Aech.boot${SGE_TASK_ID}.msmc2 ${boots}
  combineCrossCoal.py ./results/Ains.Aech.boot${SGE_TASK_ID}.msmc2.final.txt ./results/Ains.boot${SGE_TASK_ID}.msmc2.final.txt \
        ./results/Aech.boot${SGE_TASK_ID}.msmc2.final.txt > ./results/combined.boot${SGE_TASK_ID}.msmc2.final.txt


```
# Batch script to run all bootstraps
```bash
#!/bin/bash
#$ -S /bin/bash # tell the SGE to use bash
#$ -cwd         #run in current working dir
#$ -pe smp 2    # how many CPUs?
#$ -l h_vmem=2G # how much RAM?
#$ -w e         # reject jobs with error
#$ -V           # export all environment variables
#$ -N boots     # job name
#$ -e ./bootstrap.err
#$ -o ./bootstrap.out

# Usage
#qsub -v -t 1-100 batch.boots.qsub.sh

boots=$(ls ./_${SGE_TASK_ID}/*multihetsep.*.txt|tr "\n" " ")
msmc2 -t 2 -o ./results/Aech.boot${SGE_TASK_ID}.msmc2 -p 1*2+15*1+1*2 -I 0-1,2-3 ${boots}
msmc2 -t 2 -o ./results/Ains.boot${SGE_TASK_ID}.msmc2 -p 1*2+15*1+1*2 -I 4-5,6-7 ${boots}
msmc2 -t 2 -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -s -p 1*2+15*1+1*2 -o ./results/Ains.Aech.boot${SGE_TASK_ID}.msmc2 ${boots}
combineCrossCoal.py ./results/Ains.Aech.boot${SGE_TASK_ID}.msmc2.final.txt ./results/Ains.boot${SGE_TASK_ID}.msmc2.final.txt \
    ./results/Aech.boot${SGE_TASK_ID}.msmc2.final.txt > ./results/combined.boot${SGE_TASK_ID}.msmc2.final.txt


```
