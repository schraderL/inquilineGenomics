

# find mitochondrial and fungal scaffolds
#############################################

base=/usr/local/home/lschrader/data/inqGen18/QC/filterContaminants

cd $base
mkdir $base/blsResults
mkdir $base/taxonLists
mkdir $ba$base/non_insect
mkdir $ba$base/mitochondrial
mkdir $ba$base/mitochondrial/blsResults
mkdir $base/v2.1/


#######################################
# blast genomes against NT database

AchaG="~/data/genomes/inquilines/v1.1/Acha/Acromyrmex_charruanus.fa"
AheyG="~/data/genomes/inquilines/v1.1/Ahey/Acromyrmex_heyeri.fa"
AinsG="~/data/genomes/inquilines/v1.1/Ains/Acromyrmex_insinuator.fa"
PargG="~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.fa"

blastn -query $AchaG -db $base/blastDB/nt -num_threads 20 -evalue 1e-10 -num_alignments 1 -outfmt 6 > $base/blsResults/Acha.blast.tsv
blastn -query $AheyG -db $base/blastDB/nt -num_threads 20 -evalue 1e-10 -num_alignments 1 -outfmt 6 > $base/blsResults/Ahey.blast.tsv
blastn -query $AinsG -db $base/blastDB/nt -num_threads 20 -evalue 1e-10 -num_alignments 1 -outfmt 6 > $base/blsResults/Ains.blast.tsv
blastn -query $PargG -db $base/blastDB/nt -num_threads 20 -evalue 1e-10 -num_alignments 1 -outfmt 6 > $base/blsResults/Parg.blast.tsv

#######################################


#######################################
# retrieve taxon information for each NT hit.
# get gi2tax.db
#wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz

bash  scripts/retrieveTaxonForBlast.sh $base/blsResults/Acha.blast.tsv Acha &
bash  scripts/retrieveTaxonForBlast.sh $base/blsResults/Ahey.blast.tsv Ahey &
bash  scripts/retrieveTaxonForBlast.sh $base/blsResults/Ains.blast.tsv Ains &
bash  scripts/retrieveTaxonForBlast.sh $base/blsResults/Parg.blast.tsv Parg &

mv *.blast.taxon.tsv taxonLists
mv *.taxon.list taxonLists
#######################################

#######################################
# identify non-bilaterian hits using blast taxon results

# starting genomes are inquilines_1.1
# ~/data/genomes/inquilines/v1.1/


cd $base/non_insect
for species in Acromyrmex_charruanus.v1.1 Acromyrmex_heyeri.v1.1 Acromyrmex_insinuator.v1.1 Pseudoatta_argentina.v1.1
do
speciesShort=$(echo $species|perl -pe 's/(.).*_(...).*/$1$2/g')
awk -F $'\t' '{if ($14 !~/Bilateria/) print $1,$14}' ../taxonLists/$speciesShort.blast.taxon.tsv|sort|uniq > $species.contaminants.lst
done

# check if all the contaminant scaffolds are not hitting insects in blast
for species in Acromyrmex_charruanus.v1.1 Acromyrmex_heyeri.v1.1 Acromyrmex_insinuator.v1.1 Pseudoatta_argentina.v1.1
do
  grep "Bilateria" $species.contaminants.recursive.lst|perl -pe 's/.*(scaffold[0-9]*|C[0-9]*)\s.*/$1/g' |grep -f - $species.contaminants.recursive.lst;
done

# Some scaffolds also hit Vollenhovia. I assume that the Vollenhovia still contains some bacterial scaffolds
#######################################

#######################################
# identify mitochondrial scaffolds

# starting genomes are inquilines_2.0
# ~/data/genomes/inquilines_2.0

cd $ba$base/mitochondrial

#get mitochondrial genome of Wasmannia auropunctata
referenceMito="KX146469.1.Waur.mitome.fa"
esearch -db nucleotide -query "KX146469.1"| efetch -format fasta > ../$referenceMito
genomeFolder=/usr/local/home/lschrader/data/genomes/inquilines_2.0

# blast mitochndrium against inquiline genomes
for blastDB in $genomeFolder/Acromyrmex_charruanus/genome/Acromyrmex_charruanus.v2.0.fa $genomeFolder/Acromyrmex_heyeri/genome/Acromyrmex_heyeri.v2.0.fa $genomeFolder/Acromyrmex_insinuator/genome/Acromyrmex_insinuator.v2.0.fa $genomeFolder/Pseudoatta_argentina/genome/Pseudoatta_argentina.v2.0.fa
do
  species=$(echo $blastDB|perl -pe 's/.*\///g')
  blastn -db $blastDB -query $referenceMito -outfmt 6 > blsResults/$species.mito.bls
  # run R script to compute coverage across entire scaffold. Print only those that have >50 % mitochondrial hits
  Rscript $base/scripts/processHits.R blsResults/$species.mito.bls $blastDB.fai $species.mito.coverage 0.5
done

#######################################
# Create clean genomes for v2.1

# 1. remove contaminant contigs/scaffolds from fasta and create fungal and mitochondrial subset fastas
# 2. filter gffs and produce clean gene annotation gff
# 3. produce clean CDS and pep fasta files from new gffs
# 4. filter functional annotations according to v2.1 genes

# input files are inquilines_2.0
genomeFolder=/usr/local/home/lschrader/data/genomes/inquilines_2.0

for speciesName in Pseudoatta_argentina Acromyrmex_insinuator Acromyrmex_charruanus Acromyrmex_heyeri
do
  cd $base/v2.1/

#######################################
  # define files in hash
  declare -A Data=(\
  ["Base"]="$genomeFolder/$speciesName" \
  ["abbr"]=$(echo $speciesName|perl -pe 's/(.).*?_(...).*/uc($1.$2)/ge;') \
  ["GenomeFa"]="${Data[Base]}/genome/$speciesName.v2.0.fa" \
  ["BaseAnnotation"]="${Data[Base]}/annotation/gene_annotation/" \
  ["BaseFuncAnnotation"]="${Data[Base]}/annotation/function_annotation/" \
  ["pepFa"]="${Data[BaseAnnotation]}$speciesName.v2.0.pep.fa" \
  ["cdsFa"]="${Data[BaseAnnotation]}$speciesName.v2.0.cds.fa" \
  ["utrGff"]="${Data[BaseAnnotation]}$speciesName.v2.0.UTR.gff" \
  ["geneGff"]="${Data[BaseAnnotation]}$speciesName.v2.0.gff" \
  ["swiss"]="${Data[BaseFuncAnnotation]}$speciesName.v2.0.pep.fa.blast.swissprot.best" \
  ["GO"]="${Data[BaseFuncAnnotation]}$speciesName.v2.0.pep.fa.iprscan.gene.GO" \
  ["ipr"]="${Data[BaseFuncAnnotation]}$speciesName.v2.0.pep.fa.iprscan.gene.ipr" \
  ["K"]="${Data[BaseFuncAnnotation]}$speciesName.v2.0.pep.fa.KEGG.K" \
  ["ko"]="${Data[BaseFuncAnnotation]}$speciesName.v2.0.pep.fa.KEGG.ko" \
  ["repGff"]="${Data[Base]}/annotation/repeat_annotation/all.repeat.gff" \
  )
#######################################

#######################################
  # check if all folders are present
  if [ ! -d ${Data["Base"]} ]; then echo "File ${Data["Base"]} not found!"; break 1; fi
  if [ ! -d ${Data["BaseAnnotation"]} ]; then echo "File ${Data["BaseAnnotation"]} not found!"; break 1; fi
  if [ ! -d ${Data["BaseFuncAnnotation"]} ]; then echo "File ${Data["BaseFuncAnnotation"]} not found!"; break 1; fi

  # check if all input files are present
  if [ ! -f ${Data["GenomeFa"]} ]; then echo "File ${Data["GenomeFa"]} not found!"; break 1; fi
  if [ ! -f ${Data["pepFa"]} ]; then echo "File ${Data["pepFa"]} not found!"; break 1; fi
  if [ ! -f ${Data["cdsFa"]} ]; then echo "File ${Data[$"cdsFa"]} not found!"; break 1; fi
  if [ ! -f ${Data["utrGff"]} ]; then echo "File ${Data["utrGff"]} not found!"; break 1; fi
  if [ ! -f ${Data["geneGff"]} ]; then echo "File ${Data["geneGff"]} not found!"; break 1; fi
  if [ ! -f ${Data["swiss"]} ]; then echo "File ${Data["swiss"]} not found!"; break 1; fi
  if [ ! -f ${Data["GO"]} ]; then echo "File ${Data["GO"]} not found!"; break 1; fi
  if [ ! -f ${Data["ipr"]} ]; then echo "File ${Data["ipr"]} not found!"; break 1; fi
  if [ ! -f ${Data["K"]} ]; then echo "File ${Data["K"]} not found!"; break 1; fi
  if [ ! -f ${Data["ko"]} ]; then echo "File ${Data["ko"]} not found!"; break 1; fi
  if [ ! -f ${Data["repGff"]} ]; then echo "File ${Data["repGff"]} not found!"; break 1; fi

  if [ ! -f $base/non_insect/$speciesName*.contaminants.lst ]; then echo "File $base/non_insect/$speciesName*.contaminants.lst not found!"; break 1; fi
  if [ ! -f $base/mitochondrial/$speciesName*.mito.coverage ]; then echo "File $base/mitochondrial/$speciesName*.mito.coverage not found!"; break 1; fi

  echo "All good for $speciesName"
#######################################
  # prepare directories
  cd $base/v2.1/
  mkdir $speciesName.2.1
  cd $base/v2.1/$speciesName.2.1
  mkdir genome
  mkdir annotation
  mkdir annotation/function_annotation
  mkdir annotation/gene_annotation
  mkdir annotation/repeat_annotation/
  mkdir lists
#######################################

#######################################
  # select non-insect scaffolds to remove
  perl -pe 's/.*(scaffold[0-9]*|C[0-9]*)\s.*/$1/g' $base/non_insect/$speciesName*.contaminants.lst|sort|uniq > $speciesName.filter.non-bilateria.lst
  grep ";Bacteria;" $base/non_insect/$speciesName*.contaminants.lst |perl -pe 's/.*(scaffold[0-9]*|C[0-9]*)\s.*/$1/g' |sort|uniq> $speciesName.filter.bacterial.lst
  grep ";Fungi;" $base/non_insect/$speciesName*.contaminants.lst |perl -pe 's/.*(scaffold[0-9]*|C[0-9]*)\s.*/$1/g' |sort|uniq> $speciesName.filter.fungal.lst

  # select mitochondrial scaffolds to remove
  perl -pe 's/.*(scaffold[0-9]*|C[0-9]*)\s.*/$1/g' $base/mitochondrial/$speciesName*.mito.coverage|egrep -v '^#' > $speciesName.filter.mitochondrial.lst

  # cat mitochondrial and non-insect entries
  cat $speciesName.filter.non-bilateria.lst $speciesName.filter.mitochondrial.lst > $speciesName.contaminants.lst
#######################################

#######################################
  # remove contaminants from genome fasta and create mitome, fungal and bacterial fastas
  faSomeRecords -exclude ${Data[GenomeFa]} $speciesName.contaminants.lst genome/$speciesName.v2.1.fa
  # index fasta with samtools as gffread indexing is messed up
  samtools faidx genome/$speciesName.v2.1.fa
  faSomeRecords ${Data[GenomeFa]} $speciesName.filter.mitochondrial.lst genome/$speciesName.mitome.v2.1.fa
  faSomeRecords ${Data[GenomeFa]} $speciesName.filter.bacterial.lst genome/$speciesName.bacterial.v2.1.fa
  faSomeRecords ${Data[GenomeFa]} $speciesName.filter.fungal.lst genome/$speciesName.fungal.v2.1.fa
#######################################

#######################################
  # remove contaminant genes from annotation gff and prepare a proper gff3
  cat $speciesName.contaminants.lst|perl -pe 's/^(.*)$/^$1\t/g'|egrep -v -f - ${Data["geneGff"]} > $speciesName.tmp1.gff
  Rscript $base/scripts/addGene2gff.R $speciesName.tmp1.gff $speciesName.tmp.gff ${Data["abbr"]}
  # check which genes are removed
  echo "The following genes are removed from $speciesName"
  echo "###################################################"
  cat $speciesName.contaminants.lst|perl -pe 's/^(.*)$/^$1\t/g'|egrep -f - ${Data["geneGff"]}
  echo "###################################################"
  GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= $speciesName.tmp.gff

  ~/software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -addintrons $speciesName.tmp_clean.gff > annotation/gene_annotation/$speciesName.v2.1.gff3
  cat annotation/gene_annotation/$speciesName.v2.1.gff3|awk '{if ($3=="mRNA")print $1":"$4"-"$5" "$9}'|perl -pe 's/(.*?ID=)(.*?)([;.*\n|\n])/$2 $1 $3/g'|perl -pe 's/(.*) ID=.*/$1/g' > annotation/gene_annotation/$speciesName.v2.1.pos
  awk '{print ">"$1"\t>"$0}' annotation/gene_annotation/$speciesName.v2.1.pos > $speciesName.replacementTable.tsv

#######################################

#######################################
  # clean UTR gff
  cat $speciesName.contaminants.lst|perl -pe 's/^(.*)$/^$1\t/g'|egrep -v -f - ${Data["utrGff"]} > $speciesName.tmputr.gff
  # check which genes are removed
  cat $speciesName.contaminants.lst|perl -pe 's/^(.*)$/^$1\t/g'|egrep -f - ${Data["utrGff"]}
  GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= $speciesName.tmputr.gff
  ~/software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -addintrons $speciesName.tmputr_clean.gff > $speciesName.tmputr_clean2.gff
  awk '{if ($3~"UTR") print $0}' $speciesName.tmputr_clean2.gff > annotation/gene_annotation/$speciesName.v2.1.utr.gff3

#######################################


#######################################
  #prepare cds fasta file from cleaned gff
  gffread -w tmp1.fasta -g  genome/$speciesName.v2.1.fa annotation/gene_annotation/$speciesName.v2.1.gff3
  cut -f 1 -d " " tmp1.fasta|awk -F "\t" ' FNR==NR { a[$1]=$2; next } $1 in a { $1=a[$1] }1' $speciesName.replacementTable.tsv - > annotation/gene_annotation/$speciesName.v2.1.cds.fa
  rm tmp1.fasta
#######################################

#######################################
  #prepare pep fasta file
  gffread -y tmp2.fasta -g  genome/$speciesName.v2.1.fa annotation/gene_annotation/$speciesName.v2.1.gff3
  cut -f 1 -d " " tmp2.fasta|awk -F "\t" ' FNR==NR { a[$1]=$2; next } $1 in a { $1=a[$1] }1' $speciesName.replacementTable.tsv - |perl -pe 's/\./\*/g' > annotation/gene_annotation/$speciesName.v2.1.pep.fa
  rm  tmp2.fasta
#######################################

#######################################
  # clean functional annotations
  # create list of all genes to keep
  cut -f 1 annotation/gene_annotation/$speciesName.v2.1.pos -d " " |perl -pe "s/-mRNA//g" > $speciesName.geneList
  # GO annotations
  cat $speciesName.geneList | perl -pe 's/(^.*$)/^$1\t/g'| egrep -f - ${Data["GO"]}    > annotation/function_annotation/$speciesName.v2.1.pep.fa.iprscan.gene.GO
  # ipr annotations
  cat $speciesName.geneList | perl -pe 's/(^.*$)/^$1\t/g'| egrep -f - ${Data["ipr"]}   > annotation/function_annotation/$speciesName.v2.1.pep.fa.iprscan.gene.ipr
  # swissprot annotations
  head ${Data["swiss"]} -n 1 > annotation/function_annotation/$speciesName.v2.1.pep.fa.blast.swissprot.best
  cat $speciesName.geneList | perl -pe 's/(^.*$)/^$1\t/g'| egrep -f - ${Data["swiss"]} >> annotation/function_annotation/$speciesName.v2.1.pep.fa.blast.swissprot.best

  # kegg annotations
  head  ${Data["K"]} -n4 >  annotation/function_annotation/$speciesName.v2.1.pep.fa.KEGG.K
  cat $speciesName.geneList | perl -pe 's/(^.*$)/^$1\t/g'| egrep -f - ${Data["K"]} >> annotation/function_annotation/$speciesName.v2.1.pep.fa.KEGG.K
  # kegg ontology annotations
  cat $speciesName.geneList | perl -pe 's/(^.*$)/^$1\t/g'| egrep -f - ${Data["ko"]} > annotation/function_annotation/$speciesName.v2.1.pep.fa.KEGG.ko
#######################################

#######################################
  # remove repeat annotations of contaminating scaffolds.
  cat $speciesName.contaminants.lst|perl -pe 's/^(.*)$/^$1\t/g'|egrep -v -f - ${Data["repGff"]} > annotation/repeat_annotation/$speciesName.v2.1.repeats.gff3
#######################################

#######################################
  #cleanup folders
  mv *lst lists/
  mv $speciesName.geneList lists
  mv $speciesName.replacementTable.tsv lists
  rm *tmp*.gff
  tree > $speciesName.v2.1.files.tree
  tree
#######################################

#######################################
# run GAG
  python ~/software/GAG/gag.py -f genome/$speciesName.v2.1.fa -g  annotation/gene_annotation/$speciesName.v2.1.gff3
#######################################

done

# mv folder
$base/v2.1/ ~/data/inqGen18/inquilines_v2.1
