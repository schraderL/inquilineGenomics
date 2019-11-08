# 1. Annotate OR genes with Zhou et al. pipeline.
# 2. Annotate OR genes with GeMoMa
# 3. incorporate RNAseq support (GMAP, PASA)
# 4. Combine evidence with EVM

###################################################
#==================================================
# REQUIRED SOFTWARE
#==================================================
###################################################

# GeMoMa  http://www.jstacs.de/index.php/GeMoMa
# pfam_scan.pl  ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz #I think that's the right link
# BLAST # should be installed already
# GFFcleaner https://gffutils.readthedocs.io/en/latest/GFFcleaner.html
# genometools http://genometools.org/
# evidenceModeler (EVM) https://evidencemodeler.github.io/
# samtools # should be installed
# exonerate https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
# gffutils http://ccb.jhu.edu/software/stringtie/gff.shtml
# GFFUtils
###################################################
#==================================================
# GENERAL NOTES
#==================================================
###################################################

#- Used improved manual annotations of Aech ORs by Sean McKenzie
#/Users/lukas/CSE/inquilineGenomics/ORevo/ManualCuration/Aech/allAechOrs.clean.final.fa
#/Users/lukas/CSE/inquilineGenomics/ORevo/ManualCuration/Aech/allAechOrs.clean.final.gff3

#- Used improved manual annotations of Ains ORs by Lukas Schrader
#/Users/lukas/CSE/inquilineGenomics/ORevo/ManualCuration/Aech/allAechOrs.clean.final.fa
#/Users/lukas/CSE/inquilineGenomics/ORevo/ManualCuration/Aech/allAechOrs.clean.final.gff3

###################################################
#==================================================
# 0. Define reusable variables
#==================================================
###################################################


############################
#Folders
############################
# set base directory
export base=/usr/local/home/lschrader/data/inqGen18/ORs
# set directory where all scripts are located
export scripts=$base/scripts
# set directory where software is installed
export software=/usr/local/home/lschrader/software
# set number of CPUs to use
export cpu=8
# set jar excecutable of GeMoMa
export GeMoMaProgram=$software/GeMoMa-1.5.2/GeMoMa-1.5.2.jar
export EVMtools=$software/EVidenceModeler-1.1.1/EvmUtils/
export gtools=$software/genometools-1.5.9/bin

############################
# Reference OR annotations
############################
# Define 3 genomes that are used as reference for annotations

# Aech
export q1=Aech
export ORfa1=$base/data/$q1/allAechOrs.clean.final.fa
export ORgff1=$base/data/$q1/allAechOrs.clean.final.gff3
export queryGenome1=/usr/local/home/lschrader/data/genomes/Aech/Aech_2.0_scaffolds.clean.fa

# Sinv
export q2=Ains
export ORfa2=$base/data/$q2/Ains.OR.manualAnnotations.fa
export ORgff2=$base/data/$q2/Ains.OR.manualAnnotations.renamed.cleaned.gff3
export queryGenome2=/usr/local/home/lschrader/data/genomes/inquilines_2.0/Acromyrmex_insinuator/genome/Acromyrmex_insinuator.v2.0.fa

# Acep
export q3=Acep
export ORfa3=$base/data/$q3/AcepORs.McKenzie.fasta
export ORgff3=$base/data/$q3/SupplementaryFile4.AcepOrs.final.gff
export queryGenome3=/usr/local/home/lschrader/data/genomes/attines/assembly/Acep.genome.fa

# cat OR fastas
export ORfa=$base/data/$q1.$q2.$q1.ORs.fa
#cat $ORfa1 $ORfa2 $ORfa3 > $ORfa

# prepare clean and well working gff files of reference annotations if necessary
###################################################################################
# for f in $ORgff1 $ORgff2 $ORgff3
# do
# gffread $f -o tmp.gff3 --force-exons -E
# GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
# rm tmp.gff3
# $software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
# rm tmp_clean.gff
# perl $scripts/addGeneGFF.pl tmp2.gff > $f
# done
###################################################################################

###################################################################################
# Species data
###################################################################################

#set which genome to annotate. Requires a species abbreviation, a folder where the assembly fasta file is located, a link to the fasta file
#and a file with transcript data if available.

# usage: bash scripts/ORannotation.190405.sh Ahey /corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_heyeri.2.1/genome Acromyrmex_heyeri.v2.1.fa
export species=$1
export genomeFolder=$2
export genomeFa=$2"/"$3
export transcript=

echo "echo s1: $q1 s2: $q2 s3: $q3"
echo "base: $base"
echo ""
echo "$0 for $species on $2/$3"

read -p "start? y/n " confirmation

if [[ ! "$confirmation"  == "y" ]]
then
    echo "Aborted."
    exit 1
fi


###################################################
#==================================================
# 0b. prepare folders
#==================================================
###################################################

# create all directories that will be required for the pipeline
echo "create all directories that will be required for the pipeline"
echo "$base/$species"
echo "$base/$species/EVM"
echo "$base/$species/EVM/tmpGff"
echo "$base/$species/exonerate"
echo "$base/$species/exonerate/gff"
echo "$base/$species/exonerate/raw"
echo "$base/$species/exonerate/fa"
echo "$base/$species/exonerate/PROTEIN"
echo "$base/$species/exonerate/protFa/"
echo "$base/$species/GeMoMa/"
echo "$base/$species/GeMoMa/$q1"
echo "$base/$species/GeMoMa/$q2"
echo "$base/$species/GeMoMa/$q3"
echo "$base/$species/singleExon/"
echo "$base/$species/final"
echo "$base/$species/final/filtering"
echo "$base/$species/pfam"

mkdir $base/$species
mkdir $base/$species/EVM
mkdir $base/$species/EVM/tmpGff
mkdir $base/$species/exonerate
mkdir $base/$species/exonerate/gff
mkdir $base/$species/exonerate/raw
mkdir $base/$species/exonerate/fa
mkdir $base/$species/exonerate/PROTEIN
mkdir $base/$species/exonerate/protFa/
mkdir $base/$species/GeMoMa/
mkdir $base/$species/GeMoMa/$q1
mkdir $base/$species/GeMoMa/$q2
mkdir $base/$species/GeMoMa/$q3
mkdir $base/$species/singleExon/
mkdir $base/$species/final
mkdir $base/$species/final/filtering
mkdir $base/$species/pfam

# prepare config files that will be used
#cp ~/data/OR/config_raw $base
#cp ~/data/OR/config_raw2 $base
#cp ~/data/OR/EVM.weights.txt $base

###################################################
#==================================================
# 0b. Run Exonerate, GeMoMa and Blast
#==================================================
###################################################

cd $base/$species

# run in parallel
# Predict ORs with Exonerate pipeline (cf Zhou et al. 2004)
echo "Predict ORs with Exonerate pipeline (cf Zhou et al. 2004)"
nice bash $scripts/runExonerate.sh > runExonerate.output 2> runExonerate.err &
# Predict ORs with GeMoMa
echo "# Predict ORs with GeMoMa"
nice bash $scripts/runGeMoMa.pipe.sh > runGeMoMa.output 2> runGeMoMa.err &
# blast reference OR exons against target genome for further evidence
echo "blast reference OR exons against target genome for further evidence"
nice bash $scripts/runExonBlast.sh > runExonBlast.output 2> runExonBlast.err &

# wait till the 3 jobs above are finished (takes several hours)
wait;

# Run EVM to combine models from the 3 jobs
cd $base/$species
# Set weighting of different lines of evidence in EVM.weights.txt
cp $base/EVM.weights.txt $base/$species
bash $scripts/prepareEVM.sh > prepareEVM.output 2> prepareEVM.err;
bash $scripts/runEVM.sh > runEVM.output 2> runEVM.err


############################################################
# 4. Add incomplete genes from GeMoMa Aech mapping
############################################################
bash $scripts/addIncompleteGenes.sh
############################################################
# 5. Scan for domains with pfam_scan.pl
############################################################
bash $scripts/filterForDomain.sh
#################################################################
# 6. Scan proteins without 7tm6 for similarity to ORs with blastp
#################################################################
cd $base/$species/final/filtering

# blastp against OR db
makeblastdb -in $ORfa -dbtype prot
nice blastp -db $ORfa -query $base/$species/final/$species.OR.no7tm6.fa -out $species.OR.no7tm6.bls -num_threads $cpu -evalue 1e-5 -outfmt 6
# extract best hit
cat $species.OR.no7tm6.bls |sort -gk11,11 |sort -buk 1,1 >  $species.OR.no7tm6.tophit

# List containing all gene models that have a hit against an OR
# $species.OR.no7tm6.tophit

#################################################################
# 7. Filter and rename
#################################################################

cd $base/$species/final
# check fullOR annotation for pfam domains
#nice perl $software/PfamScan/pfam_scan.pl -fasta $species.fullORannotation.fa -dir ~/data/pfam/ -outfile $species.fullORannotation.domains.out -cpu 30

#Naming conventions
#- count as hit if e-value < 1e-5 [Ains.OR.no7tm6.tophit]
#- what if two separate models should in fact be one model? Ignore...
#- If M is missing > NTE
#- If * is missing > CTE
#- If both are missing > NC
#- If 6tm7 is missing > INT
#- non-canonical splice site?
#- if 6tm7 truncated (i.e. more than 50 aa missing, based on hmm alignment) > fd or df
#- H if homologous to other OR

cd $base/$species/final
mkdir finalSet
# use complex awk expression to create a table to use for relabelling
cat $species.fullORannotation.fa|tr " " "@" |tr "\t" "@"\
|awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' \
|  awk 'BEGIN {OFS="\t"}; { \
      if ($2 !~ /^M.*/ && $2 !~ /.*\*$/) print $0,"NC";
      else if ($2 !~ /^M.*/ && $2 ~ /.*\*$/) print $0,"NTE";
      else if ($2 ~ /^M.*/ && $2 !~ /.*\*$/) print $0,"CTE";
      else if ($2 ~ /^M.*\*/) print $0," "
      }'|tr "@" " "> $species.fullORannotation.qc

cat $species.fullORannotation.domains.out|perl -pe 's/ +/\t/g' > tmp.out

# Rename ORs in R
Rscript $scripts/rename.R $base"/"$species"/final" $species
# Fix genes where gff is not perfect yet.
egrep $species"Or-" renamed.gff3 > tmp.gff3
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
# Clean up the gff
$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries tmp_clean.gff > $species.OR.gff3
# validate that the gff is ok now.
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.OR.gff3
cat renamed.tsv|awk 'BEGIN {FS="\t"}; {print ">"$33,$34"\n"$3}'|fold -w 80 > $species.OR.fa
rm tmp.gff3 renamed.gff3 tmp_clean.gff
mv $species.OR.gff3 finalSet
mv $species.OR.fa finalSet
mv renamed.tsv finalSet/$species.OR.tsv

echo "Files: " > finalSet/readme.txt
echo " $species.OR.gff3 : GFF3 file of OR annotations" >> finalSet/readme.txt
echo " $species.OR.fa : protein fa file of OR annotations" >> finalSet/readme.txt
echo " $species.OR.tsv : annotation details of OR annotations" >> finalSet/readme.txt
echo "Nomenclature: " >> finalSet/readme.txt
echo " NTE = missing start, i.e. missing N-terminus." >> finalSet/readme.txt
echo " CTE = missing stop, i.e. missing C-terminus." >> finalSet/readme.txt
echo " NC = missing start and stop, i.e. missing N- and C-terminus." >> finalSet/readme.txt
echo " NC = missing start and stop, i.e. missing N- and C-terminus." >> finalSet/readme.txt
echo " fd = domain fragmented at n-terminus. Missing N-terminal piece of domain (>50 aa)." >> finalSet/readme.txt
echo " df = domain fragmented at c-terminus. Missing C-terminal piece of domain (>50 aa)." >> finalSet/readme.txt
echo " fdf = domain fragmented at n- and c-terminus. Missing pieces of domain >50 aa." >> finalSet/readme.txt
echo " dH = domain missing, but protein homologous to other ORs. (blast 1e-10)" >> finalSet/readme.txt

#################################################################
# 8. Compute some statistics
#################################################################
cd $base/$species/final/finalSet
mkdir stats
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat --genelengthdistri -o stats/$species.genelength.stats.txt
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -genescoredistri -o stats/$species.genescore.stats.txt
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -exonlengthdistri -o stats/$species.exonlength.stats.txt
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -exonnumberdistri -o stats/$species.exonnumber.stats.txt
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -intronlengthdistri -o stats/$species.intronlength.stats.txt
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -cdslengthdistri -o stats/$species.cdslength.stats.txt
python $software/GAG/gag.py -f $genomeFa -g $species.OR.gff3
