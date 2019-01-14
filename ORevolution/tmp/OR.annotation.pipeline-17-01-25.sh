# 1. Annotate OR genes with Zhou et al. pipeline.
# 2. Annotate OR genes with GeMoMa
# 3. incorporate RNAseq support (GMAP, PASA)
# 4. Combine evidence with EVM
#https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz
# 5. inspect models with IGV?

#https://www.biostars.org/p/181286/
#https://github.com/nextgenusfs/funannotate

###################################################
#==================================================
# GENERAL NOTES
#==================================================
###################################################

#- Used manual annotations of Aech ORs by Sean McKenzie

#- Legend for OR annotation

#NTE	Missing sequence at N terminus
#INT	Missing sequence in the middle of gene
#CTE	Missing sequence at C terminus
#NI	Missing N terminus and section in the middle of gene
#NC	Missing N and C terminus
#IC	Missing section in the middle of gene and C terminus
#PSE	Pseudogene
#P+N/I/C	Pseudogene and missing sequence
#(F)	Could be functional with assembly-introduced false frameshift
#(S)	Could be functional with non-canonical splice sites


#- maybe use prrn as well (http://www.genome.ist.i.kyoto-u.ac.jp/~aln_user/prrn/)


# pipeline doesn't work for genes with very long introns.
# For Aech these are
# AechOr219
# AechOr278
# AechOr424
# AechOr222
# AechOr421
# AechOr241CTE
# AechOr392
# AechOr240CTE
# AechOr277CTE



###################################################
#==================================================
# THOUGHTS
#==================================================
###################################################

#- use --maxintron 2000 4000 and no max for exonerate
#- blast individual OR exons against genomes
#- check for 6tm domain in each exonerate, cufflinks and GeMoMa model and save with different gff name tags. Then in EVM, we can use this info to assign weights to these different models. I.E. if 2 models are conflicting, the one that codes for a 6tm should be weighted much higher.
#- https://bitbucket.org/jnmaloof/gfftools
#- Remember to use wolbachia-free genome sequences!

###################################################
#==================================================
# 0. Define reusable variables
#==================================================
###################################################

#Ahey Test scf49
#export species=Ahey.test
#cd ~/data/OR/Ahey.test
#export genomeFolder=/usr/local/home/lschrader/data/genomes/inquilines/v1.1/Ahey
#export genomeFa=$genomeFolder"/Acromyrmex_heyeri.fa"
#faSomeRecords $genomeFa scf.select Ahey.scf49.fa
#export genomeFolder=/usr/local/home/lschrader/data/OR/Ahey.test
#export genomeFa=$genomeFolder"/Ahey.scf49.fa"
#export transcript=$genomeFolder/Acromyrmex_heyeri.transcripts.gtf

###################################################################################
#Folders
###################################################################################
export base=/usr/local/home/lschrader/data/OR
export scripts=$base/scripts
export software=/usr/local/home/lschrader/software
export cpu=8
#cat $ORfa1 $ORfa2 $ORfa3 > $ORfa
export ORfa=$base/data/Aech.Sinv.Acep.ORs.fa
#export GeMoMaProgram=$software/GeMoMa_1.3/GeMoMa-1.3.jar
#export GeMoMaProgram=$software/GeMoMa_1.3.2/GeMoMa-1.3.2.jar
export GeMoMaProgram=$software/GeMoMa_1.4/GeMoMa-1.4.jar

# Aech
export q1=Aech
export ORfa1=$base/data/allAechOrs.clean.final.fa
export ORgff1=$base/data/allAechOrs.clean.final.gff3
export queryGenome1=/usr/local/home/lschrader/data/genomes/Aech/Aech_2.0_scaffolds.clean.fa

# Sinv
export q2=Sinv
export ORfa2=$base/data/SinvORs.McKenzie.fasta
export ORgff2=$base/data/SupplementaryFile4.SinvOrs.final.gff
export queryGenome2=/usr/local/home/lschrader/data/genomes/Sinv/Sinv.Si_gnF.scf.fa

# Acep
export q3=Acep
export ORfa3=$base/data/AcepORs.McKenzie.fasta
export ORgff3=$base/data/SupplementaryFile4.AcepOrs.final.gff
export queryGenome3=/usr/local/home/lschrader/data/genomes/attines/assembly/Acep.genome.fa





# prepare gff if necessary
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

# add Sinv, Acep (and other) manual OR annotations for the prediction
###################################################################################

###################################################################################
# Species data
###################################################################################

#Aech
export species=Aech
export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
export genomeFa=$genomeFolder"/"$species".genome.fa"
export transcript=

#Acep
export species=Acep
export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
export genomeFa=$genomeFolder"/Acep.genome.fa"
export transcript=

#Acol
export species=Acol
export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
export genomeFa=$genomeFolder"/Acol.v1.0.mask.fa"
export transcript=
#/usr/local/home/lschrader/data/genomes/attines/cufflinks_transcripts/$species/transcripts.gtf


#Tsep
export species=Tsep
export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
export genomeFa=$genomeFolder"/$species.v1.0.filterBacteria.mask.fa"
export transcript=/usr/local/home/lschrader/data/genomes/attines/cufflinks_transcripts/$species/transcripts.gtf

#Tzet
export species=Tzet
export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
export genomeFa=$genomeFolder"/$species.v1.0.filterBacteria.mask.fa"
export transcript=/usr/local/home/lschrader/data/genomes/attines/cufflinks_transcripts/$species/transcripts.gtf

#Tcor
export species=Tcor
export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
export genomeFa=$genomeFolder"/$species.v1.0.filterBacteria.mask.fa"
export transcript=/usr/local/home/lschrader/data/genomes/attines/cufflinks_transcripts/$species/transcripts.gtf

#Ccos
export species=Ccos
export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
export genomeFa=$genomeFolder"/$species.v1.0.filterBacteria.mask.fa"
export transcript=/usr/local/home/lschrader/data/genomes/attines/cufflinks_transcripts/$species/transcripts.gtf



#Parg
export species=Parg
export genomeFolder=/usr/local/home/lschrader/data/genomes/inquilines/v1.1/Parg
export genomeFa=$genomeFolder"/Pseudoatta_argentina.fa"
export transcript=
#~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf

#Ains
export species=Ains
export genomeFolder=/usr/local/home/lschrader/data/genomes/inquilines/v1.0/Ains
export genomeFa=$genomeFolder"/Acromyrmex_insinuator.fa"
export transcript=~/data/genomes/inquilines/v1.0/Ains/Acromyrmex_insinuator.transcripts.gtf

#Ahey
export species=Ahey
export genomeFolder=/usr/local/home/lschrader/data/genomes/inquilines/v1.1/Ahey
export genomeFa=$genomeFolder"/Acromyrmex_heyeri.fa"
export transcript=$genomeFolder/Acromyrmex_heyeri.transcripts.gtf

#Acha
export species=Acha
export genomeFolder=/usr/local/home/lschrader/data/genomes/inquilines/v1.1/Acha
export genomeFa=$genomeFolder"/Acromyrmex_charruana.fa"
export transcript=
#~/data/genomes/inquilines/v1.1/Acha/Acromyrmex_charruana.transcripts.gtf



###################################################
#==================================================
# 0b. prepare folders
#==================================================
###################################################


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

cd $base/$species

# run parallel
bash $scripts/runExonerate.sh > runExonerate.output 2> runExonerate.err &
bash $scripts/runGeMoMa.pipe.sh > runGeMoMa.output 2> runGeMoMa.err &
bash $scripts/runExonBlast.sh > runExonBlast.output 2> runExonBlast.err &

#wait until finished

# currently best weights.txt settings:
# ABINITIO_PREDICTION     exonerate:protein2genome:local  1
# OTHER_PREDICTION     GeMoMa_Aech     4
# OTHER_PREDICTION     GeMoMa_Acep     4
# OTHER_PREDICTION        GeMoMa_Sinv     3
# TRANSCRIPT      transdecoder    2
# PROTEIN blast   2

cd $base/$species
cp $base/EVM.weights.txt $base/$species
bash $scripts/prepareEVM.sh > prepareEVM.output 2> prepareEVM.err
bash $scripts/runEVM.sh > runEVM.output 2> runEVM.err


############################################################
# 4. Add incomplete genes from GeMoMa Aech mapping
############################################################
#bedtools intersect -b Acep.final.OR.gff3 -a Acep.vs.Aech.GeMoMa.final.gff -wb -v|awk '{if ($3=="gene") print $0}'|perl -pe 's/.*ID=(.*?)(\n|;|\t|,)/$1\n/g'|wc -l
cd $base/$species/final/
cp $base/$species/GeMoMa/$species.vs.Aech.GeMoMa.final.gff .
bedtools intersect -b $species.final.OR.gff3 -a $species.vs.Aech.GeMoMa.final.gff -wb -v|awk '{if ($3=="gene") print $0}'|perl -pe 's/.*ID=(.*?)(\n|;|\t|,)/$1,/g'|sed -e 's/,$//' > 2add.ids.lst
rm $species.vs.Aech.GeMoMa.final.gff.db
gffutils-cli create $species.vs.Aech.GeMoMa.final.gff
#gffutils-cli children Acep.vs.Aech.GeMoMa.final.gff.db AECHOR405NTE-MRNA-1_R55-gene
sed -e s/^/'gffutils-cli children $species.vs.Aech.GeMoMa.final.gff.db '/g 2add.ids.lst >2add.sh
bash 2add.sh > 2add.dirty.gff3
awk '!a[$0]++' 2add.dirty.gff3 > 2add.gff3
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= 2add.gff3
$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -addintrons 2add_clean.gff > 2add.tmp.gff
$software/genometools-1.5.9/bin/gt merge -tidy -retainids 2add.tmp.gff   $species.final.OR.gff3 > $species.fullORannotation.gff3
rm 2add*

perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.fullORannotation.gff3

echo "Writing protein fasta"
cd $base/$species/final/
$software/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl   $species.fullORannotation.gff3 $genomeFa > $species.fullORannotation.fa


############################################################
# 5. Scan for domains with pfam_scan.pl
############################################################
#mkdir $base/$species/pfam
cd $base/$species/final
#perl $software/PfamScan/pfam_scan.pl -fasta ../GeMoMa/predicted_protein.fasta -dir ~/data/pfam/ -outfile $species.GeMoMa.domains.out -cpu 20
nice perl $software/PfamScan/pfam_scan.pl -fasta $species.fullORannotation.fa -dir ~/data/pfam/ -outfile $species.fullORannotation.domains.out -cpu 30
nice perl $software/PfamScan/pfam_scan.pl -fasta allAechOrs.clean.final.fa -dir ~/data/pfam/ -outfile allAechOrs.clean.final.pfam -cpu 30
#perl $software/PfamScan/pfam_scan.pl -fasta ../exonerate/$species.exonerate.pep -dir ~/data/pfam/ -outfile $species.exonerate.domains.out -cpu 20
grep 7tm_6 $species.fullORannotation.domains.out|cut -f 1 -d " "|sort|uniq -c|wc -l


awk '{if ($7=="7tm_6") print $1}' $species.fullORannotation.domains.out|sort|uniq > $species.ORs.7tm6.lst
wc -l $species.ORs.7tm6.lst
faSomeRecords ../final/$species.fullORannotation.fa $species.ORs.7tm6.lst $species.OR.7tm6.fa
faSomeRecords -exclude ../final/$species.fullORannotation.fa $species.ORs.7tm6.lst $species.OR.no7tm6.fa
perl $scripts/checkFastaGeneModel.pl $species.OR.7tm6.fa |grep ">"|wc -l
perl $scripts/checkFastaGeneModel.pl $species.OR.no7tm6.fa |grep ">"|wc -l
grep ">" -v ../final/$species.fullORannotation.fa|tr "\n" "-"|sed s/-//g| wc -c

# cd /usr/local/home/lschrader/data/OR/data
# awk '{if ($7=="7tm_6") print $1}' *.pfam|sort|uniq > allAechOrs.clean.final.pfam.lst
# faSomeRecords allAechOrs.clean.final.fa allAechOrs.clean.final.pfam.lst allAechOrs.clean.final.7tm_6.fa
# perl $scripts/checkFastaGeneModel.pl allAechOrs.clean.final.7tm_6.fa |grep ">"|wc -l
# grep ">" -v allAechOrs.clean.final.7tm_6.fa|tr "\n" "-"|sed s/-//g| wc -c
#

#################################################################
# 5. Scan proteins without 7tm6 for similarity to ORs with blastp
#################################################################
cd $base/$species/final/filtering

#perl $scripts/checkFastaGeneModel.pl ../$species.OR.7tm6.fa |nice blastp -db $ORfa -query - -num_threads 30 -evalue 1e-5 -outfmt 6 |sort -gk11,11 |sort -buk 1,1
#perl $scripts/checkFastaGeneModel.pl ../$species.OR.no7tm6.fa |nice blastp -db $ORfa -query - -num_threads 30 -evalue 1e-5 -outfmt 6 |sort -gk11,11 |sort -buk 1,1
nice blastp -db $ORfa -query $base/$species/final/$species.OR.no7tm6.fa -out $species.OR.no7tm6.bls -num_threads 30 -evalue 1e-5 -outfmt 6
cat $species.OR.no7tm6.bls |sort -gk11,11 |sort -buk 1,1 >  $species.OR.no7tm6.tophit
wc -l $species.OR.no7tm6.tophit

# count as one if
  # matches to same subject protein (e.g. )

#	domains	stop	Start	total	all3
#Acep 360	302		313		426		275
#Acol 396	338		385		498		314
#Parg 236 	208		208		282		183
#Acha 300 	226		275		431		198
#Ahey 361	300		323		424		275
#Ains 391	337		352		447		322
#Aech 420	410		379		432		357
#Tsep 361	305		323		453		289
#Tzet 332	286		291		374		269
#Tcor 392	331		329		446		308
#Ccos 321	274		289		364		262
#grep 7tm_6 *.final.domains.out|cut -f 1 -d " "|sort|uniq -c|wc -l



grep -E "AechOr" $species.final.domains.out |grep 7tm_6|cut -f 1 -d " " > $species.ORs.7tm6.lst
faSomeRecords AechORs.McKenzie.fasta Aech.ORs.7tm6.lst $species.OR.7tm6.fa
grep ">" Aech.OR.7tm6.fa -B1|grep -E "\*"|wc -l
grep ">" Aech.OR.7tm6.fa -A1|grep -E "^M"|wc -l
perl $scripts/checkFastaGeneModel.pl Aech.OR.7tm6.fa |grep ">"|wc -l

cd $base/$species/final
cat $species.final.OR.gff3|$software/genometools-1.5.9/bin/gt stat --genelengthdistri -o $species.genelength.stats.txt
cat $species.final.OR.gff3|$software/genometools-1.5.9/bin/gt stat -genescoredistri -o $species.genescore.stats.txt
cat $species.final.OR.gff3|$software/genometools-1.5.9/bin/gt stat -exonlengthdistri -o $species.exonlength.stats.txt
cat $species.final.OR.gff3|$software/genometools-1.5.9/bin/gt stat -exonnumberdistri -o $species.exonnumber.stats.txt
cat $species.final.OR.gff3|$software/genometools-1.5.9/bin/gt stat -intronlengthdistri -o $species.intronlength.stats.txt
cat $species.final.OR.gff3|$software/genometools-1.5.9/bin/gt stat -cdslengthdistri -o $species.cdslength.stats.txt

#	domains	complete	total # domains	total # genes	total aa
#Acep	452	397	463	742	206713
#Acol
#Parg	262	244	273	1619	330886
#Acha
#Ahey	453	392	474	631	181644
#Ains	421	410	428	495	172080
#Aech	454	419	460	558	179631
#Tsep
#Tzet
#Tcor
#Ccos
