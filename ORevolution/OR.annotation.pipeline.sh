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

#NTE	Missing sequence at N terminus#INT	Missing sequence in the middle of gene#CTE	Missing sequence at C terminus#NI	Missing N terminus and section in the middle of gene#NC	Missing N and C terminus#IC	Missing section in the middle of gene and C terminus
#PSE	Pseudogene
#P+N/I/C	Pseudogene and missing sequence
#(F)	Could be functional with assembly-introduced false frameshift
#(S)	Could be functional with non-canonical splice sites


#- maybe use prrn as well (http://www.genome.ist.i.kyoto-u.ac.jp/~aln_user/prrn/)
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
export ORfa=$base/data/Aech.Sinv.Acep.ORs.fa
export GeMoMaProgram=$software/GeMoMa_1.3/GeMoMa-1.3.jar

# Aech
export q1=Aech
export ORfa1=$base/data/AechORs.McKenzie.fasta
export ORgff1=$base/data/SupplementaryFile4.AechOrs.final.gff
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
#Acep
export species=Acep
export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
export genomeFa=$genomeFolder"/Acep.genome.fa"
export transcript=

#Acol
export species=Acol
export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
export genomeFa=$genomeFolder"/Acol.v1.0.mask.fa"
export transcript=/usr/local/home/lschrader/data/genomes/attines/cufflinks_transcripts/$species/transcripts.gtf


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
export transcript=~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf

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
export transcript=~/data/genomes/inquilines/v1.1/Acha/Acromyrmex_charruana.transcripts.gtf



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
mkdir $base/$species/exonerate/protFa/
mkdir $base/$species/GeMoMa/
mkdir $base/$species/GeMoMa/$q1
mkdir $base/$species/GeMoMa/$q2
mkdir $base/$species/GeMoMa/$q3
mkdir $base/$species/singleExon/
mkdir $base/$species/final
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



bash $scripts/prepareEVM.sh > prepareEVM.output 2> prepareEVM.err
bash $scripts/runEVM.sh > runEVM.output 2> runEVM.err


############################################################
# 4. Scan for domains with pfam_scan.pl
############################################################
#mkdir $base/$species/pfam
cd $base/$species/pfam
#perl $software/PfamScan/pfam_scan.pl -fasta ../GeMoMa/predicted_protein.fasta -dir ~/data/pfam/ -outfile $species.GeMoMa.domains.out -cpu 20
perl $software/PfamScan/pfam_scan.pl -fasta ../EVM/$species.final.OR.fa -dir ~/data/pfam/ -outfile $species.final.domains.out -cpu 20
#perl $software/PfamScan/pfam_scan.pl -fasta ../exonerate/$species.exonerate.pep -dir ~/data/pfam/ -outfile $species.exonerate.domains.out -cpu 20
grep 7tm_6 $species.final.domains.out|cut -f 1 -d " "|sort|uniq -c|wc -l

cd $base/$species/final
grep ">" $species.final.OR.fa -B1|grep -E "\*"|wc -l
grep ">" $species.final.OR.fa -A1|grep -E "^M"|wc -l


grep -E "\evm\." ../pfam/$species.final.domains.out |grep 7tm_6|cut -f 1 -d " "|sort|uniq > $species.ORs.7tm6.lst
wc -l $species.ORs.7tm6.lst
faSomeRecords $species.final.OR.fa $species.ORs.7tm6.lst $species.OR.7tm6.fa
grep ">" $species.OR.7tm6.fa -B1|grep -E "\*"|wc -l
grep ">" $species.OR.7tm6.fa -A1|grep -E "^M"|wc -l
perl $scripts/checkFastaGeneModel.pl $species.OR.7tm6.fa |grep ">"|wc -l


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



grep -E "AechOr" Aech.final.domains.out |grep 7tm_6|cut -f 1 -d " " > Aech.ORs.7tm6.lst
faSomeRecords AechORs.McKenzie.fasta Aech.ORs.7tm6.lst Aech.OR.7tm6.fa
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
