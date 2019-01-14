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



#Folders
export base=/usr/local/home/lschrader/data/OR
export scripts=$base/scripts
export ORfa=$base/AechORs.McKenzie.fasta
export ORgff=$base/SupplementaryFile4.AechOrs.final.gff
export queryGenome=/usr/local/home/lschrader/data/genomes/Aech/Aech_2.0_scaffolds.clean.fa
export software=/usr/local/home/lschrader/software


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
mkdir $base/$species/singleExon/
mkdir $base/$species/final
mkdir $base/$species/pfam
###################################################
#==================================================
# 1. Annotate OR genes with Zhou et al. pipeline.
#==================================================
###################################################

######################################################################################################
#1/ run the perl script identify_geneLS.pl on the target genome using known genes as queries
######################################################################################################

if [ ! -f $genomeFa.nhr ]; 
then
	cd $genomeFolder
	makeblastdb -in $genomeFa -dbtype nucl
fi
#mkdir $base/$species
cd $base/$species
#edited config to eval 1e-4 for Ahey.test
cat ../config_raw |envsubst > config_$species

perl $scripts/identify_geneLS.pl config_$species $species.bls > $species.out

######################################################################################################	
#2/ run the perl script "hsp2bed.pl" to get loci
######################################################################################################

# 100 is the aa cutoff
perl $scripts/hsp2bed.pl $species.out 100 > $species.bed
	

######################################################################################################	
#3/ Prepare exonerate input
######################################################################################################

cat ../config_raw2 |envsubst > config_$species"2"
perl $scripts/runExonerate.pl config_$species"2"


cd $base/$species
#mkdir exonerate
#mkdir exonerate/gff
#mkdir exonerate/raw
#mkdir exonerate/fa
#mkdir exonerate/protFa/

######################################################################################################	
#4/ Run exonerate on each locus
######################################################################################################
#AecOr291 exonerate
#f=AechFa/AechOr100PSE.fa
#f=AechFa/AechOr21.fa
for f in AechFa/*.fa
do
	gene=$(echo $f | cut -f 2 -d "/"|cut -f 1 -d ".")
	echo "working on $gene"
	target=$species".Loci/$gene.target.fa"
	nice exonerate --query $f --target $target -m protein2genome --showtargetgff yes --showalignment no --ryo '###---FASTASTART---###\n>%qi length=%tl alnlen=%tal target=%ti_(%tcb-%tce)\n%tcs###---FASTAEND---###\n' > exonerate/raw/$gene.exonerate
	perl $scripts/processExonerate2GFF.pl exonerate/raw/$gene.exonerate | sed s/cds/CDS/g > dirty.gff
	perl $scripts/cleanExonerateGFF.pl dirty.gff > tmp.gff
	rm dirty.gff
	GFFcleaner --add-missing-ids --add-exon-ids --insert-missing= tmp.gff
	rm tmp.gff
	$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > exonerate/gff/$gene.gff
	perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl exonerate/gff/$gene.gff
	perl $scripts/processExonerate2Fa.pl exonerate/raw/$gene.exonerate > exonerate/fa/$gene.exonerate.fa
	perl $scripts/translateCDS2AA.pl exonerate/fa/$gene.exonerate.fa > exonerate/protFa/$gene.exonerate.pep
	echo "$gene done. :)"
done

cat exonerate/protFa/* > exonerate/$species.exonerate.pep


###################################################
#==================================================
# 2. Annotate OR genes with GeMoMa
#==================================================
###################################################


###################################################
#1) run GeMoMa
###################################################

#mkdir $base/$species/GeMoMa/
cd $base/$species/GeMoMa/
cp ../../runGeMoMa.sh .
bash ./runGeMoMa.sh $ORgff ~/data/genomes/Aech/Aech_2.0_scaffolds.clean.fa $genomeFa $base/$species/GeMoMa/



###################################################
#==================================================
# 3. Blast single exons against target
#==================================================
###################################################


###################################################
#1) prepare single exon fa
###################################################

#mkdir $base/$species/singleExon/
cd $base/$species/singleExon/
awk '{if ($3 == "exon") print $0}' $ORgff > ORexons.gff
awk 'BEGIN {OFS=""};{ gsub(/\;Parent.*/, "")};{ gsub(/ID=/, "")};{print "samtools faidx $queryGenome ",$1,":",$4,"-",$5," |sed s/\\>.*/\\>",$9,":",$1,":",$4,"-",$5,"/g >> ORexons.fa"}' ORexons.gff > samtools.ORexons.sh
rm -f ORexons.fa 
bash samtools.ORexons.sh #AechORexons.fa
blastn -db $genomeFa -query ORexons.fa -out ORexons.bls -num_threads 16 -evalue 1e-2 -outfmt 6 
perl $scripts/bls2gff.pl ORexons.bls > $species.OR.exons.bls.gff
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.OR.exons.bls.gff

#############################################
#===========================================
#	4. COMBINE EVIDENCE WITH EVM
#===========================================
#############################################
#mkdir $base/$species/EVM/

#############################################
#1 reformat GeMoMa gff to pass EVM specs
#############################################

cd $base/$species/GeMoMa
gffread predicted_annotation.gff -o tmp.gff3 --force-exons -E 
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3 
rm tmp.gff3
$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
rm tmp_clean.gff
perl $scripts/addGeneGFF.pl tmp2.gff > $species.GeMoMa.final.gff
# could also test ~/.local/bin/PASApipeline-2.0.2/misc_utilities/gtf_to_gff3_format.pl
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.GeMoMa.final.gff
cp $species.GeMoMa.final.gff $base/$species/EVM

#############################################
#2 reformat exonerate gff to pass EVM specs
#############################################

cd $base/$species/exonerate/gff
 $software/genometools-1.5.9/bin/gt merge -tidy -retainids *.gff  > ../$species.exonerate.gff
cd ..
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.exonerate.gff
cp $species.exonerate.gff $base/$species/EVM/$species.exonerate.final.gff


#############################################
#3 reformat cufflinks gtf to pass EVM specs
#############################################
#http://transdecoder.github.io/

cd $base/$species/EVM/
bedtools intersect -wa -a $transcript -b $species.evid.gff |grep gene|uniq |grep -P "gene_id \".+?\";" -o|sed -r s/'gene_id |\;'//g|sed s/\"//g > OR.raw.transcripts.lst
grep -f OR.raw.transcripts.lst $transcript > OR.raw.transcripts.gff3
$software/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl OR.raw.transcripts.gff3 $genomeFa > $species.cufflinks.fa
$software/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl OR.raw.transcripts.gff3 >  $species.transcript.gff
TransDecoder.LongOrfs -t $species.cufflinks.fa 
#nice $software/TransDecoder-3.0.1/TransDecoder.Predict -t $species.cufflinks.fa --cpu 15
$software/TransDecoder-3.0.1/util/cdna_alignment_orf_to_genome_orf.pl $species.cufflinks.fa.transdecoder_dir/longest_orfs.gff3  $species.transcript.gff $species.cufflinks.fa > transcripts.fasta.transdecoder.genome.gff3
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl transcripts.fasta.transdecoder.genome.gff3
cp transcripts.fasta.transdecoder.genome.gff3 OR.transcripts.gff3
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl OR.transcripts.gff3

#############################################
#4 get single exon blast
#############################################
cd $base/$species/EVM/
cp ../singleExon/$species.OR.exons.bls.gff .

#############################################
#5 get relevant transcripts for OR genes
#############################################

# need to concatenate the gfffiles for exonerate and GeMoMa 
cd $base/$species/EVM
$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes $species.OR.exons.bls.gff > $species.OR.exons.bls.sorted.gff

$software/genometools-1.5.9/bin/gt merge -tidy -retainids $species.exonerate.final.gff $species.GeMoMa.final.gff  $species.OR.exons.bls.sorted.gff > $species.evid.gff
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.evid.gff

$software/genometools-1.5.9/bin/gt merge -tidy -retainids $species.exonerate.final.gff $species.GeMoMa.final.gff > $species.prot.gff
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.prot.gff

#GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= $species.prot.gff
#$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes $species.prot_clean.gff > $species.prot.gff


# get intersecting cufflinks transcripts with $species.prot.sort.gff #is done already
#cd $base/$species/EVM
#bedtools intersect -wa -a transcripts.fasta.transdecoder.genome.gff3 -b $species.evid.gff |grep gene|uniq |grep -P "ID=.+?;" -o|sed -r s/'ID=|\;'//g > OR.transcripts.lst
#grep -f OR.transcripts.lst transcripts.fasta.transdecoder.genome.gff3 > OR.transcripts.gff3
#perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl OR.transcripts.gff3

#############################################
#6 Filter Fasta for relevant scaffolds
#############################################
cat $species.exonerate.final.gff $species.GeMoMa.final.gff $species.OR.exons.bls.gff |cut -f 1 |grep -P "^#" -v|sort|uniq > OR.scf.lst
faSomeRecords $genomeFa OR.scf.lst $species.OR.scf.fa
genomeFa=$species.OR.scf.fa

#############################################
#7 Prepare EVM input files
#############################################

cd $base/$species/EVM
cp ../../EVM.weights.txt weights.txt

#weights.txt
#ABINITIO_PREDICTION	exonerate:protein2genome:local	3
#ABINITIO_PREDICTION	GeMoMa	4
#TRANSCRIPT	transdecoder	10 #excluded in Acep
#PROTEIN	blast	3



# Partitioning  
# --transcript_alignments removed for Acep
genomeFa=$species.OR.scf.fa
$software/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome $genomeFa \
     --gene_predictions $species.prot.gff \
     --transcript_alignments OR.transcripts.gff3 \
     --protein_alignments $species.OR.exons.bls.gff \
     --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out  


# write commands
# --transcript_alignments removed for Acep
$software/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome $genomeFa --weights $base/$species/EVM/weights.txt \
      --gene_predictions $species.prot.gff \
      --transcript_alignments OR.transcripts.gff3 \
      --protein_alignments $species.OR.exons.bls.gff \
      --search_long_introns 2000 \
      --re_search_intergenic 6000 \
      --output_file_name evm.out  --partitions partitions_list.out >  commands.list  



#run
#$software/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log

# run in parallel
#split --number=l/6 ${fspec} xyzzy
lines=$(wc -l commands.list|awk '{print $1/10}')
wc -l commands.list  
split --lines=$lines commands.list
for f in xa*
do
 $software/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl $f > run.$f.log &
done

#wait till all runs have finished
rm xa*

#combine
$software/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

#write gffs
$software/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome $genomeFa

rm -rf tmpGff/
mkdir tmpGff
for f in */evm.out.gff3
do
	scf=$(echo $f | cut -f 1 -d "/")
	$software/genometools-1.5.9/bin/gt gff3 -sort -tidy -retainids -checkids -fixregionboundaries -addintrons $f > tmpGff/$scf.gff
done

$software/genometools-1.5.9/bin/gt merge -tidy -retainids tmpGff/*.gff > $species.final.OR.gff3

#write protein fasta

for f in */evm.out.gff3
do
if [ -s $f ] 
	then
	echo "$software/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl  $f $genomeFa > $f.prot.fa" >> getFasta.sh
fi	
done

wc -l getFasta.sh
lines=$(wc -l getFasta.sh|awk '{print int($1/8)}')
split --lines=$lines getFasta.sh
for f in xa*
do
 bash $f &
done

cat */evm.out.gff3.prot.fa > $species.final.OR.fa




mkdir ../final
cp $species.final.OR.gff3 ../final
cp $species.final.OR.fa ../final
cp $species.prot.gff ../final
cp OR.transcripts.gff3 ../final/$species.OR.transcripts.gff3
cp $species.OR.scf.fa ../final
cp $species.GeMoMa.final.gff ../final/
cp $species.exonerate.final.gff ../final/
cp ../GeMoMa/predicted_protein.fasta ../final/$species.GeMoMa.fa
cp ../exonerate/$species.exonerate.pep ../final/
cp $species.OR.exons.bls.sorted.gff ../final/

cd ../final/

awk '{if ($3=="gene") print $0}' $species.final.OR.gff3 | wc -l
grep ">" *.final.OR.fa -B1|grep -E "\*"|wc -l
grep ">" *.final.OR.fa -A1|grep -E "^M"|wc -l


############################################################
# 4. Scan for domains with pfam_scan.pl
############################################################
#mkdir $base/$species/pfam
cd $base/$species/pfam
#perl $software/PfamScan/pfam_scan.pl -fasta ../GeMoMa/predicted_protein.fasta -dir ~/data/pfam/ -outfile $species.GeMoMa.domains.out -cpu 20
perl $software/PfamScan/pfam_scan.pl -fasta ../EVM/$species.final.OR.fa -dir ~/data/pfam/ -outfile $species.final.domains.out -cpu 20
perl $software/PfamScan/pfam_scan.pl -fasta ../exonerate/$species.exonerate.pep -dir ~/data/pfam/ -outfile $species.exonerate.domains.out -cpu 20
grep 7tm_6 $species.final.domains.out|cut -f 1 -d " "|sort|uniq -c|wc -l





#	domains	stop	Start	total
#Acep 360	302		313		426
#Acol 396	338		385		498
#Parg 236 	208		208		282
#Acha 300 	226		275		431
#Ahey 361	300		323		424
#Ains 391	337		352		447
#Aech 420	410		379		432
#Tsep 361	305		323		453
#grep 7tm_6 *.final.domains.out|cut -f 1 -d " "|sort|uniq -c|wc -l