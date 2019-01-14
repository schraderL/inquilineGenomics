
#############################################
# ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! 
#	1 reformat GeMoMa gff to pass EVM specs
# ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! 
#############################################

#############################################
#1 reformat GeMoMa gff to pass EVM specs
#############################################

species=Parg
genomeFolder=/usr/local/home/lschrader/data/genomes/inquilines/v1.1/Parg/
genomeFa=$genomeFolder"Pseudoatta_argentina.fa"
transcript=~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf


cd ~/data/OR/$species/GeMoMa
gffread predicted_annotation.gff -o tmp.gff3 --force-exons -E 
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3 
 ~/software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
perl ../../scripts/addGeneGFF.pl tmp2.gff > $species.GeMoMa.final.gff
# could also test ~/.local/bin/PASApipeline-2.0.2/misc_utilities/gtf_to_gff3_format.pl
perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.GeMoMa.final.gff

#############################################
#2 reformat exonerate gff to pass EVM specs
#############################################

cd ~/data/OR/$species/exonerate/gff
cat * > ../$species.exonerate.gff
cd ..


gffread $species.exonerate.gff -o tmp.gff3 --force-exons -E 
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3 
 ~/software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff


perl ../../scripts/addmRNAGFF.pl tmp2.gff > $species.exonerate.final.gff

perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.exonerate.final.gff


#############################################
#3 reformat cufflinks gtf to pass EVM specs
#############################################


#http://transdecoder.github.io/
mkdir ~/data/OR/$species/EVM/
cd ~/data/OR/$species/EVM/
~/software/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.fa > Parg.cufflinks.fa
~/software/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf >  Parg.transcript.gff
TransDecoder.LongOrfs -t Parg.cufflinks.fa 
 ~/software/TransDecoder-3.0.1/util/cdna_alignment_orf_to_genome_orf.pl Parg.cufflinks.fa.transdecoder_dir/longest_orfs.gff3  Parg.transcript.gff Parg.cufflinks.fa > transcripts.fasta.transdecoder.genome.gff3
perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl transcripts.fasta.transdecoder.genome.gff3

#############################################
#4 get relevant transcripts for OR genes
#############################################

# get intersecting cufflinks transcripts with Parg.prot.sort.gff 
bedtools intersect -wa -a transcripts.fasta.transdecoder.genome.gff3 -b $species.prot.sorted.gff |grep gene|uniq |grep -P "ID=.+?;" -o|sed -r s/'ID=|\;'//g > OR.transcripts.lst
grep -f OR.transcripts.lst transcripts.fasta.transdecoder.genome.gff3 > OR.transcripts.gff3

perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl OR.transcripts.gff3

#############################################
#5 Filter Fasta for relevant scaffolds
#############################################
cat $species.exonerate.final.gff $species.GeMoMa.final.gff |cut -f 1 |grep -P "^#" -v|sort|uniq > OR.scf.lst
faSomeRecords $genomeFa OR.scf.lst $species.OR.scf.fa
genomeFa=$species.OR.scf.fa



#############################################
#6 Prepare EVM input files
#############################################

# need to concatenate the gfffiles for exonerate and GeMoMa
# 
cd ~/data/OR/$species/EVM
cat $species.exonerate.final.gff $species.GeMoMa.final.gff > $species.prot.gff
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= $species.prot.gff
~/software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes $species.prot_clean.gff > $species.prot.sorted.gff


#weights.txt
#PROTEIN	exonerate_protein2genome_local	1
#PROTEIN	GeMoMa	5
#TRANSCRIPT	transdecoder	10


#~/software/EVidenceModeler-1.1.1/evidence_modeler.pl --genome $genomeFa \
 #                      --weights ./weights.txt \
 #                      --gene_predictions $species.prot.sorted.gff \
 #                      --transcript transcripts.fasta.transdecoder.genome.gff3 \
 #                    > evm.out 
                     
                     
#~/software/EVidenceModeler-1.1.1/EvmUtils/EVM_to_GFF3.pl evm.out scaffold77 > evm.out.tmp.gff3

# Partitioning  
~/software/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome $genomeFa \
     --gene_predictions $species.prot.sorted.gff \
     --transcript_alignments OR.transcripts.gff3 \
     --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out  

# write commands
~/software/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome $genomeFa --weights /usr/local/home/lschrader/data/OR/Parg/EVM/weights.txt \
      --gene_predictions $species.prot.sorted.gff \
      --transcript_alignments OR.transcripts.gff3 \
      --output_file_name evm.out  --partitions partitions_list.out >  commands.list  


#run
~/software/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log

# run in parallel
split --lines=176 commands.list
for f in xa*
do
 ~/software/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl $f > run.$f.log &
done
#rm xa*

#combine
~/software/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

#write gffs
~/software/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome $genomeFa


#write protein fasta

for f in */evm.out.gff3
do
~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl  $f $genomeFa > $f.prot.fa &
done

cat */evm.out.gff3.prot.fa > Parg.final.OR.fa

  
#############################################
#JUNK
############################################# 


https://sourceforge.net/p/pasa/mailman/message/32505345/

cd ~/data/OR/Parg
gffread $transcript -o tmp.transcript.gff3 --force-exons -E -L 
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.transcript.gff3 
 ~/software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp.transcript_clean.gff > tmp2.transcript.gff

perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl tmp.transcript.gff3
perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl tmp2.transcript.gff
perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl tmp.transcript_clean.gff
perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf 


https://groups.google.com/forum/#!topic/trinityrnaseq-users/Nc8t1ZxtuMI

#verify transcripts
~/software/PASA_Lite-0.1.0/PASA.alignmentValidator --genome ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.fa ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf 

#assemble alignments
~/software/PASA_Lite-0.1.0/PASA.alignmentAssembler pasa_lite.valid_alignments.gtf



perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl tmp3.gff


 pasa_lite.valid_alignments.gtf  > pasa_lite.valid_alignments.gff
~/.local/bin/PASApipeline-2.0.2/misc_utilities/gtf_to_gff3_format.pl ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf  > Parg.transcripts.gff 


~/software/PASA_Lite-0.1.0/PASA.alignmentAssembler pasa_lite.valid_alignments.gtf 
perl ~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl pasa_lite.pasa_assembled_alignments.gtf


#http://transdecoder.github.io/
~/data/OR/Parg/EVM/
~/software/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.fa > Parg.cufflinks.fa
~/software/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf >  Parg.transcript.gff
TransDecoder.LongOrfs -t Parg.cufflinks.fa 
 ~/software/TransDecoder-3.0.1/util/cdna_alignment_orf_to_genome_orf.pl Parg.cufflinks.fa.transdecoder_dir/longest_orfs.gff3  Parg.transcript.gff Parg.cufflinks.fa > transcripts.fasta.transdecoder.genome.gff3
