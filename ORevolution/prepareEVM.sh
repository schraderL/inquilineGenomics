#############################################
#===========================================
#	4. COMBINE EVIDENCE WITH EVM
#===========================================
#############################################
#mkdir $base/$species/EVM/
cd $base/$species/EVM/
#############################################
#1 reformat GeMoMa gff to pass EVM specs
#############################################
echo "preparing GeMoMa gff"
cp $base/$species/GeMoMa/$species.GeMoMa.final.gff $base/$species/EVM

#############################################
#1b reformat GeMoMa GAF gff to pass EVM specs
#############################################
echo "preparing GeMoMa GAF gff"
cp $base/$species/GeMoMa/abinitio/GAF.final.gff $base/$species/EVM

#############################################
#1c reformat GeMoMa OTHER gff to pass EVM specs
#############################################
echo "preparing GeMoMa OTHER gff"
cp $base/$species/GeMoMa/OTHER/OTHER.final.gff $base/$species/EVM


#############################################
#2 reformat exonerate gff to pass EVM specs
#############################################
echo "preparing exonerate gff"
cp $base/$species/exonerate/$species.exonerate.gff $base/$species/EVM/$species.exonerate.final.gff

#############################################
#2b reformat exonerate PROTEIN mappings
#############################################
echo "preparing exonerate gff"
cp $base/$species/exonerate/$species.exonerate.PROTEIN.gff $base/$species/EVM/$species.exonerate.PROTEIN.gff

#############################################
#3 get single exon blast
#############################################
echo "preparing single exon gff"
cp $base/$species/singleExon/$species.OR.exons.bls.gff $base/$species/EVM/

#############################################
#4 concatenate gffs to $species.evid.gff
#############################################
echo "Concatenating gffs for ABINITIO and OTHER"
# need to concatenate the gfffiles for exonerate and GeMoMa
cd $base/$species/EVM
$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes $species.OR.exons.bls.gff > $species.OR.exons.bls.sorted.gff
$software/genometools-1.5.9/bin/gt merge -tidy -retainids $species.exonerate.final.gff $species.GeMoMa.final.gff  $species.OR.exons.bls.sorted.gff > $species.evid.gff
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.evid.gff

$software/genometools-1.5.9/bin/gt merge -tidy -retainids $species.exonerate.final.gff $species.GeMoMa.final.gff > $species.1.gff
$software/genometools-1.5.9/bin/gt merge -tidy -retainids $species.1.gff GAF.final.gff > $species.2.gff
$software/genometools-1.5.9/bin/gt merge -tidy -retainids $species.2.gff OTHER.final.gff > $species.prot.gff

perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.prot.gff

echo "Concatenating gffs for PROTEIN"
# need to concatenate the gffs of single exon blast and exonerate
cat $base/$species/EVM/$species.OR.exons.bls.gff $base/$species/EVM/$species.exonerate.PROTEIN.gff > $base/$species/EVM/$species.PROTEIN.gff
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.PROTEIN.gff


#############################################
#5 reformat cufflinks gtf to pass EVM specs
#############################################
#http://transdecoder.github.io/
echo "Extract relevant transcripts"
cd $base/$species/EVM/

if [ -z "$transcript" ]
then
  echo "No transcript gtf file provided. Skipping transcript selection";
else
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
fi
#############################################
#6 get relevant transcripts for OR genes
#############################################
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
echo "Extract relevant scaffolds"
cat $species.exonerate.final.gff $species.GeMoMa.final.gff $species.OR.exons.bls.gff |cut -f 1 |grep -P "^#" -v|sort|uniq > OR.scf.lst
faSomeRecords $genomeFa OR.scf.lst $species.OR.scf.fa
genomeFa=$species.OR.scf.fa
