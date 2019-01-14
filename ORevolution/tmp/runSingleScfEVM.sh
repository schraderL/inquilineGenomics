#mv EVM EVM1
#mkdir EVM
#cd EVM
#cp ../EVM1/*.*  .
cat $species.exonerate.final.gff $species.GeMoMa.final.gff $species.OR.exons.bls.gff |cut -f 1 |grep -P "^#" -v|sort|uniq|awk '{if ($0 == "scaffold02788") print $0}' > OR.scf.lst

faSomeRecords $genomeFolder"/"$species".genome.fa" OR.scf.lst $species.OR.scf.fa
genomeFa=$species.OR.scf.fa


awk '{if ($1 == "scaffold02788") print $0}' $species.prot.gff > tmp.prot.gff
awk '{if ($1 == "scaffold02788") print $0}' $species.PROTEIN.gff > tmp.PROTEIN.gff

$software/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome $genomeFa \
     --gene_predictions tmp.prot.gff \
     --protein_alignments tmp.PROTEIN.gff \
     --segmentSize 100000 --overlapSize 40000 --partition_listing partitions_list.out
#      --transcript_alignments OR.transcripts.gff3 \

 awk '{if ($1=="scaffold02788") print $0}' partitions_list.out > part.tmp

$software/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome $genomeFa --weights $base/$species/EVM/weights.edit.txt \
           --gene_predictions tmp.prot.gff \
           --protein_alignments tmp.PROTEIN.gff \
           --search_long_introns 1000 \
           --re_search_intergenic 4000 \
           --output_file_name evm.out  --partitions part.tmp >  commands.scaffold02788.list
           #--search_long_introns 2000 \
           #--re_search_intergenic 6000 \#
           #--transcript_alignments OR.transcripts.gff3 \

$software/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.scaffold02788.list > tmp.evm.out
$software/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions part.tmp --output_file_name evm.out
$software/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome $genomeFa
mv scaffold02788/evm.out.gff3 scaffold02788/scaffold02788.gff3


# ABINITIO_PREDICTION     exonerate:protein2genome:local  10
# ABINITIO_PREDICTION     GAF     1
# OTHER_PREDICTION        complete        5
# OTHER_PREDICTION        GeMoMa_Aech     15
# OTHER_PREDICTION        GeMoMa_Acep     4
# OTHER_PREDICTION        GeMoMa_Sinv     1
# TRANSCRIPT      transdecoder    2
# PROTEIN blast   3
# PROTEIN exonerate       5

ABINITIO_PREDICTION     exonerate:protein2genome:local  10
ABINITIO_PREDICTION     GAF     1
OTHER_PREDICTION        complete        5
OTHER_PREDICTION        GeMoMa_Aech     40
OTHER_PREDICTION        GeMoMa_Acep     4
OTHER_PREDICTION        GeMoMa_Sinv     1
TRANSCRIPT      transdecoder    2
PROTEIN blast   3
PROTEIN exonerate       20
