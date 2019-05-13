#!/usr/bin/bash
#############################################
#7 Prepare EVM input files
#############################################

cd $base/$species/EVM
cp $base/EVM.weights.txt weights.txt

# Partitioning
echo "Partitioning"

# Ignoring transcript data in general appears to be a good idea... At least the files created by transdecoder are not useful and overinterpret the data. The raw cufflinks transcriots could maybe be used instead.
genomeFa=$species.OR.scf.fa

  $EVMtools/partition_EVM_inputs.pl --genome $genomeFa \
       --gene_predictions $species.prot.gff \
       --protein_alignments $species.PROTEIN.gff \
       --segmentSize 100000 --overlapSize 40000 --partition_listing partitions_list.out
#      --transcript_alignments OR.transcripts.gff3 \

       echo "Writing commands"
       # --transcript_alignments removed for Acep
       $EVMtools/write_EVM_commands.pl --genome $genomeFa --weights $base/$species/EVM/weights.txt \
             --gene_predictions $species.prot.gff \
             --protein_alignments $species.PROTEIN.gff \
             --search_long_introns 1000 \
             --re_search_intergenic 4000 \
             --output_file_name evm.out  --partitions partitions_list.out >  commands.list
             #--search_long_introns 2000 \
             #--re_search_intergenic 6000 \#
             #--transcript_alignments OR.transcripts.gff3 \

# run in parallel
if [ -z ${cpu+x} ]; then echo "cpu is unset, setting it to 8."; cpu=8; else echo "cpu is set to '$cpu'"; fi
echo "Running commands with $cpu cores"

lines=$(wc -l commands.list|awk '{print int($1/'$cpu')}')

split --lines=$lines commands.list

ls xa* | parallel --no-notice $EVMtools/execute_EVM_commands.pl {} '>' run.{.}.log
rm xa*

#combine
echo "Combining gffs per scaffold"
$EVMtools/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

#write gffs
echo "writing gffs"
$EVMtools/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome $genomeFa

rm -rf tmpGff/
mkdir tmpGff
for f in */evm.out.gff3
do
	scf=$(echo $f | cut -f 1 -d "/")
	gt gff3 -sort -tidy -retainids -checkids -fixregionboundaries -addintrons $f > tmpGff/$scf.gff
done

echo "Merging gffs for all scaffolds"
gt merge -tidy -retainids tmpGff/*.gff > $species.final.OR.gff3


#write protein fasta
echo "Writing protein fasta"
cd $base/$species/EVM/
for f in */evm.out.gff3
do
if [ -s $f ]
	then
	echo "$EVMtools/gff3_file_to_proteins.pl  $f $genomeFa > $f.prot.fa" >> getFasta.sh
fi
done

wc -l getFasta.sh
lines=$(wc -l getFasta.sh|awk '{print int($1/'$cpu')}')
split --lines=$lines getFasta.sh getFasta.split.

ls getFasta.split.a* | parallel --no-notice bash {}

cat */evm.out.gff3.prot.fa > $species.final.OR.fa

cd $base/$species/EVM/
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
