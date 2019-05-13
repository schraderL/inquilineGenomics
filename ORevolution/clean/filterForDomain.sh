
cd $base/$species/final
# check fullOR annotation for pfam domains
nice pfam_scan.pl -fasta $species.fullORannotation.fa -dir $pfam -outfile $species.fullORannotation.domains.out -cpu $cpu

# extract only those models that have a 7tm_6 domain
awk '{if ($7=="7tm_6") print $1}' $species.fullORannotation.domains.out|sort|uniq > $species.ORs.7tm6.lst
faSomeRecords ../final/$species.fullORannotation.fa $species.ORs.7tm6.lst $species.OR.7tm6.fa
faSomeRecords -exclude ../final/$species.fullORannotation.fa $species.ORs.7tm6.lst $species.OR.no7tm6.fa

# create gff for OR genes with 7tm domain
###################
sed s/evm.model.//g $species.ORs.7tm6.lst | grep -f - $species.fullORannotation.gff3 > tmp.gff
GFFcleaner --clean-replace-attributes  --add-missing-ids --add-exon-ids tmp.gff
gt gff3 -sort -tidy -addids -checkids yes -retainids yes -fixregionboundaries -addintrons tmp_clean.gff > $species.OR.7tm6.gff3
rm tmp*
perl $EVMtools/gff3_gene_prediction_file_validator.pl $species.OR.7tm6.gff3
