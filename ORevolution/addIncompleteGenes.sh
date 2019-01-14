############################################################
# 4. Add incomplete genes from GeMoMa Aech mapping
############################################################

cd $base/$species/final/
#cp $base/$species/GeMoMa/$species.vs.Aech.GeMoMa.final.gff .
cp $base/$species/GeMoMa/$species.vs.$q1.GeMoMa.final.gff .
# extract those models that are not overlapping with proper ORs
bedtools intersect -b $species.final.OR.gff3 -a $species.vs.$q1.GeMoMa.final.gff -wb -v|awk '{if ($3=="gene") print $0}'|perl -pe 's/.*ID=(.*?)(\n|;|\t|,)/$1,/g'|sed -e 's/,$//' > 2add.ids.lst
rm $species.vs.$q1.GeMoMa.final.gff.db
gffutils-cli create $species.vs.$q1.GeMoMa.final.gff

# add gff entries for the non-overlapping models
sed -e s/^/'gffutils-cli children $species.vs.$q1.GeMoMa.final.gff.db '/g 2add.ids.lst >2add.sh
bash 2add.sh > 2add.dirty.gff3
awk '!a[$0]++' 2add.dirty.gff3 > 2add.gff3

# clean gff
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= 2add.gff3
$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -addintrons 2add_clean.gff > 2add.tmp.gff
$software/genometools-1.5.9/bin/gt merge -tidy -retainids 2add.tmp.gff   $species.final.OR.gff3 > $species.fullORannotation.gff3
rm 2add*
perl $EVMtools/gff3_gene_prediction_file_validator.pl $species.fullORannotation.gff3

echo "Writing protein fasta"
cd $base/$species/final/
$EVMtools/gff3_file_to_proteins.pl   $species.fullORannotation.gff3 $genomeFa > $species.fullORannotation.fa

# gff & fa containing full and partial models
# GFF:    $base/$species/final/$species.fullORannotation.gff3
# Fasta:  $base/$species/final/$species.fullORannotation.fa
