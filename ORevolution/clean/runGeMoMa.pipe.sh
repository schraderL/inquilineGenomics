###################################################
#==================================================
# 2. Annotate OR genes with GeMoMa
#==================================================

###################################################
#1) run GeMoMa
###################################################

#mkdir $base/$species/GeMoMa/
cd $base/$species/GeMoMa/
cp $scripts/runGeMoMa.sh .

# GeMoMa filters all genes with premature stop codons.
# From ~/data/OR/Acep/GeMoMa/Aech/protocol_Extractor.txt
# reasons for discarding transcripts:
# ambigious nucleotide    5
# start phase not zero    0
# missing start   0
# missing stop    0
# premature stop  17
# no DNA  0
# wrong phase     1
# conflicting phase       0

bash ./runGeMoMa.sh $ORgff1 $queryGenome1 $genomeFa $base/$species/GeMoMa/$q1 $GeMoMaProgram
bash ./runGeMoMa.sh $ORgff2 $queryGenome2 $genomeFa $base/$species/GeMoMa/$q2 $GeMoMaProgram
bash ./runGeMoMa.sh $ORgff3 $queryGenome3 $genomeFa $base/$species/GeMoMa/$q3 $GeMoMaProgram

####################################################################################################
# Prepare GAF gff
####################################################################################################
mkdir $base/$species/GeMoMa/abinitio
cd $base/$species/GeMoMa/abinitio
rm filtered_predictions*.gff
rm protocol_GAF*
java -jar $GeMoMaProgram CLI GAF g=../$q1/predicted_annotation.gff g=../$q2/predicted_annotation.gff g=../$q3/predicted_annotation.gff c=false r=0 m=true cbf=0
grep -v -E ".*GAF.*" ./filtered_predictions.gff |sed s/prediction/mRNA/g> tmp.gff
gffread tmp.gff -o tmp.gff3 --force-exons -E
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
rm tmp.gff3
gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
rm tmp_clean.gff
sed s/';geneID=.*'//g tmp2.gff| perl $scripts/addGeneGFF.pl - | sed s/GeMoMa/GAF/g > GAF.final.gff
perl $EVMtools/gff3_gene_prediction_file_validator.pl GAF.final.gff
####################################################################################################

####################################################################################################
# Prepare OTHER  gff
####################################################################################################
mkdir $base/$species/GeMoMa/OTHER
cd $base/$species/GeMoMa/OTHER
rm filtered_predictions*.gff
rm protocol_GAF*
java -jar $GeMoMaProgram CLI GAF g=../$q1/predicted_annotation.gff g=../$q2/predicted_annotation.gff g=../$q3/predicted_annotation.gff c=true
mv filtered_predictions.gff filtered_predictions_complete.gff
grep -v -E ".*GAF.*" filtered_predictions_complete.gff|sed s/prediction/mRNA/g> tmp.gff
gffread tmp.gff -o tmp.gff3 --force-exons -E
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
rm tmp.gff3
gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
rm tmp_clean.gff
sed s/';geneID=.*'//g tmp2.gff| perl $scripts/addGeneGFF.pl - | sed s/GeMoMa/complete/g > OTHER.final.gff
perl $EVMtools/gff3_gene_prediction_file_validator.pl OTHER.final.gff
####################################################################################################


####################################################################################################
# Prepare UNIFIED  gff
####################################################################################################
mkdir $base/$species/GeMoMa/OTHER
cd $base/$species/GeMoMa/OTHER
rm filtered_predictions*.gff
rm protocol_GAF*
java -jar $GeMoMaProgram CLI GAF g=../$q1/predicted_annotation.gff g=../$q2/predicted_annotation.gff g=../$q3/predicted_annotation.gff r=0 c=false outdir=unified
mv unified/filtered_predictions.gff filtered_predictions_unified.gff
grep -v -E ".*GAF.*" filtered_predictions_unified.gff|sed s/prediction/mRNA/g> tmp.gff
gffread tmp.gff -o tmp.gff3 --force-exons -E
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
rm tmp.gff3
gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
rm tmp_clean.gff
sed s/';geneID=.*'//g tmp2.gff| perl $scripts/addGeneGFF.pl - | sed s/GeMoMa/complete/g > unified.final.gff
perl $EVMtools/gff3_gene_prediction_file_validator.pl unified.final.gff
####################################################################################################

cd $base/$species/GeMoMa
for i in $q1 $q2 $q3
do
  cd $base/$species/GeMoMa/$i
  grep -v -E ".*GAF.*" filtered_predictions.gff |sed s/prediction/mRNA/g> tmp.gff
  gffread tmp.gff -o tmp.gff3 --force-exons -E
  GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
  rm tmp.gff3
  gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
  rm tmp_clean.gff
  sed s/';geneID=.*'//g tmp2.gff| perl $scripts/addGeneGFF.pl - | sed s/GeMoMa/GeMoMa_$i/g > $species.vs.$i.GeMoMa.final.gff
  perl $EVMtools/gff3_gene_prediction_file_validator.pl $species.vs.$i.GeMoMa.final.gff
  cp $species.vs.$i.GeMoMa.final.gff $base/$species/GeMoMa
done

cd $base/$species/GeMoMa
gt merge -tidy -retainids *.vs.*.GeMoMa.final.gff  > $species.GeMoMa.final.gff
