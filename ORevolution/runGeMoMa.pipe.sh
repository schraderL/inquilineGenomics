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

export EVMutils=$software/EVidenceModeler-1.1.1/EvmUtils
export gtools=$gtools
####################################################################################################
# Prepare GAF  gff
####################################################################################################
mkdir $base/$species/GeMoMa/abinitio
cd $base/$species/GeMoMa/abinitio
java -jar $GeMoMaProgram CLI GAF g=../$q1/predicted_annotation.gff g=../$q2/predicted_annotation.gff g=../$q3/predicted_annotation.gff c=false e=0.1 r=0 e=0 m=true cbf=0
#java -jar $GeMoMaProgram CLI GAF g=../$q1/predicted_annotation.gff g=../$q2/predicted_annotation.gff g=../$q3/predicted_annotation.gff c=false e=0.1 r=-1000 e=0 m=true cbf=0
grep -v -E ".*GAF.*" ./filtered_predictions.gff |sed s/prediction/mRNA/g> tmp.gff
gffread tmp.gff -o tmp.gff3 --force-exons -E
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
rm tmp.gff3
$gtools/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
rm tmp_clean.gff
sed s/';geneID=.*'//g tmp2.gff| perl $scripts/addGeneGFF.pl - | sed s/GeMoMa/GAF/g > GAF.final.gff
perl $EVMutils/gff3_gene_prediction_file_validator.pl GAF.final.gff

#cp $base/$species/GeMoMa/abinitio/GAF.final.gff $base/$species/EVM
####################################################################################################

####################################################################################################
# Prepare OTHER  gff
####################################################################################################
mkdir $base/$species/GeMoMa/OTHER
cd $base/$species/GeMoMa/OTHER

java -jar $GeMoMaProgram CLI GAF g=../$q1/predicted_annotation.gff g=../$q2/predicted_annotation.gff g=../$q3/predicted_annotation.gff c=true

mv filtered_predictions.gff filtered_predictions_complete.gff
grep -v -E ".*GAF.*" filtered_predictions_complete.gff|sed s/prediction/mRNA/g> tmp.gff
gffread tmp.gff -o tmp.gff3 --force-exons -E
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
rm tmp.gff3
$gtools/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
rm tmp_clean.gff
sed s/';geneID=.*'//g tmp2.gff| perl $scripts/addGeneGFF.pl - | sed s/GeMoMa/complete/g > OTHER.final.gff
# could also test ~/.local/bin/PASApipeline-2.0.2/misc_utilities/gtf_to_gff3_format.pl
perl $EVMutils/gff3_gene_prediction_file_validator.pl OTHER.final.gff
####################################################################################################


####################################################################################################
# Prepare UNIFIED  gff
####################################################################################################
mkdir $base/$species/GeMoMa/OTHER
cd $base/$species/GeMoMa/OTHER

java -jar $GeMoMaProgram CLI GAF g=../$q1/predicted_annotation.gff g=../$q2/predicted_annotation.gff g=../$q3/predicted_annotation.gff r=0 c=false
mv filtered_predictions.gff filtered_predictions_unified.gff
grep -v -E ".*GAF.*" filtered_predictions_unified.gff|sed s/prediction/mRNA/g> tmp.gff
gffread tmp.gff -o tmp.gff3 --force-exons -E
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
rm tmp.gff3
$gtools/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
rm tmp_clean.gff
sed s/';geneID=.*'//g tmp2.gff| perl $scripts/addGeneGFF.pl - | sed s/GeMoMa/complete/g > unified.final.gff
# could also test ~/.local/bin/PASApipeline-2.0.2/misc_utilities/gtf_to_gff3_format.pl
perl $EVMutils/gff3_gene_prediction_file_validator.pl unified.final.gff
####################################################################################################


cd $base/$species/GeMoMa
for i in $q1 $q2 $q3
do
cd $base/$species/GeMoMa/$i
grep -v -E ".*GAF.*" filtered_predictions.gff |sed s/prediction/mRNA/g> tmp.gff
gffread tmp.gff -o tmp.gff3 --force-exons -E
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
rm tmp.gff3
$gtools/gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > tmp2.gff
rm tmp_clean.gff

sed s/';geneID=.*'//g tmp2.gff| perl $scripts/addGeneGFF.pl - | sed s/GeMoMa/GeMoMa_$i/g > $species.vs.$i.GeMoMa.final.gff
# could also test ~/.local/bin/PASApipeline-2.0.2/misc_utilities/gtf_to_gff3_format.pl
perl $EVMutils/gff3_gene_prediction_file_validator.pl $species.vs.$i.GeMoMa.final.gff
cp $species.vs.$i.GeMoMa.final.gff $base/$species/GeMoMa
done

cd $base/$species/GeMoMa
$gtools/gt merge -tidy -retainids *.vs.*.GeMoMa.final.gff  > $species.GeMoMa.final.gff

# export base=/usr/local/home/lschrader/data/OR
# export scripts=$base/scripts
# export ORfa=$base/AechORs.McKenzie.fasta
# export ORgff=$base/SupplementaryFile4.AechOrs.final.gff
# export queryGenome=/usr/local/home/lschrader/data/genomes/Aech/Aech_2.0_scaffolds.clean.fa
# export software=/usr/local/home/lschrader/software
# export cpu=8
# export target=/usr/local/home/lschrader/data/OR/Acep/GeMoMa2
# export target=/usr/local/home/lschrader/data/OR/Acep/GeMoMa3
# export species=Acep
# export genomeFolder=/usr/local/home/lschrader/data/genomes/attines/assembly
# export genomeFa=$genomeFolder"/Acep.genome.fa"
#
# #java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$scf.fa  outdir=$target p=100 i=50 m=1500 timeout=1000
# #java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=100 i=50 ct=0.2
# java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=100 ct=0.2
# java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=100 ct=0.1 r=0.4 i=75
# java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=10 ct=0.05
# java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=10 ct=0.01 r=0.4 i=75
# java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=10 ct=0.01 r=0.4 i=100 m=5000
# java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=10 ct=0.1 r=0.4 i=100 m=5000
#
# # pretty good but produces too many models
# java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=100 ct=0.1 r=0.1 i=75 m=6000 #Acep.GeMoMa.ct0.1.p100.r0.1.i75.m4000.gff
#java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=100 ct=0.1 r=0.1 i=75 m=6000 #Acep.GeMoMa.ct0.1.p100.r0.1.i75.m4000.gff

#GeMoMa might have a problem with tandem duplicated genes, as it could use exons of both genes to make one prediction. In this case you could play around with the region threshold “rt”, formerly “r”.

#java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=$target/tblastn.txt c=$target/cds-parts.fasta a=$target/assignment.tabular tg=$genomeFa  outdir=$target p=50 ct=0.1 r=0.9 i=75 m=6000

#Acep.GeMoMa.ct0.1.p100.r0.1.i75.m4000.gff


#If you like to have more predictions you should set the parameter contig threshold (ct) and predictions (p).
#For each reference transcript/CDS, GeMoMa initially makes simple prediction and uses this prediction to determine
#whether a contig is promising and should be used to determine the final predictions. You have to decrease ct and
#increase p to have more contigs in the final prediction. Increasing the number of predictions allows GeMoMa to
#output more predictions that have been computed. Decreasing the contig threshold allows to increase the number
#of predictions that are (internally) computed. Increasing p to a very large number without decreasing ct does not help.
#We are using ct=0.4 and p=10 for our studies right now.
