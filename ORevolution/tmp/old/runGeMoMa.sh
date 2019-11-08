# This script allows to run GeMoMa with minimal work from command line. A simple example is
# ./run.sh test_data/annotation.gff test_data/reference.fasta test_data/contig.fasta
# ./run.sh /usr/local/home/lschrader/data/OR/SupplementaryFile4.AechOrs.final.gff ~/data/genomes/Aech/Aech_2.0_scaffolds.clean.fa ~/data/OR/genomes/NCBI/Tsep/GCF_001594115.1_Tsep1.0_genomic.fna results/Tsep/
# In this case, the prediction is located in ./predicted_annotation.gff

#parameters
annotation=$1
#e.g. annotation=/usr/local/home/lschrader/data/OR/SupplementaryFile4.AechOrs.final.gff
reference=$2
#e.g. reference=~/data/genomes/Aech/Aech_2.0_scaffolds.clean.fa
target=$3
#e.g. target=~/data/OR/genomes/NCBI/Tsep/GCF_001594115.1_Tsep1.0_genomic.fna
#outfolder
out=$4
#/usr/local/home/lschrader/software/GeMoMa_1.3/GeMoMa-1.3.jar
GeMoMaJar=$5

echo "============================================================================="
echo ""
echo "check versions:"
echo ""

java -version

echo ""

tblastn -version

echo ""
echo "============================================================================="
echo ""

java -jar $GeMoMaJar CLI Extractor a=${annotation} g=${reference} outdir=${out} f=false r=true Ambiguity=AMBIGUOUS
#java -jar $GeMoMaJar CLI Extractor a=${annotation} g=${reference} outdir=${out} f=true

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo "makeblastdb:"
echo ""

makeblastdb -out ${out}/blastdb -hash_index -in ${target} -title "target" -dbtype nucl

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo "tblastn:"
echo ""

#tblastn -query ${out}/cds-parts.fasta -db ${out}/blastdb -evalue 100.0 -out ${out}/tblastn.txt -outfmt "6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles" -db_gencode 1 -matrix BLOSUM62 -seg no -word_size 3 -comp_based_stats F -gapopen 11 -gapextend 1 -max_hsps 0
tblastn -query ${out}/cds-parts.fasta -db ${out}/blastdb -evalue 100.0 -out ${out}/tblastn.txt -outfmt "6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles" -db_gencode 1 -matrix BLOSUM62 -seg no -word_size 3 -comp_based_stats F -gapopen 11 -gapextend 1 -num_threads 20


echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""

#java -jar /usr/local/home/lschrader/software/GeMoMa_Official/GeMoMa-1.1.3.jar CLI GeMoMa t=${out}/tblastn.txt c=${out}/cds-parts.fasta a=${out}/assignment.tabular tg=${target}  outdir=${out} p=100
#java -jar $GeMoMaJar CLI GeMoMa t=${out}/tblastn.txt c=${out}/cds-parts.fasta a=${out}/assignment.tabular tg=${target}  outdir=${out} p=30 ct=0.1 rt=0.1 intron-loss-gain-penalty=100 m=4000
#java -jar $GeMoMaJar CLI GeMoMa t=${out}/tblastn.txt c=${out}/cds-parts.fasta a=${out}/assignment.tabular tg=${target}  outdir=${out} p=30 ct=0.1 rt=0.1  m=8000
java -jar $GeMoMaJar CLI GeMoMa t=${out}/tblastn.txt c=${out}/cds-parts.fasta a=${out}/assignment.tabular tg=${target}  outdir=${out} p=60 ct=0.1 rt=0.1 m=20000
#java -jar $GeMoMaJar CLI GeMoMa t=${out}/tblastn.txt c=${out}/cds-parts.fasta a=${out}/assignment.tabular tg=${target}  outdir=${out} p=60 ct=0.1 rt=0.1 m=100000

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""

cd ${out}
#java -jar /usr/local/home/lschrader/software/GeMoMa_1.3/GeMoMa-1.3.jar CLI GAF g=${out}/predicted_annotation.gff r=-10000 c=F
java -jar $GeMoMaJar CLI GAF g=${out}/predicted_annotation.gff r=-10000 c=F
