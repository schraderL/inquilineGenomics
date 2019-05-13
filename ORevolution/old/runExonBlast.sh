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

export AqGenome=($queryGenome1 $queryGenome2 $queryGenome3)
export AqGff=($ORgff1 $ORgff2 $ORgff3)
export qspecies=($q1 $q2 $q3)

export EVMutils=$software/EVidenceModeler-1.1.1/EvmUtils

for i in 0 1 2
do
mkdir $base/$species/singleExon/${qspecies[$i]}
cd $base/$species/singleExon/${qspecies[$i]}
export qGenome=${AqGenome[$i]}
export qORgff=${AqGff[$i]}
awk '{if ($3 == "exon") print $0}' $qORgff > ORexons.gff
awk 'BEGIN {OFS=""};{ gsub(/\;Parent.*/, "")};{ gsub(/ID=/, "")};{print "samtools faidx $qGenome ",$1,":",$4,"-",$5," |sed s/\\>.*/\\>",$9,":",$1,":",$4,"-",$5,"/g >> ORexons.fa"}' ORexons.gff > samtools.ORexons.sh
rm -f ORexons.fa
bash samtools.ORexons.sh #AechORexons.fa
#blastn -db $genomeFa -query ORexons.fa -out ORexons.bls -num_threads 16 -evalue 1e-2 -outfmt 6
tblastx -db $genomeFa -query ORexons.fa -out ORexons.bls -num_threads 16 -evalue 1e-2 -outfmt 6
#tblastn -db $genomeFa -query tmp.fa -out ORexons.bls -num_threads 16 -evalue 1 -outfmt 6
perl $scripts/bls2gff.pl ORexons.bls > $species.OR.exons.${qspecies[$i]}.bls.gff3
perl $EVMutils/gff3_gene_prediction_file_validator.pl $species.OR.exons.${qspecies[$i]}.bls.gff3
cp $species.OR.exons.${qspecies[$i]}.bls.gff3 ..
done

cd $base/$species/singleExon/
rm $species.OR.exons.bls.gff
cat *.OR.exons.*.bls.gff3  |sort -k1 -k 4 -n> $species.OR.exons.bls.gff
