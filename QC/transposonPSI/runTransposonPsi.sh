#################################################
# transposonPSI
#################################################

#############################
# Run transposon psi to identify potential TE-derived or TE-contaminated genes.

#edited the TransposonPSI.pl  to use 20 threads in blastall (-a 20). Uses blastgpg though unless nuc database
base=/usr/local/home/lschrader/data/inqGen18/QC/TransposonPSI/
cd $base/data

# pep.fa file for new genomes
find /usr/local/home/lschrader/data/genomes/inquilines_2.0/  -name *.pep.fa > pepFiles.new.sh
# pep.fa file for attine genomes
find /usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/  -name *.pep.fa > pepFiles.attines.sh
#############################


#############################
# run TransposonPSI on each genome in parallel
cd $base/output
ls $base/data/*.pep.fa|parallel --nice 10 "perl ~/software/TransposonPSI_08222010/transposonPSI.pl {} prot 1>{}.PSI.out 2> {}.blast.err"
############################

#############################
# transition from 2.0 to 2.1
mkdir $base/inquilines_2.1
cd $base/inquilines_2.1

genomes=/usr/local/home/lschrader/data/inqGen18/inquilines_v2.1
ACHA=$genomes/Acromyrmex_charruanus.2.1/annotation/gene_annotation/Acromyrmex_charruanus.v2.1.pep.fa
AHEY=$genomes/Acromyrmex_heyeri.2.1/annotation/gene_annotation/Acromyrmex_heyeri.v2.1.pep.fa
AINS=$genomes/Acromyrmex_insinuator.2.1/annotation/gene_annotation/Acromyrmex_insinuator.v2.1.pep.fa
PARG=$genomes/Pseudoatta_argentina.2.1/annotation/gene_annotation/Pseudoatta_argentina.v2.1.pep.fa


#Acha
grep ">" $ACHA|cut -f 1 -d "-"|cut -f 2 -d ">" > Acha.v2.1.prot.lst
cat $base/output/Acha.pep.fa.TPSI.allHits|perl -pe 's/Acha_([0-9]{5})/ACHA$1/g' |awk -F'\t' 'NR==FNR{c[$1]++;next};c[$6]' Acha.v2.1.prot.lst - > Acha.v2.1.pep.fa.TPSI.allHits
cat $base/output/Acha.pep.fa.TPSI.topHits|perl -pe 's/Acha_([0-9]{5})/ACHA$1/g' |awk -F'\t' 'NR==FNR{c[$1]++;next};c[$6]' Acha.v2.1.prot.lst - > Acha.v2.1.pep.fa.TPSI.topHits

#Ahey
grep ">" $AHEY|cut -f 1 -d "-"|cut -f 2 -d ">" > Ahey.v2.1.prot.lst
cat $base/output/Ahey.pep.fa.TPSI.allHits|perl -pe 's/Ahey_([0-9]{5})/AHEY$1/g' |awk -F'\t' 'NR==FNR{c[$1]++;next};c[$6]' Ahey.v2.1.prot.lst - > Ahey.v2.1.pep.fa.TPSI.allHits
cat $base/output/Ahey.pep.fa.TPSI.topHits|perl -pe 's/Ahey_([0-9]{5})/AHEY$1/g' |awk -F'\t' 'NR==FNR{c[$1]++;next};c[$6]' Ahey.v2.1.prot.lst - > Ahey.v2.1.pep.fa.TPSI.topHits

#Ains
grep ">" $AINS|cut -f 1 -d "-"|cut -f 2 -d ">" > Ains.v2.1.prot.lst
cat $base/output/Ains.pep.fa.TPSI.allHits|perl -pe 's/Ains_([0-9]{5})/AINS$1/g' |awk -F'\t' 'NR==FNR{c[$1]++;next};c[$6]' Ains.v2.1.prot.lst - > Ains.v2.1.pep.fa.TPSI.allHits
cat $base/output/Ains.pep.fa.TPSI.topHits|perl -pe 's/Ains_([0-9]{5})/AINS$1/g' |awk -F'\t' 'NR==FNR{c[$1]++;next};c[$6]' Ains.v2.1.prot.lst - > Ains.v2.1.pep.fa.TPSI.topHits

#Parg
grep ">" $PARG|cut -f 1 -d "-"|cut -f 2 -d ">" > Parg.v2.1.prot.lst
cat $base/output/Parg.pep.fa.TPSI.allHits|perl -pe 's/Parg_([0-9]{5})/PARG$1/g' |awk -F'\t' 'NR==FNR{c[$1]++;next};c[$6]' Parg.v2.1.prot.lst - > Parg.v2.1.pep.fa.TPSI.allHits
cat $base/output/Parg.pep.fa.TPSI.topHits|perl -pe 's/Parg_([0-9]{5})/PARG$1/g' |awk -F'\t' 'NR==FNR{c[$1]++;next};c[$6]' Parg.v2.1.prot.lst - > Parg.v2.1.pep.fa.TPSI.topHits

rm *.prot.lst
############################

cd $base/output
ls *pep.fa.TPSI.topHits| parallel --nice 10 'echo {}|cut -f 1 -d "."' > ../species.tsv
ls *pep.fa.TPSI.topHits| parallel --nice 10 'egrep "^//" {} -v|cut -f 6|wc -l' > ../topHits.tsv
paste ../species.tsv ../topHits.tsv > ../transposonPsiResults.out
