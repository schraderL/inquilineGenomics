#https://github.com/alexharkess/PBIO_6550/wiki/Lab-4-:-Alignment-and-Phylogeny

base=~/data/inqGen18/ORphylogeny/
cd $base
mkdir data
cd $base/data

for i in $(readlink -f ../../ORs/*/final/finalSet/*.fa)
do
ln -s $i .
done

rm Aaem.* Agra.* Alob* Pcal* Plob*

cat *.OR.fa > allORs.fa
cd $base
unset MAFFT_BINARIES
nice prank -d=$base/data/allORs.fa
cp ../results/tree/output.best.fas prank.aln.fa


~/software/trimal/source/trimal -in prank.aln.fa  -out prank.aln.trimmed.fa  -strict -fasta
#~/software/alan/alan ~/data/inqGen18/ORevo/phylogeny/prank.aln.trimmed.fa
~/software/FastTree/FastTreeMP $base/prank.aln.trimmed.fa > $base/prank.aln.trimmed.fa
