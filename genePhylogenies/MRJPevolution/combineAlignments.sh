# Combine Ant and Amel alignments
base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/MRJPs/aln/mergeAln
cd $base
prank -d1=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/MRJPs/aln/output.best.fas -d2=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/MRJPs/Amel.reference/aln/output.best.fas
 cat output.fas |perl -pe 's/(>.*?-.*?-.*?)-.*/$1/g'|perl -pe 's/.*gene=\"(.*?)\";.*/>$1/g' > output.aln.fas
FastTree output.aln.fas > output.aln.tre

# RaxML on aa aln
nice raxmlHPC-PTHREADS -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE -soutput.aln.fas -n r1 -T 20
