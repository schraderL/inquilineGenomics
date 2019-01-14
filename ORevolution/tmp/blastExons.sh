db=~/data/genomes/inquilines/v1.1/Ahey/Acromyrmex_heyeri.fa


cd /usr/local/home/lschrader/data/OR/singleExons
awk '{if ($3 == "exon") print $0}' ../SupplementaryFile4.AechOrs.final.gff > Aech.ORexons.gff
awk  'BEGIN {OFS=""};{ gsub(/\;Parent.*/, "")};{ gsub(/ID=/, "")};{print "samtools faidx ~/data/genomes/Aech/Aech_2.0_scaffolds.clean.fa ",$1,":",$4,"-",$5," >> AechORexons.fa"}' Aech.ORexons.gff > samtools.Aech.ORexons.sh
rm AechORexons.fa
bash samtools.Aech.ORexons.sh #AechORexons.fa
blastn -db $db -query Aech.ORexons.fa -out Ahey.ORexons.bls -num_threads 16 -evalue 1e-2 -outfmt 6 
perl $scripts/bls2gff.pl Ahey.ORexons.bls > Ahey.OR.exons.gff

