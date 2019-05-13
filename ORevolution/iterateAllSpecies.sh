bash scripts/ORannotation.190405.sh Ahey /corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_heyeri.2.1/genome/ Acromyrmex_heyeri.v2.1.fa
bash scripts/ORannotation.190405.sh Acha /corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_charruanus.2.1/genome Acromyrmex_charruanus.v2.1.fa

bash scripts/ORannotation.190405.sh Aech /corefac/cse/lukas/genomes/NCBI/Aech GCF_000204515.1_Aech_3.9_genomic.fna

bash scripts/ORannotation.190405.sh Aech /corefac/cse/lukas/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/genome Acromyrmex_echinatior.v2.0.fa
bash scripts/ORannotation.190405.sh Acep /corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_cephalotes/genome Atta_cephalotes.v2.0.fa
bash scripts/ORannotation.190405.sh Acol /corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_colombica/genome Atta_colombica.v2.0.fa
bash scripts/ORannotation.190405.sh Ccos /corefac/cse/lukas/genomes/reannotations_of_7_ants/Cyphomyrmex_costatus/genome Cyphomyrmex_costatus.v2.0.fa
bash scripts/ORannotation.190405.sh Tcor /corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_cornetzi/genome Trachymyrmex_cornetzi.v2.0.fa
bash scripts/ORannotation.190405.sh Tsep /corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_septentrionalis/genome Trachymyrmex_septentrionalis.v2.0.fa
bash scripts/ORannotation.190405.sh Tzet /corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_zeteki/genome Trachymyrmex_zeteki.v2.0.fa

bash scripts/ORannotation.190405.sh Tzet2 /corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_zeteki/genome Trachymyrmex_zeteki.v2.0.fa


database=/corefac/cse/lukas/genomes/reannotations_of_7_ants
for i in Acromyrmex_echinatior/genome/Acromyrmex_echinatior.v2.0.fa Atta_cephalotes/genome/Atta_cephalotes.v2.0.fa Atta_colombica/genome/Atta_colombica.v2.0.fa Cyphomyrmex_costatus/genome/Cyphomyrmex_costatus.v2.0.fa Trachymyrmex_cornetzi/genome/Trachymyrmex_cornetzi.v2.0.fa Trachymyrmex_septentrionalis/genome/Trachymyrmex_septentrionalis.v2.0.fa Trachymyrmex_zeteki/genome/Trachymyrmex_zeteki.v2.0.fa
do
 species=$(echo $i |cut -f 1 -d "/"|perl -pe 's/(.).*_(...).*/$1$2/g')
 echo "bash scripts/ORannotation.190405.sh $species $database/$i"|sed 's|genome/|genome/ |g'
done > runGenomes.sh
