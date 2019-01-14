base=/corefac/cse/lukas/inqGen18/QC/quast

cd $base

genomeBase=~/data/inqGen18/inquilines_v2.1/
ACHAG=$genomeBase/Acromyrmex_charruanus.2.1/genome/Acromyrmex_charruanus.v2.1.fa
AHEYG=$genomeBase/Acromyrmex_heyeri.2.1/genome/Acromyrmex_heyeri.v2.1.fa
AINSG=$genomeBase/Acromyrmex_insinuator.2.1/genome/Acromyrmex_insinuator.v2.1.fa
PARGG=$genomeBase/Pseudoatta_argentina.2.1/genome/Pseudoatta_argentina.v2.1.fa

ACHA=$genomeBase/Acromyrmex_charruanus.2.1/annotation/gene_annotation/Acromyrmex_charruanus.v2.1.gff3
AHEY=$genomeBase/Acromyrmex_heyeri.2.1/annotation/gene_annotation/Acromyrmex_heyeri.v2.1.gff3
AINS=$genomeBase/Acromyrmex_insinuator.2.1/annotation/gene_annotation/Acromyrmex_insinuator.v2.1.gff3
PARG=$genomeBase/Pseudoatta_argentina.2.1/annotation/gene_annotation/Pseudoatta_argentina.v2.1.gff3

nice quast.py $ACHAG -g $ACHA -o Acha.quast.out --split-scaffolds --eukaryote
nice quast.py $PARGG -g $PARG -o Parg.quast.out --split-scaffolds --eukaryote
nice quast.py $AINSG -g $AINS -o Ains.quast.out --split-scaffolds --eukaryote
nice quast.py $AHEYG -g $AHEY -o Ahey.quast.out --split-scaffolds --eukaryote


AECHG=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/genome/Acromyrmex_echinatior.v2.0.fa
ACEPG=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_cephalotes/genome/Atta_cephalotes.v2.0.fa
ACOLG=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_colombica/genome/Atta_colombica.v2.0.fa
CCOSG=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Cyphomyrmex_costatus/genome/Cyphomyrmex_costatus.v2.0.fa
TCORG=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_cornetzi/genome/Trachymyrmex_cornetzi.v2.0.fa
TSEPG=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_septentrionalis/genome/Trachymyrmex_septentrionalis.v2.0.fa
TZETG=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_zeteki/genome/Trachymyrmex_zeteki.v2.0.fa

AECH=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/annotation/gene_annotation/Acromyrmex_echinatior.v2.0.gff
ACEP=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_cephalotes/annotation/gene_annotation/Atta_cephalotes.v2.0.gff
ACOL=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_colombica/annotation/gene_annotation/Atta_colombica.v2.0.gff
CCOS=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Cyphomyrmex_costatus/annotation/gene_annotation/Cyphomyrmex_costatus.v2.0.gff
TCOR=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_cornetzi/annotation/gene_annotation/Trachymyrmex_cornetzi.v2.0.gff
TSEP=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_septentrionalis/annotation/gene_annotation/Trachymyrmex_septentrionalis.v2.0.gff
TZET=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Trachymyrmex_zeteki/annotation/gene_annotation/Trachymyrmex_zeteki.v2.0.gff


nice quast.py $AECHG -g $AECH -o Aech.quast.out --split-scaffolds --eukaryote
nice quast.py $ACEPG -g $ACEP -o Acep.quast.out --split-scaffolds --eukaryote
nice quast.py $ACOLG -g $ACOL -o Acol.quast.out --split-scaffolds --eukaryote
nice quast.py $CCOSG -g $CCOS -o Ccos.quast.out --split-scaffolds --eukaryote
nice quast.py $TCORG -g $TCOR -o Tcor.quast.out --split-scaffolds --eukaryote
