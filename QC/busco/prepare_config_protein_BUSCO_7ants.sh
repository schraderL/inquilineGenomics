# generate config files for BUSCO protein runs for 7 ant genomes

# Acep protein config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Atta_cephalotes/annotation/gene_annotation/Atta_cephalotes.v2.0.pep.fa\
|g" tmp.protein.config.ini |\
sed "s|out = outRun|out = Acep_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Acep/protein/|g"\
> Acep.protein.config.ini


# Acol protein config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Atta_colombica/annotation/gene_annotation/Atta_colombica.v2.0.pep.fa.fa\
|g" tmp.protein.config.ini |\
sed "s|out = outRun|out = Acol_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Acol/protein/|g"\
> Acol.protein.config.ini


# Tzet protein config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Trachymyrmex_zeteki/annotation/gene_annotation/Trachymyrmex_zeteki.v2.0.pep.fa\
|g" tmp.protein.config.ini |\
sed "s|out = outRun|out = Tzet_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Tzet/protein/|g"\
> Tzet.protein.config.ini


# Tsep protein config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Trachymyrmex_septentrionalis/annotation/gene_annotation/Trachymyrmex_septentrionalis.v2.0.pep.fa\
|g" tmp.protein.config.ini |\
sed "s|out = outRun|out = Tsep_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Tsep/protein/|g"\
> Tsep.protein.config.ini


# Tcor protein config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Trachymyrmex_cornetzi/annotation/gene_annotation/Trachymyrmex_cornetzi.v2.0.pep.fa\
|g" tmp.protein.config.ini |\
sed "s|out = outRun|out = Tcor_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Tcor/protein/|g"\
> Tcor.protein.config.ini

# Ccos protein config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Cyphomyrmex_costatus/annotation/gene_annotation/Cyphomyrmex_costatus.v2.0.pep.fa\
|g" tmp.protein.config.ini |\
sed "s|out = outRun|out = Ccos_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Ccos/protein/|g"\
> Ccos.protein.config.ini


# Aech protein config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/annotation/gene_annotation/Acromyrmex_echinatior.v2.0.pep.fa\
|g" tmp.protein.config.ini |\
sed "s|out = outRun|out = Aech_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Aech/protein/|g"\
> Aech.protein.config.ini
