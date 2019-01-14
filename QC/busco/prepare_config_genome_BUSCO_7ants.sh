# generate config files for BUSCO genome runs for 7 ant genomes

# Acep genome config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Atta_cephalotes/genome/Atta_cephalotes.v2.0.fa\
|g" tmp.genome.config.ini |\
sed "s|out = outRun|out = Acep_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Acep/genome/|g"\
> Acep.genome.config.ini


# Acol genome config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Atta_colombica/genome/Atta_colombica.v2.0.fa\
|g" tmp.genome.config.ini |\
sed "s|out = outRun|out = Acol_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Acol/genome/|g"\
> Acol.genome.config.ini


# Tzet genome config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Trachymyrmex_zeteki/genome/Trachymyrmex_zeteki.v2.0.fa\
|g" tmp.genome.config.ini |\
sed "s|out = outRun|out = Tzet_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Tzet/genome/|g"\
> Tzet.genome.config.ini


# Tsep genome config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Trachymyrmex_septentrionalis/genome/Trachymyrmex_septentrionalis.v2.0.fa\
|g" tmp.genome.config.ini |\
sed "s|out = outRun|out = Tsep_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Tsep/genome/|g"\
> Tsep.genome.config.ini


# Tcor genome config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Trachymyrmex_cornetzi/genome/Trachymyrmex_cornetzi.v2.0.fa\
|g" tmp.genome.config.ini |\
sed "s|out = outRun|out = Tcor_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Tcor/genome/|g"\
> Tcor.genome.config.ini

# Ccos genome config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Cyphomyrmex_costatus/genome/Cyphomyrmex_costatus.v2.0.fa\
|g" tmp.genome.config.ini |\
sed "s|out = outRun|out = Ccos_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Ccos/genome/|g"\
> Ccos.genome.config.ini


# Aech genome config
# replace genome fasta, then output name, then output folder
sed "s|\
inputFasta|\
/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/genome/Acromyrmex_echinatior.v2.0.fa\
|g" tmp.genome.config.ini |\
sed "s|out = outRun|out = Aech_genome|g"|\
sed "s|outPath|/usr/local/home/lschrader/data/inqGen18/QC/BUSCO_7ants/Aech/genome/|g"\
> Aech.genome.config.ini
