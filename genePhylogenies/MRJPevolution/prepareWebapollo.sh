# ACOL
species="Acol"
folder="/global/homes/jg/schradel/data/genomes/reannotations_of_7_ants/Atta_colombica/"
genome=$(readlink -f $folder"/genome/Atta_colombica.v2.0.fa")
geneGFF=$folder"/annotation/gene_annotation/Atta_colombica.v2.0.gff"
ls $geneGFF
software=/home/s/schradel/software/JBrowse-1.15.0/bin
perl $software/prepare-refseqs.pl --fasta $genome --out ./$species
perl $software/flatfile-to-json.pl --gff $geneGFF --type mRNA --trackLabel GLEAN --out ./$species
#perl $software/flatfile-to-json.pl --gff $geneGFFncbi --type mRNA --trackLabel ncbi --out ./$species
#perl $software/flatfile-to-json.pl --gff $species"_data/"$species".exonerate.final.gff" --type mRNA --trackLabel exonerate --out $species
#perl $software/flatfile-to-json.pl --gff $species"_data/"*.GeMoMa.final.gff --type mRNA --trackLabel GeMoMa --out ./$species
#perl $software/flatfile-to-json.pl --gff $species"_data/"$species.OR.gff3 --type mRNA --trackLabel ORs --out ./$species
#perl $software/flatfile-to-json.pl --gff $species"_data/"*.OR.transcripts.gff3 --type mRNA --trackLabel prot --out ./$species
#perl $software/flatfile-to-json.pl --gff $geneGFF2 --type mRNA --trackLabel prot --out ./$species
#perl $software/flatfile-to-json.pl --gff $species"_data/"*.vs.Aech.GeMoMa.final.gff --type mRNA --trackLabel $species.vs.Aech.GeMoMa.final.gff --out ./$species
perl $software/generate-names.pl -v --out ./$species
# tar folder and copy
ls
tar -zcvf $species.tar.gz $species
scp $species.tar.gz schradel@iebapollo:/opt/apollo
ssh iebapollo

# at iebapollo
cd /opt/apollo
tar xvf Acol.tar.gz
chmod -R 755 Acol

# ACOL
species="Acep"
folder="/global/homes/jg/schradel/data/genomes/reannotations_of_7_ants/Atta_cephalotes/"
genome=$(readlink -f $folder"/genome/Atta_cephalotes.v2.0.fa")
geneGFF=$folder"/annotation/gene_annotation/Atta_cephalotes.v2.0.gff"
ls $geneGFF
software=/home/s/schradel/software/JBrowse-1.15.0/bin
perl $software/prepare-refseqs.pl --fasta $genome --out ./$species
perl $software/flatfile-to-json.pl --gff $geneGFF --type mRNA --trackLabel GLEAN --out ./$species
#perl $software/flatfile-to-json.pl --gff $geneGFFncbi --type mRNA --trackLabel ncbi --out ./$species
#perl $software/flatfile-to-json.pl --gff $species"_data/"$species".exonerate.final.gff" --type mRNA --trackLabel exonerate --out $species
#perl $software/flatfile-to-json.pl --gff $species"_data/"*.GeMoMa.final.gff --type mRNA --trackLabel GeMoMa --out ./$species
#perl $software/flatfile-to-json.pl --gff $species"_data/"$species.OR.gff3 --type mRNA --trackLabel ORs --out ./$species
#perl $software/flatfile-to-json.pl --gff $species"_data/"*.OR.transcripts.gff3 --type mRNA --trackLabel prot --out ./$species
#perl $software/flatfile-to-json.pl --gff $geneGFF2 --type mRNA --trackLabel prot --out ./$species
#perl $software/flatfile-to-json.pl --gff $species"_data/"*.vs.Aech.GeMoMa.final.gff --type mRNA --trackLabel $species.vs.Aech.GeMoMa.final.gff --out ./$species
perl $software/generate-names.pl -v --out ./$species
# tar folder and copy
#ls
tar -zcvf $species.tar.gz $species
scp $species.tar.gz schradel@iebapollo:/opt/apollo
ssh iebapollo

# at iebapollo
cd /opt/apollo
tar xvf Acep.tar.gz
chmod -R 755 Acep
