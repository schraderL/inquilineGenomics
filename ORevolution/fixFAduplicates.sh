export EVMtools=$software/EVidenceModeler-1.1.1/EvmUtils/
############################
#Folders
############################
# set base directory
export base=/usr/local/home/lschrader/data/inqGen18/ORs
# set directory where all scripts are located
export scripts=$base/scripts
# set directory where software is installed
export software=/usr/local/home/lschrader/software
# set number of CPUs to use
export cpu=8
# set jar excecutable of GeMoMa
export GeMoMaProgram=$software/GeMoMa-1.5.2/GeMoMa-1.5.2.jar

#Aech Acromyrmex echinatior
export species=Aech


# With the following 2 commands you can check whether the fasta contains genes not present in the gff
#awk '{if ($3=="gene") print $9}' finalSet/$species.OR.gff3 |cut -f 3 -d "=" > finalSet/tmp.lst
#grep ">" finalSet/$species.OR.fa |grep -f finalSet/tmp.lst -v

#Parg Pseudoatta_argentina
export species=Parg
Tsep
#for species in Ains Aech Acep Acha Acol Agra Ahey Alob Ccos Pcal Plob Tcor Tsep Tzet
Ahey
for species in  Pcal Tsep
do

cd $base/$species/final

# Rename ORs in R
Rscript $scripts/rename.R $base"/"$species"/final" $species
# Fix genes where gff is not perfect yet.
egrep $species"Or-" renamed.gff3 > tmp.gff3
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= tmp.gff3
# Clean up the gff
$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries tmp_clean.gff > $species.OR.gff3
# validate that the gff is ok now.
perl $software/EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl $species.OR.gff3
cat renamed.tsv|awk 'BEGIN {FS="\t"}; {print ">"$33,$34"\n"$3}'|fold -w 80 > $species.OR.fa
rm tmp.gff3 renamed.gff3 tmp_clean.gff
mv $species.OR.gff3 finalSet
mv $species.OR.fa finalSet
mv renamed.tsv finalSet/$species.OR.tsv

echo "Files: " > finalSet/readme.txt
echo " $species.OR.gff3 : GFF3 file of OR annotations" >> finalSet/readme.txt
echo " $species.OR.fa : protein fa file of OR annotations" >> finalSet/readme.txt
echo " $species.OR.tsv : annotation details of OR annotations" >> finalSet/readme.txt
echo "Nomenclature: " >> finalSet/readme.txt
echo " NTE = missing start, i.e. missing N-terminus." >> finalSet/readme.txt
echo " CTE = missing stop, i.e. missing C-terminus." >> finalSet/readme.txt
echo " NC = missing start and stop, i.e. missing N- and C-terminus." >> finalSet/readme.txt
echo " NC = missing start and stop, i.e. missing N- and C-terminus." >> finalSet/readme.txt
echo " fd = domain fragmented at n-terminus. Missing N-terminal piece of domain (>50 aa)." >> finalSet/readme.txt
echo " df = domain fragmented at c-terminus. Missing C-terminal piece of domain (>50 aa)." >> finalSet/readme.txt
echo " fdf = domain fragmented at n- and c-terminus. Missing pieces of domain >50 aa." >> finalSet/readme.txt
echo " dH = domain missing, but protein homologous to other ORs. (blast 1e-10)" >> finalSet/readme.txt
done
#################################################################
# 8. Compute some statistics
#################################################################

#for species in Parg Aech Acep Acha Acol Agra Ahey Alob Ccos Pcal Plob Tcor Tsep Tzet
for species in Ahey Pcal Tsep
do
cd $base/$species/final/finalSet
mkdir stats
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat --genelengthdistri -o stats/$species.genelength.stats.txt --force
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -genescoredistri -o stats/$species.genescore.stats.txt --force
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -exonlengthdistri -o stats/$species.exonlength.stats.txt --force
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -exonnumberdistri -o stats/$species.exonnumber.stats.txt --force
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -intronlengthdistri -o stats/$species.intronlength.stats.txt --force
cat $species.OR.gff3|$software/genometools-1.5.9/bin/gt stat -cdslengthdistri -o stats/$species.cdslength.stats.txt --force
python $software/GAG/gag.py -f $genomeFa -g $species.OR.gff3

done


#for species in Parg Aech Acep Acha Acol Agra Ahey Alob Ccos Pcal Plob Tcor Tsep Tzet
for species in Ains Ains2.0 Parg Aech Acep Acha Acol Agra Ahey Alob Ccos Pcal Plob Tcor Tsep Tzet
do
  cd $base/$species/final/finalSet
  #check again whether gff and fastas are ok now.
  echo $species
  awk '{if ($3=="gene") print $9}' *.OR.gff3 |cut -f 3 -d "=" > tmp.lst
  wc -l tmp.lst
  grep ">" *.OR.fa |grep -f tmp.lst -v
  cut -f 1,3,4,5 *.OR.gff3|grep "#" -v |sort|uniq -d
  cd $base
done


for species in Ains Ains2.0
do
echo $species":"
cut -f 1 -d " " $base/$species/final/finalSet/*.OR.fa|grep ">"|cut -f 3 -d "-"|sort|uniq -c
done
