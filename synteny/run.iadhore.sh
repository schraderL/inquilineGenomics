
# Synteny with i-adhore
https://isugenomics.github.io/bioinformatics-workbook/dataAnalysis/ComparativeGenomics/OrthofinderSynteny.html

less ../Results_Orthofinder_Jun19/Orthogroups.txt |tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix ../Results_Orthofinder_Jun19/Orthogroups_UnassignedGenes.csv |awk 'NR>1 {print $2,$1}') |tr " " "\t" >Orthologues.list

base=~/data/inqGen18/synteny
unassigned=/corefac/cse/lukas/inqGen18/orthologs/orthofinder/Results_Sep03/Orthogroups_UnassignedGenes.csv
orthogroups=/corefac/cse/lukas/inqGen18/orthologs/orthofinder/Results_Sep03/Orthogroups.txt


cd $base
mkdir $base/data
cat $orthogroups|tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix $unassigned | awk 'NR>1 {print $2,$1}') |tr " " "\t" > $base/data/Orthologues.full.list


####### test case with 2 species
cat $base/data/Orthologues.full.list|egrep "AINS|AECH" > $base/data/Orthologues.list

#This is the output that is needed family of orthologues (FYI: iadhore calls this a blast table). Essentially it is getting EVERY gene from your "Orthogroups.txt" and "Orthogroups_UnassignedGenes.csv", and creating a table assigning genes to a family.


# Create lists of files named by scaffold name
#Essentially here, you are extracting: scaffold name, protein name, strand orientation.

species1=AINS
gff1=/usr/local/home/lschrader/data/inqGen18//inquilines_v2.1/Acromyrmex_insinuator.2.1/annotation/gene_annotation/Acromyrmex_insinuator.v2.1.gff3
species2=AECH
gff2=/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/annotation/gene_annotation/Acromyrmex_echinatior.v2.0.gff

less $gff1 |awk '$3=="CDS"' |sed 's/;/\t&/g' |awk '{print $9}' |sort|uniq|sed 's/Parent=//g' |wc

#now make all the query files for iadhore
### This will make a two columned file with (protein_nameOrientation scaffold_name).  The scaffold name will be used to name the files in the last bit of code below, but we will get rid of it later.
mkdir $base/data/$species1
cd $base/data/$species1
rm $base/data/$species1/*
cat $gff1 |awk '{if ($3=="CDS") print $1,$7,$9}'|perl -p -e 's/(.*) (.) .*?Parent=(.*?)-.*/$3$2 $1/g'|sort|uniq|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst

mkdir $base/data/$species2
cd $base/data/$species2
cat $gff2 |awk '$3=="CDS"' |sed 's/;/\t&/g' |awk '{print $1,$7,$9}' |sed 's/Parent=//g'|awk '{print $3$2,$1}'  |sort|uniq|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
#head $base/data/AECH/C1964625.lst

cd $base/data/

ls $species1/*lst >input.txt
paste <(cut -f 1 -d "." input.txt|cut -f 2 -d "/") <(awk '{print "data/"$1}' input.txt)|sed "s/\t/ /g">$species1.ini

ls $species2/*lst >input.txt
paste <(cut -f 1 -d "." input.txt|cut -f 2 -d "/") <(awk '{print "data/"$1}' input.txt)|sed "s/\t/ /g">$species2.ini

cd $base
cat data/query.ini data/subject.ini |tr "\t" " " >iadhore.ini

echo -e "\nblast_table=data/Orthologues.list
table_type=family
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
write_stats=true
output_path= output3
alignment_method=gg2
gap_size=55
cluster_gap=60
level_2_only=true
q_value=.05" > blast.settings.ini

echo "genome=$species1" > $species1.tag
echo "genome=$species2" > $species2.tag

cat $species1.tag data/$species1.ini $species2.tag data/$species2.ini blast.settings.ini > iadhore.ini

cd $base
i-adhore iadhore.ini
