
# Synteny with i-adhore
https://isugenomics.github.io/bioinformatics-workbook/dataAnalysis/ComparativeGenomics/OrthofinderSynteny.html

base=~/data/inqGen18/synteny/OR
unassigned=/corefac/cse/lukas/inqGen18/ORphylogeny/orthofinder/Results_Apr19/Orthogroups_UnassignedGenes.csv
orthogroups=/corefac/cse/lukas/inqGen18/ORphylogeny/orthofinder/Results_Apr19/Orthogroups.txt
bd=/usr/local/home/lschrader/data/inqGen18/ORs/
cd $base
mkdir $base/data
cat $orthogroups|tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix $unassigned | awk 'NR>1 {print $2,$1}') |tr " " "\t" > $base/data/Orthologues.full.list

cat $base/data/Orthologues.full.list|sed 's/-//g' > $base/data/Orthologues.list
#This is the output that is needed family of orthologues (FYI: iadhore calls this a blast table). Essentially it is getting EVERY gene from your "Orthogroups.txt" and "Orthogroups_UnassignedGenes.csv", and creating a table assigning genes to a family.
# Create lists of files named by scaffold name
#Essentially here, you are extracting: scaffold name, protein name, strand orientation.

for file in $(find $bd/*/final/finalSet/  -name "*OR.gff3")
do
ln -s $file $base/data/.
done


species1=AINS
gff1=$base/data/Ains.OR.gff3
species2=AECH
gff2=$base/data/Aech.OR.gff3

# more species
species3=ACHA
gff3=$base/data/Acha.OR.gff3
species4=AHEY
gff4=$base/data/Ahey.OR.gff3
species5=PARG
gff5=$base/data/Parg.OR.gff3

species6=ACOL
gff6=$base/data/Acol.OR.gff3

#now make all the query files for iadhore
### This will make a two columned file with (protein_nameOrientation scaffold_name).  The scaffold name will be used to name the files in the last bit of code below, but we will get rid of it later.
mkdir $base/data/$species1
cd $base/data/$species1
rm $base/data/$species1/*
cat $gff1 |awk -v FS='\t' '{if ($3=="gene") gsub("-","",$9)}$1' |awk '{if ($3=="gene") print $1,$7,$9}'|perl -p -e 's/(.*) (.) .*?ID=(.*?);.*/$3$2 $1/g'|sort|uniq|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst

mkdir $base/data/$species2
cd $base/data/$species2
rm $base/data/$species2/*
cat $gff2 |awk -v FS='\t' '{if ($3=="gene") gsub("-","",$9)}$1' |awk '{if ($3=="gene") print $1,$7,$9}'|perl -p -e 's/(.*) (.) .*?ID=(.*?);.*/$3$2 $1/g'|sort|uniq|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst

mkdir $base/data/$species3
cd $base/data/$species3
rm $base/data/$species3/*
cat $gff3 |awk -v FS='\t' '{if ($3=="gene") gsub("-","",$9)}$1' |awk '{if ($3=="gene") print $1,$7,$9}'|perl -p -e 's/(.*) (.) .*?ID=(.*?);.*/$3$2 $1/g'|sort|uniq|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst

mkdir $base/data/$species4
cd $base/data/$species4
rm $base/data/$species4/*
cat $gff4 |awk -v FS='\t' '{if ($3=="gene") gsub("-","",$9)}$1' |awk '{if ($3=="gene") print $1,$7,$9}'|perl -p -e 's/(.*) (.) .*?ID=(.*?);.*/$3$2 $1/g'|sort|uniq|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst

mkdir $base/data/$species5
cd $base/data/$species5
rm $base/data/$species5/*
cat $gff5 |awk -v FS='\t' '{if ($3=="gene") gsub("-","",$9)}$1' |awk '{if ($3=="gene") print $1,$7,$9}'|perl -p -e 's/(.*) (.) .*?ID=(.*?);.*/$3$2 $1/g'|sort|uniq|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst

mkdir $base/data/$species6
cd $base/data/$species6
rm $base/data/$species6/*
cat $gff6 |awk -v FS='\t' '{if ($3=="gene") gsub("-","",$9)}$1' |awk '{if ($3=="gene") print $1,$7,$9}'|perl -p -e 's/(.*) (.) .*?ID=(.*?);.*/$3$2 $1/g'|sort|uniq|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst

cd $base/data/

ls $species1/*lst >input.txt
paste <(cut -f 1 -d "." input.txt|cut -f 2 -d "/") <(awk '{print "data/"$1}' input.txt)|sed "s/\t/ /g">$species1.ini
ls $species2/*lst >input.txt
paste <(cut -f 1 -d "." input.txt|cut -f 2 -d "/") <(awk '{print "data/"$1}' input.txt)|sed "s/\t/ /g">$species2.ini
ls $species3/*lst >input.txt
paste <(cut -f 1 -d "." input.txt|cut -f 2 -d "/") <(awk '{print "data/"$1}' input.txt)|sed "s/\t/ /g">$species3.ini
ls $species4/*lst >input.txt
paste <(cut -f 1 -d "." input.txt|cut -f 2 -d "/") <(awk '{print "data/"$1}' input.txt)|sed "s/\t/ /g">$species4.ini
ls $species5/*lst >input.txt
paste <(cut -f 1 -d "." input.txt|cut -f 2 -d "/") <(awk '{print "data/"$1}' input.txt)|sed "s/\t/ /g">$species5.ini
ls $species6/*lst >input.txt
paste <(cut -f 1 -d "." input.txt|cut -f 2 -d "/") <(awk '{print "data/"$1}' input.txt)|sed "s/\t/ /g">$species6.ini
cd $base

## Simple run
echo -e "\nblast_table=data/Orthologues.list
table_type=family
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
write_stats=true
output_path= output1
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=true
q_value=.05" > blast.settings.ini

## Simple run 2
echo -e "\nblast_table=data/Orthologues.list
table_type=family
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
write_stats=true
output_path= output2
alignment_method=gg2
gap_size=40
cluster_gap=55
level_2_only=true
q_value=.05" > blast.settings.r2.ini


echo -e "\nblast_table=data/Orthologues.list
table_type= family
output_path= output.multiple
alignment_method=gg2
gap_size=20
cluster_gap=35
q_value=0.75
prob_cutoff=0.01
anchor_points=3
level_2_only=false
number_of_threads=4
visualizeAlignment=true
" > blast.settings.cloud.ini

echo -e "\nblast_table=data/Orthologues.list
table_type= family
output_path= output.multiple2
alignment_method=gg2
gap_size=20
cluster_gap=35
q_value=0.99
prob_cutoff=0.1
anchor_points=4
level_2_only=false
number_of_threads=4
visualizeAlignment=true
" > blast.settings.cloud2.ini

# Cloud mode and hybrid modes fail
echo -e "\nblast_table=data/Orthologues.list
table_type= family
output_path= output.multiple3
alignment_method=gg2
gap_size=40
cluster_gap=40
q_value=0.4
prob_cutoff=0.01
anchor_points=3
level_2_only=false
number_of_threads=4
visualizeAlignment=true
write_stats=true
verbose_output=true
" > blast.settings.cloud3.ini

# Cloud mode and hybrid modes fail
echo -e "\nblast_table=data/Orthologues.list
table_type= family
output_path= output.multiple4
alignment_method=gg2
gap_size=20
cluster_gap=20
q_value=0.8
prob_cutoff=0.001
anchor_points=3
level_2_only=false
number_of_threads=4
visualizeAlignment=true
write_stats=true
verbose_output=true
" > blast.settings.cloud4.ini


echo "genome=$species1" > $species1.tag
echo "genome=$species2" > $species2.tag
echo "genome=$species3" > $species3.tag
echo "genome=$species4" > $species4.tag
echo "genome=$species5" > $species5.tag
echo "genome=$species6" > $species6.tag

cd $base/
cat $species1.tag data/$species1.ini $species2.tag data/$species2.ini $species3.tag data/$species3.ini $species4.tag data/$species4.ini $species5.tag data/$species5.ini $species6.tag data/$species6.ini blast.settings.ini > iadhore.ini
cd $base
i-adhore  iadhore.ini

# Best multiple
cd $base
cat $species1.tag data/$species1.ini $species2.tag data/$species2.ini $species3.tag data/$species3.ini $species4.tag data/$species4.ini $species5.tag data/$species5.ini $species6.tag data/$species6.ini blast.settings.cloud.ini > iadhore.cloud.ini
cd $base
i-adhore  iadhore.cloud.ini

cd $base
cat $species1.tag data/$species1.ini $species2.tag data/$species2.ini $species3.tag data/$species3.ini $species4.tag data/$species4.ini $species5.tag data/$species5.ini $species6.tag data/$species6.ini blast.settings.cloud2.ini > iadhore.cloud2.ini
cd $base
i-adhore  iadhore.cloud2.ini

cd $base
cat $species1.tag data/$species1.ini $species2.tag data/$species2.ini $species3.tag data/$species3.ini $species4.tag data/$species4.ini $species5.tag data/$species5.ini $species6.tag data/$species6.ini blast.settings.cloud3.ini > iadhore.cloud3.ini
i-adhore  iadhore.cloud3.ini

cd $base
cat $species1.tag data/$species1.ini $species2.tag data/$species2.ini $species3.tag data/$species3.ini $species4.tag data/$species4.ini $species5.tag data/$species5.ini $species6.tag data/$species6.ini blast.settings.cloud4.ini > iadhore.cloud4.ini
i-adhore  iadhore.cloud4.ini


grep "of the genome" output.multiple3/collinear_portions.txt
