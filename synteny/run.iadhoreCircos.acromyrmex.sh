
#I create a new folder above the 02_iadhore folder. Again, some of my commands depend on this.
#/work/GIF/remkv6/USDA/15_OrthoFinderSynteny/03_circos

#Softlink all relevant files:(genome, GFF, segments.txt)
base=~/data/inqGen18/synteny.acromyrmex
mkdir $base/circos
cd $base/circos


species1=AINS
gff1=/usr/local/home/lschrader/data/inqGen18//inquilines_v2.1/Acromyrmex_insinuator.2.1/annotation/gene_annotation/Acromyrmex_insinuator.v2.1.gff3
species2=AECH
gff2=/usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/annotation/gene_annotation/Acromyrmex_echinatior.v2.0.gff
# more species
species3=ACHA
gff3=/usr/local/home/lschrader/data/inqGen18//inquilines_v2.1/Acromyrmex_charruanus.2.1/annotation/gene_annotation/Acromyrmex_charruanus.v2.1.gff3
species4=AHEY
gff4=/usr/local/home/lschrader/data/inqGen18//inquilines_v2.1/Acromyrmex_heyeri.2.1/annotation/gene_annotation/Acromyrmex_heyeri.v2.1.gff3
species5=PARG
gff5=/usr/local/home/lschrader/data/inqGen18//inquilines_v2.1/Pseudoatta_argentina.2.1/annotation/gene_annotation/Pseudoatta_argentina.v2.1.gff3


ln -s /usr/local/home/lschrader/data/inqGen18//inquilines_v2.1/Acromyrmex_insinuator.2.1/genome/Acromyrmex_insinuator.v2.1.fa ./$species1.fna
ln -s /usr/local/home/lschrader/data/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/genome/Acromyrmex_echinatior.v2.0.fa ./$species2.fna
ln -s /usr/local/home/lschrader/data/inqGen18//inquilines_v2.1/Acromyrmex_charruanus.2.1/genome/Acromyrmex_charruanus.v2.1.fa ./$species3.fna
ln -s /usr/local/home/lschrader/data/inqGen18//inquilines_v2.1/Acromyrmex_heyeri.2.1/genome/Acromyrmex_heyeri.v2.1.fa ./$species4.fna
ln -s /usr/local/home/lschrader/data/inqGen18//inquilines_v2.1/Pseudoatta_argentina.2.1/genome/Pseudoatta_argentina.v2.1.fa ./$species5.fna

ln -s $gff1 ./$species1.gff
ln -s $gff2 ./$species2.gff
ln -s $gff3 ./$species3.gff
ln -s $gff4 ./$species4.gff
ln -s $gff5 ./$species5.gff
ln -s $base/output/segments.txt

#All I am doing here is creating a column in the gff that has only the protein name.  This way I can use grep -w to get exact matches later on.

sed 's/;/\t/g' $species1.gff |sed 's/Parent=//g' |sed -e 's/-mRNA$//g' |awk '$3=="mRNA"'|cut -f 1-8,10 > $species1.GrepMod.gff
sed 's/;/\t/g' $species2.gff |sed 's/Parent=//g' |sed -e 's/-mRNA$//g' |awk '$3=="mRNA"'|cut -f 1-9|sed "s/ID=//g"  > $species2.GrepMod.gff
sed 's/;/\t/g' $species3.gff |sed 's/Parent=//g' |sed -e 's/-mRNA$//g' |awk '$3=="mRNA"'|cut -f 1-8,10 > $species3.GrepMod.gff
sed 's/;/\t/g' $species4.gff |sed 's/Parent=//g' |sed -e 's/-mRNA$//g' |awk '$3=="mRNA"'|cut -f 1-8,10 > $species4.GrepMod.gff
sed 's/;/\t/g' $species5.gff |sed 's/Parent=//g' |sed -e 's/-mRNA$//g' |awk '$3=="mRNA"'|cut -f 1-8,10 > $species5.GrepMod.gff

# The five scripts below will just work if you change 3 things.
# 1.  change "pathenogenic" to whatever you named your second genome in your iadhore.ini file.
# 2.  change 24891GrepMod.gff to the grepMod gff that is associated with the second genome name in your iadhore.ini file.
# 3. change the 5845GrepMod.gff to the grepMod gff you created that is associated with the first genome name in your iadhore.ini file.
#Essentially what is happening below is that you swapping columns in segments.txt until you get pathogenic all on one side.  Then I extract the 5' position for the 5' syntenic gene and the 3' position for the 3' syntenic gene for each genome

# AINS vs AECH
cat *.GrepMod.gff > GrepMod.gff
cat ../output/segments.txt|awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk -v var=$species2 '{if($5==var){print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line GrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
cat ../output/segments.txt|awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk -v var=$species2 '{if($5==var) {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line GrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8
cat ../output/segments.txt|awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk -v var=$species2 '{if($5==var) {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line GrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
cat ../output/segments.txt|awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk -v var=$species2 '{if($5==var){print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line GrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4

#This last step adds the scaffold names to the gene positions extracted above.
cat ../output/segments.txt|awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk -v var=$species2 '{if($5==var) {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $1"."$2,$5"."$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf

# add species tags to scaffoldIDs

#Here is the SyntenicRibbons.conf file #scaffold position position,scaffold, position, position
#These two commands are essentially extracting the scaffold lengths in your genome and putting them in the proper format.
bioawk -c fastx '{print $name,length($seq)}' $species1.fna  |awk -v var=$species1 '{print "chr","-",var"."$1,$1,"0",$2,"blue"}'  > $species1.conf
bioawk -c fastx '{print $name,length($seq)}' $species2.fna  |awk -v var=$species2 '{print "chr","-",var"."$1,$1,"0",$2,"green"}' > $species2.conf

#The next six scripts below are essentially extracting the scaffolds that have some synteny. You dont want to display those scaffolds that do not have any information, right?.  Make sure you have the proper column for each extraction.  Remember column 1 is one species' scaffolds, and column 4 is the other species' scaffolds
awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' $species1.conf >>$species1.conf1";done >$species1.sh
rm $species1.conf1
sh $species1.sh

awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' $species2.conf >>$species2.conf1";done >$species2.sh
rm $species2.conf1
sh $species2.sh

cat <(sort $species1.conf1 |uniq) <(sort $species2.conf1 |uniq) >karyotype.conf


#Now lets reduce the number of times the circos synteny plot lines overlap, so it is more pleasing to the eye.
#I just download this tool everytime because it is small and easier than finding the original circos installation directory

 #We will use the tmpKaryotype.conf1 file to get the scaffold names that we want grouped together.  You can also use tmpKaryotype.conf2 to do this.  I would suggest using the file that is the smallest.
 #the below script generates the command.
 sort $species2.conf1 |uniq|awk '{print $3}' |tr "\n" "," |sed 's/.$//' |awk '{print "~/software/circos-tools-0.23/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - "$0" -static_rx "$0 }' |less
#the last bit is what we want chromosomes_order = .....

chromosomes_order = AINS.scaffold49,AECH.scaffold207,AECH.scaffold52,AECH.scaffold178,AECH.scaffold341,AECH.scaffold85,AECH.scaffold193,AECH.scaffold25,AECH.scaffold335,AECH.scaffold201,AECH.scaffold952,AECH.scaffold536,AECH.scaffold44,AECH.scaffold202,AECH.scaffold682,AECH.scaffold48,AECH.scaffold5,AECH.scaffold95,AINS.scaffold275,AECH.scaffold568,AECH.scaffold10,AECH.scaffold642,AECH.scaffold16,AECH.scaffold185,AECH.scaffold412,AECH.scaffold108,AECH.scaffold408,AECH.scaffold258,AECH.scaffold138,AECH.scaffold326,AINS.scaffold507,AINS.scaffold12,AINS.scaffold408,AINS.scaffold556,AINS.scaffold29,AINS.scaffold153,AINS.scaffold18,AINS.scaffold307,AINS.scaffold179,AECH.scaffold26,AECH.scaffold227,AECH.scaffold109,AECH.scaffold358,AECH.scaffold513,AECH.scaffold419,AECH.scaffold65,AINS.scaffold21,AECH.scaffold70,AECH.scaffold569,AECH.scaffold621,AECH.scaffold120,AECH.scaffold171,AECH.scaffold650,AECH.scaffold136,AINS.scaffold402,AECH.scaffold40,AECH.scaffold27,AECH.scaffold57,AECH.scaffold673,AINS.scaffold73,AECH.scaffold571,AECH.scaffold22,AINS.scaffold26,AECH.scaffold438,AECH.scaffold483,AECH.scaffold248,AECH.scaffold250,AECH.scaffold179,AECH.scaffold338,AINS.scaffold66,AECH.scaffold300,AECH.scaffold47,AECH.scaffold121,AINS.scaffold397,AINS.scaffold762,AINS.scaffold348,AINS.scaffold576,AECH.scaffold156,AINS.scaffold521,AINS.scaffold260,AINS.scaffold787,AINS.scaffold844,AINS.scaffold25,AINS.scaffold559,AINS.scaffold638,AINS.scaffold102,AINS.scaffold301,AINS.scaffold288,AINS.scaffold41,AINS.scaffold229,AECH.scaffold274,AINS.scaffold580,AINS.scaffold567,AINS.scaffold242,AINS.scaffold388,AINS.scaffold212,AINS.scaffold345,AECH.scaffold456,AECH.scaffold390,AECH.scaffold18,AECH.scaffold303,AECH.scaffold517,AINS.scaffold184,AECH.scaffold142,AECH.scaffold501,AECH.scaffold69,AECH.scaffold135,AECH.scaffold1539,AINS.scaffold45,AECH.scaffold343,AECH.scaffold88,AECH.scaffold50,AECH.scaffold442,AINS.scaffold63,AECH.scaffold283,AECH.scaffold99,AECH.scaffold182,AECH.scaffold627,AECH.scaffold817,AINS.scaffold700,AINS.scaffold51,AINS.scaffold392,AINS.scaffold35,AINS.scaffold138,AECH.scaffold357,AECH.scaffold458,AECH.scaffold144,AECH.scaffold107,AINS.scaffold238,AECH.scaffold541,AECH.scaffold319,AECH.scaffold480,AECH.scaffold591,AECH.scaffold152,AECH.scaffold91,AECH.scaffold612,AECH.scaffold342,AECH.scaffold192,AINS.scaffold349,AECH.scaffold140,AECH.scaffold377,AECH.scaffold110,AECH.scaffold555,AECH.scaffold400,AECH.scaffold23,AINS.scaffold396,AINS.scaffold304,AINS.scaffold34,AECH.scaffold206,AINS.scaffold158,AINS.scaffold327,AINS.scaffold857,AINS.scaffold230,AINS.scaffold600,AINS.scaffold190,AECH.scaffold232,AINS.scaffold684,AINS.scaffold644,AINS.scaffold671,AINS.scaffold78,AECH.scaffold244,AINS.scaffold116,AINS.scaffold163,AINS.scaffold316,AINS.scaffold169,AINS.scaffold92,AINS.scaffold303,AINS.scaffold19,AINS.scaffold159,AECH.scaffold7,AINS.scaffold338,AINS.scaffold389,AINS.scaffold96,AINS.scaffold155,AECH.scaffold459,AECH.scaffold471,AECH.scaffold151,AECH.scaffold511,AECH.scaffold553,AECH.scaffold970,AECH.scaffold162,AINS.scaffold156,AINS.scaffold126,AINS.scaffold198,AINS.scaffold448,AINS.scaffold196,AECH.scaffold490,AECH.scaffold21,AINS.scaffold337,AECH.scaffold316,AECH.scaffold385,AECH.scaffold184,AECH.scaffold143,AECH.scaffold340,AINS.scaffold434,AECH.scaffold174,AECH.scaffold299,AINS.scaffold263,AINS.scaffold406,AINS.scaffold333,AINS.scaffold622,AECH.scaffold153,AINS.scaffold226,AINS.scaffold4,AECH.scaffold196,AINS.scaffold149,AINS.scaffold634,AINS.scaffold218,AECH.scaffold215,AINS.scaffold512,AINS.scaffold204,AINS.scaffold462,AECH.scaffold28,AINS.scaffold114,AINS.scaffold312,AINS.scaffold145,AECH.scaffold29,AINS.scaffold385,AINS.scaffold501,AINS.scaffold426,AECH.scaffold317,AINS.scaffold357,AINS.scaffold645,AINS.scaffold605,AECH.scaffold331,AINS.scaffold117,AINS.scaffold160,AINS.scaffold463,AECH.scaffold43,AINS.scaffold761,AINS.scaffold7,AINS.scaffold703,AECH.scaffold527,AINS.scaffold589,AINS.scaffold560,AINS.scaffold67,AINS.scaffold447,AINS.scaffold550,AECH.scaffold585,AINS.scaffold530,AINS.scaffold248,AINS.scaffold428,AECH.scaffold67,AINS.scaffold805,AINS.scaffold245,AINS.scaffold336,AINS.scaffold289,AINS.scaffold140,AINS.scaffold712,AECH.scaffold176,AECH.scaffold322,AINS.scaffold144,AECH.scaffold379,AECH.scaffold116,AECH.scaffold424,AINS.scaffold261,AECH.scaffold288,AECH.scaffold397,AECH.scaffold239,AECH.scaffold646,AECH.scaffold127,AECH.scaffold350,AECH.scaffold487,AECH.scaffold375,AINS.scaffold410,AECH.scaffold294,AECH.scaffold519,AECH.scaffold807,AINS.scaffold467,AECH.scaffold166,AECH.scaffold577,AECH.scaffold1043,AINS.scaffold538,AECH.scaffold262,AECH.scaffold658,AECH.scaffold1054,AINS.scaffold232,AECH.scaffold633,AINS.scaffold57,AECH.scaffold175,AECH.scaffold545,AINS.scaffold588,AECH.scaffold708,AECH.scaffold159,AECH.scaffold808,AINS.scaffold80,AECH.scaffold289,AECH.scaffold411,AECH.scaffold719,AECH.scaffold373,AINS.scaffold99,AECH.scaffold181,AECH.scaffold1,AECH.scaffold255,AECH.scaffold61,AINS.scaffold119,AINS.scaffold235,AINS.scaffold109,AECH.scaffold686,AECH.scaffold360,AECH.scaffold237,AINS.scaffold143,AECH.scaffold520,AECH.scaffold671,AECH.scaffold586,AINS.scaffold150,AECH.scaffold39,AECH.scaffold220,AINS.scaffold166,AECH.scaffold56,AECH.scaffold614,AECH.scaffold309,AINS.scaffold191,AECH.scaffold392,AECH.scaffold445,AECH.scaffold533,AECH.scaffold947,AINS.scaffold455,AECH.scaffold307,AECH.scaffold540,AECH.scaffold242,AINS.scaffold221,AECH.scaffold190,AECH.scaffold183,AINS.scaffold241,AECH.scaffold59,AECH.scaffold702,AINS.scaffold250,AECH.scaffold380,AECH.scaffold315,AINS.scaffold255,AECH.scaffold298,AECH.scaffold141,AINS.scaffold269,AECH.scaffold285,AECH.scaffold81,AINS.scaffold272,AECH.scaffold172,AINS.scaffold273,AECH.scaffold64,AINS.scaffold274,AECH.scaffold940,AECH.scaffold189,AINS.scaffold319,AECH.scaffold784,AECH.scaffold751,AINS.scaffold325,AECH.scaffold267,AECH.scaffold763,AECH.scaffold93,AINS.scaffold360,AECH.scaffold507,AECH.scaffold132,AINS.scaffold362,AECH.scaffold296,AINS.scaffold377,AECH.scaffold12,AECH.scaffold327,AECH.scaffold767,AECH.scaffold474,AINS.scaffold401,AECH.scaffold678,AINS.scaffold412,AECH.scaffold90,AECH.scaffold635,AINS.scaffold424,AECH.scaffold19,AINS.scaffold425,AECH.scaffold104,AECH.scaffold649,AECH.scaffold293,AINS.scaffold211,AECH.scaffold750,AINS.scaffold174,AECH.scaffold532,AINS.scaffold481,AECH.scaffold131,AECH.scaffold518,AINS.scaffold522,AECH.scaffold351,AECH.scaffold503,AINS.scaffold55,AECH.scaffold177,AECH.scaffold106,AECH.scaffold544,AINS.scaffold592,AECH.scaffold123,AECH.scaffold450,AINS.scaffold606,AECH.scaffold663,AECH.scaffold96,AECH.scaffold195,AECH.scaffold86,AINS.scaffold649,AECH.scaffold356,AINS.scaffold127,AINS.scaffold98,AECH.scaffold584,AECH.scaffold66,AINS.scaffold135,AINS.scaffold172,AINS.scaffold543,AINS.scaffold749,AINS.scaffold178,AINS.scaffold1057,AINS.scaffold50,AINS.scaffold564,AECH.scaffold186,AINS.scaffold210,AINS.scaffold104,AINS.scaffold710,AINS.scaffold148,AECH.scaffold221,AINS.scaffold491,AINS.scaffold195,AINS.scaffold110,AINS.scaffold755,AINS.scaffold972,AINS.scaffold61,AINS.scaffold652,AECH.scaffold324,AINS.scaffold299,AINS.scaffold122,AECH.scaffold346,AINS.scaffold499,AINS.scaffold452,AINS.scaffold873,AECH.scaffold35,AINS.scaffold220,AINS.scaffold514,AECH.scaffold367,AINS.scaffold285,AINS.scaffold774,AINS.scaffold38,AECH.scaffold482,AINS.scaffold361,AINS.scaffold133,AINS.scaffold165,AINS.scaffold390,AINS.scaffold39,AECH.scaffold531,AINS.scaffold440,AINS.scaffold356,AINS.scaffold777,AINS.scaffold702,AINS.scaffold214,AECH.scaffold580,AINS.scaffold369,AINS.scaffold189,AINS.scaffold394,AINS.scaffold31,AINS.scaffold458,AINS.scaffold271,AINS.scaffold85,AINS.scaffold436,AINS.scaffold42,AECH.scaffold8,AINS.scaffold1,AINS.scaffold404,AINS.scaffold37,AINS.scaffold311,AECH.scaffold92,AINS.scaffold411,AINS.scaffold3,AECH.scaffold604,AINS.scaffold141,AECH.scaffold226,AECH.scaffold439,AECH.scaffold669,AINS.scaffold162,AECH.scaffold89,AECH.scaffold429,AECH.scaffold473,AINS.scaffold27,AECH.scaffold508,AECH.scaffold810,AECH.scaffold229,AINS.scaffold32,AECH.scaffold740,AECH.scaffold657,AINS.scaffold366,AECH.scaffold344,AECH.scaffold765,AINS.scaffold379,AECH.scaffold271,AECH.scaffold526,AECH.scaffold489,AINS.scaffold427,AECH.scaffold495,AECH.scaffold818,AECH.scaffold510,AECH.scaffold825,AINS.scaffold468,AECH.scaffold433,AECH.scaffold685,AINS.scaffold48,AECH.scaffold477,AECH.scaffold864,AINS.scaffold46,AECH.scaffold75,AECH.scaffold210,AECH.scaffold692,AINS.scaffold5,AECH.scaffold874,AECH.scaffold829,AINS.scaffold546,AECH.scaffold546,AECH.scaffold594,AECH.scaffold834,AINS.scaffold62,AECH.scaffold45,AINS.scaffold810,AECH.scaffold1499,AECH.scaffold1010,AINS.scaffold82,AECH.scaffold188,AINS.scaffold30,AINS.scaffold329,AINS.scaffold215,AINS.scaffold137,AINS.scaffold81,AINS.scaffold152,AINS.scaffold111,AECH.scaffold405,AINS.scaffold103,AECH.scaffold485,AINS.scaffold120,AECH.scaffold15,AINS.scaffold306,AECH.scaffold102,AINS.scaffold504,AECH.scaffold17,AINS.scaffold13,AECH.scaffold447,AINS.scaffold16,AECH.scaffold730,AINS.scaffold339,AECH.scaffold55,AINS.scaffold185,AECH.scaffold320,AINS.scaffold187,AECH.scaffold31,AINS.scaffold173,AECH.scaffold754,AINS.scaffold217,AECH.scaffold609,AINS.scaffold228,AECH.scaffold204,AINS.scaffold22,AECH.scaffold366,AINS.scaffold23,AECH.scaffold115,AINS.scaffold194,AECH.scaffold180,AINS.scaffold247,AECH.scaffold943,AINS.scaffold298,AECH.scaffold788,AINS.scaffold308,AECH.scaffold425,AECH.scaffold746,AINS.scaffold322,AECH.scaffold252,AINS.scaffold328,AECH.scaffold154,AINS.scaffold498,AECH.scaffold946,AINS.scaffold444,AECH.scaffold261,AINS.scaffold359,AECH.scaffold491,AINS.scaffold355,AECH.scaffold386,AINS.scaffold371,AECH.scaffold233,AINS.scaffold378,AECH.scaffold465,AINS.scaffold399,AECH.scaffold436,AINS.scaffold347,AECH.scaffold398,AINS.scaffold421,AECH.scaffold369,AINS.scaffold422,AECH.scaffold20,AINS.scaffold442,AECH.scaffold209,AINS.scaffold418,AECH.scaffold199,AINS.scaffold17,AECH.scaffold378,AINS.scaffold456,AECH.scaffold610,AINS.scaffold476,AECH.scaffold611,AINS.scaffold164,AECH.scaffold434,AINS.scaffold97,AECH.scaffold191,AINS.scaffold701,AECH.scaffold130,AINS.scaffold557,AECH.scaffold681,AINS.scaffold691,AECH.scaffold827,AINS.scaffold544,AECH.scaffold655,AINS.scaffold729,AECH.scaffold758,AINS.scaffold735,AECH.scaffold111,AINS.scaffold823,AECH.scaffold574,AINS.scaffold171,AECH.scaffold101,AINS.scaffold93,AECH.scaffold155,AINS.scaffold585,AECH.scaffold329,AINS.scaffold188,AECH.scaffold311,AINS.scaffold295,AECH.scaffold1055,AINS.scaffold365,AECH.scaffold2,AINS.scaffold387,AECH.scaffold14,AINS.scaffold56,AECH.scaffold33

# testruns




 bioawk -c fastx '{print $name,length($seq)}' ../$species2.fna|sort -k2 -nr|head -n 20|cut -f 1|egrep "^$" -v|sed -e "s/^/$species2./g" > Aech.top10.lst
 bioawk -c fastx '{print $name,length($seq)}' ../$species1.fna|sort -k2 -nr|head -n 30|cut -f 1|egrep "^$" -v|sed -e "s/^/$species1./g" > Ains.top10.lst

 grep -wf Aech.top10.lst ../SyntenicRibbons.conf|grep -wf Ains.top10.lst -> SyntenicRibbons.conf

 grep "AECH" ../karyotype.conf| grep -wf Aech.top10.lst - > karyotype.tmp1
 grep "AINS" ../karyotype.conf|grep -wf Ains.top10.lst - > karyotype.tmp2

 sort karyotype.tmp1 |uniq|awk '{print $3}' |tr "\n" "," |sed 's/.$//' |awk '{print "~/software/circos-tools-0.23/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - "$0" -static_rx "$0 }' |less

chromosomes_order = AECH.scaffold121,AINS.scaffold39,AINS.scaffold190,AECH.scaffold159,AINS.scaffold178,AINS.scaffold588,AECH.scaffold232,AINS.scaffold684,AINS.scaffold406,AINS.scaffold301,AECH.scaffold326,AECH.scaffold192,AINS.scaffold45,AECH.scaffold88,AECH.scaffold50,AECH.scaffold140,AINS.scaffold349,AECH.scaffold142,AINS.scaffold345,AECH.scaffold220,AINS.scaffold150,AECH.scaffold244,AINS.scaffold116,AECH.scaffold316,AINS.scaffold261,AECH.scaffold35,AINS.scaffold220,AECH.scaffold585,AINS.scaffold248,AECH.scaffold89,AINS.scaffold162,AECH.scaffold99,AINS.scaffold63

 cat karyotype.tmp1 karyotype.tmp2 > karyotype.conf
 circos
