################################################################################
# RUN DNDS ANALYSES ON OR CLADES
################################################################################
base=/usr/local/home/lschrader/data/inqGen18/OR.SOS
software=/usr/local/home/lschrader/software
################################################################################
# TAG Run ABSREL:
################################################################################
cd $base
mkdir $base/cdsAln
mkdir $base/absREL
mkdir $base/RELAX

cat ../ORs/*/final/finalSet/gag_output/genome.proteins.fasta|perl -pe 's/\>.*Parent=(.*)\|.*/>$1/g' > $base/fullFasta/pep.fa
cat ../ORs/*/final/finalSet/gag_output/genome.mrna.fasta|perl -pe 's/\>.*Parent=(.*)\|.*/>$1/g' > $base/fullFasta/cds.fa

################################################################################
# TAG Prepare alignments and trees for All
################################################################################
target=all
datafolder=$base/data/$target
mkdir $base/cdsAln/$target
mkdir $base/absREL/$target


# get subsetted fastas for each clade, only remove dH, NC, NTE, and CTE models
ls $datafolder/*.tsv|parallel --nice 10 "cut -f 1 {}|sed 1d |egrep -v 'dH|NTE|CTE|NC' > {.}.lst"
ls $datafolder/*.lst|parallel --nice 10 "faSomeRecords $base/fullFasta/cds.fa {} {.}.cds.fa"
ls $datafolder/*.lst|parallel --nice 10 "faSomeRecords $base/fullFasta/pep.fa {} {.}.pep.fa"
# align complete gene models with prank
ls $datafolder/*.pep.fa|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "prank  -d={}.pep.fa -t={}.tre -o={.}.pep.aln -prunetree -showtree -once"
# create cdsAln folders for each clade
ls $datafolder/*.pep.fa|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "mkdir $base/cdsAln/$target/{/.}"
# reverse translate to CDS alignments
ls $datafolder/*.pep.fa|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "perl ~/software/pal2nal.v14/pal2nal.pl  {}.pep.aln.best.fas {.}.cds.fa -nogap -nomismatch -output fasta > $base/cdsAln/$target/{/.}/{/.}.aln"
# A-Clade-033 and XA-Clade-011 only contain fragments
################################################################################
# TAG Run absREL on all
################################################################################
# create absREL folders for each clade
ls $datafolder/*.pep.fa|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "mkdir $base/absREL/$target/{/.}"
# remove empty files
find $datafolder/ -size 0 -delete
# softlink alignment to absREL folder
ls $datafolder/*.pep.fa|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "ln -s $base/cdsAln/$target/{/.}/{/.}.aln $base/absREL/$target/{/.}/{/.}.aln"

HMdir=/corefac/cse/lukas/software/hyphy-2.3.14
bf=$HMdir/res/TemplateBatchFiles/SelectionAnalyses/aBSREL.bf
ls $base/absREL/$target/*/*.aln|parallel -j2 --nice 10 "HYPHYMP CPU=20 LIBPATH=$HMdir/res/ $bf 'Universal' '{}' '$datafolder/{/.}.pep.aln.best.dnd' 'All' '{}.absREL.out'"
################################################################################
# TAG Run RELAX on all
################################################################################
target=all
datafolder=$base/data/$target
mkdir $base/RELAX/$target
# copy alignments to RELAX folder
ls $datafolder/*.pep.fa|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "mkdir $base/RELAX/$target/{/.}"
ls $datafolder/*.pep.fa|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "ln -s $base/cdsAln/$target/{/.}/{/.}.aln $base/RELAX/$target/{/.}/{/.}.aln"
# Create Tree with {T} and {R} tags. T=All parasite leaves, R=All leaf-cutter leaves
ls $datafolder/*pep.aln.best.dnd|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "cat {}.pep.aln.best.dnd| perl -pe 's/(PargOr-[0-9][0-9][0-9]\:[0-9]\.[0-9]*)|(AinsOr-[0-9][0-9][0-9]\:[0-9]\.[0-9]*)|(AchaOr-[0-9][0-9][0-9]\:[0-9]\.[0-9]*)/\$1\$2\$3\{T\}/g'| perl -pe 's/(AcolOr-[0-9][0-9][0-9]\:[0-9]\.[0-9]*)|(AcepOr-[0-9][0-9][0-9]\:[0-9]\.[0-9]*)|(AheyOr-[0-9][0-9][0-9]\:[0-9]\.[0-9]*)|(AechOr-[0-9][0-9][0-9]\:[0-9]\.[0-9]*)/\$1\$2\$3\$4\{R\}/g' > $base/RELAX/$target/{/.}/{/.}.relax.tre"

HMdir=/corefac/cse/lukas/software/hyphy-2.3.14
bf=$HMdir/res/TemplateBatchFiles/SelectionAnalyses/RELAX.bf
ls $base/RELAX/$target/*/*.aln|parallel --nice 10 -j2 "HYPHYMP CPU=20 LIBPATH=$HMdir/res/ $bf 'Universal' '{}' '{.}.relax.tre' 'T' 'R' 'Minimal'"

################################################################################
# TAG Run RELAX (full model) on all
################################################################################
target=all
datafolder=$base/data/$target
mkdir $base/RELAX.full/
mkdir $base/RELAX.full/$target
# copy alignments to RELAX folder
ls $datafolder/*.pep.fa|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "mkdir $base/RELAX.full/$target/{/.}"
ls $datafolder/*.pep.fa|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "ln -s $base/cdsAln/$target/{/.}/{/.}.aln $base/RELAX.full/$target/{/.}/{/.}.aln"
# Create Tree with {T} and {R} tags. T=All parasite leaves, R=All leaf-cutter leaves
ls $base/RELAX/$target/*/*.relax.tre|perl -pe 's/(.*\/.*?)\..*/$1/g'|parallel --nice 10 "ln -s {}.relax.tre $base/RELAX.full/$target/{/.}/{/.}.relax.tre"

HMdir=/corefac/cse/lukas/software/hyphy-2.3.14
bf=$HMdir/res/TemplateBatchFiles/SelectionAnalyses/RELAX.bf
ls $base/RELAX/$target/*/*.aln|parallel --nice 10 -j2 "HYPHYMP CPU=20 LIBPATH=$HMdir/res/ $bf 'Universal' '{}' '{.}.relax.tre' 'T' 'R' 'All'"
# results strored e.g. here ~/data/inqGen18/OR.SOS/RELAX/all/E-Clade-030
