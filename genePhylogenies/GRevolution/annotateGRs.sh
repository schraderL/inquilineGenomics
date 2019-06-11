######################
#
# gustatory receptors
#
######################

##################################################################################################
# Retrieve all annotated GRs in ants from NCBI
#################################################
export software=/usr/local/home/lschrader/software
export base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/grs
mkdir $base
cd $base

# retrieve all ant genes annotated as GRs
esearch -db protein -query "(IPR013604 OR PF08395 OR (gustatory AND receptor)) AND txid34695[Organism:exp] AND refseq[filter]" | elink -target gene|efetch -format tabular > ant.GRs.tsv
# get details about completeness of gene models and gene lengths
esearch -db protein -query "(IPR013604 OR PF08395 OR (gustatory AND receptor)) AND txid34695[Organism:exp] AND refseq[filter]" | efetch -format gb|egrep "\/gene=.*|COMPLETENESS.*|LOCUS.*|   ORGANISM.*" |tr "\n" "\t" |perl -pe 's/(gene\=.*?)\t/$1\n/g' |perl -pe 's/  +/\t/g'|cut -f 2,3,7,9|perl -pe "s/ aa//g"|perl -pe 's/\/gene=\"(.+?)\"/$1/g'> ant.GRs.details.tsv

# select only complete  GRs between 50-500 aa
cat ant.GRs.details.tsv|grep "full length"|awk '{if ($2<=500 && $2>=50) print $0}'|cut -f 4 > ant.GRs.selection.lst
grep -f ant.GRs.selection.lst ant.GRs.tsv > ant.GRs.selection.tsv


# Download all ant genome gffs from NCBI
https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=NCID_1_5317497_130.14.22.76_5555_1544272010_2902454397_0MetA0_S_HStore&QueryKey=13&ReleaseType=RefSeq&FileType=genomic.gff.gz&Flat=true
# Download all ant genome fastas from NCBI
https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=NCID_1_5317497_130.14.22.76_5555_1544272010_2902454397_0MetA0_S_HStore&QueryKey=13&ReleaseType=RefSeq&FileType=genomic.fna.gz&Flat=true

data=/usr/local/home/lschrader/data/genomes/NCBI/ncbi-genomes-2018-12-08
cd $base

# retrieve all accessions for ants
for i in $(ls $data|cut -f 1-2 -d "_" |grep GCF|uniq)
do
  echo -ne "$i\t"
  esearch -query $i -db assembly|efetch -format docsum|egrep Organism|perl -pe 's/.*\>(.*?)\<.*/$1/g'
done > ant.gcf.tsv

# prepare subsetted gffs for all species with GRs annotated
mkdir $base/data/
#get species: cat ant.gcf.tsv |perl -pe 's/.*\t(.*? .*?) .*/$1/g'
for i in $( cut -f 1 ant.gcf.tsv)
do
  genus=$(grep $i ant.gcf.tsv |perl -pe 's/.*\t(.*?) .*? .*/$1/g')
  species=$(grep $i ant.gcf.tsv |perl -pe 's/.*\t.*? (.*?) .*/$1/g')
  genes=$(grep "$genus $species" ant.GRs.tsv|cut -f 6|tr "\n" "|"|sed 's/.$//')
  if ((${#genes} > 0))
  then
    zcat $data"/"$i*gff*| egrep "$genes " > $base/data/$i.$genus.$species.GRs.gff
    echo "$genes are for $genus $species"
  fi
done

# A. charruanus
target_genome1=/corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_charruanus.2.1/genome/Acromyrmex_charruanus.v2.1.fa
target_genome2=/corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_heyeri.2.1/genome/Acromyrmex_heyeri.v2.1.fa
target_genome3=/corefac/cse/lukas/inqGen18/inquilines_v2.1/Pseudoatta_argentina.2.1/genome/Pseudoatta_argentina.v2.1.fa
target_genome4=/corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_insinuator.2.1/genome/Acromyrmex_insinuator.v2.1.fa
target_genome5=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/genome/Acromyrmex_echinatior.v2.0.fa
target_genome6=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_colombica/genome/Atta_colombica.v2.0.fa
target_genome7=/corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_cephalotes/genome/Atta_cephalotes.v2.0.fa

declare -a arr=("/corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_charruanus.2.1/genome/Acromyrmex_charruanus.v2.1.fa" \
                "/corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_heyeri.2.1/genome/Acromyrmex_heyeri.v2.1.fa" \
                "/corefac/cse/lukas/inqGen18/inquilines_v2.1/Pseudoatta_argentina.2.1/genome/Pseudoatta_argentina.v2.1.fa" \
                "/corefac/cse/lukas/inqGen18/inquilines_v2.1/Acromyrmex_insinuator.2.1/genome/Acromyrmex_insinuator.v2.1.fa" \
                "/corefac/cse/lukas/genomes/reannotations_of_7_ants/Acromyrmex_echinatior/genome/Acromyrmex_echinatior.v2.0.fa") #\
                #"/corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_colombica/genome/Atta_colombica.v2.0.fa" \
                #"/corefac/cse/lukas/genomes/reannotations_of_7_ants/Atta_cephalotes/genome/Atta_cephalotes.v2.0.fa")


export software=/usr/local/home/lschrader/software
export base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/grs
export EVMutils=$software/EVidenceModeler-1.1.1/EvmUtils
cd $base
for target_genome in "${arr[@]}"
do
  cd $base
  out=$(echo $target_genome |perl -pe 's/.*\/(.*?)\..*/$1.GRs/g')
  refGenomes=()
  refAnnot=()

  mkdir ./$out
  command=
  for i in $(ls $base/data|cut -f 1-2 -d "."|head -n 10); do
      var1=$(readlink -f $base/data/$i*)
      var2=$(readlink -f $data/$i*genomic.fna)
      refAnnot+=("$var1")
      refGenomes+=("$var2")
      command=$(echo $command s=own a=$var1 g=$var2)
  done

  nice java -jar ~/software/GeMoMa-1.5.3/GeMoMa-1.5.3.jar CLI GeMoMaPipeline threads=10 Extractor.p=T GAF.c=F GeMoMa.prefix=GR t=$target_genome outdir=$out $command 1> $out/gemoma.out 2> $out/gemoma.err

  cd $base/$out

  GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --report-duplicates --insert-missing= filtered_predictions.gff 2>/dev/null
  $software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries filtered_predictions_clean.gff| sed s/prediction/mRNA/g|sed s/GAF/GeMoMa/g> tmp.gff 2>/dev/null
  gt gff3 -sort -tidy  -addids -fixregionboundaries -checkids yes -retainids yes tmp.gff > tmp2.gff 2>/dev/null
  awk '{if ( $3=="gene") print $0 }' tmp.gff > tmp.genes.gff 2>/dev/nullhsa
  gffread tmp.gff -o tmp3.gff --force-exons -E 2>/dev/null
  cat tmp.genes.gff tmp3.gff |$software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -|egrep "^###" -v|sed s/geneID/Parent/g > tmp4.gff
  GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --report-duplicates --insert-missing= tmp4.gff -o GRs.GeMoMa.gff 2>/dev/null
  perl $EVMutils/gff3_gene_prediction_file_validator.pl GRs.GeMoMa.gff
  perl $EVMutils/gff3_file_to_proteins.pl   GRs.GeMoMa.gff $target_genome > GRs.GeMoMa.fa 2>/dev/null
  rm tmp*
  nice perl $software/PfamScan/pfam_scan.pl -fasta GRs.GeMoMa.fa -dir ~/data/pfam/ -outfile GRs.GeMoMa.pfam -cpu 10

  cut -f 1 GRs.GeMoMa.pfam |grep PF08395.10 |cut -f 1 -d " "|sort |uniq > GRs.GeMoMa.PF08395.10.lst
  faSomeRecords GRs.GeMoMa.fa GRs.GeMoMa.PF08395.10.lst GRs.GeMoMa.PF08395.10.fa
done

for folder in $(readlink -f *.GRs)
do
  echo $folder
  cd $folder
  nice perl $software/PfamScan/pfam_scan.pl -fasta GRs.GeMoMa.fa -dir ~/data/pfam/ -outfile GRs.GeMoMa.lowE.pfam -cpu 20  -e_seq 0.01
  cd $base
done

# cleanup and rename

for i in $(readlink -f *.GRs/)
do
  cd $i
  species=$(echo $i|rev|cut -f 1 -d "/"|rev)
  # create replacement table to rename GRs
  cat GRs.GeMoMa.gff |grep "ID=GR"|perl -pe 's/.*ID=(GRRNA.*?)\;.*/$1/g'|sort|uniq| awk '{printf("%s\tGR-%02d\n",$1,NR)}' > replacement.tsv
  seqkit fx2tab GRs.GeMoMa.fa |sed -r 's/GeMoMa.+gene/GeMoMa\tgene/g'|sed -r 's/ +/\t/g'> GRs.GeMoMa.tsv
  while read p; do
    gene=$(echo "$p"|cut -f 1)
    newgene=$(echo -n "$p"|cut -f 2)
    echo -ne "$newgene\t"
    if grep -Fxq $gene GRs.GeMoMa.PF08395.10.lst
    then
      echo -ne "7_TM7\t"
    else
      echo -ne "NONE\t"
    fi
    STARTSTOP=$(grep $gene GRs.GeMoMa.tsv|cut -f 5|perl -pe 's/(.).+(.)/$1$2/g')
    SP=$(grep  $gene GRs.GeMoMa.tsv|cut -f 10)
    echo -ne "$STARTSTOP\t$SP\t"
    grep $gene GRs.GeMoMa.tsv
  done <replacement.tsv > $species.GeMoMa.details.tsv

  sed -r 's|(\S+)\s*(\S+).*|s/\1\/\2/g|' replacement.tsv | sed -f - GRs.GeMoMa.gff > $species.GeMoMa.gff
  sed -r 's|(\S+)\s*(\S+).*|s/\1\/\2/g|' replacement.tsv | sed -f - GRs.GeMoMa.fa > $species.GeMoMa.fa
cd $base
done

for i in $(readlink -f *.GRs/)
do
  echo -ne "$i\t"
  awk '{if($2=="7_TM7" && $3=="M*") print $0}' $i/*.GRs.GeMoMa.details.tsv|wc -l
done
## see https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002930
for i in $(readlink -f *.GRs/)
do
  echo -ne "$i\t"
  awk '{if($3=="M*") print $0}' $i/*.GRs.GeMoMa.details.tsv|wc -l
done
