
#######################################
# run GAG for attines

#######################################
base=/corefac/cse/lukas/inqGen18/QC/7ants_GAG
base=/usr/local/home/lschrader/data/inqGen18/QC/filterContaminants


export genomeBase=~/data/genomes/reannotations_of_7_ants/
for gffFile in $(readlink -f  $genomeBase/*/annotation/gene_annotation/*.v2.0.gff)
do
  species=$(echo $gffFile|perl -pe 's/.*\/(.).*\_(...).*/$1$2/g')
  basepath=$(echo $gffFile|perl -pe 's/(.*)\/annotation\/.*/$1/g')
  genome=$(readlink -f $basepath//genome/*.v2.0.fa)
  out=$(echo $gffFile|rev|cut -f 1 -d "/" |rev)
  # clean up gff
  SPECIES=$(echo $species| tr a-z A-Z)
  Rscript addGene2gff.R $gffFile $species.tmp.gff $SPECIES
  # check which genes are removed
  GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --insert-missing= $species.tmp.gff
  ~/software/genometools-1.5.9/bin/gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -addintrons $species.tmp_clean.gff > $species.v2.0.gff3
  python ~/software/GAG/gag.py -f $genome -g $species.v2.0.gff3
  mkdir $species
  mv gag_output $species/$species.gag
  mv $species.* $species
done


for gffFile in $(readlink -f  $genomeBase/*/annotation/gene_annotation/*.v2.0.gff)
do
  species=$(echo $gffFile|perl -pe 's/.*\/(.).*\_(...).*/$1$2/g')
  egrep "Total sequence length|Number of genes|Number of mRNAs|Number of CDS|Overlapping genes|Contained genes|Total gene length|Total CDS length|Shortest gene|Longest gene|mean gene length|% of genome covered by genes" $species/$species.gag/genome.stats > $species.gag.stats
  echo -e " \t $species" > $species.gag.tmp
  cat $species.gag.stats |cut -c 34- >> $species.gag.tmp
done

echo "Species|Total sequence length|Number of genes|Number of mRNAs|Number of CDS|Overlapping genes|Contained genes|Total gene length|Total CDS length|Shortest gene|Longest gene|mean gene length|% of genome covered by genes"|tr "\|" "\n" > gag.tmp
paste gag.tmp *.gag.tmp > 7attines.gag.stats
