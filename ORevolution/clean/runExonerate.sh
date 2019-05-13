###################################################
#==================================================
# 1. Annotate OR genes with Zhou et al. pipeline.
#==================================================
###################################################

######################################################################################################
#1/ run the perl script identify_geneLS.pl on the target genome using known genes as queries
######################################################################################################

if [ ! -f $genomeFa.nhr ];
then
	cd $genomeFolder
	makeblastdb -in $genomeFa -dbtype nucl
fi
cd $base/$species
#edited config to eval 1e-4 for Ahey.test
cat ../config_raw |envsubst > config_$species

perl $scripts/identify_geneLS.pl config_$species $species.bls > $species.out

######################################################################################################
#2/ run the perl script "hsp2bed.pl" to get loci
######################################################################################################

# 150 is the aa cutoff
# bioawk -c fastx '{ print $name, length($seq) }' < $ORfa | sort -k 2 -nr|tail -n1
perl $scripts/hsp2bed.pl $species.out 100 > $species.bed


######################################################################################################
#3/ Prepare exonerate input
######################################################################################################

cat ../config_raw2 |envsubst > config_$species"2"
perl $scripts/runExonerate.pl config_$species"2"

cd $base/$species

######################################################################################################
#4/ Run exonerate on each locus
######################################################################################################
#AecOr291 exonerate
#f=AechFa/AechOr100PSE.fa
#f=AechFa/AechOr21.fa
for f in $species.Loci/*.fa
do
	gene=$(echo $f | cut -f 2 -d "/"|perl -pe 's/.target.fa//g')
	echo "working on $gene"
	query=$(echo "ORFa/"$gene|perl -pe 's/\.\d+$/.fa/g')

	#target=$species".Loci/$gene.target.fa"
	target=$f
	nice exonerate --query $query --target $target -m protein2genome --showtargetgff yes --showalignment no --ryo '###---FASTASTART---###\n>%qi length=%tl alnlen=%tal target=%ti_(%tcb-%tce)\n%tcs###---FASTAEND---###\n' > exonerate/raw/$gene.exonerate
	perl $scripts/processExonerate2GFF.pl exonerate/raw/$gene.exonerate | sed s/cds/CDS/g > dirty.gff
	perl $scripts/cleanExonerateGFF.pl dirty.gff > tmp.gff
	rm dirty.gff
	GFFcleaner --add-missing-ids --add-exon-ids --insert-missing= tmp.gff
	rm tmp.gff
	gt gff3 -sort -tidy  -checkids yes -retainids yes tmp_clean.gff > exonerate/gff/$gene.gff
	perl $EVMtools/gff3_gene_prediction_file_validator.pl exonerate/gff/$gene.gff
	perl $scripts/processExonerate2Fa.pl exonerate/raw/$gene.exonerate > exonerate/fa/$gene.exonerate.fa
	perl $scripts/translateCDS2AA.pl exonerate/fa/$gene.exonerate.fa > exonerate/protFa/$gene.exonerate.pep
	perl $EVMtools/misc/exonerate_gff_to_alignment_gff3.pl exonerate/raw/$gene.exonerate > exonerate/PROTEIN/$gene.PROTEIN.gff
	echo "$gene done. :)"
done

cat exonerate/protFa/* > exonerate/$species.exonerate.pep

gt merge -tidy -retainids $base/$species/exonerate/gff/*.gff  > $base/$species/exonerate/$species.exonerate.gff
cat $base/$species/exonerate/PROTEIN/*.PROTEIN.gff  > $base/$species/exonerate/$species.exonerate.PROTEIN.gff
perl $EVMtools/gff3_gene_prediction_file_validator.pl $base/$species/exonerate/$species.exonerate.gff
perl $EVMtools/gff3_gene_prediction_file_validator.pl $base/$species/exonerate/$species.exonerate.PROTEIN.gff
