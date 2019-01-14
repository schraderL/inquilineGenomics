###################################################
#==================================================
# 0. Define reusable variables
#==================================================
###################################################



#Parg
export species=Parg
export genomeFolder=/usr/local/home/lschrader/data/genomes/inquilines/v1.1/Parg
export genomeFa=$genomeFolder"/Pseudoatta_argentina.fa"
export transcript=~/data/genomes/inquilines/v1.1/Parg/Pseudoatta_argentina.transcripts.gtf

#Folders
export base=/usr/local/home/lschrader/data/OR
export scripts=$base/scripts
export ORfa=$base/AechORs.McKenzie.fasta
export ORgff=$base/SupplementaryFile4.AechOrs.final.gff
export queryGenome=/usr/local/home/lschrader/data/genomes/Aech/Aech_2.0_scaffolds.clean.fa


###################################################
#==================================================
# 1. Annotate OR genes with Zhou et al. pipeline.
#==================================================
###################################################

bash $scripts/runExonerate.sh


###################################################
#==================================================
# 2. Annotate OR genes with GeMoMa
#==================================================
###################################################

bash $scripts/runGeMoMa.sh

###################################################
#==================================================
# 3. Blast single exons against target
#==================================================
###################################################

bash $scripts/ExonBlast.sh

#############################################
#===========================================
#	4. COMBINE EVIDENCE WITH EVM
#===========================================
#############################################

bash $scripts/runEVM.sh



#run
#~/software/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log

# run in parallel
#split --number=l/6 ${fspec} xyzzy
split --lines=130 commands.list
for f in xa*
do
 ~/software/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl $f > run.$f.log &
done

#wait till all runs have finished
rm xa*

#combine
~/software/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

#write gffs
~/software/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome $genomeFa


#write protein fasta

for f in */evm.out.gff3
do
if [ -s $f ] 
	then
	echo "~/software/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl  $f $genomeFa > $f.prot.fa" >> getFasta.sh
fi	
done

split --lines=10 getFasta.sh
for f in xa*
do
 bash $f &
done

cat */evm.out.gff3.prot.fa > $species.final.OR.fa

rm -rf tmpGff/
mkdir tmpGff
for f in */evm.out.gff3
do
	scf=$(echo $f | cut -f 1 -d "/")
	~/software/genometools-1.5.9/bin/gt gff3 -sort -tidy -retainids -checkids -fixregionboundaries -addintrons $f > tmpGff/$scf.gff
done


~/software/genometools-1.5.9/bin/gt merge -tidy -retainids tmpGff/*.gff > $species.final.OR.gff3



mkdir ../final
cp $species.final.OR.gff3 ../final
cp $species.final.OR.fa ../final
cp $species.prot.gff ../final
cp OR.transcripts.gff3 ../final/$species.OR.transcripts.gff3
cp $species.OR.scf.fa ../final
cp $species.GeMoMa.final.gff ../final/
cp $species.exonerate.final.gff ../final/
cp ../GeMoMa/predicted_protein.fasta ../final/$species.GeMoMa.fa
cp ../exonerate/$species.exonerate.pep ../final/
cp Ahey.OR.exons.bls.sorted.gff ../final/




############################################################
# 4. Scan for domains with pfam_scan.pl
############################################################
mkdir $base/$species/pfam
cd $base/$species/pfam
#perl ~/software/PfamScan/pfam_scan.pl -fasta ../GeMoMa/predicted_protein.fasta -dir ~/data/pfam/ -outfile $species.GeMoMa.domains.out -cpu 20
perl ~/software/PfamScan/pfam_scan.pl -fasta ../EVM/$species.final.OR.fa -dir ~/data/pfam/ -outfile $species.final.domains.out -cpu 20
perl ~/software/PfamScan/pfam_scan.pl -fasta ../exonerate/$species.exonerate.pep -dir ~/data/pfam/ -outfile $species.exonerate.domains.out -cpu 20
grep 7tm_6 $species.final.domains.out |wc -l