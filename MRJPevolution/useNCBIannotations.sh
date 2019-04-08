cd /Users/lukas/sciebo/inquilineGenomics18/geneFamilies/specificFamilies/mrjps/Atta_cephalotes.MRJPs
# download sequence report
#https://www.ncbi.nlm.nih.gov/assembly#/def
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/143/395/GCF_000143395.1_Attacep1.0/GCF_000143395.1_Attacep1.0_assembly_report.txt
# download gff3
#https://www.ncbi.nlm.nih.gov/genome/?term=txid12957[Organism:noexp]
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/143/395/GCF_000143395.1_Attacep1.0/GCF_000143395.1_Attacep1.0_genomic.gff.gz

# rewrite gff
perl ~/sciebo/inquilineGenomics18/code/MRJPevolution/NCBI2BGI.scfID.pl GCF_000143395.1_Attacep1.0_assembly_report.txt GCF_000143395.1_Attacep1.0_genomic.gff > GCF_000143395.1_Attacep1.0_genomic.BGIids.gff3
