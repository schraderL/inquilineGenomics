cd /Users/lukas/CSE/inquilineGenomics/ORevo/Ae.antennae.RNAseq/ORsubset

awk -F $'\t' 'BEGIN {OFS = FS} {if ($3=="gene")
{
  if ($4>40000){
    $4=$4-40000;
  }else{
    $4=1}
  $5=$5+40000;
  print $0;
}
}' ../../data/McKenzie/SupplementaryFile4.AechOrs.final.gff > OR.loci.gff

bedtools intersect -wa -b OR.loci.gff -a ../cuffcmp.combined.gtf > OR.cuffcmp.gtf
cd /Users/lukas/CSE/inquilineGenomics/ORevo/data/McKenzie
awk '{if ($3=="gene") print $1":"$4-2000".."$5+2000,$9}' SupplementaryFile4.AechOrs.final.gff > ~/CSE/inquilineGenomics/ORevo/Ae.antennae.RNAseq/ORsubset/SupplementaryFile4.AechOrs.final.list


# Manually corrected the Aech annotations
# load OR.cuffcmp.gtf as webapollo feature track
open /Users/lukas/CSE/inquilineGenomics/ORevo/Ae.antennae.RNAseq/ORsubset/SupplementaryFile4.AechOrs.final.list
#http://hymenopteragenome.org/Apollo2/jbrowse/index.html?&organism=1683838

Canonical splice site:
3' CG-GT
5' AG-GA


#- Legend for OR annotation

#NTE	Missing sequence at N terminus
#INT	Missing sequence in the middle of gene
#CTE	Missing sequence at C terminus
#NI	Missing N terminus and section in the middle of gene
#NC	Missing N and C terminus
#IC	Missing section in the middle of gene and C terminus
#PSE	Pseudogene
#P+N/I/C	Pseudogene and missing sequence
#(F)	Could be functional with assembly-introduced false frameshift
#(S)	Could be functional with non-canonical splice sites
