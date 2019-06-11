base=/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/MRJPs/Amel.reference
cd $base
mkdir gb
esearch -db protein -query "(yellow OR MRJP) AND (Apis mellifera[organism]) AND srcdb+refseq[prop]" < /dev/null|elink -target gene|efetch -format tabular > Amel.y.tsv

for i in $(cut -f 3 Amel.y.tsv|sed 1d); do
  esearch -db protein -query "($i) AND (Apis mellifera[organism]) AND srcdb+refseq[prop]" < /dev/null|efetch -format gb > $base/gb/Amel.$i.gb

for i in $(cut -f 3 Amel.y.tsv|sed 1d); do
# convert gb entry to fasta with gene="<>" tag and all splice forms
cat gb/Amel.$i.gb|\
awk '/^ACCESSION   / {printf(">%s\t",$2);next;} \
     /[[:space:]]*\/gene\=/ {{match($0,"gene=.*",a)}{print a[0];next;}} \
     /^ORIGIN/ {inseq=1;next;} /^\/\// {inseq=0;} {if(inseq==0) next;\
      gsub(/[0-9 ]/,"",$0); printf("%s\n",toupper($0));}' > fa/Amel.$i.fa

# select longest isoform from all splice forms
cat fa/Amel.$i.fa|\
awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' |\
awk -F '\t'  '{printf("%s\t%d\n",$0,length($3));}' |\
sort -k2,2 -k4,4nr |\
sort -k2,2 -u -s |\
awk '{print $1";symbol="$2";length="$4"\n"$3}'> lfa/Amel.$i.longestIsoform.fa
done
cd $base/aln

cat ../lfa/* > Amel.fa
prank Amel.fa

#/usr/local/home/lschrader/data/inqGen18/geneFamilies/specificFamilies/MRJPs/Amel.reference/aln/output.best.fas
