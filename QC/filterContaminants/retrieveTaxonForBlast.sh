#!/usr/bin/bash
#usage filterContaminants.sh <blast.output.tab> <species name>
# setup sqlite3 db
# sqlite3 gi2tax.db < import.sql
# index database

# iterate over each blast hit
while read p; do
   # retrieve gi from blast result
   gi=$(echo $p|grep -E "gi\|[0-9]+"  -o|cut -f 2 -d "|")
   # retrieve taxonID from gi
   taxid=$(sqlite3 gi2tax.db  "select blah from gi2tax INDEXED BY gi2tax_idx where name='${gi}' LIMIT 1 ";)
   # retrieve lineage from taxid
   lineage=$(curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=${taxid}&retmode=xml" |\
   grep Lineage |\
   cut -d '>' -f 2 |\
   cut -d '<' -f 1 |\
   sed s/'; '/';'/g|sed s/' '/'^'/g)
   # print blast hit and lineage info
   echo -ne $p $taxid $lineage "\n"|awk -v OFS="\t" '$1=$1'|sed s/'\^'/' '/g
done < $1 > results/$2.blast.taxon.tsv

#compile list of taxons found
grep -o '[^;]*$' results/$2.blast.taxon.tsv -o|sort|uniq -c > $2.taxon.list

#filter out all bacterial scaffolds
#awk -F $'\t' '{if ($14 ~/Bacteria/) print $1,$14}' results/$2.blast.taxon.tsv |sort -u -k1,1 > results/$2.contaminants.list
#awk -F $'\t' '{if ($14 ~/Bacteria/) print $1,$14}' results/$2.blast.taxon.tsv |sort > results/$2.contaminants.detailed.list

#cross-check scaffolds hitting both bacterial and others
#grep -f results/$2.contaminants.list results/$2.blast.taxon.tsv |grep -v Bacteria
