# cat all ORs that contain the domain or a fragment of the domain
cat \
 ../Aech/20180727/finalSet/Aech.OR.fa \
 ../Ahey/180727/finalSet/Ahey.OR.fa \
 ../Ains/20180718/finalSet/Ains.OR.fa \
 ../Acep/180729/Acep.OR.fa \
 ../Acha/180729/Acha.OR.fa \
 ../Acol/180730/Acol.OR.fa \
 ../Parg/180728/Parg.OR.fa \
 |tr "\n" "\t" |tr ">" "\n"|grep "dH" -v|tr "\n" ">"|tr "\t" "\n" |sed "s/^>$//g" > OR.allwith.6tm7.fa

 prank  -d=OR.allwith.6tm7.fa -support

 #/usr/local/raxml/latest/raxmlHPC-PTHREADS -T 2 -f a -x 12345 -m PROTGAMMAGTR -s pastajob.marker001.interesting_seqs.pep.aln -# 100 -p 2 -n COX2alignment
