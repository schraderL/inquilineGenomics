pV=30
cV=0.1  #contig threshold
rV=0.1  #region threshold
iV=100  #intron gain/loss penalty
mV=500 #intron length
rG=0 #relative score filter
cbfG=0.75  #common border filter
#ct0.1 p10 r0.1 i100 m500
#ct0.1 p10 r0.1 i100 m6000


java -jar /usr/local/home/lschrader/software/GeMoMa_1.3/GeMoMa-1.3.jar CLI GeMoMa t=tblastn.txt c=./cds-parts.fasta a=./assignment.tabular tg=/usr/local/home/lschrader/data/genomes/attines/assembly/Acep.genome.fa  outdir=./ p=$pV ct=$cV rt=$rV intron-loss-gain-penalty=$iV m=$mV
mv predicted_annotation.gff GeMomA.p.$pV.ct.$cV.rt.$rV.i.$iV.m.$mV.gff3
java -jar /usr/local/home/lschrader/software/GeMoMa_1.3/GeMoMa-1.3.jar CLI GAF g=GeMomA.p.$pV.ct.$cV.rt.$rV.i.$iV.m.$mV.gff3 r=$rG cbf=$cbfG
mv filtered_predictions.gff GAF.r.$rG.cbf.$cbfG-p.$pV.ct.$cV.rt.$rV.i.$iV.m.$mV.gff3

refilter GeMomA.p.30.ct.0.1.rt.0.4.i.100.m.5000.gff3. It has basically all the good models. Changes from GAF.r.0.5.cbf.1-p.30.ct.0.1.rt.0.4.i.100.m.5000.gff3: r=0.3,cbf=default

#mv filtered_predictions.gff GAF.r0.1-p.$pV.ct.$cV.rt.$rV.i.$iV.m.$mV.gff3
GAF.r0.2-p.$pV.ct.$cV.rt.$rV.i.$iV.m.$mV.gff3

#looks great!!
GAF.r0.2-p.30.ct.0.2.rt.0.1.i.100.m.5000.gff3
  GeMomA.p.30.ct.0.1.rt.0.4.i.100.m.5000.gff3
GAF.r.0-p.30.ct.0.1.rt.0.1.i.100.m.4000.gff3
GAF.r.0.incomplete-p.30.ct.0.1.rt.0.1.i.100.m.4000.gff3
# now try with shorter introns
# some small models are picked up that are missed with m=4000 at m=2000

combining m=2000 and m=4000 results doesn't help.
java -jar /usr/local/home/lschrader/software/GeMoMa_1.3/GeMoMa-1.3.jar CLI GAF g=GeMomA.p.30.ct.0.1.rt.0.1.i.100.m.4000.gff3 r=-1000
g=GeMomA.p.30.ct.0.1.rt.0.1.i.100.m.2000.gff3  cbf=$cbfG
# gene missed at
scaffold00043:1,758,682-1,774,960
# still missing entire scf00058

# almost perfect results (i.e. overlap with manual annotation)
# The following commands give pretty good results for Aech annotating Acep.
java -jar $GeMoMaJar CLI GeMoMa t=${out}/tblastn.txt c=${out}/cds-parts.fasta a=${out}/assignment.tabular tg=${target}  outdir=${out} p=30 ct=0.1 rt=0.1 intron-loss-gain-penalty=100 m=4000
java -jar /usr/local/home/lschrader/software/GeMoMa_1.3/GeMoMa-1.3.jar CLI GAF g=${out}/predicted_annotation.gff r=-10000 c=F

# I am currently running the same command for Sinv and Acol manual annotations to annotate Acep ORs. I will then have to check whether the models are similarly good.
# The challenge will be to run EVM on all the evidence that we have. Now, since the GeMoMa models are filtered and of such high quality, it is likely that I can give them a very high weight compared to exonerate. The blast evidence is also a nice addition. RNAseq data has proven to be difficult because it can lead to wrong predictions by forcing the incorporation of an unrelated transcript at the same locus (e.g. from the opposite strand)

# almost perfect results (i.e. overlap with manual annotation)
java -jar /usr/local/home/lschrader/software/GeMoMa_1.3/GeMoMa-1.3.jar CLI GAF g=GeMomA.p.30.ct.0.1.rt.0.1.i.100.m.4000.gff3 r=-10000 c=F
