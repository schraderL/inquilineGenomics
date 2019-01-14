#!/usr/bin/env Rscript
# usage: addGene2gff.R infile.gff outfile.gff  <gene abbreviation (e.g. AINS)>
args = commandArgs(trailingOnly=TRUE)
a<-read.table(args[1])

#a<-read.table("/Users/lukas/sciebo/inquilineGenomics18/QC/contaminants/tmp.gff")
b<-split(a,a$V3)
gene<-b$mRNA
gene$V3<-"gene"
gene$V6<-"."
abbr<-args[3]
genPattern<-paste("(",abbr,"[0-9]{5})",sep="")
b<-lapply(b, function(x){
  x$V9<-gsub(genPattern,"\\1-mRNA",x$V9,perl=T)
  return(x)
})

b$mRNA$V9<-gsub("((....[0-9]{5})-mRNA)","\\1;Parent=\\2",b$mRNA$V9,perl=T)
b$gene<-gene
gff<-do.call("rbind", b)
gff<-gff[with(gff, order(gff$V1, gff$V4)),]
write.table(gff,args[2],sep='\t',row.names=F,quote=F,col.names=F)

