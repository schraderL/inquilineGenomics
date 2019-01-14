#!/usr/bin/env Rscript
# Rscript analyseAln.R ~/data/inqGen18/phylogeny/SCO/4d/OGs/ "*.4d$" test
args = commandArgs(trailingOnly=TRUE)
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#install.packages("msaR")
#install.packages("AlignStat")
#install.packages("strataG")
library(msaR)
library(strataG)
#path<-"~/data/inqGen18/phylogeny/SCO/4d/OGs/"
path<-args[1]
#file.names<-list.files(path,pattern="*.4d$")
file.names<-list.files(args[1],args[2])

aln<-list()
aln2<-list()
titvs<-list()
dists<-list()
tits<-NA
mdists<-NA

for(i in 1:length(file.names)){
  file <- paste(path,file.names[i],sep="")
  aln[[i]]<-read.fasta(file)
  titvs[[i]]<-TiTvRatio(aln[[i]])
  tits[i]<-titvs[[i]][3]
  #aln2[[i]]<-read.alignment(alignmentfile,format="fasta")
  dists[[i]]<-dist.dna(aln[[i]])
  mdists[i]<-mean(dist.dna(aln[[i]]))
}
save(dists,titvs,aln,file=paste(args[3],"dists.titvs.aln.Rfile",sep=""))

draw<-as.data.frame(cbind(tits,mdists))
draw$name<-file.names

pdf(paste(args[3],".pdf",sep=""))
  plot(log(draw$tits,2),draw$mdists,cex=.2,col=rgb(0,0,0,.5),xlab="ti/tv",ylab="mean(K80 distance)",xaxt="n")
  axis(1,at=-1:5,labels=c(0.5,1,2,4,8,16,32))
  write.table(draw,paste(args[3],".raw.tsv",sep=""),quote=F,sep="\t",row.names=F)

  t<-log(tits,10)
  m<-log(mdists,10)
  d<-as.data.frame(cbind(t,m))
  d$name<-file.names
  d<-d[is.finite(d[,1]) & is.finite(d[,2]),]
  lm1<-lm(d$t~d$m)
  print(summary(lm1))
  abline(lm1,lty=2,col="red",lwd=2)
dev.off()


