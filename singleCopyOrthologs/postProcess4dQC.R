#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#load("/Users/lukas/sciebo/inquilineGenomics18/phylogeny/SCOs/QC4d/QC4ddists.titvs.aln.Rfile")
#draw<-read.table("/Users/lukas/sciebo/inquilineGenomics18/phylogeny/SCOs/QC4d/QC4d.raw.tsv",T)
#Rscript path/to/QC4ddists.titvs.aln.Rfile path/to/QC4d.raw.tsv outname
load(args[1])
draw<-read.table(args[2],T)

library(msaR)
library(strataG)

#par(mfrow=c(5,5))
#par(mai=c(0.25,0.2,0.1,0.1))
#par(mgp=c(0,0,0))
mods<-list()
coefs<-NA
dist<-list()
dist.corrected<-list()
for (i in 1:length(aln)){
  dist[[i]]<-dist.dna(aln[[i]], model="raw")
  #titv<-
  dist.corrected[[i]]<-dist.dna(aln[[i]], model="TN93")
  ###Make plot###
  if(is.finite(sum(dist.corrected[[i]]))){
      lm_coef<-coef(lm(dist[[i]]~dist.corrected[[i]]))
      mods[[i]]<-lm_coef
      coefs[i]<-lm_coef[2]
      }
}

bp<-boxplot(coefs,plot=F)
saturated<-NA

for (i in 1:length(aln)){
  if(is.finite(sum(dist.corrected[[i]]))){
      if((mods[[i]][2])<bp$stats[1,1]){
        #print(summary(lm(dist~dist.corrected)))
        saturated<-rbind(saturated,cbind(draw[i,],coefs[i]))
#       print(draw$name[i])
#       plot(dist[[i]]~dist.corrected[[i]], pch=20, col="red", xlab="TrN model distance", ylab="Uncorr. gen. distance", main="",xaxt="n",yaxt="n")
#       abline(0,1, lty=2)
#       abline(mods[[i]], lwd=3,lty=2)
#       text(0.2,0.05,bquote(y == .(mods[[i]][2])*x))
#       mtext(text = gsub("(OG.*?)\\..*","\\1",draw$name[i]),3,padj = 2,adj=.2,cex=.6)
      }
      }
    }

saturated<-saturated[-1,]

write.table(saturated,paste(args[3],"saturated.tsv",sep=""),sep="\t",quote=F,row.names = F)

