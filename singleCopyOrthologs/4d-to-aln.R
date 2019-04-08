# calculate alignment lengths per ortholog for inquiline phylogeny
a<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/phylogeny/SCOs/QC4d/CDSaln4d.statistics.tsv",sep="\t",F)
dim(a)
boxplot(a$V3)
pdf("/Users/lukas/sciebo/inquilineGenomics18/phylogeny/SCOs/QC4d/4d-to-aln.pdf",width=5,height=4)
plot(a$V3,a$V5,pch=19,cex=.6,col=rgb(0,0,0,.1),log="xy",xlab="alignment length",ylab="number of 4-fold degenerate sites")
dev.off()
lm(a$V3~a$V5)


       