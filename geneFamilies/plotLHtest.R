#lh<-read.table("~/CSE/inquilineGenomics/0006-geneFamilyEvolution/CAFE/lhtest")
#lh<-read.table("/Users/lukas/sciebo/inquilineGenomics18/geneFamilies/CAFE/results/Acro.Atta/LRT/lhtest_result.txt")
lh<-read.table("/Users/lukas/sciebo/inquilineGenomics18/geneFamilies/CAFE/lhtest_result_testing.txt")

colnames(lh)<-c("LHglobal","globalLambda","LHmulti","estLambda1","estLambda2","estLambda3","estLambda4","estLambda5")

#LR = 2*(score of global lambda model ??? score of multi- lambda model)
lh$LR<-2*(lh$LHglobal-lh$LHmulti)
lh<-lh[order(lh$LR),]
lh<-lh[-nrow(lh),]

#calculate observed likelihood ratio
# Score global Lambda: 27947.2
# Score 5 Lambda: 26365.1
LRobs<-2*(-27947.2 - (-26365.1))
hist(lh$LR,xlab="likelihood Ratio",main="",10,col=rgb(.7,.3,.2,.8))
sum(lh$LR<LRobs)/length(lh$LR)

plot(head(lh$LR,73),type="l")

#calculate quantiles
qs<-NA
p<-1
seqs<-seq(.01,1,length=10)
for(i in seqs){
	qs[p]<-quantile(lh$LR,i)
	
	p<-p+1
	
}

plot(qs~seqs,type="l",xlab="quantile",ylab="simulation LR")

abline(h=LRobs,col="red",lty=2)
abline(v=seqs[max(which(qs<LRobs))],col="red",lty=2)
par(mgp = c(3, 2, 0))
axis(1,at=seqs[max(which(qs<LRobs))],labels=trunc(seqs[max(which(qs<LRobs))]*10000)/10000,col="red",col.axis="red",cex=.7)

boxplot(lh$estLambda1,lh$estLambda2,outline=F)
abline(h=0.0028,col="red",lty=2)

#p-value
sum(lh$LR<LRobs)
