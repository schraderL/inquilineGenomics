library(VGAM)
p<-0.6
hist(rskewnorm(p*1000000, location = 17.67667, scale = 2, shape = -1),xlim=c(5,22),200,col="gray80",border="gray80",main="",xlab="S2N(0.6,17.68,2,-1,12.64,2,1)",ylab="",yaxt="n")
hist(rskewnorm((1-p)*1000000, location = 12.64, scale = 2, shape = 1),add=T,200,col="gray80",border="gray80")
axis(2,at=c(0,(15000)),labels=c(0,0.015))
dev.print(pdf,"/Users/lukas/sciebo/inquilineGenomics18/phylogeny/TreeDating/MCMCtree/results/S2N.pdf")
