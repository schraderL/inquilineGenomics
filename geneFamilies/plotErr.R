setwd("/Users/lukas/sciebo/inquilineGenomics18/geneFamilies/CAFE/results/Acro.Atta/errorRates/")
results<-c("/Users/lukas/sciebo/inquilineGenomics18/geneFamilies/CAFE/results/Acro.Atta/errorRates/errorTest.smallSet.filter100","/Users/lukas/sciebo/inquilineGenomics18/geneFamilies/CAFE/results/Acro.Atta/errorRates/errorTest.smallSet.filter100.RUN2")

#Define function to plot the error model likelihoods
plotErr <- function(species){
	
  a<-read.table(paste(results[1],"/err/",species,".err",sep=""))
  b<-read.table(paste(results[2],"/err/",species,".err",sep=""))
	colnames(a)<-c("ErrorModel","Score")
	colnames(b)<-c("ErrorModel","Score")
	plot(a$ErrorModel,a$Score,type = "l",main=species,xlab="ErrorModel",ylab="Score")
	points(a$ErrorModel,a$Score,col="red")
	min<- a$ErrorModel[a$Score==min(a$Score)]
	segments(min,0,min,1000000,col="red")
	
	points(b$ErrorModel,b$Score,col="red",pch=2)
	min<- b$ErrorModel[b$Score==min(b$Score)]
	segments(min,0,min,1000000,col="red")
	par(xpd=NA)
	text(min,mean(b$Score),round(min,5),cex=1,pos=4,col="red",font=2)
	par(xpd=F)

	
	
	
	
	
	
	}

par(mfrow=c(2,3))
plotErr("PARG")
plotErr("ACOL")
plotErr("AINS")
plotErr("ACHA")
plotErr("AECH")
plotErr("AHEY")
