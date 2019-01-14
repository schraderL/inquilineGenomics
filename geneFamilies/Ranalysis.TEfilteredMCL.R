########################################################
# R script to remove MCL gene clusters with TE contamination
# and plot hierachical clustering and some more QC plots
########################################################

# load libraries
#install.packages("pvclust")
library(pvclust)
source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("limma")
library(limma)

# Define the coefficient of variance function
CV <- function(values){
      (sd(values)/mean(values))*100
      }

#set the path
#local
#path<-"/Users/lukas/sciebo/inquilineGenomics18/geneFamilies/MCLclustering/mcl.TEfiltered"
#outpath<-"/Users/lukas/sciebo/inquilineGenomics18/geneFamilies/MCLclustering/mcl.clean/"

#pallas
#path<-"~/data/inqGen18/geneFamilies/results/mcl.TEfiltered"
#outpath<-"~/data/inqGen18/geneFamilies/results/mcl.clean/"

# prepare environment
a<-list()
clean<-list()
TEonly<-list()
setwd(path)
file.names <- dir(path, pattern ="largeQC")


# loop over each inflation factor (here 1.0, 1.5, and 3.0)
totalCol<-8
range_of_genomes<-2:7
idCol<-9
DescriptionCol<-14
for(i in 1:length(file.names)){
	#load the large QC file in a list
	a[[i]] <- read.csv(file.names[i],T,sep="\t")

	#create TE-free cluster file
	# select all cluster that don't contain any TEs
	clean[[i]]<-a[[i]][a[[i]]$TEcount==0,]
	# select all cluster that do contain TEs
	TEonly[[i]]<-a[[i]][a[[i]]$TEcount>1,]

	names(a)[i]<-file.names[i]
	names(clean)[i]<-file.names[i]

	#Calculate some QC measures
	clusters<-clean[[i]]
	cv<-NA
	maxdiff<-NA
	medianSize<-NA
	deviation<-NA


	for(count in 1:length(clusters[,1])){
		# calculate coefficient of variance for each cluster (i.e. how strong does the gene family size differ between species?)
	  cv[count]<-CV(as.numeric(as.numeric(clusters[count,range_of_genomes])))
	  # calculate maximal difference between species
		maxdiff[count]<-max(as.numeric(clusters[count,range_of_genomes]))-min(as.numeric(clusters[count,range_of_genomes]))
		# calculate median size of a given cluster
		medianSize[count]<-median(as.numeric(clusters[count,range_of_genomes]))
		# calculate maximal deviation from median size
		deviation[count]<-max(as.numeric(clusters[count,range_of_genomes]))-medianSize[count]
		}
	qc<-deviation/(medianSize+1)
	clusters<-cbind(clusters,cv,maxdiff,medianSize,deviation,qc)
	clean[[i]]<-clusters
}

#save fullQC for I15 (clean[[2]]), as it will be used later
file.names[1]
write.table(clean[[1]],paste(outpath,"/finalQC.I01.5.tsv",sep=""),quote=F,row.names=F,sep="\t")

#Plot hierachical clustering results for the TE-infested clusters with more than one member

pdf(paste(outpath,"/hierachicalClustering.TE-infested.pdf",sep=""))
par(mfrow=c(1,3))

for (i in 1:length(a)){
#Define subset of Clusters
  # all TE-infested clusters with more than 2 members
	sub<-subset(TEonly[[i]],rowSums(TEonly[[i]][,range_of_genomes]>0,na.rm=T)>=2)
	# Ward Hierarchical Clustering with Bootstrapped p values
	data<-(as.matrix(sub[,range_of_genomes]/sub[,totalCol]))
	fit <- pvclust(data, method.hclust="ward.D2",method.dist="euclidean",nboot=1000,parallel=T)
	plot(fit,main=names(clean)[i]) # dendogram with p values
	pvrect(fit, alpha=.9,max.only=F)
}
dev.off()

#Plot hierachical clustering results for all TE free clusters with more than one member
pdf(paste(outpath,"/hierachicalClustering.TE-free.pdf",sep=""))
par(mfrow=c(1,3))
for (i in 1:length(clean)){
#Define subset of Clusters
  # all clean (i.e. TE-free) clusters with more than 2 members
	sub<-subset(clean[[i]],rowSums(clean[[i]][,range_of_genomes]>0,na.rm=T)>=2)
	# Ward Hierarchical Clustering with Bootstrapped p values
	data<-(as.matrix(sub[,range_of_genomes]/sub[,totalCol]))
	fit <- pvclust(data, method.hclust="ward.D2",method.dist="euclidean",nboot=1000,parallel=T)
	plot(fit,main=names(clean)[i]) # dendogram with p values
	pvrect(fit, alpha=.9,max.only=F)
}
dev.off()


#plot MDS
pdf(paste(outpath,"/mds.TE-free.pdf",sep=""))
par(mfrow=c(1,3))
	for (i in 1:length(clean)){
	sub<-subset(clean[[i]], rowSums(clean[[i]][,range_of_genomes]>0,na.rm=T)>=2)
	plotMDS(scale(sub[,range_of_genomes]/sub[,totalCol]),main=names(clean)[i])
}
dev.off()

#plot MDS
pdf(paste(outpath,"/mds.TE-infested.pdf",sep=""))
par(mfrow=c(1,3))
	for (i in 1:length(clean)){
	sub<-subset(TEonly[[i]], rowSums(TEonly[[i]][,range_of_genomes]>0,na.rm=T)>=2)
	plotMDS(scale(sub[,range_of_genomes]/sub[,totalCol]),main=names(clean)[i])
}
dev.off()


#Plot Basic Stats
pdf(paste(outpath,"/BasicStats.pdf",sep=""))
par(mfrow=c(4,1))
#plot FamilySizes for each inflation parameter
plot(clean[[1]]$FamSize,type="l",col="red",xlim=c(0,100),main="Family size",ylab="Cluster Size")
for (i in 2:length(clean)){
	points(clean[[i]]$FamSize,type="l")
}
#plot Density distribution for the qc parameter for each inflation parameter
for (i in 1:length(clean)){
	plot(density(clean[[i]]$qc),xlim=c(0,3),ylim=c(0,2.5),main=names(clean)[i])
}
dev.off()

#Prepare CAFE files
cafe<-list()
tmp<-list()

for (i in 1:length(clean)){
	tmp[[i]]<-clean[[i]][rowSums(clean[[i]][,range_of_genomes]>0,na.rm=T)>=2,]
	cafe[[i]]<-tmp[[i]][,c(DescriptionCol,idCol,range_of_genomes)]
	colnames(cafe[[i]])[1:2]<-c("Description","ID")
	names(cafe)[i]<-file.names[i]
	fileOut<-paste(outpath,gsub("tsv","cafe.tsv",file.names[i]),sep="")
  write.table(cafe[[i]],fileOut,quote=F,sep="\t",row.names=F)
}
