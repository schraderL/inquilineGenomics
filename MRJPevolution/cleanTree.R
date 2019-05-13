
######################
# Install packages
######################
#chooseCRANmirror()
#install.packages("BiocManager")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Biostrings", version = "3.8")
#BiocManager::install("ggtree", version = "3.8")
#BiocManager::install("treeman")

######################
# Load packages
######################
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("phytools")
library("treeman")


############################################
# Load rerooted RaxML tree
############################################
t1<-read.tree("~/sciebo/inquilineGenomics18/geneFamilies/specificFamilies/mrjps/tree/mergeAln/FigtreeRooted.bipartitions.nwk")


############################################
# Rename Amel genes where necessary
############################################
c("LOC413379","LOC413380","LOC724293","LOC726274","Mrjp3","Mrjp5","Mrjp6","Mrjp4","Mrjp7","Mrjp1","LOC102654393","Mrjp2","Mrjp8","Mrjp9","Y-e3","LOC727110","Y-y","Y-h","LOC113218568","Y-f")
c("Y-g1"     ,"Y-g2"     ,"Y-x1"     ,"Y-e1"     ,"Mrjp3","Mrjp5","Mrjp6","Mrjp4","Mrjp7","Mrjp1","Mrjp2-like"  ,"Mrjp2","Mrjp8","Mrjp9","Y-e3","Y-x2"     ,"Y-y","Y-h","Y-like","Y-f")
for (i in 1:length(old)){
  t1$tip.label[t1$tip.label==old[i]]<-new[i]
}

############################################
# Rename yellow genes in ants
############################################

overview<-list()
for (query in c("Y-h","Y-y","Y-g1","Y-g2","Y-f","Y-e1","Y-x2","Y-e3","Y-like","Y-x1")){
  t2 <- as(t1, 'TreeMan')
  sis<-getNdSstr(t2,query)
  childs<-getNdKids(t2,sis)
  #gsub(t1$tip.label[t1$tip.label %in% childs])
  t1$tip.label[t1$tip.label %in% childs]<-gsub(x=childs,pattern = "(....).*",paste("\\1",query,sep="-"),perl=T)
  overview[[query]]<-data.frame(childs,gsub(x=childs,pattern = "(....).*",paste("\\1",query,sep="-"),perl=T))
  colnames(overview[[query]])<-c("original","new")
}

############################################
# Plot
############################################

plot(t1,cex=.3,type = "fan")
dev.print(pdf,"~/sciebo/inquilineGenomics18/geneFamilies/specificFamilies/mrjps/tree/mergeAln/Ants.MRJP-yellow.pdf",width=7,height=7)
overviewTable<-do.call(rbind.data.frame, overview)
write.table(overviewTable,"~/sciebo/inquilineGenomics18/geneFamilies/specificFamilies/mrjps/tree/mergeAln/Ants.MRJP-yellow.overviewTable.tsv",sep="\t",quote = F,row.names = F)
write.tree(t1,"~/sciebo/inquilineGenomics18/geneFamilies/specificFamilies/mrjps/tree/mergeAln/Ants.MRJP-yellow.txt")

#https://itol.embl.de/tree/12817610684345251554995787