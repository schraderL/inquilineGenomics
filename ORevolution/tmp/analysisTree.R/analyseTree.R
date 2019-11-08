
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
library(diversitree)
library(phytools)
library(geiger)

############################################
# Load rerooted RaxML tree
############################################
t1<-read.tree("/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/tree/ant.OR.rooted.tre")
tree <- as(t1, 'TreeMan')

############################################
# Define Clades
############################################

#https://www.r-phylo.org/wiki/HowTo/DataTreeManipulation
#How do I calculate the patristic distance between two taxa?
distances<-cophenetic(t1)

par(mfrow=c(5,5))
par(mai=c(0,0,0,0))
clade<-list()
element<-1
distances<-distances[grep("Ahey|Aech|Acol|Acep|Parg|Ains|Acha",rownames(distances),perl=T),]
for(p in c(1:dim(distances)[1])){
  distance<-sort(distances[p,])
  #plot(distance[1:18])
  innerdistances<-NA
  for(i in 1:18){
    innerdistances[i]<-abs(distance[i+2]-distance[i+3])
  }
  biggestStep<-which(innerdistances==max(innerdistances))+4
  abline(v=biggestStep+0.5,lty=2,col="red")
  clade[[element]]<-sort(names(distance[1:biggestStep]))
  element<-element+1
}

par(mfrow=c(5,5))
par(mai=c(0.1,0.1,0.3,0))
clade2<-unique(clade)
#clade3<-duplicated(clade)
cladeTrees<-list()
p<-1
for(i in 1:length(clade2)){
#for(i in 26:100){
  if(length(clade2[[i]])<4){
    next
  }
  else{
    parentNode<-getPrnt(tree, clade2[[i]])
    tn<-getSubtree(tree,parentNode)
    tnApe <- as(tn, 'phylo')
    cladeTrees[[p]]<-tnApe
    p<-p+1
    #plot(tnApe,main=i,cex=.4)
  }
}

cladesTree<-getCladesofSize(t1,11)
plot(cladeSizes)
abline(h=11,col="red",lwd=2)

cladeSizes<-NA
traits<-rep(NA,length(t1$tip.label))
traits2<-rep(NA,length(t1$tip.label))
traits3<-rep(NA,length(t1$tip.label))
names(traits)<-t1$tip.label
names(traits2)<-t1$tip.label
names(traits3)<-t1$tip.label
for(i in 1:length(cladesTree)){
  cladeSizes[i]<-(length(cladesTree[[i]]$tip.label))
  # check if clade contains Parg
  if (length(grep("Ahey",cladesTree[[i]]$tip.label))>length(grep("Aech",cladesTree[[i]]$tip.label))){
    traits[cladesTree[[i]]$tip.label]<-1
  }else{
    traits[cladesTree[[i]]$tip.label]<-0
  }
  if (length(grep("Ahey",cladesTree[[i]]$tip.label))>length(grep("Parg",cladesTree[[i]]$tip.label))){
    traits2[cladesTree[[i]]$tip.label]<-1
  }else{
    traits2[cladesTree[[i]]$tip.label]<-0
  }
  if (length(grep("Aech",cladesTree[[i]]$tip.label))>length(grep("Ains",cladesTree[[i]]$tip.label))){
    traits3[cladesTree[[i]]$tip.label]<-1
  }else{
    traits3[cladesTree[[i]]$tip.label]<-0
  }
  
}

traitsAll<-data.frame(cbind(traits,traits2,traits3))
trait.plot(t1,traitsAll ,list(traits = c("gray30", rgb(1,0,0,1)),traits2 = c("gray20", rgb(1,.1,0,1)),traits3 = c("gray10", rgb(1,.1,.1,1))),cex.lab = 0.00001,w=1)
############################################
# Remove all suboptimal clades
############################################
# check clade for duplicates and remove accordingly (keep the one clade whose correspondng cladeTree is closest to 11 members)
"TzetOr-236"

# loop over all tips
for (p in 1:length(ttmp$tip.label)){
  filterClade<-grep(ttmp$tip.label[p],clade2)
  filter<-NA
  for (i in 1:length(filterClade)) {
    filter[i]<-(length(cladeTrees[[filterClade[i]]]$tip.label))
  }
  # keep clade 1
  print(which(filter==min(filter)))
}
############################################
# Analyse clades
############################################

checkParg<-NA
checkAhey<-NA
checkTsep<-NA
checkAech<-NA
checkAins<-NA
for(i in 1:length(cladeTrees)){
  #checkParg[i]<-any(grepl("Parg",cladeTrees[[i]]$tip.label))
  #checkAhey[i]<-any(grepl("Ahey",cladeTrees[[i]]$tip.label,perl=T))
  #checkTsep[i]<-any(grepl("Tsep",cladeTrees[[i]]$tip.label,perl=T))
  #checkAech[i]<-any(grepl("Aech",cladeTrees[[i]]$tip.label,perl=T))
  #checkAins[i]<-any(grepl("Ains",cladeTrees[[i]]$tip.label,perl=T))
  checkParg[i]<-sum(grepl("Parg",clade2[[i]],perl=T))
  checkAhey[i]<-sum(grepl("Ahey",clade2[[i]],perl=T))
  checkTsep[i]<-sum(grepl("Tsep",clade2[[i]],perl=T))
  checkAech[i]<-sum(grepl("Aech",clade2[[i]],perl=T))
  checkAins[i]<-sum(grepl("Ains",clade2[[i]],perl=T))
}

############################################
# Plot Clades
############################################

# plot red if no parasite gene in clade
# plot green if all 3 parasites in clade
# for each tip, find clade and colour according to above condition
ttmp<-t1
traits<-rep(NA,length(ttmp$tip.label))
#traits2<-rep(0,length(ttmp$tip.label))
#traits3<-rep(0,length(ttmp$tip.label))
names(traits)<-ttmp$tip.label
#names(traits2)<-ttmp$tip.label
#names(traits3)<-ttmp$tip.label
traits[grep("Parg|Acha|Ains",names(traits),perl=T)]<-1
traits[grep("Acep|Acol|Ccos|Tzet|Ahey|Aech|Tcor|Tsep",names(traits),perl=T)]<-0
traits3[grep("Acha",names(traits3))]<-1
traitsAll<-data.frame(cbind(traits,traits2,traits3))
#trait.plot(ttmp,traitsAll ,list(traits = c("lightgray", "red"),traits2 = c("lightgray", "blue"),traits3 = c("lightgray", "green")),cex.lab = 0.00001,w=.1)
trait.plot(ttmp,traitsAll ,list(traits = c("lightgray", "red")),cex.lab = 0.00001,w=1)

########################################################################################
# Find ORCO

orcoNode<-getNdSstr(tree,"CcosOr-001")
getCnnctdNds(tree,c("CcosOr-001","TzetOr-001","AcolOr-001","AechOr-001"))
tOrco<-getSubtree(tree,orcoNode)
tOrco2 <- as(tOrco, 'phylo')
getOtgrp(tree,id="CcosOr-001")
getPath(tree,"CcosOr-001","PargOr-001")
parent<-getPrnt(tree, c("CcosOr-001","PargOr-001"))
tOr2<-getSubtree(tree,parent)
tOr2ape <- as(tOr2, 'phylo')
plot(tOr2ape)

n<-getPrnt(tree,c("CcosOr-100","AechOr-100"))
tn<-getSubtree(tree,n)
tnApe <- as(tn, 'phylo')
getSpnsAge(tree,c("CcosOr-001","AechOr-001"),1)
tids<-c("CcosOr-100","AechOr-100")
ids<-c("CcosOr-100","AechOr-100")
calcFrPrp(tn, tids, progress = "none")
calcNdsBlnc(tree, ids, parallel = FALSE, progress = "none")
plot(tnApe,cex=.2)

library(diversitree)
library(phytools)
library(geiger)
#https://www.r-phylo.org/wiki/HowTo/DataTreeManipulation
tips(tnApe, "AechOr-100")
#How do I calculate the patristic distance between two taxa?
cophenetic(t1)["TzetOr-002", "AechOr-002"]
tmp<-cophenetic(t1)

par(mfrow=c(5,5))
par(mai=c(0,0,0,0))
clade<-list()
element<-1
for(p in c(2000:2099)){
  test<-sort(tmp[p,])
  #plot(test[1:20])
  bla<-NA
  for(i in 1:20){
    bla[i]<-abs(test[i]-test[i+1])
  }
  biggestStep<-which(bla==max(bla))
  #abline(v=biggestStep+0.5,lty=2,col="red")
  clade[[element]]<-sort(names(test[1:biggestStep]))
  
  element<-element+1
}
par(mfrow=c(5,5))
clade2<-unique(clade)
for(i in 1:length(clade2)){
  if(length(clade2[[i]])==1){
    next
    }
  else{
    parentNode<-getPrnt(tree, clade2[[i]])
    tn<-getSubtree(tree,parentNode)
    tnApe <- as(tn, 'phylo')
    plot(tnApe)
  }
}

cophenetic(t1)["TzetOr-001", "AechOr-001"]
parent<-getPrnt(tree, c("CcosOr-001","CcosOr-293"))
#sis<-getNdSstr(tree,"CcosOr-005")
#parent<-getPrnt(tree, sis)
tOr2<-getSubtree(tree,parent)
tOr2ape <- as(tOr2, 'phylo')
#plot(tOr2ape)


ttmp<-t1
traits<-rep(NA,length(ttmp$tip.label))
traits2<-rep(0,length(ttmp$tip.label))
traits3<-rep(0,length(ttmp$tip.label))
names(traits)<-ttmp$tip.label
names(traits2)<-ttmp$tip.label
names(traits3)<-ttmp$tip.label
traits[grep("Parg|Acha|Ains",names(traits),perl=T)]<-1
traits[grep("Acep|Acol|Ccos|Tzet|Ahey|Aech|Tcor|Tsep",names(traits),perl=T)]<-0
traits3[grep("Acha",names(traits3))]<-1
traitsAll<-data.frame(cbind(traits,traits2,traits3))
#trait.plot(ttmp,traitsAll ,list(traits = c("lightgray", "red"),traits2 = c("lightgray", "blue"),traits3 = c("lightgray", "green")),cex.lab = 0.00001,w=.1)
trait.plot(ttmp,traitsAll ,list(traits = c("lightgray", "red")),cex.lab = 0.00001,w=1)

getNdPtids(tree, sis)
test<-getPrnt(tree, c(sis,"PargOr-002"))
childs<-getNdKids(tree,test)

