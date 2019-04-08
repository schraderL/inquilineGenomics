syn<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/synteny/results/SyntenicRibbons.acro.conf",sep=" ",F)


hist(abs(syn$V5-syn$V6),100)
hist(abs(syn$V2-syn$V3),100,add=T,col=rgb(1,0,0,.2))



syn$lengthAins<-(abs(syn$V5-syn$V6))
syn<-syn[order(syn$lengthAins,decreasing=T),]
head(syn)

syn$s1<-substr(syn$V1,1,4)
syn$s2<-substr(syn$V4,1,4)
syn$pair<-as.factor(paste(syn$s1,syn$s2,sep="."))
a<-split(syn,syn$pair)
boxplot(a[[1]]$lengthAins,a[[2]]$lengthAins,a[[3]]$lengthAins,a[[4]]$lengthAins,a[[5]]$lengthAins,outline=F)
length(a)
par(mai=c(1,2,0,0))
lengths<-lapply(a, '[[', 7)
boxplot(lengths,las=2,outline=F,cex.names=.5,horizontal=T)
lengths

head(syn)
#mean(a[[1]]$lengthAins)
#mean(a[[2]]$lengthAins)
#mean(a[[3]]$lengthAins)
(a[[4]]$pair[1])
mean(a[[4]]$lengthAins)
(a[[5]]$pair[1])
mean(a[[5]]$lengthAins)
(a[[6]]$pair[1])
mean(a[[6]]$lengthAins)


head(syn)
bla<-NA
for (i in 1:length(a)){
  bla[i]<-(sum(a[[i]]$lengthAins<10000)/dim(a[[i]])[1])
  names(bla)[i]<-(as.character(a[[i]]$pair[1]))
}
barplot(bla,cex.names=.6,las=2,ylim=c(0,1))

