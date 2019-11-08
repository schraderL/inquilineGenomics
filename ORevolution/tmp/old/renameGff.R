library(rtracklayer)
gff<-readGFF("/Users/lukas/Dropbox/inquilineGenomics18/ORevolution/manualCuratedSets/Ains/LukasImproved/Ains.OR.manualAnnotations.gff3")
e<-data.frame(cbind(unlist(gff$ID),paste("AinsPred",1:length(unlist(gff$ID)),sep="-")))
colnames(e)<-c("name","newname")
#ID
unlisted <- unlist(gff$ID)
unlisted2<-e$newname[match( unlisted,c(as.character(e$name)))]
relisted <- relist(unlisted2, gff$ID)
gff$ID<-relisted

#Name
unlisted <- unlist(gff$Name)
unlisted2<-e$newname[match( unlisted,c(as.character(e$name)))]
relisted <- relist(unlisted2, gff$Name)
gff$Name<-relisted

#Parent
unlisted <- unlist(gff$Parent)
unlisted2<-e$newname[match( unlisted,c(as.character(e$name)))]
relisted <- relist(unlisted2, gff$Parent)
gff$Parent<-relisted


export(gff, "/Users/lukas/Dropbox/inquilineGenomics18/ORevolution/manualCuratedSets/Ains/LukasImproved/Ains.OR.manualAnnotations.renamed.gff3", format = "GFF3")
