#R

args <- commandArgs()
# args[6] = $base
# args[7] = $species
print(args)
setwd(args[6])
library(stringr)
print(paste(args[7],".fullORannotation.qc",sep=""))
a<-read.csv(paste(args[7],".fullORannotation.qc",sep=""),sep = "\t",F)
#a<-read.csv(paste("/Users/lukas/Dropbox/inquilineGenomics18/ORevolution/results/Ains/20180718/Ains",".fullORannotation.qc",sep=""),sep = "\t",F)
#a<-read.csv(paste("/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/Ains/20180718/Ains",".fullORannotation.qc",sep=""),sep = "\t",F)
#a<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/Ahey/190404/Ahey.fullORannotation.qc",sep = "\t",F)
#head(a)
b<-read.csv("tmp.out",sep = "\t",F,comment.char="#")
#b<-read.csv("/Users/lukas/Dropbox/inquilineGenomics18/ORevolution/results/Ains/20180718/tmp.out",sep = "\t",F,comment.char="#")
#b<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/Ains/20180718/tmp.out",sep = "\t",F,comment.char="#")
#b<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/Ahey/190404/tmp.out",sep = "\t",F,comment.char="#")
a$name<-gsub(pattern = ">(.*?)\\s+.*",replacement = "\\1",x = a$V1,perl = T)
colnames(b)<-c("seqid","alignmentstart","alignmentend","envelopestart","envelopeend","hmmacc","hmmname","type","hmmstart","hmmend","hmmlength","bitscore","hmmEvalue","significance","clan")
m1<-merge(a,b,by.x="name",by.y="seqid",all.x=T)

m1$d[m1$hmmstart>50 & m1$hmmend < m1$hmmlength-50]<-"fdf"
m1$d[m1$hmmstart>50 & m1$hmmend >= m1$hmmlength-50]<-"fd"
m1$d[m1$hmmstart <= 50 & m1$hmmend < m1$hmmlength-50]<-"df"
m1$d[m1$hmmstart <= 50 & m1$hmmend >= m1$hmmlength-50]<-""
m1$d[m1$hmmname != "7tm_6" | is.na(m1$hmmname)]<-"d"

blast<-read.csv(paste("filtering/",args[7],".OR.no7tm6.tophit",sep=""),sep="\t",F)
#blast<-read.csv("/Users/lukas/Dropbox/inquilineGenomics18/ORevolution/results/Ains/20180718/filtering/Ains.OR.no7tm6.tophit",sep="\t",F)
#blast<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/Ains/20180718/filtering/Ains.OR.no7tm6.tophit",sep="\t",F)
#blast<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/Ahey/190404/filtering/Ahey.OR.no7tm6.tophit",sep="\t",F)

colnames(blast)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
e<-merge(m1,blast,by.x="name",by.y="qseqid",all.x=T)
# remove some entries that are duplicate due to pfam scan issues
e<-subset(e, !duplicated(e$name))

e$blast<-" "
e$blast[e$evalue < 1e-10 | !is.na(e$evalue)]<-"H"
e<-e[order(e$hmmEvalue),]
#remove entries that have no domain and no homology
# "select those that have homology or have the domain."
e<-e[e$blast=="H" | e$d != "d",]
e$newname<-paste(args[7],"Or-",str_pad(1:length(e$name), 3, pad = "0"),"-",e$V3,e$d,e$blast,sep="")
#e$newname<-paste("Ains","Or-",str_pad(1:length(e$name), 3, pad = "0"),"-",e$V3,e$d,e$blast,sep="")
#e$newname<-paste("Ahey","Or-",str_pad(1:length(e$name), 3, pad = "0"),"-",e$V3,e$d,e$blast,sep="")
e$pos<-gsub(x=e$V1,pattern=".*\\s",perl=T,replacement="")
e$newname<-gsub(pattern="\\s",replacement="",perl=T,x=e$newname)
e$newname<-gsub(pattern="-$",replacement="",perl=T,x=e$newname)
#head(e)
e$Parent<-gsub(pattern=".*?\034(.*?)\\s+.*",replacement="\\1",perl=T,x=e$V1)
e<-subset(e, !duplicated(paste(e$pos,e$V2)))
write.table(e,"renamed.tsv",quote=F,sep="\t",col.names=F,row.names=F)

library(rtracklayer)
gff<-readGFF(paste(args[7],".fullORannotation.gff3",sep=""))
#gff<-readGFF("/Users/lukas/Dropbox/inquilineGenomics18/ORevolution/results/Ains/20180718/Ains.fullORannotation.gff3")
#gff<-readGFF("/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/Ains/20180718/Ains.fullORannotation.gff3")
#gff<-readGFF("/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/Ahey/190404/Ahey.fullORannotation.gff3")
gff<-subset(gff, !duplicated(paste(gff$seqid,gff$type,gff$start,gff$end,gff$score)))
unlisted <- unlist(gff$Parent)
unlisted2 <- paste(e$newname[(match( unlisted,c(e$name)))],"-mRNA",sep="")
unlisted2[unlisted2 == "NA-mRNA"]<-NA
unlisted3 <- e$newname[(match( unlisted,c(e$Parent)))]
unlisted2[is.na(unlisted2)]<-unlisted3[is.na(unlisted2)]
#head(unlisted2)
relisted <- relist(unlisted2, gff$Parent)
gff$Parent<-relisted

unlisted <- unlist(gff$ID)
unlisted2 <- paste(e$newname[(match( unlisted,c(e$name)))],"-mRNA",sep="")
unlisted2[unlisted2 == "NA-mRNA"]<-NA
unlisted3 <- e$newname[(match( unlisted,c(e$Parent)))]
unlisted2[is.na(unlisted2)]<-unlisted3[is.na(unlisted2)]
#head(unlisted2)
relisted <- relist(unlisted2, gff$ID)
gff$ID<-relisted
gff$Name<-gff$ID

export(gff, "renamed.gff3", format = "GFF3")
#export(gff, "/Users/lukas/sciebo/inquilineGenomics18/ORevolution/results/tmp/renamed.gff3", format = "GFF3")
