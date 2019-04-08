source("~/sciebo/librarySchrader.R")


PARGi<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/PARG.ins.bed",sep="\t",comment.char = "#",F)
ACHAi<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/ACHA.ins.bed",sep="\t",comment.char = "#",F)
AHEYi<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/AHEY.ins.bed",sep="\t",comment.char = "#",F)
Anc09i<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/Anc09.ins.bed",sep="\t",comment.char = "#",F)
AINSi<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/AINS.ins.bed",sep="\t",comment.char = "#",F)
AECHi<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/AECH.ins.bed",sep="\t",comment.char = "#",F)

colnames(ACHAi)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(AHEYi)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(PARGi)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(AECHi)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(AINSi)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(Anc09i)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")

PARGd<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/PARG.del.bed",sep="\t",comment.char = "#",F)
ACHAd<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/ACHA.del.bed",sep="\t",comment.char = "#",F)
AHEYd<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/AHEY.del.bed",sep="\t",comment.char = "#",F)
AECHd<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/AECH.del.bed",sep="\t",comment.char = "#",F)
AINSd<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/AINS.del.bed",sep="\t",comment.char = "#",F)
Anc09d<-read.csv("/Users/lukas/sciebo/inquilineGenomics18/WGA/results/Anc09.del.bed",sep="\t",comment.char = "#",F)

colnames(ACHAd)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(AHEYd)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(PARGd)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(AECHd)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(AINSd)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")
colnames(Anc09d)<-c("Sequence", "Start", "End","MutationID","ParentGenome","ChildGenome")

data<-list()
data[["ACHA"]]<-rbind(ACHAd,ACHAi)
data[["PARG"]]<-rbind(PARGd,PARGi)
data[["AECH"]]<-rbind(AECHd,AECHi)
data[["AINS"]]<-rbind(AINSd,AINSi)
data[["AHEY"]]<-rbind(AHEYd,AHEYi)
data[["Anc09"]]<-rbind(Anc09d,Anc09i)

data2<-list()
for (species in c("ACHA","PARG","AHEY","AINS","AECH","Anc09")){
  data2[[species]]<-split(data[[species]],f=data[[species]]$MutationID)
  }

# P=Transposition
# I=Insertion 
# D=Deletion 
# V=Inversion 
# GI(D)=GapInsertion(GapDeletion) 
# U=Duplication 
# DB=Deletion Breakpoint GDB=Gap Deletion Breakpoint

P<-list()
D<-list()
V<-list()
U<-list()
Psummary<-list()
Dsummary<-list()
Vsummary<-list()
Usummary<-list()

PRanges<-c(0,100,1000,2000,5000,10000,20000,30000,40000,500000,1000000,2000000,2000000000)
DRanges<-c(0,50,100,500,1000,1500,2000,2500,5000,10000,20000,30000,40000)
VRanges<-c(0,50,100,500,1000,1500,2000,2500,5000,10000,20000,2000000000)
URanges<-c(0,50,100,500,1000,1500,2000,2500,5000,10000,20000,30000,40000)

for (species in c("ACHA","PARG","AHEY","AINS","AECH","Anc09")){
  P[[species]]<-data2[[species]]$P
  P[[species]]$length<-P[[species]]$End-P[[species]]$Start
  D[[species]]<-data2[[species]]$D
  D[[species]]$length<-D[[species]]$End-D[[species]]$Start
  V[[species]]<-data2[[species]]$V
  V[[species]]$length<-V[[species]]$End-V[[species]]$Start
  U[[species]]<-data2[[species]]$U
  U[[species]]$length<-U[[species]]$End-U[[species]]$Start
  

  P[[species]]$group<-cut(P[[species]]$length, breaks = PRanges, right = TRUE)
  Psummary[[species]]<-P[[species]] %>% 
    group_by(group) %>%
    count(group)
  Psummary[[species]]$species<-species
  
  D[[species]]$group<-cut(D[[species]]$length, breaks = DRanges, right = TRUE)
  Dsummary[[species]]<-D[[species]] %>% 
    group_by(group) %>%
    count(group)
  Dsummary[[species]]$species<-species
  
  V[[species]]$group<-cut(V[[species]]$length, breaks = VRanges, right = TRUE)
  Vsummary[[species]]<-V[[species]] %>% 
    group_by(group) %>%
    count(group)
  Vsummary[[species]]$species<-species
  
  
  U[[species]]$group<-cut(U[[species]]$length, breaks = URanges, right = TRUE)
  Usummary[[species]]<-U[[species]] %>% 
    group_by(group) %>%
    count(group)
  Usummary[[species]]$species<-species
}


##################################################################
# Deletions
##################################################################

Vall<-do.call("rbind", Vsummary)
Vall$group<-gsub("]",")",as.character(Vall$group))

# rename labels for ranges
Vall$group[Vall$group=="(0,50)"]<-" 0.05 kb"
Vall$group[Vall$group=="(50,100)"]<-" 0.05-0.1 kb"
Vall$group[Vall$group=="(100,500)"]<-" 0.1-0.5 kb"
Vall$group[Vall$group=="(500,1e+03)"]<-" 0.5-1.0 kb"
Vall$group[Vall$group=="(1e+03,1.5e+03)"]<-" 1.0-1.5 kb"
Vall$group[Vall$group=="(1.5e+03,2e+03)"]<-" 1.5-2.0 kb"
Vall$group[Vall$group=="(2e+03,2.5e+03)"]<-" 2.0-2.5 kb"
Vall$group[Vall$group=="(2.5e+03,5e+03)"]<-" 2.5-5.0 kb"
Vall$group[Vall$group=="(5e+03,1e+04)"]<-" 5.0-10 kb"
Vall$group[Vall$group=="(1e+04,2e+04)"]<-" 10-20 kb"
Vall$group[Vall$group=="(2e+04,2e+09)"]<-" 20 kb"

Vall$group<-as.factor(Vall$group)
#levels(Vall$group)<-levels(Vall$group)[c(9,1:8)]
#levels(Vall)<-rev(levels(Vall))

library(dplyr) 
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

theme_set(theme_pubclean())
theme_set()
cols <- colorRampPalette(brewer.pal(8,"Blues"))(length(levels(Vall$group)))
ggplot(Vall, aes(x = species, y = n)) +
  geom_bar(aes(color = Vall$group, fill = Vall$group), stat = "identity") +
  scale_color_manual(values = rep("black",10))  + theme(legend.position="right") + 
  theme_minimal() + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = cols) +
  ylab("Count") +
  xlab("Branch")
  

##################################################################
# Deletions
##################################################################

Dall<-do.call("rbind", Dsummary)
Dall$group<-gsub("]",")",as.character(Dall$group))

#c(0,50,100,500,1000,1500,2000,2500,5000,10000,20000,30000,40000)
# rename labels for ranges
Dall$group[Dall$group=="(0,50)"]<-" 0.05 kb"
Dall$group[Dall$group=="(50,100)"]<-" 0.05-0.1 kb"
Dall$group[Dall$group=="(100,500)"]<-" 0.1-0.5 kb"
Dall$group[Dall$group=="(500,1e+03)"]<-" 0.5-1.0 kb"
Dall$group[Dall$group=="(1e+03,1.5e+03)"]<-" 1.0-1.5 kb"
Dall$group[Dall$group=="(1.5e+03,2e+03)"]<-" 1.5-2.0 kb"
Dall$group[Dall$group=="(2e+03,2.5e+03)"]<-" >2.0 kb"

Dall$group<-as.factor(Dall$group)
#levels(Dall$group)<-levels(Dall$group)[c(9,1:8)]
#levels(Dall)<-rev(levels(Dall))

library(dplyr) 
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

theme_set(theme_pubclean())
theme_set()
cols <- colorRampPalette(brewer.pal(8,"Blues"))(length(levels(Dall$group)))
ggplot(Dall, aes(x = species, y = n)) +
  geom_bar(aes(color = Dall$group, fill = Dall$group), stat = "identity") +
  scale_color_manual(values = rep("black",10))  + theme(legend.position="right") + 
  theme_minimal() + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = cols) +
  ylab("Count") +
  xlab("Branch")

##################################################################
# Transpositions
##################################################################

Pall<-do.call("rbind", Psummary)
Pall$group<-gsub("]",")",as.character(Pall$group))


# rename labels for ranges
Pall$group[Pall$group=="(0,100)"]<-" 0.1 kb"
Pall$group[Pall$group=="(100,1e+03)"]<-" 0.1-1.0 kb"
Pall$group[Pall$group=="(1e+03,2e+03)"]<-" 1.0-2.0 kb"
Pall$group[Pall$group=="(2e+03,5e+03)"]<-" 2.0-5.0 kb"
Pall$group[Pall$group=="(5e+03,1e+04)"]<-" 5.0-10 kb"
Pall$group[Pall$group=="(1e+04,2e+04)"]<-" 10-20 kb"
Pall$group[Pall$group=="(2e+04,3e+04)"]<-" 20-30 kb"
Pall$group[Pall$group=="(3e+04,4e+04)"]<-" 30-40 kb"
Pall$group[Pall$group=="(4e+04,5e+05)"]<-" >40 kb"



Pall$group<-as.factor(Pall$group)
#levels(Pall$group)<-levels(Pall$group)[c(9,1:8)]
#levels(Pall)<-rev(levels(Pall))

library(dplyr) 
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

theme_set(theme_pubclean())
theme_set()
cols <- colorRampPalette(brewer.pal(8,"Blues"))(length(levels(Pall$group)))
ggplot(Pall, aes(x = species, y = n)) +
  geom_bar(aes(color = Pall$group, fill = Pall$group), stat = "identity") +
  scale_color_manual(values = rep("black",10))  + theme(legend.position="right") + 
  theme_minimal() + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = cols) +
  ylab("Count") +
  xlab("Branch")

##################################################################
# U=Duplication
##################################################################

Uall<-do.call("rbind", Usummary)
Uall$group<-gsub("]",")",as.character(Uall$group))


# rename labels for ranges

Uall$group[Uall$group=="(0,50)"]<-" 0.05 kb"
Uall$group[Uall$group=="(50,100)"]<-" 0.05-0.1 kb"
Uall$group[Uall$group=="(100,500)"]<-" 0.1-0.5 kb"
Uall$group[Uall$group=="(500,1e+03)"]<-" 0.5-1.0 kb"
Uall$group[Uall$group=="(1e+03,1.5e+03)"]<-" 1.0-1.5 kb"
Uall$group[Uall$group=="(1.5e+03,2e+03)"]<-" 1.5-2.0 kb"
Uall$group[Uall$group=="(2e+03,2.5e+03)"]<-" 2.0-2.5 kb"
Uall$group[Uall$group=="(2.5e+03,5e+03)"]<-" 2.5-5.0 kb"

Uall$group[Uall$group=="(5e+03,1e+04)"]<-" 5.0-10 kb"
Uall$group[Uall$group=="(1e+04,2e+04)"]<-" 10-20 kb"
Uall$group[Uall$group=="(2e+04,3e+04)"]<-" 20-30 kb"
Uall$group[Uall$group=="(3e+04,4e+04)"]<-" 30-40 kb"
Uall$group[Uall$group=="(4e+04,5e+05)"]<-" >40 kb"



Uall$group<-as.factor(Uall$group)
#levels(Uall$group)<-levels(Uall$group)[c(9,1:8)]
#levels(Uall)<-rev(levels(Uall))

library(dplyr) 
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

theme_set(theme_pubclean())

cols <- colorRampPalette(brewer.pal(8,"Blues"))(length(levels(Uall$group)))
ggplot(Uall, aes(x = species, y = n)) +
  geom_bar(aes(color = Uall$group, fill = Uall$group), stat = "identity") +
  scale_color_manual(values = rep("black",10))  + theme(legend.position="right") + 
  theme_minimal() + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = cols) +
  ylab("Count") +
  xlab("Branch")
