

##### bash
perl ../scripts/GOannotation2topGO.pl /usr/local/home/lschrader/data/inqGen18/inquilines_v2.1/Acromyrmex_heyeri.2.1/annotation/function_annotation/Acromyrmex_heyeri.v2.1.pep.fa.iprscan.gene.GO > gene2GO.Ahey.tsv



#source("https://bioconductor.org/biocLite.R")
#biocLite("GO.db")
library("topGO")

#setwd("/Users/lukas/CSE/inquilineGenomics/SoS/absREL2/GOenrichment")

##### get Data
#all
files <- list.files(path="/Users/lukas/CSE/inquilineGenomics/SoS/absREL2/absREL.results", pattern=".out", all.files=T, full.names=T)
l<-read.table("/Users/lukas/CSE/inquilineGenomics/SoS/absREL2/trimmedAln.length",F)

colnames(l)<-c("OGid","length")

data<-list()
i<-1
for (file in files) {

########################
#SELECT RIGHT DATA
########################
# all set
data[[i]]<-read.csv(file)
names(data)[i]<-substr(strsplit(file,"/")[[1]][9],4,13)
data[[i]]$tag<-sub("_*[0-9]{5}","",data[[i]]$Branch,perl=T)
i<-i+1
}


parg<-NA
p<-1
sign<-list()
for (i in 1:length(data)){
	data[[i]]$FDR<-p.adjust(data[[i]]$p_Holm,method="fdr",n=length(data))
	if (sum(data[[i]]$tag %in% "Parg")==0){
		next
		}
	tmp<-(data[[i]][data[[i]]$tag=="Parg",])
	if (tmp$FDR<0.05){
		sign[[p]]<-data[[i]]
		parg[p]<-as.character(tmp$Branch)
		p<-p+1
}
}





# read in list of GO terms:
###########################
geneID2GO<-readMappings(file="gene2GO.Parg.tsv")
#geneID2GO<-readMappings(file="gene2GO.Ains.tsv")
#geneID2GO<-readMappings(file="gene2GO.Acha.tsv")
geneNames<-names(geneID2GO)


# read in list of interesting genes:
####################################
#BUSTED<-read.table("BUSTED.results.tsv")
#colnames(BUSTED)<-c("OG","Acep","Aech","Acha","Ahey","Ains","Parg","LR","p")
#BUSTED$FDR<-p.adjust(BUSTED$p,method="bonferroni")

#Plot histogram of p-values
######################################


#select significant genes
######################################
#posSel<-BUSTED[BUSTED$p<0.05,]

#FDR correction
#posSel<-BUSTED[BUSTED$FDR<0.05,]

#restrict gene space to tested 1-to-1 orthologs
######################################
#geneNames<-geneNames[geneNames %in% parg]


# create a 0/1 matrix (interesting/not interesting)
######################################

geneList <- factor(as.integer(geneNames %in% parg))

table(geneList)

names(geneList) <- geneNames

# create object for Biological Process, Molecular Function and Cellular component:
##################################################################################

GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)


# Molecular Function
####################

resultFis_MF 	<- runTest(GOdata_MF, algorithm = "parentchild", statistic = "fisher")

table_MF <- GenTable(GOdata_MF, parentChild = resultFis_MF, orderBy = "parentChild", ranksOf = "parentChild", topNodes = 200)


# Biological Process
####################

resultFis_BP <- runTest(GOdata_BP, algorithm = "parentchild", statistic = "fisher")

table_BP <- GenTable(GOdata_BP, parentChild = resultFis_BP, orderBy = "parentChild", ranksOf = "parentChild", topNodes = 200)


# Cellular Component
####################

resultFis_CC <- runTest(GOdata_CC, algorithm = "parentchild", statistic = "fisher")

table_CC <- GenTable(GOdata_CC, parentChild = resultFis_CC, orderBy = "parentChild", ranksOf = "parentChild", topNodes = 200)


# filter out significant GO terms:
###################################

MF <- subset(table_MF, parentChild < 0.05)
BP <- subset(table_BP, parentChild < 0.05)
CC <- subset(table_CC, parentChild < 0.05)

MF$ontology<-"MF"
BP$ontology<-"BP"
CC$ontology<-"CC"

sig <- rbind(MF, BP, CC)


# write out tables:
###################
#write.table(sig, "./Parg.sig.GOterms.fullGeneSpace.txt", col.names=T, row.names=F, quote=F, sep="\t")
#write.table(sig, "./Ains.sig.GOterms.fullGeneSpace.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(sig, "./Acha.sig.GOterms.fullGeneSpace.txt", col.names=T, row.names=F, quote=F, sep="\t")
#write.table(sig, "./Parg.sig.GOterms.BUSTED.smallGeneSpace.FDR.txt", col.names=T, row.names=F, quote=F, sep="\t")

# in bash:
############
#cd /Users/lukas/CSE/inquilineGenomics/SoS/absREL2/GOenrichment
#cat Parg.sig.GOterms.fullGeneSpace.txt | awk '{print $1" "$(NF-1)}'

# => insert this lists incl p-values in REVIGO!
