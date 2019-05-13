#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require("rphast")
library(RColorBrewer)
library(methods)


setwd("/Users/lukas/CSE/inquilineGenomics/ORevo/processing")
GFF<-"/Users/lukas/CSE/inquilineGenomics/ORevo/results/Acep.final.OR.gff3"
genome<-"/Users/lukas/CSE/attine_genomes/Aech/Aech_2.0_scaffolds.fa"
TEgff<-"/Users/lukas/CSE/attine_genomes/Aech/Aech_v1.0_repbase.gff"
##########################################
# promotor annotation function
##########################################
annotate.promotor <- function(feats, promotorsize) {
#	starts<-feats[ feats$feature=="start_codon", ]
######################################
# test if this works !!!!!!!!!!!!!!!!
######################################
# UTR5 not part of promotor
#http://biology.stackexchange.com/questions/38783/are-eukaroytic-promoters-located-in-the-5-utr-region
	starts<-feats[ feats$feature=="mRNA", ]
	size<-promotorsize
	label<-paste(promotorsize/1000,"k.promotor",sep="")
	starts$feature<-as.character(starts$feature)
	for (i in 1:length(starts$strand)){
		if (starts$strand[i]=="+"){
			starts$end[i]<-starts$start[i]-1
			starts$start[i]<-starts$start[i]-size
			starts$feature[i]<-label
		}else{
			starts$start[i]<-starts$end[i]+1
			starts$end[i]<-starts$end[i]+size
			starts$feature[i]<-label
		}
		if (starts$end[i] < 0){starts$end[i]<-0}
		if (starts$start[i] < 0){starts$start[i]<-0}
	}
	starts$feature<-as.factor(starts$feature)
	return(starts)
}


################################################################################
# perfect feature object
################################################################################
##########################################

	feats <- read.feat(GFF)
	tmp<-feats[feats$feature=="mRNA",]
	feats<-add.signals.feat(feats)

	feats<-feats[order(feats$feature,decreasing=T),]
	feats<-feats[order(feats$start),]
	feats<-feats[order(feats$seq),]
	promotor<-annotate.promotor(feats,2000)
	feats <- rbind.feat(feats, promotor)
	feats<-feats[order(feats$feature,decreasing=T),]
	feats<-feats[order(feats$start),]
	feats<-feats[order(feats$seq),]


	##########################################
	# read TEs
	##########################################
	tes <- read.feat(TEgff)
	featsTE <- rbind.feat(feats, tes)
	featsTE<-feats[order(feats$feature,decreasing=T),]
	featsTE<-feats[order(feats$start),]
	featsTE<-feats[order(feats$seq),]


	##########################################
	# write feat file
	##########################################

	table(feats$feature)
	write.feat(feats,paste(folders[1],"Aech.",chromosome,".features.gff",sep=""))

}else{
	feats <- read.feat(paste(folders[1],"Aech.",chromosome,".features.gff",sep=""))
}


library(GenomicFeatures)
txdb <- makeTranscriptDbFromGFF(file = GFF, format = "gff")
