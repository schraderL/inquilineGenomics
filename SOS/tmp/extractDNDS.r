#install.packages("rjson")
library("rjson")

folder<-"/Users/lukas/sciebo/inquilineGenomics18/SOS/results/"
files<-dir(folder,pattern =".json")

omegas<-list()
branchNames<-list()
pvalues<-list()
others<-list()
dNs<-list()
dSs<-list()

for (q in 1:length(files)){
  OG<-gsub(pattern="(OG[0-9]{7})\\..*","\\1",perl=T,x=files[q])
  result <- fromJSON(file = paste(folder,files[q],sep=""))
  
  omega<-NA
  pvals<-NA
  other<-NA
  dN<-NA
  dS<-NA
  for (i in 1:length(result$`branch attributes`$`0`)){
    omega[i]<-result$`branch attributes`$`0`[[i]]$`Baseline MG94xREV omega ratio`
    pvals[i]<-result$`branch attributes`$`0`[[i]]$`Corrected P-value`
    other[i]<-result$`branch attributes`$`0`[[i]]$`Baseline MG94xREV`
    dN[i]<-result$`branch attributes`$`0`[[i]]$`Rate Distributions`[[1]][1]
    dS[i]<-result$`branch attributes`$`0`[[i]]$`Rate Distributions`[[1]][2]
    names(omega)[i]<-names(result$`branch attributes`$`0`)[i]
    names(pvals)[i]<-names(result$`branch attributes`$`0`)[i]
    omegas[[OG]]<-omega
    pvalues[[OG]]<-pvals
    others[[OG]]<-other
    dNs[[OG]]<-dN
    dSs[[OG]]<-dS
  }
  names(pvalues[[OG]])<-gsub(x=names(pvals),pattern="([A-Z]{4})[0-9]{5}","\\1",perl=T)
  names(others[[OG]])<-gsub(x=names(other),pattern="([A-Z]{4})[0-9]{5}","\\1",perl=T)
  branchNames[[OG]]<-names(omega)
  
  names(omegas[[OG]])<-gsub(x=names(omega),pattern="([A-Z]{4})[0-9]{5}","\\1",perl=T)
  
  reorder<-order(names(omegas[[OG]]))
  
  omegas[[OG]]<-omegas[[OG]][reorder]
  pvalues[[OG]]<-pvalues[[OG]][reorder]
  branchNames[[OG]]<-branchNames[[OG]][reorder]
  others[[OG]]<-others[[OG]][reorder]
}

om <- data.frame(do.call("rbind", omegas))
bn <- data.frame(do.call("rbind", branchNames))
pv <- data.frame(do.call("rbind", pvalues))
ot <- data.frame(do.call("rbind", others))
dn <- data.frame(do.call("rbind", dNs))
ds <- data.frame(do.call("rbind", dSs))
boxplot(dn,outline=F)
boxplot(ds,outline=T)
boxplot(om,outline=F,las=2)
boxplot(ot,outline=F,las=2)
boxplot(subset(om,PARG<100 & AHEY<100 & P1 <100 & AINS<100 & ACHA <100 & ACOL<100 & AECH<100 & A1<100 & A2<100),outline=F,las=2)

