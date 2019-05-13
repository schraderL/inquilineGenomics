
library(readxl)
library(dplyr)



#my_data <- read_excel("/Users/lukas/sciebo/inquilineGenomics18/geneFamilies/specificFamilies/mrjps/Atta_colombica.MRJPs/Atta_colombica.MRJPs.GeMoMa.details.xlsx")
folder.path<-"/Users/lukas/sciebo/inquilineGenomics18/geneFamilies/specificFamilies/mrjps/"
file.names<-dir(folder.path,pattern="*.xlsx",recursive=T)
data<-list()
for(i in 1:length(file.names)){
  data[[i]] <- read_excel(paste(folder.path,file.names[i],sep=""))
  data[[i]]$species<-gsub(pattern = "(.*?)\\..*","\\1",perl=T,file.names[i])
}


allData<-do.call("rbind", data)
s<-subset(allData,TYPE!="NA",c(species,TYPE))
s2<-apply(s,1,paste,collapse=" ")
st<-s %>% group_by(species,TYPE) %>% tally()
write.table(st,paste(folder.path,"Overview.tsv",sep=""),sep="\t",quote=F,row.names=F)

