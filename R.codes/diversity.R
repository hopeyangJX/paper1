OTU<-read.csv("Pathothoph.csv",row.names = 1)
head(OTU)
library(vegan)
library(permute)
library(lattice)
Shannon<-diversity(OTU)
Shannon
Inv_Simpson<-diversity(OTU,"inv")
Inv_Simpson
S<-specnumber(OTU)
S
Pielou_evenness<-Shannon/log(S)
Pielou_evenness
Simpson_evenness<-Inv_Simpson/S
Simpson_evenness
report=cbind(Shannon, Inv_Simpson,Pielou_evenness, Simpson_evenness,S)
head(report)
write.table(report,file="Pathothoph.txt",sep="\t",col.names=NA)
bray.dist<-vegdist(OTU)
head(bray.dist)
write.table(as.matrix(bray.dist),file="beta1.txt",sep="\t",col.names=NA)