data1<-read.csv("mydatan+relation.csv")
library(ggplot2)
library(RColorBrewer)
library(reshape2)
mat<-round(cor(data1),1)
mydata<-melt(mat)
colnames(mydata)<-c("Var1","Var2","value")
mydata$AbsValue<-abs(mydata$value)
ggplot(mydata,aes(x=Var1,y=Var2))+
  geom_point(aes(size=AbsValue,fill=value),shape=21,colour="black")+
  scale_fill_gradientn(colours = c(brewer.pal(7,"Set1")[2],"white",brewer.pal(7,"Set1")[1]),na.value = NA)

