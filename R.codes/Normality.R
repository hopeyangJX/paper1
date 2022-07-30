mydata<-read.csv("Fungidata.csv",header = TRUE)
head(mydata)
x<-mydata[,9]
head(x)
shapiro.test(x[11:21])
#y<-x[1:10]
#head(y)
bartlett.test(mydata$Frichness~mydata$Nlevel)
fligner.test(mydata$Frichness~mydata$Nlevel)
