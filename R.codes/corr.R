data1<-read.csv("Fungidatacorr.csv",header = T)
#View(data1)
data1<-na.omit(data1)
head(data1)
mydata<-data1[,c(24,25,26,27,28,29,30,31)]
head(mydata,6)
res<-cor(mydata)##计算相关矩阵并赋值给res
round(res,2)##保留两位小数
library(Hmisc)##计算相关系数的显著性
res2<-rcorr(as.matrix(mydata),type="pearson")##选择“pearson"或者”spearman"
res2
res2$r
res2$P
library(corrplot)##安装corrplot包后，载入
corrplot(res2$r,type="upper",order="hclust",p.mat=res2$P,sig.level=0.05,insig="blank")
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
chart.Correlation(mydata, histogram=TRUE, pch=12)

