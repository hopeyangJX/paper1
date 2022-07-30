##yjx
newdata<-read.csv("bacteria20172019.csv",header = TRUE)
fix(newdata)
newdata<-na.omit(newdata)
#fix(newdata)
#install.packages("nlme")
head(newdata)
library(nlme)
head(newdata)
m1.nlme=lme( lpd ~Nadd*lNlevel*Year,random = ~1|Block,data = newdata)
anova(m1.nlme)##提取方差分析的结果
summary(m1.nlme)##提取混合线性模型的结果
shapiro.test(m1.nlme$residuals)##对模型的残差的正态性检验
library(sjPlot)
##py
mydata<-read.csv("2019-enzyme-fixed.csv",header = TRUE)
head(mydata)
newdata<-mydata[,c(2,4,12)]
head(newdata)
newdata<-na.omit(newdata)
head(newdata)
#fix(newdata)
library(nlme)
m1.nlme=lme(PPO~Factor,random = ~1|Block,data = newdata)
#anova(m1.nlme)##提取方差分析的结果
summary(m1.nlme)##提取混合线性模型的结果


