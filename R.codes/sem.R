###Ncess
data1<-read.csv("Fungidata-N+.csv")
#fix(data1)
data1<-na.omit(data1)
library(lavaan)
head(data1)
sem1<-'LnpH~Lnlevel
NPC1~Lnlevel+LnpH
Lbiomass~LnpH+NPC1+Lnlevel
PatPC1~LnpH+NPC1+Lbiomass'
fitsem1<-sem(sem1,data=data1)
summary(fitsem1,rsq=T,fit.measures=TRUE,standardized=T)
resid(fitsem1,type="standardized")
modindices(fitsem1)
summary(fitsem1,modindices=T)##根据mi是否大于3.84补丢失的路径
summary(fitsem1,standardized=T,rsq=T)
library(semPlot)
semPaths(fitsem1,what="stand",layout="tree2")
###Ncont
data1<-read.csv("mydatan+ 2017.csv")
fix(data1)
library(lavaan)
sem1<-'LNH4~LN
LCEC~LNH4
LNO3~LN+LNH4+LCEC
LB~LN+LCEC
LchiA~LCEC+LNH4
LAOA~LCEC+LB+LNH4+LNO3+LchiA
LAOB~LNH4+LCEC+LN
LnirK~LCEC
LnirS~LCEC+LNH4+LnirK
LnosZ~LCEC+LNH4+LNO3+LnirS'
fitsem1<-sem(sem1,data=data1)
summary(fitsem1,rsq=T,fit.measures=TRUE)
summary(fitsem1,standardized=T,rsq=T)
library(semPlot)
semPaths(fitsem1,what="stand",layout="tree2")
?semPaths
?semPlot

