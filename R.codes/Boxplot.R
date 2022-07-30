##YJX
mydata<-read.csv("NP-70-fixed.csv",header = TRUE)
head(mydata)
mydata$Nlevel<-factor(mydata$Nlevel,levels = c("N0","N2","N5","N10","N15","N20","N50"))
library(ggplot2)
newdata<-mydata[,c(3,4,9)]
head(newdata)
newdata<-na.omit(newdata)
ggplot(newdata,aes(x=Nlevel,y=NP))+geom_boxplot(outlier.size = 1,aes(fill=factor(Nadd)),
    position=position_dodge(0.8),size=0.5)+guides(fill=guide_legend(title="Nadd"))
+theme_light()##双数据系列的箱型图
##PY 21个处理
mydata<-read.csv("2019-enzyme-fixed.csv",header = TRUE)
head(mydata)
mydata$Tre<-factor(mydata$Tre,levels = c("1","2","3","11","12","13","21","22","23","31","32","33","41","42","43","51","52","53","61","62","63"))
library(ggplot2)
library(RColorBrewer)
newdata<-mydata[,c(3,12)]
head(newdata)
ggplot(newdata,aes(Tre,PPO))+
  geom_boxplot(aes(fill=Tre),notch=FALSE)+
  #geom_jitter(binaxis="y",position=position_jitter(0.3),stackdir="center",dotsize=0.4)#+
 # scale_fill_manual(values= c(brewer.pal(15,"Set2")[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]))#+
  theme_light()+
  labs(title="PY",
  x="Number of added element",y ="The activities of PPO (nmol/g h)",fill ="Factor")# fill为修改图例标题

##PY 8个处理
mydata<-read.csv("2019-enzyme-fixed.csv",header = TRUE)
head(mydata)
mydata$Factor<-factor(mydata$Factor,levels = c("0","2","3","4","5","6","7","8"))
library(ggplot2)
library(RColorBrewer)
newdata<-mydata[,c(4,12)]
head(newdata)
newdata<-na.omit(newdata)
ggplot(newdata,aes(Factor,PPO))+
  geom_boxplot(aes(fill=Factor),notch=FALSE)+
  geom_jitter(binaxis="y",position=position_jitter(0.05),stackdir="center",dotsize=0.4)+
  scale_fill_manual(values= c(brewer.pal(15,"Set2")[c(1,2,3,4,5,6,7,8)]))+
  theme_light()+
  labs(title="PY",
       x="Number of added element",y ="The activities of PPO [nmol/(g.h)]",fill ="Factor")# fill为修改图例标题

