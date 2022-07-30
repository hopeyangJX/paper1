library(vegan)
library(ggplot2)
sp<-read.csv("Fungi20172019.csv",header=TRUE,sep=",")
head(sp)
sp$NY<-factor(sp$NY)
sp$Nlevel<-factor(sp$Nlevel)
t<-factor(sp$NY,levels=c("2017Ncont","2017Ncess","2019Ncont","2019Ncess",labels=c("2017Ncont","2017Ncess","2019Ncont","2019Ncess")))
shape_NY<-c("2017Ncont"="0","2017Ncess"="22","2019Ncont"="2","2019Ncess"="24")
shape<-shape_NY[as.numeric(factor(t))]
ggplot(sp, aes(x=NMDS1, y=NMDS2,colour=Nlevel,shape=NY))+
  geom_point(alpha = 2.0,size=5)+
  scale_color_manual(values=c("blue", "red","green","pink","yellow","orange","purple"))+
    scale_shape_manual(values=c(15,0,19,1))+
  stat_ellipse(level = 0.95, show.legend = F,size=0.5,linetype=2)+theme_bw()+
  theme(panel.grid=element_blank())






## 置换多元方差分析 ##
library(permute)
library(vegan)
#abund_table <- read.csv("E:\\R\\R_Microbes\\Bacteria_OTU_20172018.CSV.csv",row.names=1)
abund_table <- read.csv("bacteria20172019Otu.csv",row.names=1)
head(abund_table )
#group<-read.csv("E:\\R\\R_Microbes\\Group.csv")
group<-read.csv("bacteria20172019Group.csv")
head(group)
group<-data.frame((group[,c(2:4)]),row.names=rownames(abund_table))
adonis(abund_table ~ Nadd*Nlevel*Year, data=group, permutations=9999)

ggsave("1-Pathotroph NMDS.tiff", scale = 1, width = 145, height = 120, units="mm",dpi = 300)
