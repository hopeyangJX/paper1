library(vegan)
library(ggplot2)
library(ggpubr)

#otu <- read.delim('otu.txt',  row.names = 1, sep = '\t')
otu<-read.csv("OTUF20172019.csv",header=TRUE,sep=",")
#otu <- data.frame(t(otu))
otu<-data.frame(otu)
head(otu)

#ngroup <- read.table('group.txt', header = T)
ngroup<-read.csv("OTUF20172019group.csv")
names(ngroup)
attach(ngroup)

##method 1###
bray_dis <- vegdist(otu, method = 'bray')
nmds_dis <- metaMDS(bray_dis, k = 3)
nmds_dis_site <- data.frame(nmds_dis$points)
nmds_dis_site$name <- rownames(nmds_dis_site)
write.table(nmds_dis$points, 'NMDS3.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')

otunmds <- ggplot(data = nmds_dis_site , aes(MDS1,MDS2,label = ngroup.N_level), colour = ngroup.N_level) + 
  theme_bw() +
  geom_point(data = nmds_dis_site, aes(x=MDS1,y=MDS2,colour=ngroup.N_level), size = 6) +
  scale_color_viridis_d(option = "viridis") +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, "bold"), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
  labs(colour="N_level")
otunmds

??plot_model
### method 2#
model_nmds<-metaMDS(otu,distance="bray")
 NMDS = data.frame(MDS1 = model_nmds$points[,1], MDS2 = model_nmds$points[,2], ngroup$N_level)
attach(NMDS)
 print(model_nmds$stress)
head(model_nmds$points)
NMDS$site <- rownames(NMDS)
write.table(model_nmds$points, 'NMDS1.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
# otunmds <- ggplot(data = NMDS , aes(MDS1,MDS2,label = ngroup.N_level), colour = ngroup.N_level) + 
#   theme_bw() +
#   geom_point(data = NMDS, aes(x=MDS1,y=MDS2,colour=ngroup.N_level), size = 6) +
#   scale_color_viridis_d(option = "viridis") +
#   theme(panel.grid = element_blank()) +
#   theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, "bold"), 
#         legend.text = element_text(size = 12), legend.title = element_text(size = 14))+
#   labs(colour="N_level")
# otunmds

