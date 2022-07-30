library(gcookbook)#加载包
library(ggplot2)
setwd("C:\\Users\\liuzh\\Desktop\\博士实验\\2018年10月土壤种子库\\PCOI")
otu<- read.delim('1.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
library(vegan)
distance <- vegdist(otu, method = 'bray')
pcoa <- cmdscale(distance, k = (nrow(otu) - 1), eig = TRUE)
ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')
summary(pcoa)
pcoa$eig
point <- data.frame(pcoa$point)
species <- wascores(pcoa$points[,1:2], otu)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_site <- data.frame({pcoa$point})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
#导入分组数据
mydata=read.table("clipboard",header = T)
sample_site <- merge(sample_site, mydata, by = 'names', all.x = TRUE)
sample_site$deal <- factor(sample_site$deal, levels = c('0', '1',"2","3","5","10","15","20","50"))
library(plyr)
group_border <- ddply(sample_site, 'deal', function(df) df[chull(df[[2]], df[[3]]), ])#这里的deal是分组里的处理那个列的名称
library(ggplot2)

pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2, group = deal)) +
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.key = element_rect(fill = 'transparent')) + #去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_polygon(data = group_border, aes(fill = deal),alpha = 0.3) + #绘制多边形区域
  geom_point(aes(color = deal), size = 1.5, alpha = 0.8) + #可在这里修改点的透明度、大小
  #geom_point( size = 1.5, alpha = 0.8,color="yellow") + #可在这里修改点的透明度、大小,没有加点的形状，将点的颜色统一设为黄色了
  scale_shape_manual(values = c(17, 16,15,14,13,12,11,10,9)) + #可在这里修改点的形状
  scale_color_manual(values = c('yellow', 'orange', 'red', 'red4','#C673FF2E', '#73D5FF2E', '#49C35A2E', '#FF985C2E','#C673FF')) + #可在这里修改点的颜色
  scale_fill_manual(values = c('#C673FF2E', '#73D5FF2E', '#49C35A2E', '#FF985C2E','yellow', 'orange', 'red', 'red4','#C673FF')) + #可在这里修改区块的颜色
  guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) + #设置图例展示顺序
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%'))+
  #可通过修改下面四句中的点坐标、大小、颜色等，修改“A、B、C、D”标签
  annotate('text', label = '0', x = -0.31, y = -0.15, size = 5, colour = '#C673FF') +
  annotate('text', label = '1', x = -0.1, y = 0.3, size = 5, colour = '#73D5FF') +
  annotate('text', label = '2', x = 0.1, y = 0.15, size = 5, colour = '#49C35A') +
  annotate('text', label = '3', x = 0.35, y = 0, size = 5, colour = '#FF985C')+
  annotate('text', label = '5', x = -0.31, y = -0.15, size = 5, colour = '#C673FF') +
  annotate('text', label = '10', x = -0.1, y = 0.3, size = 5, colour = '#73D5FF') +
  annotate('text', label = '15', x = 0.1, y = 0.15, size = 5, colour = '#49C35A') +
  annotate('text', label = '20', x = 0.35, y = 0, size = 5, colour = '#FF985C')+
  annotate('text', label = '50', x = 0.35, y = 0, size = 5, colour = '#FF985C')+
  theme_bw()
pcoa_plot
library(permute)
library(vegan)
setwd("C:\\Users\\liuzh\\Desktop\\博士实验\\2018年10月土壤种子库\\PCOI")
A<-read.delim('1.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
A<- data.frame(t(A))
head(A)
B=read.table("clipboard",header = T)#导入分组数据
adonis_result_otu <- adonis(A~deal, B, permutations = 999, distance = 'bray') #整体分析处理是否有差异
adonis_result_otu
summary(adonis_result_otu)
otuput <- data.frame(adonis_result_otu$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'PERMANOVA.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
#进行两两比较
group_name <- unique(B$deal)

adonis_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(B, deal %in% c(group_name[i], group_name[j]))
    otu_ij <- A[group_ij$names,]
    adonis_result_otu_ij <- adonis(otu_ij~deal, group_ij, permutations = 999, distance = 'bray')     #随机置换检验 999 次
    adonis_result_two <- rbind(adonis_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', unlist(data.frame(adonis_result_otu_ij$aov.tab, check.names = FALSE)[1, ])))
  }
}
adonis_result_two <- data.frame(adonis_result_two, stringsAsFactors = FALSE)
names(adonis_result_two) <- c('B', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
#可选标记 “*” 显著性
for (i in 1:nrow(adonis_result_two)) {
  if (adonis_result_two[i, 'Pr (>F)'] <= 0.001) adonis_result_two[i, 'Sig'] <- '***'
  else if (adonis_result_two[i, 'Pr (>F)'] <= 0.01) adonis_result_two[i, 'Sig'] <- '**'
  else if (adonis_result_two[i, 'Pr (>F)'] <= 0.05) adonis_result_two[i, 'Sig'] <- '*'
}
#输出
write.table(adonis_result_two, 'PERMANOVA.result_two.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
多元置换方差分析
