library(vegan)
phylum <- read.delim('phylum_table.txt', row.names = 1, sep = '\t',
                     stringsAsFactors = FALSE, check.names = FALSE)
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE)
#head(phylum)
#head(env)
#物种数据 Hellinger 预转化（处理包含很多 0 值的群落物种数据时，推荐使用）
phylum_hel <- decostand(phylum, method = 'hellinger')#采用total标准化以后再取平方根
#head(phylum_hel)
#使用全部的环境数据
rda_tb <- rda(phylum_hel~., env, scale = FALSE)
summary(rda_tb)
#例如控制土壤 pH 影响后（pH 作为协变量），观测其它环境因素的影响；
#物种数据 Hellinger 预转化
#rda_part <- rda(phylum_hel~TC+DOC+SOM+TN+NO3+NH4+AP+AK+Condition(pH), data = env, scale = FALSE)
#summary(rda_part)
##db-RDA（以 Bray-curtis 距离为例，注意这里直接使用了原始的物种丰度矩阵）
#计算距离
#dis_bray <- vegdist(phylum, method = 'bray')
#PCoA 排序
#pcoa <- cmdscale(dis_bray, k = nrow(phylum) - 1, eig = TRUE, add = TRUE)
#提取 PCoA 样方得分（坐标）
#pcoa_site <- pcoa$point
#db-RDA，使用全部的环境数据
#rda_db <- rda(pcoa_site, env, scale = FALSE)
#或者，capscale() 提供了直接运行的方法
rda_db <- capscale(phylum~., env, distance = 'bray', add = TRUE)
summary(rda_db)
#若基于欧氏距离，则和常规 RDA 的结果一致
#物种数据 Hellinger 转化后，分别执行 tb-RDA 与使用欧氏距离的 db-RDA
rda_tb_test <- rda(phylum_hel~., env)
summary(rda_tb_test)
rda_db_test <- capscale(phylum_hel~., env, distance = 'euclidean')
summary(rda_db_test)
par(mfrow = c(1, 2))
plot(rda_tb_test, scaling = 1)
plot(rda_db_test, scaling = 1)
#查看统计结果信息，以 I 型标尺为例
rda_tb.scaling1 <- summary(rda_tb, scaling = 1)
rda_tb.scaling1   
#作图查看排序结果，三序图，包含 I 型标尺和 II 型标尺
par(mfrow = c(1, 2))
plot(rda_tb, scaling = 1, main = 'I 型标尺', display = c('wa', 'sp', 'cn'))
rda_sp.scaling1 <- scores(rda_tb, choices = 1:2, scaling = 1, display = 'sp')
arrows(0, 0, rda_sp.scaling1[ ,1], rda_sp.scaling1[ ,2], length =  0, lty = 1, col = 'red')
plot(rda_tb, scaling = 2, main = 'II 型标尺', display = c('wa', 'sp', 'cn'))
rda_sp.scaling2 <- scores(rda_tb, choices = 1:2, scaling = 2, display = 'sp')
arrows(0, 0, rda_sp.scaling2[ ,1], rda_sp.scaling2[ ,2], length =  0, lty = 1, col = 'red')     
#隐藏物种信息，以 I 型标尺为例展示双序图，并查看分别使用物种加权计算的样方坐标以及拟合的样方坐标的差异

par(mfrow = c(1, 2))
plot(rda_tb, scaling = 1, main = 'I 型标尺，加权', display = c('wa', 'cn'))
plot(rda_tb, scaling = 1, main = 'I 型标尺，拟合', display = c('lc', 'cn'))
##RDA 结果提取
#scores() 提取排序得分（坐标），以 I 型标尺为例，分别提取前两个约束轴中的样方（排序对象）、物种（响应变量）、环境因子（解释变量）的排序坐标
rda_site.scaling1 <- scores(rda_tb, choices = 1:2, scaling = 1, display = 'wa')	#使用物种加权和计算的样方坐标
rda_sp.scaling1 <- scores(rda_tb, choices = 1:2, scaling = 1, display = 'sp')
rda_env.scaling1 <- scores(rda_tb, choices = 1:2, scaling = 1, display = 'cn')
#若需要输出在本地，以样方坐标为例，以 csv 格式为例
write.csv(data.frame(rda_site.scaling1), 'rda_site.scaling1.csv')
#coef() 提取典范系数
rda_coef <- coef(rda_tb)
#其它提取方式，首先查看结果包含的所有信息
names(rda_tb.scaling1)
#然后提取相关内容，例如我们想提取前两轴的“样方得分”，可如此做
rda_site.scaling1 <- rda_tb.scaling1$sites[ ,1:2]
#该结果和上述“scores(rda_tb, choices = 1:2, scaling = 1, display = 'wa')”的结果是一致的
#同理，可如此输出在本地
write.csv(data.frame(rda_site.scaling1), 'rda_site.scaling1.csv')
##R2 校正
#RsquareAdj() 提取 R2
r2 <- RsquareAdj(rda_tb)
rda_noadj <- r2$r.squared	#原始 R2
rda_adj <- r2$adj.r.squared	#校正后的 R2
rda_noadj
rda_adj
##置换检验
#所有约束轴的置换检验，以 999 次为例
rda_tb_test <- anova(rda_tb, permutations = 999)
summary(rda_tb_test)
#或者使用
rda_tb_test <- anova.cca(rda_tb, step = 1000)
#各约束轴逐一检验，以 999 次为例
rda_tb_test_axis <- anova(rda_tb, by = 'axis', permutations = 999)
#或者使用
rda_tb_test_axis <- anova.cca(rda_tb, by = 'axis', step = 1000)

#p 值校正（Bonferroni 为例）
rda_tb_test_axis$`Pr(>F)` <- p.adjust(rda_tb_test_axis$`Pr(>F)`, method = 'bonferroni')
##断棍模型和 Kaiser-Guttman 准则确定残差轴
#提取残差特征值
pca_eig <- rda_tb$CA$eig

#Kaiser-Guttman 准则
pca_eig[pca_eig > mean(pca_eig)]
#断棍模型
n <- length(pca_eig)
bsm <- data.frame(j=seq(1:n), p = 0)
bsm$p[1] <- 1/n
for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
bsm$p <- 100*bsm$p/n
bsm
# 绘制每轴的特征根和方差百分比 
par(mfrow = c(2, 1))
barplot(pca_eig, main = '特征根', col = 'bisque', las = 2)
abline(h = mean(pca_eig), col = 'red')
legend('topright', '平均特征根', lwd = 1, col = 2, bty = 'n')
barplot(t(cbind(100 * pca_eig/sum(pca_eig), bsm$p[n:1])), beside = TRUE, main = '% 变差', col = c('bisque', 2), las = 2)
legend('topright', c('% 特征根', '断棍模型'), pch = 15, col = c('bisque', 2), bty = 'n')
##变量选择
#计算方差膨胀因子
vif.cca(rda_tb)
#vegan 包 ordistep() 前向选择，基于 999 次置换检验
rda_tb_forward_p <- ordistep(rda(phylum_hel~1, env, scale = FALSE), scope = formula(rda_tb), direction = 'forward', permutations = 999)
#vegan 包 ordiR2step() 前向选择，基于 999 次置换检验
rda_tb_forward_r <- ordiR2step(rda(phylum_hel~1, env, scale = FALSE), scope = formula(rda_tb), R2scope = rda_adj, direction = 'forward', permutations = 999)
#以 rda_tb 和 rda_tb_forward_r 为例，简要绘制双序图比较变量选择前后结果
par(mfrow = c(1, 2))
plot(rda_tb, scaling = 1, main = '原始模型，I 型标尺', display = c('wa', 'cn'))
plot(rda_tb_forward_r, scaling = 1, main = '前向选择后，I 型标尺', display = c('wa', 'cn'))
#细节部分查看
summary(rda_tb_forward_r, scaling = 1)
#比较选择前后校正后 R2 的差异
RsquareAdj(rda_tb)$adj.r.squared
RsquareAdj(rda_tb_forward_r)$adj.r.squared
#packfor 包 forward.sel() 前向选择
install.packages("adespatial")
library(adespatial)
forward.sel(phylum_hel, env, adjR2thresh = rda_adj)
##变差分解 varpart()，以前向选择后的简约模型 rda_tb_forward_r 为例（包含 6 个环境解释变量）
#以两组环境变量为例，运行变差分解
rda_tb_forward_vp <- varpart(phylum_hel, env['pH'], env[c('DOC', 'SOM', 'AP', 'AK', 'NH4')])
rda_tb_forward_vp
plot(rda_tb_forward_vp, digits = 2, Xnames = c('pH', 'CNPK'), bg = c('blue', 'red'))
#查看前向选择中被剔除的环境变量“TC”，与这 6 个被保留的环境变量之间解释变差的“共享程度”
rda_tb_forward_vp <- varpart(phylum_hel, env['TC'], env[c('pH', 'DOC', 'SOM', 'AP', 'AK', 'NH4')])
plot(rda_tb_forward_vp, digits = 2, Xnames = c('TC', 'forward_env'), bg = c('blue', 'red'))
#解释变差的置换检验，以 pH 所能解释的全部变差为例；999 次置换
anova(rda(phylum_hel, env['pH']), permutations = 999)
#若考虑 pH 单独解释的变差部分，需将其它变量作为协变量；999 次置换
anova(rda(phylum_hel, env['pH'], env[c('DOC', 'SOM', 'AP', 'AK', 'NH4')]), permutations = 999)

##plot() 作图示例，以前向选择后的简约模型 rda_tb_forward_r 为例，展示前两轴，II 型标尺，双序图，默认使用物种加权和计算的样方坐标
pdf('rda_test1.pdf', width = 5, height = 5)
plot(rda_tb_forward_r, choices = c(1, 2), scaling = 2, type = 'n')
text(rda_tb_forward_r, choices = c(1, 2), scaling = 2, dis = 'cn', col = 'blue', cex = 0.8)
points(rda_tb_forward_r, choices = c(1, 2), scaling = 2, pch = 21, bg = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), col = NA, cex = 1.2)
dev.off()
##ggplot2 作图，以前向选择后的简约模型 rda_tb_forward_r 为例，展示前两轴，II 型标尺，双序图，默认使用物种加权和计算的样方坐标
#提取样方和环境因子排序坐标，前两轴，II 型标尺
rda_tb_forward_r.scaling2 <- summary(rda_tb_forward_r, scaling = 2)
rda_tb_forward_r.site <- data.frame(rda_tb_forward_r.scaling2$sites)[1:2]
rda_tb_forward_r.env <- data.frame(rda_tb_forward_r.scaling2$biplot)[1:2]
#读取样本分组数据（附件“group.txt”）
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#合并样本分组信息，构建 ggplot2 作图数据集
rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
rda_tb_forward_r.site <- merge(rda_tb_forward_r.site, group, by = 'sample')
rda_tb_forward_r.env$sample <- NA
rda_tb_forward_r.env$group <- rownames(rda_tb_forward_r.env)
#ggplot2 作图
library(ggplot2)

p <- ggplot(rda_tb_forward_r.site, aes(RDA1, RDA2)) +
  geom_point(aes(color = group)) +
  scale_color_manual(values = c('red', 'orange', 'green3')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) + 
  labs(x = 'RDA1 (42.91%)', y = 'RDA2 (9.80%)') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_tb_forward_r.env, aes(x = 0,y = 0, xend = RDA1,yend = RDA2), arrow = arrow(length = unit(0.1, 'cm')), size = 0.3, color = 'blue') +
  geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = group), color = 'blue', size = 3)

ggsave('rda_test2.pdf', p, width = 5, height = 4)
ggsave('rda_test2.png', p, width = 5, height = 4)
plot(p)

