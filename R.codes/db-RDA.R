# Redundancy Analysis (RDA)
library(adespatial);
packageVersion("adespatial")
library(vegan)
library(ggplot2)
#setwd("D:/HST.AMF.2017/results/AMF/RDA")

# Read data
spe.amf <- read.csv(file = "OTUF20172019.csv", header = TRUE, row.names = 1)
#spe.amf <- t(spe.amf)
head(spe.amf)

env <- read.csv(file = "env.csv", header = TRUE, row.names = 1)
head(env)
View(env)
#envplant <- read.csv(file = "envplant.csv", header = TRUE, row.names = 1)
#head(envplant)
#envchem <- read.csv(file = "envchem.csv", header = TRUE, row.names = 1)
#head(envchem)
#envtreament <- read.csv(file = "envtreatment.csv", header = TRUE, row.names = 1)
#head(envtreament)

# db-RDA
(spe.amf.dbRDA <- capscale(spe.amf ~  Biomass + DOC + NH4 + NO3 + AP + pH + Richness, 
                           env, distance = "bray", comm = spe.amf))

#help(capscale)
#View(spe.amf.dbRDA)
spe.amf.dbRDA
plot(spe.amf.dbRDA)
anova(spe.amf.dbRDA)
anova(spe.amf.dbRDA, permutations = how(nperm = 9999))
anova(spe.amf.dbRDA, permutations = how(nperm = 9999), by = "axis")
summary(spe.amf.dbRDA)

(R2a.amf.dbRDA <- RsquareAdj(spe.amf.dbRDA)$adj.r.squared)
(R2 <- RsquareAdj(spe.amf.dbRDA)$r.squared)
vif.cca(spe.amf.dbRDA)

# using ggplot2 make figures
spe.amf.dbRDA.scaling1 <- summary(spe.amf.dbRDA, scaling = 1) # extract scaling 1
View(spe.amf.dbRDA.scaling1)
spe.amf.dbRDA.site <- data.frame(spe.amf.dbRDA.scaling1$sites)[, 1:2]
View(spe.amf.dbRDA.site)
# spe.amf.dbRDA.site.constraints <- data.frame(spe.amf.dbRDA.scaling1$constraints)[, 1:2]
spe.amf.dbRDA.env <- data.frame(spe.amf.dbRDA.scaling1$biplot)[, 1:2]
View(spe.amf.dbRDA.env)

# read group data
group <- read.csv('group.csv', header = TRUE, row.names = 1)
head(group)
#group$Treatment <- factor(group$Treatment, levels =c("N0", "N2", "N5", "N10", "N15", "N20", "N50"))
head(group)
View(group)
head(group$Treatment)
# merge group information
spe.amf.dbRDA.site$sample <- row.names(spe.amf.dbRDA.site)
View(spe.amf.dbRDA.site$sample)
View(spe.amf.dbRDA.site)
head(group)
spe.amf.dbRDA.site <- merge(spe.amf.dbRDA.site, group, by = 'sample')
head(spe.amf.dbRDA.site)

spe.amf.dbRDA.env$sample <- NA
spe.amf.dbRDA.env$group <- row.names(spe.amf.dbRDA.env)
head(spe.amf.dbRDA.env)
View(spe.amf.dbRDA.env)
spe.amf.dbRDA.env <- spe.amf.dbRDA.env[1:2, ]
head(spe.amf.dbRDA.env)
head(spe.amf.dbRDA.scaling1)
spe.amf.dbRDA.env.tre <- data.frame(spe.amf.dbRDA.scaling1$centroids)[, 1:2]
head(spe.amf.dbRDA.scaling1)
colnames(spe.amf.dbRDA.env.tre) <- c("T1", "T2")
head(spe.amf.dbRDA.env.tre)

p <- ggplot(spe.amf.dbRDA.site, aes(CAP1, CAP2)) +
  geom_point(aes(color = Treatment), size = 5) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  scale_color_manual(values = c("#8DC63F", "#FBB040", "#F06B22", "#39B54A", "#009444")) +
  geom_vline(xintercept = 0, color = "gray", size = 0.5) +
  geom_hline(yintercept = 0, color = "gray", size = 0.5) +
  geom_point(data = spe.amf.dbRDA.env.tre, aes(x = T1, y = T2), color ="blue", size = 5) +
  geom_text(data = spe.amf.dbRDA.env.tre, aes(x = T1, y = T2, label = rownames(spe.amf.dbRDA.env.tre)), color ="blue", size = 3) +
  geom_segment(data = spe.amf.dbRDA.env, 
               aes(x = 0,y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.1, "cm")), 
               size = 0.3, color = "blue") +
  geom_text(data = spe.amf.dbRDA.env, aes(CAP1 * 1.1, CAP2 * 1.1, label = group), color ="blue", size = 3) +
  labs(x = 'db-RDA1 (18.81%)', y = 'db-RDA2 (3.50%)') + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = "black", fill ="transparent"), 
        legend.title = element_text(), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.position = c(0.85, 0.85))
p

p1 <- p + stat_ellipse(aes(color = Treatment), level = 0.65, linetype = 2)
p1

