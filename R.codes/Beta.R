library(vegan)
dis <- read.delim('bray.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
dis <- as.dist(dis)	
otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE)
adonis_result_dis <- adonis(dis~site, group, permutations = 999)
adonis_result_otu <- adonis(otu~site, group, permutations = 999, distance = 'bray') 
dis1 <- vegdist(otu, method = 'bray')
adonis_result_dis1 <- adonis(dis1~site, group, permutations = 999)
adonis_result_dis
summary(adonis_result_dis)
adonis_result_dis$aov.tab 
otuput <- data.frame(adonis_result_dis$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'PERMANOVA.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
group_name <- unique(group$site)

adonis_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
	for (j in (i + 1):length(group_name)) {
		group_ij <- subset(group, site %in% c(group_name[i], group_name[j]))
		otu_ij <- otu[group_ij$names, ]
		adonis_result_otu_ij <- adonis(otu_ij~site, group_ij, permutations = 999, distance = 'bray')	
		adonis_result_two <- rbind(adonis_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', unlist(data.frame(adonis_result_otu_ij$aov.tab, check.names = FALSE)[1, ])))
	}
}
adonis_result_two <- data.frame(adonis_result_two, stringsAsFactors = FALSE)
names(adonis_result_two) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')

for (i in 1:nrow(adonis_result_two)) {
	if (adonis_result_two[i, 'Pr (>F)'] <= 0.001) adonis_result_two[i, 'Sig'] <- '***'
	else if (adonis_result_two[i, 'Pr (>F)'] <= 0.01) adonis_result_two[i, 'Sig'] <- '**'
	else if (adonis_result_two[i, 'Pr (>F)'] <= 0.05) adonis_result_two[i, 'Sig'] <- '*'
}

write.table(adonis_result_two, 'PERMANOVA.result_two.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
