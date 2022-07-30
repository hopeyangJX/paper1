# Network construction
library(igraph)
library(psych)
library(vegan)
library(Matrix)
library(RColorBrewer)

setwd("D:/Study/PhD/Writing/Extreme drought/HST.AMF.2017/Data analyses/results/AMF/networks/species.net")

otu.table <- read.csv(file = "AMF_plant_ASVtab.csv", header = T, row.names = 1)
head(t(otu.table))

# Filtering low abundacne OTUs, maintaining OTUs appear more than 6 times of all the samples
otu.table.filter <- otu.table[, specnumber(t(otu.table)) >= 3]
print(c(ncol(otu.table), "versus", ncol(otu.table.filter))) # OTUs left

# caculate pairwise 
occor <- corr.test(otu.table.filter, use = "pairwise", method = "spearman", adjust = "BH", alpha = .05)
occor.r <- occor$r
occor.p <- occor$p

#  Filter the association based on p-values and level of correlations
occor.r[occor.p > 0.01 | abs(occor.r) < 0.3] = 0

# Create igraph object
net.graph <- graph_from_adjacency_matrix(occor.r, mode = "undirected", weighted = TRUE, diag = FALSE)
net.graph

# Creating a vector to remove the isolated nodes (nodes with no interactions)
bad.vs <- V(net.graph)[degree(net.graph) == 0] 

# Removing the isolated nodes from the graph object using the function delete.vertices()
net.graph <- delete.vertices(net.graph, bad.vs)
net.graph

# Assign the igraph weight attribute to net.graph.weight
net.graph.weight <- E(net.graph)$weight

# remove the weight before plotting (may influence the layouts)
E(net.graph)$weight <-NA

edge.list <- as_edgelist(net.graph)
write.csv(edge.list, file = "species.net.edgelist.filter3.csv")
# Plot the graph object
set.seed(123)
plot(net.graph, main = "AMF-PLANT NETWORK FILTER5", edge.lty = 1, margin=c(0,0,0,0),
     vertex.size = 5, vertex.frame.color = NA, edge.curved = T, edge.width = 1, 
     layout = layout_with_fr(net.graph), vertex.label = NA)

# compute node degrees (links) and use that to set node size:
deg <- degree(net.graph, mode = "all")
V(net.graph)$size <- deg*2
plot(net.graph, layout = layout_with_fr(net.graph), main = "AMF-PLANT NETWORK FILTER5", vertex.label = NA)

mycolors <- brewer.pal(12, "Paired") #library(RColorBrewer)
mycolors
# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"

node.attr <- read.csv(file = "groups.csv", header = T, row.names = 1)
# set vertices color
net.graph.col <- node.attr[V(net.graph)$name, ]
net.graph.col$Genus <- factor(net.graph.col$Genus, levels = c("Ambispora", "Archaeospora", "Claroideoglomus", "Diversispora", "Dominikia", "Glomus", 
                                                              "Kamienskia", "Paraglomus", "Funneliformis", "Rhizophagus", "Septoglomus", "Unclassified", "plant"))
levels(net.graph.col$Genus) <- c("#A6CEE3", "#FFFF99", "gray80", "gray30", "#1F78B4", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FB9A99", "#B15928", "#33A02C")
V(net.graph)$color <- as.character(net.graph.col$Genus)

plot(net.graph, edge.curved = T, vertex.frame.color = NA, layout = layout_with_fr(net.graph), 
     main = "Plant-AM fungal co-occurrence network", vertex.label = NA, edge.width	= net.graph.weight*3, 
     edge.color = "gray80", vertex.size = deg*1.5)

legend(x = -1.5, y = -0.5, c("Ambispora", "Archaeospora", "Claroideoglomus", "Diversispora", "Dominikia", "Glomus", "Kamienskia", "Paraglomus", "Funneliformis", "Rhizophagus", "Septoglomus", "Unclassified", "Plant"), pch=21, col="#777777", pt.bg = levels(net.graph.col$Genus), pt.cex=2, cex=.8, bty="n", ncol=1)

# network attributes
# 1. number of edges
num.edges <- length(E(net.graph))
num.edges

positive.edges <- sum(net.graph.weight > 0) # positive edges
positive.edges

negative.edges <- sum(net.graph.weight < 0) # negative edges
negative.edges

# 2. number of vertices
num.vertices <- length(V(net.graph))
num.vertices

# 3. number of connectance
connectance <- edge_density(net.graph, loops = FALSE) # loops = TRUE means self-connectance i.e A-A, B-B
connectance

# 4. average degree
average.degree <- mean(igraph::degree(net.graph))
average.degree

# 5. average path length
average.path.length <- average.path.length(net.graph) # or mean_distance(net.graph)
average.path.length

# 6.diameter
diameter <- diameter(net.graph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter

# 7. adge connectivity / group adhesion
edge.connectivity <- edge_connectivity(net.graph)
edge.connectivity

# 8. clustering coefficient
clustering.coefficient <- transitivity(net.graph)
clustering.coefficient

no.clusters <- no.clusters(net.graph)
no.clusters

# 9. betweenness centralization
centralization.betweenness <- centralization.betweenness(net.graph)$centralization
centralization.betweenness

# 10. degree centralization
centralization.degree <- centralization.degree(net.graph)$centralization
centralization.degree