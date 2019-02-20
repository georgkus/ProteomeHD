## This script plots a set of correlation networks using the igraph and ggnet2 packages

# Load the necessary libraries
library(data.table); library(ggplot2); library(igraph); library(GGally); library(intergraph)


#### Prepare the data ####

# Load the data
DT <- fread("coregulation_scores.csv")

# Which co-regulation score cut-off corresponds to the top 0.5% of protein pairs
score_cut_off_0.5 <- DT[ order(-coregulation_score) ][ floor( DT[,.N] * 0.005 ) , coregulation_score ]

# Expand into a full "correlation" matrix
DT <- rbind(DT, DT[, .(Protein_1 = Protein_2, Protein_2 = Protein_1, coregulation_score)])
DT <- dcast(DT, Protein_1 ~ Protein_2, value.var = "coregulation_score")
DT <- as.matrix( data.frame(DT, row.names = "Protein_1"))

# Remove values below cut-off by setting them to zero
DT[ DT < score_cut_off_0.5 ] <- 0


#### Create the network (igraph) ####

# Create a weighted, undirected network using igraph
net <- graph_from_adjacency_matrix( DT, mode = "undirected", weighted = TRUE, diag = FALSE )

# Inspect network properties (should be U for undirected, N for named, W for weighted and have 5,013 nodes and 62,812 edges)
net


#### Plot the network with different layouts (ggnet2) ####

# Default layout (Fruchtermanreingold)
pFR <- ggnet2( net,  mode = "fruchtermanreingold" ,
             node.shape = 16, node.size = 0.3, node.color = "steelblue", node.alpha = 0.7,
             edge.size = 0.15, edge.color = "black" , edge.alpha = 0.7)

ggsave("Net_Fruchtermanreingold.png", pFR, width = 5, height = 5, units = "cm" , dpi = 900 )

# Random (Gaussian donut)
pcircrand <- ggnet2( net,  mode = "circrand" ,
                   node.shape = 16, node.size = 0.3, node.color = "steelblue", node.alpha = 0.7,
                   edge.size = 0.15, edge.color = "black" , edge.alpha = 0.7)

ggsave("Net_circrand.png", pcircrand, width = 5, height = 5, units = "cm" , dpi = 900 )

# Kamadakawai
pkamadakawai <- ggnet2( net,  mode = "kamadakawai" ,
                     node.shape = 16, node.size = 0.3, node.color = "steelblue", node.alpha = 0.7,
                     edge.size = 0.15, edge.color = "black" , edge.alpha = 0.7)

ggsave("Net_kamadakawai.png", pkamadakawai, width = 5, height = 5, units = "cm" , dpi = 900 )

# Spring
pspring <- ggnet2( net,  mode = "spring" ,
                node.shape = 16, node.size = 0.3, node.color = "steelblue", node.alpha = 0.7,
                edge.size = 0.15, edge.color = "black" , edge.alpha = 0.7)

ggsave("Net_spring.png", pspring, width = 5, height = 5, units = "cm" , dpi = 900 )









