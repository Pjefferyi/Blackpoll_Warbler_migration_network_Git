# Network analysis for the blackpoll warbler migration network
library(anytime)
library(dplyr)
library(ggplot2)
library(sf)
library(rgdal)
library(terra)
library(RColorBrewer)

# Network metrics and plotting
library(igraph)
library(ggnetwork)
library(intergraph)
library(tnet)
library(networktools)
library(qgraph)

# analyses of network communities 
library(clustAnalytics)
library(robin)

#source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Network_construction.R")

# NOTE: THIS SCRIPT IS OPTIMIZED FOR DATA FROM THE FILE: Geolocator_analysis_V3_BI_mask

# Prepare for fall network analysis ############################################

# Load the fall graph
fall.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.graph.edge.list.txt", directed = TRUE)

# Load fall graph node metadata 
meta.fall.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.csv")

# Load fall graph edge weights
fall.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.edge.weights.csv")

# Import polygon of the region affected by the equinox
equinox_region <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Mapping_components/Data/Fall_Equinox_affected_regionV5.shp")

# Import polygons of BPW range
bpw_range <-  read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Birdlife_international_species_distribution/SppDataRequest.shp")

# Plot the fall graph ###########################################################
data("wrld_simpl")
type.palette <- c("#440154FF", "#FDE725FF", "#21908CFF")

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)

plot(as_Spatial(equinox_region), add = T, col = "#D9D5B2", lwd = 0.000000001)
#plot(as_Spatial(bpw_range), col = "#F4F0D3", lwd = 0.0001, add = T)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -35), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)

# Edge colours for spring and fall edges (spring edges should not appear in this plot)
edge.cols <- fall.con.ab %>% mutate(col = case_when(
  edge.type == "fall" ~ adjustcolor("darkgray", alpha.f = 0.9),
  edge.type == "spring" ~ NA))

plot(fall.graph, vertex.size = 300, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*45, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = type.palette[meta.fall.ab$node.type.num], vertex.label.dist = 30,
     add = T, edge.color = edge.cols$col, vertex.label = NA)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.fall.ab$node.type.num)],
       pch = 19, title = "Node type",cex = 0.8)

# Plot the fall graph with node use pie charts #################################

# Create a palette for site use by range region  
use.palette <- list(c("#440154FF", "#FDE725FF", "#21908CFF"))

#create a column that can be converted to a numeric vector
meta.fall.ab <- transform(meta.fall.ab, use.ab.vector = asplit(cbind(use.breeding.ab,  use.stopover.ab, use.nonbreeding.ab), 1))

# Create vector of vertex shapes 
meta.fall.ab <- meta.fall.ab %>% rowwise() %>%
  mutate(shape_single.use = length(which(c(use.breeding.ab, 
                                           use.nonbreeding.ab, 
                                           use.stopover.ab)==0))) %>%
  ungroup() %>%
  mutate(shape_single.use = ifelse(shape_single.use == 2, "circle", "none")) %>%
  mutate(shape_multiple.use = ifelse(shape_single.use== "none", "pie", "none")) %>%
  mutate(shape_colour_multiple.use = case_when(shape_single.use != "none" & use.breeding.ab != 0 ~ type.palette[1] ,
                                         shape_single != "none" & use.stopover.ab != 0 ~ type.palette[2] ,
                                         shape_single != "none" & use.nonbreeding.ab != 0 ~ type.palette[3] ,
                                         .default = NA)) %>%
  mutate(shape_colour_single.use = case_when(shape_single.use!= "circle" & use.breeding.ab != 0 ~ type.palette[1],
                                             shape_single.use != "circle" & use.stopover.ab != 0 ~ type.palette[2],
                                             shape_single.use != "circle" & use.nonbreeding.ab!= 0 ~ type.palette[3],
                                                  .default = NA))


plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)

plot(as_Spatial(equinox_region), add = T, col = "#D9D5B2", lwd = 0.000000001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -35), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)

plot(fall.graph, vertex.size = 300, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*40, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = type.palette[meta.fall.ab$node.type.num], vertex.label.dist = 30,
     add = T, vertex.label = NA, edge.color = adjustcolor("darkgray", alpha.f = 0.6),
     vertex.shape = meta.fall.ab$shape_single.use)

plot(fall.graph, vertex.size = 300, vertex.size2 = 200,
     edge.width = 0, edge.arrow.size = 0, edge.arrow.width = 0,
     vertex.shape = meta.fall.ab$shape_multiple.use, vertex.pie = meta.fall.ab$ use.ab.vector,
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.pie.color = use.palette, vertex.color = NA,
     add = T, vertex.label = NA, edge.color = adjustcolor("darkgray", alpha.f = 0.6))

legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.fall.ab$node.type.num)],
       pch = 19, title = "Node type",cex = 0.8)

## Fall network analysis ###########################################################

## Fall betweenness centrality (exclude breeding sites)
fall.betw.c <- betweenness(fall.graph, directed = T, weights = 1/fall.con.ab$weight)

# plot of the betweenness centrality of each node 
betw.c.palette <- hcl.colors(n = length(seq(0, max(fall.betw.c) + 1)), palette = "Reds 3", rev = T) 
names(betw.c.palette) <- seq(0, max(fall.betw.c) +1)

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(fall.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = betw.c.palette[as.character(round(fall.betw.c))], add = T, vertex.label = NA)

## Fall degree centrality ---
fall.degree.c <- degree(fall.graph, mode = "all")

# plot of the degree centrality of each node 
deg.c.palette <- hcl.colors(n = length(seq(1, max(fall.degree.c ))), palette = "Reds 3", rev = T) 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(fall.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = deg.c.palette[fall.degree.c],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = NA)

# Fall weighed degree centrality ----
fall.strength.c <- strength(fall.graph, mode = "in", weight = fall.con.ab$weight)

# plot of the strength centrality of each node 
deg.c.palette <- hcl.colors(n = length(seq(0, max(fall.strength.c + 0.01), 0.001)), palette = "Reds 3", rev = T) 
names(deg.c.palette) <- seq(0, max(fall.strength.c+ 0.01), 0.001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = deg.c.palette[as.character(round(fall.strength.c, digits = 3))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = round(fall.strength.c, digits = 2))

# Fall degree centrality using the method developed by Opshal et al. ---
E(fall.graph)$weight <- fall.con.ab$weight
fall.edge.list <- cbind(get.edgelist(fall.graph) , round( E(fall.graph)$weight, 3 ))
fall.net <- as.tnet(fall.edge.list, type = "weighted one-mode tnet")

fall.degree.TO <- degree_w(fall.edge.list, measure=c("degree", "output", "alpha"), type = "in", alpha = 0.5)

deg.TO.palette <- hcl.colors(n = length(seq(0, max(fall.degree.TO[,"alpha"]+ 0.0001),0.0001)), palette = "Reds 3", rev = T) 
names(deg.TO.palette) <- seq(0, max(fall.degree.TO[,"alpha"] + 0.0001),0.0001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
     #main = "Degree centrality")

plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = deg.TO.palette[as.character(round(fall.degree.TO[,"alpha"],4))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = round(fall.degree.TO , digits = 4))
 
legend_image <- as.raster(matrix(hcl.colors(n = length(seq(0, max(fall.degree.TO[,"alpha"]), max(fall.degree.TO[,"alpha"])/20)), palette = "Reds 3", rev = F) , ncol=1))
plot(c(0,2),c(-0.01,1),type = 'n', axes = F,xlab = '', ylab = '')
title(main = 'In-degree centrality', cex.main = 0.9, line  = 0.8)
text(x= 1.25, y = seq(0,1,l=5), labels = round(seq(0,max(fall.degree.TO[,"alpha"]),l=5), digits = 1), cex = 1)
rasterImage(legend_image, 0.75, 0, 1,1)

# Fall betweenness centrality using the method developed by Opshal et al. ---
fall.betweenness.TO <- betweenness_w(fall.edge.list, directed = T, alpha = 0.5)

betweenness.TO.palette <- hcl.colors(n = length(seq(0, max(fall.betweenness.TO[,"betweenness"]),0.5)), palette = "Reds 3", rev = T) 
names(betweenness.TO.palette) <- seq(0, max(fall.betweenness.TO[,"betweenness"]),0.5)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
     #main = "Betweenness centrality")

plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = betweenness.TO.palette[as.character(fall.betweenness.TO[,"betweenness"])],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = round(fall.betweenness.TO, digits = 4))

legend_image <- as.raster(matrix(hcl.colors(n = length(seq(0, max(fall.betweenness.TO[,"betweenness"]), 1)), palette = "Reds 3", rev = F) , ncol=1))
plot(c(0,2),c(-0.01,1),type = 'n', axes = F,xlab = '', ylab = '')
title(main = 'Betweenness centrality', cex.main = 0.9, line  = 0.8)
text(x= 1.25, y = seq(0,1,l=5), labels = seq(0,max(fall.betweenness.TO[,"betweenness"]),l=5), cex = 1)
rasterImage(legend_image, 0.75, 0, 1,1)

### community detection using propagating labels ----------  
E(fall.graph)$weight <- fall.con.ab$weight

undirected.fall.graph <- as.undirected(fall.graph, mode = "collapse",
                                      edge.attr.comb = "sum")

fall.communities.lab <- cluster_label_prop(undirected.fall.graph, weights = E(undirected.fall.graph)$weight)

# plot communities
plot(fall.communities.lab, fall.graph)

# plot communities on map
fall.comm.pal <- rainbow(length(seq(1, max(fall.communities.lab$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.lab$membership], add = T)
  
### community detection using the walktrap algorithm -----------
E(fall.graph)$weight <- fall.con.ab$weight

undirected.fall.graph <- as.undirected(fall.graph, mode = "collapse",
                                       edge.attr.comb = "sum")

fall.communities.lab <- cluster_walktrap(undirected.fall.graph, weights = E(undirected.fall.graph)$weight,
                                         steps = 10, modularity = T)

# plot communities
plot(fall.communities.lab, fall.graph)

# plot communities on map
fall.comm.pal <- rainbow(length(seq(1, max(fall.communities.lab$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.lab$membership], add = T)
  
### community detection using edge betweenness -----------
fall.communities.bet <- cluster_edge_betweenness(fall.graph, fall.con.ab$weight, directed  = F)

# plot communities
plot(fall.communities.bet, fall.graph)

# plot communities on map
fall.comm.pal <- rainbow(length(seq(1, max(fall.communities.bet $membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.bet$membership], add = T)

### community detection using infomap -----------
fall.communities.info <- cluster_infomap(graph = undirected.fall.graph,
                                                  modularity = T,
                                         nb.trials = 5)
# plot communities
plot(fall.communities.info, fall.graph)

# plot communities on map
fall.comm.pal <-  rainbow(length(seq(1, max(fall.communities.info$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.info$membership], add = T)

### community detection using multilevel optimization of modularity -----------
fall.communities.louvain <- cluster_louvain(graph = undirected.fall.graph, resolution = 1)

meta.fall.ab$community <- fall.communities.louvain$membership

# plot communities
plot(fall.communities.louvain, fall.graph)

# plot communities on map
fall.comm.pal <- rainbow(length(seq(1, max(fall.communities.louvain$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.louvain$membership], add = T)

### community detection using the leiden algorithm -----------
fall.communities.leiden <- cluster_leiden(graph = undirected.fall.graph, resolution = 0.005)

meta.fall.ab$community <- fall.communities.leiden$membership

# plot communities
plot(fall.communities.leiden, undirected.fall.graph)

# plot communities on map
fall.comm.pal <- rainbow(length(seq(1, max(fall.communities.leiden$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.leiden$membership], add = T)

### Community detection using the algorithm developed by lancichinetti and Fortunato (https://www.nature.com/articles/srep00336) ----
E(fall.graph)$weight <- fall.con.ab$weight
undirected.fall.graph <- as.undirected(fall.graph, mode = "collapse",
                                         edge.attr.comb = "sum")
# Run concensusCluster function 
cluster_output <- concensusCluster(graph = undirected.fall.graph, thresh = 0.5, algiter = 3000)
comms <- cluster_output$`community structure`

# plot concensus graph
fall.comm.pal <- rainbow(length(seq(1, max(comms$membership))))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -35), ylim = c(-10, 65), lwd = 0.5, col = "#F7F7F7")
plot(fall.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[comms$membership], 
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

### calculate bridging centrality using the networktools package -----
fall.graph.el <- cbind( get.edgelist(fall.graph), round( E(fall.graph)$weight, 3 ))
qgraph(fall.graph.el)

bridge.c <- bridge(fall.graph, nodes =as.character(V(fall.graph)), communities = comms$membership, useCommunities = "all", directed = T)

bridge.strengh.palette  <- hcl.colors(n = length(seq(0, max(bridge.c$`Bridge Strength`),0.01)), palette = "Reds 3", rev = T) 
names(bridge.strengh.palette ) <- seq(0, max(bridge.c$`Bridge Strength`),0.01)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = bridge.strengh.palette [as.character(round(bridge.c$`Bridge Strength`, digits = 2))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

legend_image <- as.raster(matrix(hcl.colors(n = length(seq(0, max(bridge.c$`Bridge Strength`), 0.01)), palette = "Reds 3", rev = F) , ncol=1))
plot(c(0,2),c(-0.01,1),type = 'n', axes = F,xlab = '', ylab = '')
title(main = 'Bridge strength', cex.main = 0.9, line  = 0.8)
text(x= 1.25, y = seq(0,1,l=5), labels = round(seq(0,max(bridge.c$`Bridge Strength`),l=5), digits = 1), cex = 1)
rasterImage(legend_image, 0.75, 0, 1,1)

### Assess cluster significance with robin ----

# # Create random graph
# graphRandom <- random(graph = undirected.fall.graph)
# 
# # Use concensusClusterMod helper function
# proc <- robinRobust(graph=undirected.fall.graph, graphRandom=graphRandom, measure="vi",
#                     method="other", FUN = concensusClusterMod, type="independent", weights = E(undirected.fall.graph)$weight)
# 
# plotRobin(graph=undirected.fall.graph, model1=proc$Mean, model2=proc$MeanRandom,
#           legend=c("real data", "null model"))

# Create random graph
graphRandom <- random(graph = x)

# Use concensusClusterMod helper function
proc <- robinRobust(graph=x, graphRandom=graphRandom, measure="vi",
                    method= "labelProp", type="independent", weights = E(x)$weight)

plotRobin(graph=x, model1=proc$Mean, model2=proc$MeanRandom,
          legend=c("real data", "null model"))

### Assess cluster significance with clustanalytics ----

# create multiple rewired versions of the network, apply the clustering method, then calculate scoring functions 

iter = 100

rand.data = NULL

for (i in seq(1,iter)){
  
  # Create randomized graph
  rewire.fall <- rewireCpp(g = undirected.fall.graph, weight_sel="max_weight", Q = 1)
  
  # cluster analysis
  rewire.cluster <- concensusCluster(graph = undirected.fall.graph, thresh = 0.5, algiter = 3000)
  rewire.comms <- rewire.cluster$`community structure`$membership
  
  # calculate scoring function
  rewire.score <- scoring_functions(rewire.fall, com = rewire.comms, weighted = T, type = "global")
  
  # Build data.frame with results 
  if (is.null(rand.data)){
    rand.data <- as.data.frame(rewire.score)
  }else{
    rand.data <-rbind(rand.data,  as.data.frame(rewire.score))
  }
  
  print(paste0("progress: ", as.character(i), " to ", as.character(iter)))
}

# test whether scoring functions differ between the observed and randomized networks
score.test <- as.data.frame(NULL)

# scoring functions for the observed graph
ob.score <- as.data.frame(scoring_functions(undirected.fall.graph, comms$membership, type = "global"))

for (i in seq(1, length(colnames(rand.data)))){
  
  if (i != 15){
  test.results <- t.test(rand.data[,i], mu = ob.score[,i])
  
  score.test[i,"score.function"] <- colnames(rand.data)[i]
  score.test[i,"observed.score"] <- ob.score[,i]
  score.test[i,"random.function.mean"] <- mean(rand.data[,i], na.rm = T)
  score.test[i, "random.function.se"] <- sd(rand.data[,i], na.rm = T)/sqrt(length(rand.data[,i][!is.na(rand.data[,i])]))
  score.test[i,"statistic"] <- test.results$statistic
  score.test[i,"parameter"] <- test.results$parameter
  score.test[i,"p-value"] <- test.results$p.value
  }else{
    
    score.test[i,"score.function"] <- colnames(rand.data)[i]
    score.test[i,"observed.score"] <- ob.score[,i]
    score.test[i,"random.function.mean"] <- mean (rand.data[,i])
    score.test[i, "random.function.se"] <- sd(rand.data[,i], na.rm = T)/sqrt(length(rand.data[,i][!is.na(rand.data[,i])]))
    score.test[i,"statistic"] <- NA
    score.test[i,"parameter"] <- NA
    score.test[i,"p-value"] <- NA
  }
}

boot_alg_list(g = undirected.fall.graph)

# export test results 
write.csv(score.test, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Analysis_outuputs/fall_network_community_scoring_functions_test.csv")

### Time spent in each node during the fall ---

# Load edge metadata
fall.stat.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stationary.data.csv")
fall.stat.data$StartTime <- anytime(fall.stat.data$StartTime, asUTC = T)
fall.stat.data$EndTime <- anytime(fall.stat.data$EndTime, asUTC = T)

fall.stat <- fall.stat %>% mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# calculate the number of days spent at each stopover divided by the number of bird that were present (ONLY FOR STOPS CONSIDERED AS STOPOVERS)
fall.node.time <- fall.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                              Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                              Lat.2.5., Lat.97.5., EndTime,
                                              sitenum, duration, period,study.site,
                                              Range_region, NB_count, period,
                                              site_type) %>%
  mutate(duration = ifelse(site_type == "Stopover", duration, 0)) %>%
  group_by(geo_id, cluster) %>% summarise(time.cluster.occupied = sum(duration)) %>%
  group_by(cluster) %>% summarise(mean.time.cluster.occupied = mean(time.cluster.occupied))
  #group_by(cluster) %>% summarise(sum.time.cluster.occupied = sum(time.cluster.occupied))
  
# add the times calculated to the node metadata
meta.fall.ab$time.occupied <- fall.node.time$mean.time.cluster.occupied

# plot of the average time spent in each node 
time.palette <- hcl.colors(n = length(seq(0, max(meta.fall.ab$time.occupied +1))), palette = "Blues3", rev = T) 
names(time.palette) <- seq(0, max(meta.fall.ab$time.occupied +1))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(fall.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = time.palette[as.character(round(meta.fall.ab$time.occupied))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = NA)

### estimated proportion of blackpoll warblers in each node ----
r.ab.palette <- hcl.colors(n = length(seq(0, max(meta.fall.ab$r.abundance.at.cluster + 0.01), 0.001)), palette = "Reds 3", rev = T) 
names(r.ab.palette) <- seq(0, max(meta.fall.ab$r.abundance.at.cluster + 0.01), 0.001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(fall.graph, vertex.size = 300, vertex.size2 = 200, vertex.label= NA,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = r.ab.palette[as.character(round(meta.fall.ab$r.abundance.at.cluster, digits = 3))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

### Number of breeding populations in each node ----
fall.n.br <- fall.stat %>% group_by(cluster) %>% summarize(n.br.pops = n_distinct(study.site)) 
meta.fall.ab <- merge(meta.fall.ab, fall.n.br[, c("cluster", "n.br.pops")], by.x = "vertex", by.y = "cluster")

r.n.br.palette <- hcl.colors(n = length(seq(1, max(meta.fall.ab$n.br.pops), 1)), palette = "Reds 3", rev = T) 
names(r.n.br.palette) <- seq(1, max(meta.fall.ab$n.br.pops), 1)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = r.n.br.palette [as.character(meta.fall.ab$n.br.pops)],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

# Prepare for spring network analysis ############################################

# Load the spring graph
spring.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.graph.edge.list.txt", directed = TRUE)

# Load spring graph node metadata 
meta.spring.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.node.metadata.csv")

# Load spring graph edge weights
spring.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.edge.weights.csv")

# Plot the spring graph 
data("wrld_simpl")
type.palette <- c("#440154FF", "#FDE725FF", "#21908CFF")

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
    xlim = c(-165, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)

# Edge colours for spring and fall edges (fall edges should not appear in this plot)
edge.cols <- spring.con.ab %>% mutate(col = case_when(
  edge.type == "spring" ~ adjustcolor("darkgray", alpha.f = 0.9),
  edge.type == "fall" ~ NA))

plot(spring.graph, vertex.label = NA, vertex.size = 300, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*20, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = type.palette[meta.spring.ab$node.type.num],
     edge.color= edge.cols$col, add = T)
# legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
#        col = type.palette[unique(meta.spring.ab$node.type.num)],
#        pch = 19, )

# Plot the spring graph with node use pie charts #################################

# Create a palette for site use by range region  
use.palette <- list(c("#440154FF", "#FDE725FF", "#21908CFF"))

#create a column that can be converted to a numeric vector
meta.spring.ab <- transform(meta.spring.ab, use.ab.vector = asplit(cbind(use.breeding.ab,  use.stopover.ab, use.nonbreeding.ab), 1))

# Create vector of vertex shapes 
meta.spring.ab <- meta.spring.ab %>% rowwise() %>%
  mutate(shape_single.use = length(which(c(use.breeding.ab, 
                                           use.nonbreeding.ab, 
                                           use.stopover.ab)==0))) %>%
  ungroup() %>%
  mutate(shape_single.use = ifelse(shape_single.use == 2, "circle", "none")) %>%
  mutate(shape_multiple.use = ifelse(shape_single.use== "none", "pie", "none")) %>%
  mutate(shape_colour_multiple.use = case_when(shape_single.use != "none" & use.breeding.ab != 0 ~ type.palette[1] ,
                                               shape_single != "none" & use.stopover.ab != 0 ~ type.palette[2] ,
                                               shape_single != "none" & use.nonbreeding.ab != 0 ~ type.palette[3] ,
                                               .default = NA)) %>%
  mutate(shape_colour_single.use = case_when(shape_single.use!= "circle" & use.breeding.ab != 0 ~ type.palette[1],
                                             shape_single.use != "circle" & use.stopover.ab != 0 ~ type.palette[2],
                                             shape_single.use != "circle" & use.nonbreeding.ab!= 0 ~ type.palette[3],
                                             .default = NA))


plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)

plot(as_Spatial(equinox_region), add = T, col = "#D9D5B2", lwd = 0.000000001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -35), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)

plot(spring.graph, vertex.size = 300, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*40, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = type.palette[meta.spring.ab$node.type.num], vertex.label.dist = 30,
     add = T, vertex.label = NA, edge.color = adjustcolor("darkgray", alpha.f = 0.6),
     vertex.shape = meta.spring.ab$shape_single.use)

plot(spring.graph, vertex.size = 300, vertex.size2 = 200,
     edge.width = 0, edge.arrow.size = 0, edge.arrow.width = 0,
     vertex.shape = meta.spring.ab$shape_multiple.use, vertex.pie = meta.spring.ab$ use.ab.vector,
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.pie.color = use.palette, vertex.color = NA,
     add = T, vertex.label = NA, edge.color = adjustcolor("darkgray", alpha.f = 0.6))

legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.spring.ab$node.type.num)],
       pch = 19, title = "Node type",cex = 0.8)

##Spring network analysis ###########################################################

## Spring betweenness centrality (exclude breeding sites)
spring.betw.c <- betweenness(spring.graph, directed = T, weights = 1/spring.con.ab$weight)

# plot of the betweenness centrality of each node 
betw.c.palette <- hcl.colors(n = length(seq(0, max(spring.betw.c)+1)), palette = "Reds 3", rev = T) 
names(betw.c.palette) <- seq(0, max(spring.betw.c)+1)

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = betw.c.palette[as.character(round(spring.betw.c))], add = T, vertex.label= NA)

## spring degree centrality
spring.degree.c <- degree(spring.graph, mode = "in")

# plot of the degree centrality of each node 
spring.deg.c.palette <- hcl.colors(n = length(seq(0, max(spring.degree.c))), palette = "Reds 3", rev = T) 
names(spring.deg.c.palette) <- seq(0, max(spring.degree.c))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.deg.c.palette[as.character(spring.degree.c)],
     edge.color= adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = NA)

## spring weighed degree centrality ----
spring.strength.c <- strength(spring.graph, mode = "all", weight = spring.con.ab$weight)

# plot of the weighed centrality of each node 
deg.c.palette <- hcl.colors(n = length(seq(0, max(spring.strength.c + 0.01), 0.01)), palette = "Reds 3", rev = T) 
names(deg.c.palette) <- seq(0, max(spring.strength.c+ 0.01), 0.01)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = deg.c.palette[as.character(round(spring.strength.c, digits = 2))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = NA)

# spring degree centrality using the method developed by Opshal et al. ---
E(spring.graph)$weight <- spring.con.ab$weight
spring.edge.list <- cbind(get.edgelist(spring.graph), round(E(spring.graph)$weight, 3))
spring.net <- as.tnet(spring.edge.list, type = "weighted one-mode tnet")

spring.degree.TO <- degree_w(spring.edge.list, measure=c("degree", "output", "alpha"), type = "in", alpha = 0.5)

deg.TO.palette <- hcl.colors(n = length(seq(0, max(spring.degree.TO[,"alpha"]+ 0.0001),0.0001)), palette = "Reds 3", rev = T) 
names(deg.TO.palette) <- seq(0, max(spring.degree.TO[,"alpha"] + 0.0001),0.0001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = deg.TO.palette[as.character(round(spring.degree.TO[,"alpha"],4))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = round(spring.strength.c, digits = 4))

# spring degree centrality using the method developed by Opshal et al. ---
E(spring.graph)$weight <- spring.con.ab$weight
spring.edge.list <- cbind(get.edgelist(spring.graph) , round( E(spring.graph)$weight, 3 ))
spring.net <- as.tnet(spring.edge.list, type = "weighted one-mode tnet")

spring.degree.TO <- degree_w(spring.edge.list, measure=c("degree", "output", "alpha"), type = "out", alpha = 0.5)

deg.TO.palette <- hcl.colors(n = length(seq(0, max(spring.degree.TO[,"alpha"]+ 0.0001),0.0001)), palette = "Reds 3", rev = T) 
names(deg.TO.palette) <- seq(0, max(spring.degree.TO[,"alpha"] + 0.0001),0.0001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
#main = "Degree centrality")

plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = deg.TO.palette[as.character(round(spring.degree.TO[,"alpha"],4))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = round(spring.strength.c, digits = 4))

legend_image <- as.raster(matrix(hcl.colors(n = length(seq(0, max(spring.degree.TO[,"alpha"]), max(spring.degree.TO[,"alpha"])/20)), palette = "Reds 3", rev = F) , ncol=1))
plot(c(0,2),c(-0.01,1),type = 'n', axes = F,xlab = '', ylab = '')
title(main = 'Out-degree centrality', cex.main = 0.9, line  = 0.8)
text(x= 1.25, y = seq(0,1,l=5), labels = round(seq(0,max(spring.degree.TO[,"alpha"]),l=5), digits = 1), cex = 1)
rasterImage(legend_image, 0.75, 0, 1,1)

# spring betweenness centrality using the method developed by Opshal et al. ---
spring.betweenness.TO <- betweenness_w(spring.edge.list, directed = T, alpha = 0.5)

betweenness.TO.palette <- hcl.colors(n = length(seq(0, max(spring.betweenness.TO[,"betweenness"]),1)), palette = "Reds 3", rev = T) 
names(betweenness.TO.palette) <- seq(0, max(spring.betweenness.TO[,"betweenness"]),1)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
#main = "Betweenness centrality")

plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = betweenness.TO.palette[as.character(spring.betweenness.TO[,"betweenness"])],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = round(spring.strength.c, digits = 4))

legend_image <- as.raster(matrix(hcl.colors(n = length(seq(0, max(spring.betweenness.TO[,"betweenness"]), 1)), palette = "Reds 3", rev = F) , ncol=1))
plot(c(0,2),c(-0.01,1),type = 'n', axes = F,xlab = '', ylab = '')
title(main = 'Betweenness centrality', cex.main = 0.9, line  = 0.8)
text(x= 1.25, y = seq(0,1,l=5), labels = seq(0,max(spring.betweenness.TO[,"betweenness"]),l=5), cex = 1)
rasterImage(legend_image, 0.75, 0, 1,1)

# community detection using propagating labels ----------
E(spring.graph)$weight <- spring.con.ab$weight

undirected.spring.graph <- as.undirected(spring.graph, mode = "collapse",
                                       edge.attr.comb = "sum")

spring.communities.lab <- cluster_label_prop(undirected.spring.graph, weights = E(undirected.spring.graph)$weight)

#spring.communities.lab <- cluster_label_prop(spring.graph, weights = spring.con.ab$weight)
meta.spring.ab$community <- spring.communities.lab$membership

# plot communities
plot(spring.communities.lab, spring.graph)

# plot communities on map
spring.comm.pal <- rainbow(length(seq(1, max(spring.communities.lab$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2")
plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal[spring.communities.lab$membership], add = T,
     vertex.label = NA)

# community detection using edge betweenness -----------
spring.communities.bet <- cluster_edge_betweenness(spring.graph, weights = spring.con.ab$weight, directed = F)
meta.spring.ab$community <- spring.communities.bet$membership

# plot communities
plot(spring.communities.bet, spring.graph)

# plot communities on map
spring.comm.pal <- rainbow(length(seq(1, max(spring.communities.bet$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(spring.graph.weighed.ab, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal[spring.communities.bet$membership], add = T,
     vertex.label = NA)

### Community detection using the algorithm developed by lancichinetti and Fortunato (https://www.nature.com/articles/srep00336) ----
E(spring.graph)$weight <- spring.con.ab$weight
undirected.spring.graph <- as.undirected(spring.graph, mode = "collapse",
                                       edge.attr.comb = "sum")
# Run concensusCluster function 
cluster_output <- concensusCluster(graph = undirected.spring.graph, thresh = 0.6, algiter = 3000)
comms <- cluster_output$`community structure`

# plot concensus graph
spring.comm.pal <- rainbow(length(seq(1, max(comms$membership))))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -35), ylim = c(-10, 65), lwd = 0.5, col = "#F7F7F7")
plot(spring.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal[comms$membership], 
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

### calculate bridging centrality using the networktools package -----
spring.graph.el <- cbind( get.edgelist(spring.graph), round( E(spring.graph)$weight, 3 ))

bridge.c <- bridge(spring.graph, nodes =as.character(V(spring.graph)), communities = comms$membership, useCommunities = "all", directed = T)

bridge.strengh.palette  <- hcl.colors(n = length(seq(0, max(bridge.c$`Bridge Strength`),0.01)), palette = "Reds 3", rev = T) 
names(bridge.strengh.palette ) <- seq(0, max(bridge.c$`Bridge Strength`),0.01)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = bridge.strengh.palette [as.character(round(bridge.c$`Bridge Strength`, digits = 2))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

legend_image <- as.raster(matrix(hcl.colors(n = length(seq(0, max(bridge.c$`Bridge Strength`), 0.01)), palette = "Reds 3", rev = F) , ncol=1))
plot(c(0,2),c(-0.01,1),type = 'n', axes = F,xlab = '', ylab = '')
title(main = 'Bridge strength', cex.main = 0.9, line  = 0.8)
text(x= 1.25, y = seq(0,1,l=5), labels = round(seq(0,max(bridge.c$`Bridge Strength`),l=5), digits = 1), cex = 1)
rasterImage(legend_image, 0.75, 0, 1,1)

### calculate bridge betweennness using the networktools package -----
bridge.betw.palette  <- hcl.colors(n = length(seq(0, max(bridge.c$`Bridge Betweenness`),0.01)), palette = "Reds 3", rev = T) 
names(bridge.betw.palette ) <- seq(0, max(bridge.c$`Bridge Betweenness`),0.01)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = bridge.betw.palette[as.character(round(bridge.c$`Bridge Betweenness`, digits = 2))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

legend_image <- as.raster(matrix(hcl.colors(n = length(seq(0, max(bridge.c$`Bridge Betweenness`), 0.01)), palette = "Reds 3", rev = F) , ncol=1))
plot(c(0,2),c(-0.01,1),type = 'n', axes = F,xlab = '', ylab = '')
title(main = 'Bridge betweenness', cex.main = 0.9, line  = 0.8)
text(x= 1.25, y = seq(0,1,l=5), labels = round(seq(0,max(bridge.c$`Bridge Betweenness`),l=5), digits = 1), cex = 1)
rasterImage(legend_image, 0.75, 0, 1,1)

### Assess cluster significance with clustanalytics ----

# create multiple rewired versions of the network, apply the clustering method, then calculate scoring functions 

iter = 100

rand.data = NULL

for (i in seq(1,iter)){
  
  # Create randomized graph
  rewire.spring <- rewireCpp(g = undirected.spring.graph, weight_sel="max_weight", Q = 1)
  
  # cluster analysis
  rewire.cluster <- concensusCluster(graph = undirected.spring.graph, thresh = 0.5, algiter = 3000)
  rewire.comms <- rewire.cluster$`community structure`$membership
  
  # calculate scoring function
  rewire.score <- scoring_functions(rewire.spring, com = rewire.comms, weighted = T, type = "global")
  
  # Build data.frame with results 
  if (is.null(rand.data)){
    rand.data <- as.data.frame(rewire.score)
  }else{
    rand.data <-rbind(rand.data,  as.data.frame(rewire.score))
  }
  
  print(paste0("progress: ", as.character(i), " to ", as.character(iter)))
}

# test whether scoring functions differ between the observed and randomized networks
score.test <- as.data.frame(NULL)

# scoring functions for the observed graph
ob.score <- as.data.frame(scoring_functions(undirected.spring.graph, comms$membership, type = "global"))

for (i in seq(1, length(colnames(rand.data)))){
  
  if (!(i %in% c(13, 15))){
    test.results <- t.test(rand.data[,i], mu = ob.score[,i])
    
    score.test[i,"score.function"] <- colnames(rand.data)[i]
    score.test[i,"observed.score"] <- ob.score[,i]
    score.test[i,"random.function.mean"] <- mean (rand.data[,i], na.rm = T)
    score.test[i,"statistic"] <- test.results$statistic
    score.test[i,"parameter"] <- test.results$parameter
    score.test[i,"p-value"] <- test.results$p.value
  }else{
    
    score.test[i,"score.function"] <- colnames(rand.data)[i]
    score.test[i,"observed.score"] <- ob.score[,i]
    score.test[i,"random.function.mean"] <- mean (rand.data[,i])
    score.test[i,"statistic"] <- NA
    score.test[i,"parameter"] <- NA
    score.test[i,"p-value"] <- NA
  }
}

# export test results 
write.csv(score.test, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Analysis_outuputs/spring_network_community_scoring_functions_test.csv")

### Time spent in each node during the spring ---

# Load edge metadata
spring.stat.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/spring.stationary.data.csv")
spring.stat.data$StartTime <- anytime(spring.stat.data$StartTime, asUTC = T)
spring.stat.data$EndTime <- anytime(spring.stat.data$EndTime, asUTC = T)

spring.stat <- spring.stat %>% mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# calculate the number of days spent at each stopover divided by the number of bird that were present (ONLY FOR STOPS CONSIDERED AS STOPOVERS)
spring.node.time <- spring.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                              Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                              Lat.2.5., Lat.97.5., EndTime,
                                              sitenum, duration, period,study.site,
                                              Range_region, NB_count, period,
                                              site_type) %>%
  mutate(duration = ifelse(site_type == "Stopover", duration, 0)) %>%
  group_by(geo_id, cluster) %>% summarise(time.cluster.occupied = sum(duration)) %>%
  group_by(cluster) %>% summarise(mean.time.cluster.occupied = mean(time.cluster.occupied))
  #group_by(cluster) %>% summarise(sum.time.cluster.occupied = sum(time.cluster.occupied))

# add the times calculated to the node metadata
meta.spring.ab$time.occupied <- spring.node.time$mean.time.cluster.occupied

# plot of the average time spent in each node 
time.palette <- hcl.colors(n = length(seq(0, max(meta.spring.ab$time.occupied +1))), palette = "Blues3", rev = T) 
names(time.palette) <- seq(0, max(meta.spring.ab$time.occupied +1))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(spring.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = time.palette[as.character(round(meta.spring.ab$time.occupied))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = NA) # round(meta.spring.ab$time.occupied))

### estimated proportion of blackpoll warblers in each node ----
r.ab.palette <- hcl.colors(n = length(seq(0, max(meta.spring.ab$r.abundance.at.cluster + 0.01), 0.001)), palette = "Reds 3", rev = T) 
names(r.ab.palette) <- seq(0, max(meta.spring.ab$r.abundance.at.cluster + 0.01), 0.001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = r.ab.palette[as.character(round(meta.spring.ab$r.abundance.at.cluster, digits = 3))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

### Number of breeding populations in each node ----
spring.n.br <- spring.stat %>% group_by(cluster) %>% summarize(n.br.pops = n_distinct(study.site)) 
meta.spring.ab <- merge(meta.spring.ab, spring.n.br[, c("cluster", "n.br.pops")], by.x = "vertex", by.y = "cluster")

r.n.br.palette <- hcl.colors(n = length(seq(1, max(meta.spring.ab$n.br.pops), 1)), palette = "Reds 3", rev = T) 
names(r.n.br.palette) <- seq(1, max(meta.spring.ab$n.br.pops), 1)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label= NA,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = r.n.br.palette [as.character(meta.spring.ab$n.br.pops)],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

# Code snippets ################################################################

# code to display a color vector COL by stack exchange user G5W (https://stackoverflow.com/questions/48641985/plot-a-vector-of-colors-in-r)
# plot(NULL, xlim=c(0,length(COL)), ylim=c(0,1), 
#      xlab="", ylab="", xaxt="n", yaxt="n")
# rect(0:(length(COL)-1), 0, 1:length(COL), 1, col=COL)

# mods = c()
# 
# for (i in seq(1,100)){
# 
# x <- rewire(undirected.fall.graph, with = keeping_degseq(niter = 100))
# E(x)$weight <- sample(E(undirected.fall.graph)$weight, replace = FALSE)
# 
# y <- rewireCpp(g = undirected.fall.graph, Q = 1, weight_sel = "max_weight")
# E(y)
# 
# x.cluster <- concensusCluster(graph = x, thresh = 0.5, algiter = 3000)
# y.cluster <- concensusCluster(graph = y, thresh = 0.5, algiter = 3000)
# undirected.fall.cluster <- concensusCluster(graph = undirected.fall.graph, thresh = 0.5, algiter = 3000)
# 
# plot(x, vertex.color = x.cluster$`community structure`$membership, layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), edge.width = E(x)$weight*15)
# plot(y, vertex.color = y.cluster$`community structure`$membership, layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), edge.width = E(y)$weight*15)
# plot(undirected.fall.graph, vertex.color = undirected.fall.cluster $`community structure`$membership, layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]))
# 
# strength(x)
# strength(y)
# strength(undirected.fall.graph)
# 
# mods[i] <- modularity(x, x.cluster$`community structure`$membership)
# scoring_functions(x, x.cluster$`community structure`$membership, type = "global")
# scoring_functions(y, y.cluster$`community structure`$membership, type = "global")
# }
# 
# hist(mods, breaks = 10)
# abline(v = ob.score$modularity, col = "blue")
# abline (v = mean(mods), col = "red")
# t.test(x = mods, mu = ob.score$modularity)
