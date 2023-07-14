# Network analysis for the blackpoll warbler migration network
library(igraph)
library(anytime)
library(dplyr)
library(ggplot2)
library(sf)
library(rgdal)
library(terra)
library(RColorBrewer)

# Prepare for fall network analysis ############################################

# Load the fall graph
fall.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.graph.edge.list.txt", directed = TRUE)

# Load fall graph node metadata 
meta.fall <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.csv")

# Load fall graph edge weights
fall.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.edge.weights.csv")

# For fall nodes where latitudinal accuracy is low, set location close to the coast
meta.fall[c(10, 24, 18),]$Lat.50. <- c(34, 42, 44.4)

# Import polygon of the region affected by the equinox
equinox_region <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Mapping_components/Data/Fall_Equinox_affected_regionV3.shp")

# Plot the fall graph 
data("wrld_simpl")
type.palette <- rainbow(3)  

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(as_Spatial(equinox_region), add = T, col = "#FCF2CF", lwd = 0.5)
plot(fall.graph, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*40, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = type.palette[meta.fall$node.type.num], vertex.label.dist = 30,
     add = T, vertex.label = NA)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.fall$node.type.num)],
       pch = 16)

## Fall network analysis ###########################################################

## Fall betweenness centrality (exclude breeding sites)
fall.betw.c <- betweenness(fall.graph, directed = T, weights = 1/fall.con.ab$weight)

# plot of the betweenness centrality of each node 
betw.c.palette <- hcl.colors(n = length(seq(0, max(fall.betw.c))), palette = "Reds 3", rev = T) 
names(betw.c.palette) <- seq(0, max(fall.betw.c))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = betw.c.palette[as.character(round(fall.betw.c))], add = T, vertex.label = NA)

## Fall degree centrality
fall.degree.c <- degree(fall.graph, mode = "all")

# plot of the degree centrality of each node 
deg.c.palette <- hcl.colors(n = length(seq(1, max(fall.degree.c ))), palette = "Reds 3", rev = T) 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = deg.c.palette[fall.degree.c], add = T, vertex.label = NA)

### community detection using propogating labels ----------
fall.communities.lab <- cluster_label_prop(fall.graph, weights = fall.con.ab$weight)

# plot communities
plot(fall.communities.lab, fall.graph)

# plot communities on map
fall.comm.pal <- rainbow(length(seq(1, max(fall.communities.lab$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.lab$membership], add = T)

### community detection using edge betweenness -----------
fall.communities.bet <- cluster_edge_betweenness(fall.graph, fall.con.ab$weight, directed  = T)

# plot communities
plot(fall.communities.bet, fall.graph)

# plot communities on map
fall.comm.pal <- rainbow(length(seq(1, max(fall.communities.bet $membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.bet$membership], add = T)

### community detection using infomap -----------
fall.communities.info <- cluster_infomap(graph = fall.graph,
                                                  modularity = T,
                                         nb.trials = 100)
# plot communities
plot(fall.communities.info, fall.graph)

# plot communities on map
fall.comm.pal <-  rainbow(length(seq(1, max(fall.communities.info$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.info$membership], add = T)

### community detection using multilevel optimization of modularity -----------
fall.communities.louvain <- cluster_louvain(graph = as.undirected(fall.graph))

meta.fall$community <- fall.communities.louvain$membership

# plot communities
plot(fall.communities.louvain, fall.graph)

# plot communities on map
fall.comm.pal <- rainbow(length(seq(1, max(fall.communities.louvain$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.louvain$membership], add = T)

### Time spent in each node during the fall ---

# Load edge metadata
fall.stat.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stationary.data.csv")
fall.stat.data$StartTime <- anytime(fall.stat.data$StartTime, asUTC = T)
fall.stat.data$EndTime <- anytime(fall.stat.data$EndTime, asUTC = T)

fall.stat <- fall.stat %>% mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# calculate the average time that individuals spent stationary in each node (ONLY FOR STOPS CONSIDERED AS STOPOVERS)
fall.node.time <- fall.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                              Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                              Lat.2.5., Lat.97.5., EndTime,
                                              sitenum, duration, period,study.site,
                                              Range_region, NB_count, period,
                                              site_type) %>%
  mutate(duration = ifelse(site_type == "Stopover", duration, 0)) %>%
  group_by(geo_id, cluster) %>% summarise(time.cluster.occupied = sum(duration)) %>%
  group_by(cluster) %>% summarise(mean.time.cluster.occupied = mean(time.cluster.occupied))

# add the times calculated to the node metadata
meta.fall$time.occupied <- fall.node.time$mean.time.cluster.occupied

# plot of the average time spent in each node 
time.palette <- hcl.colors(n = length(seq(0, max(meta.fall$time.occupied +1))), palette = "Blues3", rev = T) 
names(time.palette) <- seq(0, max(meta.fall$time.occupied +1))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = time.palette[as.character(round(meta.fall$time.occupied))], add = T, vertex.label = NA)
  
# Prepare for spring network analysis ############################################

# Load the spring graph
spring.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.graph.edge.list.txt", directed = TRUE)

# Load spring graph node metadata 
meta.spring <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.node.metadata.csv")

# Load spring graph edge weights
spring.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.edge.weights.csv")

# Plot the spring graph 
data("wrld_simpl")
type.palette <- rainbow(3)  

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(spring.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*20, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = type.palette[meta.spring$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.spring$node.type.num)],
       pch = 16)

##Spring network analysis ###########################################################

## Spring betweenness centrality (exclude breeding sites)
spring.betw.c <- betweenness(spring.graph, directed = T, weights = 1/spring.con.ab$weight)

# plot of the betweenness centrality of each node 
betw.c.palette <- hcl.colors(n = length(seq(0, max(spring.betw.c))), palette = "Reds 3", rev = T) 
names(betw.c.palette) <- seq(0, max(spring.betw.c))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = betw.c.palette[as.character(round(spring.betw.c))], add = T, vertex.label = NA)

## spring degree centrality
spring.degree.c <- degree(spring.graph, mode = "all")

# plot of the degree centrality of each node 
spring.deg.c.palette <- hcl.colors(n = length(seq(1, max(spring.degree.c ))), palette = "Reds 3", rev = T) 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.deg.c.palette[spring.degree.c], add = T, vertex.label = NA)

# community detection using propagating labels ----------
spring.communities.lab <- cluster_label_prop(spring.graph, weights = spring.con.ab$weight)
meta.spring$community <- spring.communities.lab$membership

# plot communities
plot(spring.communities.lab, spring.graph)

# plot communities on map
spring.comm.pal <- rainbow(length(seq(1, max(spring.communities.lab$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2")
plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal[spring.communities.lab$membership], add = T,
     vertex.label = NA)

# community detection using edge betweenness -----------
spring.communities.bet <- cluster_edge_betweenness(spring.graph, weights = spring.con.ab$weight, directed = T)
meta.spring$community <- spring.communities.bet$membership

# plot communities
plot(spring.communities.bet, spring.graph)

# plot communities on map
spring.comm.pal <- rainbow(length(seq(1, max(spring.communities.bet$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(spring.graph.weighed.ab, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal[spring.communities.bet$membership], add = T,
     vertex.label = NA)

### Time spent in each node during the spring ---

# Load edge metadata
spring.stat.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/spring.stationary.data.csv")
spring.stat.data$StartTime <- anytime(spring.stat.data$StartTime, asUTC = T)
spring.stat.data$EndTime <- anytime(spring.stat.data$EndTime, asUTC = T)

spring.stat <- spring.stat %>% mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# calculate the average time that individuals spent stationary in each node (ONLY FOR STOPS CONSIDERED AS STOPOVERS)
spring.node.time <- spring.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                              Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                              Lat.2.5., Lat.97.5., EndTime,
                                              sitenum, duration, period,study.site,
                                              Range_region, NB_count, period,
                                              site_type) %>%
  mutate(duration = ifelse(site_type == "Stopover", duration, 0)) %>%
  group_by(geo_id, cluster) %>% summarise(time.cluster.occupied = sum(duration)) %>%
  group_by(cluster) %>% summarise(mean.time.cluster.occupied = mean(time.cluster.occupied))

# add the times calculated to the node metadata
meta.spring$time.occupied <- spring.node.time$mean.time.cluster.occupied

# plot of the average time spent in each node 
time.palette <- hcl.colors(n = length(seq(0, max(meta.spring$time.occupied +1))), palette = "Blues3", rev = T) 
names(time.palette) <- seq(0, max(meta.spring$time.occupied +1))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(spring.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = time.palette[as.character(round(meta.spring$time.occupied))], add = T, vertex.label = NA)

# Code snippets ################################################################

# code to display a color vector COL by stack exchange user G5W (https://stackoverflow.com/questions/48641985/plot-a-vector-of-colors-in-r)
# plot(NULL, xlim=c(0,length(COL)), ylim=c(0,1), 
#      xlab="", ylab="", xaxt="n", yaxt="n")
# rect(0:(length(COL)-1), 0, 1:length(COL), 1, col=COL)
