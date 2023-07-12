# Network analysis for the blackpoll warbler migration network
library(igraph)
library(anytime)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(RColorBrewer)

# Prepare for fall network analysis ############################################

# Load the fall graph
fall.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.graph.edge.list.txt", directed = TRUE)

# Load fall graph node metadata 
meta.fall <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.csv")

# Load fall graph edge weights
fall.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.edge.weights.csv")

# Plot the fall graph 
data("wrld_simpl")
type.palette <- rainbow(3)  

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*10, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = type.palette[meta.fall$node.type.num], vertex.label.dist = 30,
     add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.fall$node.type.num)],
       pch = 16)

## Fall network analysis ###########################################################

## Fall betweenness centrality (exclude breeding sites)
fall.betw.c <- betweenness(fall.graph, directed = F, weights = NA)

# plot of the betweenness centrality of each node 
betw.c.palette <- hcl.colors(n = length(seq(1, max(fall.betw.c))), palette = "Reds 3", rev = T) 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = betw.c.palette[round(fall.betw.c)], add = T)

## Fall degree centrality
fall.degree.c <- degree(fall.graph, mode = "all")

# plot of the betweenness centrality of each node 
deg.c.palette <- hcl.colors(n = length(seq(1, max(fall.degree.c ))), palette = "Reds 3", rev = T) 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.size = 200, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = deg.c.palette[fall.degree.c], add = T)

### community detection using propogating labels ----------
fall.communities.lab <- cluster_label_prop(fall.graph, weights = fall.con.ab$weight)
meta.fall$community <- fall.communities.lab$membership

# plot communities
plot(fall.communities.lab, fall.graph)

# plot communities on map
fall.comm.pal <- rainbow(n_distinct(fall.communities.lab$membership))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.lab$membership], add = T)

### community detection using edge betweenness -----------
fall.communities.bet <- cluster_edge_betweenness(fall.graph, fall.con.ab$weight, directed  = T)
meta.fall$community <- fall.communities.bet$membership

# plot communities
plot(fall.communities.bet, fall.graph)

# plot communities on map
fall.comm.pal <- brewer.pal(n = n_distinct(fall.communities.bet $membership), name = "Paired") 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.comm.pal[fall.communities.bet $membership], add = T)

### community detection using infomap -----------
fall.communities.info <- cluster_infomap(graph = as.undirected(fall.graph),
                                                  modularity = T,
                                         nb.trials = 100)

meta.fall$community <- fall.communities.info$membership

# plot communities
plot(fall.communities.info, fall.graph)

# plot communities on map
fall.comm.pal <- brewer.pal(n = n_distinct(fall.communities$membership), name = "Paired") 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.communities.info$membership, add = T)


### community detection using multilevel optimization of modularity -----------
fall.communities.louvain <- cluster_louvain(graph = as.undirected(fall.graph))

meta.fall$community <- fall.communities.info$membership

# plot communities
plot(fall.communities.info, fall.graph)

# plot communities on map
fall.comm.pal <- brewer.pal(n = n_distinct(fall.communities.louvain$membership), name = "Paired") 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = fall.communities.louvain$membership, add = T)

# Prepare for spring network analysis ############################################

# Load the spring graph
spring.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.graph.edge.list.txt", directed = TRUE)

# Load spring graph node metadata 
meta.spring <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.node.metadata.csv")

# Load spring graph edge weights
spring.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.edge.weights.csv")

# Plot the springgraph 
data("wrld_simpl")
type.palette <- rainbow(3)  

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(spring.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = type.palette[meta.spring$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.spring$node.type.num)],
       pch = 16)

##Spring network analysis ###########################################################

# community detection using propogating labels ----------
spring.communities.lab <- cluster_label_prop(spring.graph, weights = spring.con.ab$weight)
meta.spring$community <- spring.communities.lab$membership

# plot communities
plot(spring.communities.lab, spring.graph)

# plot communities on map
spring.comm.pal <- brewer.pal(n = n_distinct(spring.communities.lab$membership), name = "Paired") 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(spring.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal[ spring.communities.lab$membership], add = T)

# community detection using edge betweenness -----------
spring.communities.bet <- cluster_edge_betweenness(spring.graph, spring.con.ab$weight, directed  = T)
meta.spring$community <- spring.communities.bet$membership

# plot communities
plot(spring.communities.bet, spring.graph)

# plot communities on map
spring.comm.pal <- brewer.pal(n = n_distinct(spring.communities.bet$membership), name = "Paired") 

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(spring.graph.weighed.ab, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal [spring.communities.bet$membership], add = T)
