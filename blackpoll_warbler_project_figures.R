# Thesis figures 

# Load libraries ----
library(tidyverse)
library(igraph)
library(lubridate)
library(anytime)
library(ggplot2)
library(geosphere)
library(terra)
library(sf)
library(maptools)
library(ebirdst)
library(scatterpie)
library(ggnewscale)
library(cowplot)

# network specific libraries ----
library(igraph)
library(ggnetwork)
library(intergraph)
library(networktools)

# Will need to run the network analysis and construction scripts ----
#source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/Network_Analysis.R")

# Load required data for the fall
fall.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.graph.edge.list.txt", directed = TRUE)
meta.fall.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.csv")
fall.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.edge.weights.csv")
equinox_region <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Mapping_components/Data/Fall_Equinox_affected_regionV5.shp")
bpw_range <-  read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Birdlife_international_species_distribution/SppDataRequest.shp")

# Load required data for the spring 
spring.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.graph.edge.list.txt", directed = TRUE)
meta.spring.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.node.metadata.csv")
spring.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.edge.weights.csv")

# add weights to fall and spring graph  --- 
E(fall.graph)$weight <- fall.con.ab$weight
E(spring.graph)$weight <- spring.con.ab$weight

# Other types of edge data ----
edge.cols.fall <- fall.con.ab %>% mutate(col = case_when(
  edge.type == "fall" ~ adjustcolor("darkgray", alpha.f = 0.9),
  edge.type == "spring" ~ NA)) # edge types for the fall

edge.cols.spring <- spring.con.ab %>% mutate(col = case_when(
  edge.type == "spring" ~ adjustcolor("darkgray", alpha.f = 0.9),
  edge.type == "fall" ~ NA)) # edge types for the spring 

# Edge types
E(fall.graph)$edge.type <- edge.cols.fall$edge.type
E(spring.graph)$edge.type <- edge.cols.spring$edge.type

# Node degree weights 
V(fall.graph)$node.weight <- meta.fall.ab$r.abundance.at.cluster
V(spring.graph)$node.weight <- meta.spring.ab$r.abundance.at.cluster

# Node type (breeding, stopover or nonbreeding)
V(fall.graph)$node.type <- meta.fall.ab$node.type
V(spring.graph)$node.type <- meta.spring.ab$node.type

# Node population composition 
V(fall.graph)$West <- meta.fall.ab$prop.ab.western
V(fall.graph)$Central <- meta.fall.ab$prop.ab.central
V(fall.graph)$East <- meta.fall.ab$prop.ab.eastern

V(spring.graph)$West <- meta.spring.ab$prop.ab.western
V(spring.graph)$Central <- meta.spring.ab$prop.ab.central
V(spring.graph)$East <- meta.spring.ab$prop.ab.eastern

# Node positions 
V(fall.graph)$long <- meta.fall.ab$Lon.50.
V(fall.graph)$lat <- meta.fall.ab$Lat.50.
 
V(spring.graph)$long <- meta.spring.ab$Lon.50.
V(spring.graph)$lat <- meta.spring.ab$Lat.50.

#Node composition(only individuals from one region? or more?)
fall.comp <- meta.fall.ab %>% rowwise() %>%
  mutate(comp = length(which(c(prop.ab.central, 
                               prop.ab.eastern, 
                               prop.ab.western)==0)),
         single.reg = case_when(comp == 2 & prop.ab.western != 0 ~ "West",
                                comp == 2 & prop.ab.central != 0 ~ "Central",
                                comp == 2 & prop.ab.eastern != 0 ~ "East",
                                .default = NA)) 

spring.comp <- meta.spring.ab %>% rowwise() %>%
  mutate(comp = length(which(c(prop.ab.central, 
                               prop.ab.eastern, 
                               prop.ab.western)==0)),
         single.reg = case_when(comp == 2 & prop.ab.western != 0 ~ "West",
                                comp == 2 & prop.ab.central != 0 ~ "Central",
                                comp == 2 & prop.ab.eastern != 0 ~ "East",
                                .default = NA)) 

V(fall.graph)$node.comp <- fall.comp$comp
V(spring.graph)$node.comp <- spring.comp$comp

V(fall.graph)$single.reg <- fall.comp$single.reg
V(spring.graph)$single.reg <- spring.comp$single.reg
  
# Fall network structure
undirected.fall.graph <- as.undirected(fall.graph, mode = "collapse",
                                       edge.attr.comb = list(weight = "sum", edge.type = "ignore"))

cluster_label_prop(undirected.fall.graph)

fall.label.prop <- concensusCluster(graph = undirected.fall.graph, thresh = 0.5, algiter = 1000)
fall.infomap <- cluster_infomap(fall.graph)
fall.walktrap <- cluster_walktrap(fall.graph)

modularity(fall.graph, fall.label.prop$`community structure`$membership)
modularity(fall.graph, fall.infomap$membership)
modularity(fall.graph, fall.walktrap$membership)

V(fall.graph)$label.prop.comm <- fall.label.prop$`community structure`$membership
V(fall.graph)$info.map.comm <- fall.infomap$membership
V(fall.graph)$walktrap.comm  <- fall.walktrap$membership

# spring network structure
undirected.spring.graph <- as.undirected(spring.graph, mode = "collapse",
                                         edge.attr.comb = list(weight = "sum", edge.type = "ignore"))

spring.label.prop <- concensusCluster(graph = undirected.spring.graph, thresh = 0.5, algiter = 1000)
spring.infomap <- cluster_infomap(spring.graph)
spring.walktrap <- cluster_walktrap(spring.graph)

modularity(spring.graph, spring.label.prop$`community structure`$membership)
modularity(spring.graph, spring.infomap$membership)
modularity(spring.graph, spring.walktrap$membership)

V(spring.graph)$label.prop.comm <- spring.label.prop$`community structure`$membership
V(spring.graph)$infomap.comm <- spring.infomap$membership
V(spring.graph)$walktrap.comm <- spring.walktrap$membership

# Fall and spring betweenness centrality 

# Fall migratory network without spring edges 
fall.e <- which(E(fall.graph)$edge.type == "spring")
fall.graph.disc <- fall.graph - edge(fall.e)

V(fall.graph)$betweenness <- betweenness(fall.graph.disc, directed = T, weights = 1/E(fall.graph.disc)$weight) 

# spring migratory network without fall edges  
spring.e <- which(E(spring.graph)$edge.type == "fall")
spring.graph.disc <- spring.graph - edge(spring.e)

V(spring.graph)$betweenness <- betweenness(spring.graph.disc, directed = T, weights = 1/E(spring.graph.disc)$weight) 

# Fall and spring bridge betweenness
fall.graph.brd <- fall.graph.disc 
spring.graph.brd <- spring.graph.disc 

#E(fall.graph.brd)$weight  <- 1/E(fall.graph.disc)$weight
#E(spring.graph.brd)$weight <- 1/E(spring.graph.disc)$weight

V(fall.graph)$bridge.betweenness <- bridge(fall.graph.brd,  nodes =as.character(V(fall.graph.brd)), communities = V(fall.graph)$label.prop.comm , directed = T)$`Bridge Strength`
V(spring.graph)$bridge.betweenness <- bridge(spring.graph.brd, nodes =as.character(V(spring.graph.brd)), communities = V(spring.graph)$infomap.comm, directed = T)$`Bridge Strength`

# Figure 1: Fall and spring migratory network node types and stationary location clusters ----
America <- wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),]

## fall node types ---- 
fall.ggnet <- ggnetwork(fall.graph, layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), scale = F)
fall.gplot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = equinox_region, fill = "#D9D5B2", lwd = 0.1, alpha = 1)+
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = node.type), shape=21)+
  scale_size(range = c(1, 6), guide = "none")+
  scale_fill_manual(values=c("Breeding"  = "#440154FF", "Stopover" = "#FDE725FF", "Nonbreeding" = "#21908CFF"), name = "Node type")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.4), text = element_text(size = 12), legend.key = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm"))+
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## fall stationary location clusters 
fall.clustplot<- ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = equinox_region, fill = "#D9D5B2", lwd = 0.1, alpha = 1)+
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  #geom_errorbar(data = fall.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  #geom_errorbar(data = fall.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  #geom_path(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = fall.stat[(fall.stat$site_type!= "Breeding"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster)), cex = 1) +
  geom_text(data = meta.fall.ab[meta.fall.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex))+
  labs(colour = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "None",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm"))

## Spring node types ----
spring.ggnet <- ggnetwork(spring.graph, layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), scale = F)
spring.gplot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = node.type), shape=21)+
  scale_size(range = c(1, 6), name = "Node weight")+
  scale_fill_manual(values=c("Breeding"  = "#440154FF", "Stopover" = "#FDE725FF", "Nonbreeding" = "#21908CFF"), guide = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.4), text = element_text(size = 12), legend.key = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm"))

## spring stationary location clusters 
spring.clustplot<- ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  #geom_errorbar(data = spring.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  #geom_errorbar(data = spring.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  #geom_path(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = spring.stat[(spring.stat$site_type!= "Breeding"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster)), cex = 1) +
  geom_text(data = meta.spring.ab[meta.spring.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex))+
  labs(colour = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "None",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "mm"))

## Panel ----
plot_grid(fall.gplot,  spring.gplot,
          fall.clustplot, spring.clustplot) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Figure 2: Node population composition ---- 

## fall population composition ----
fall.data <- as_data_frame(fall.graph, "vertices")
fall.gplot.comp <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_scatterpie(cols = c("West", "Central", "East"), data = fall.data[fall.data$node.comp < 2,], mapping = aes(x = long, y = lat, r = 2)) +
  geom_point(data = fall.data[fall.data$node.comp == 2,], mapping = aes(x = long, y = lat, fill = single.reg), shape= 21, cex = 5, colour = "black",  show.legend = F)+
  scale_fill_manual(values = c("West" = "#D55E00", "Central" = "#009E73", "East" = "#0072B2"), name = "Region") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.4), text = element_text(size = 12), legend.key = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank())+
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## Spring population compositon ----
spring.data <- as_data_frame(spring.graph, "vertices")
spring.gplot.comp <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  geom_scatterpie(cols = c("West", "Central", "East"), data = spring.data[spring.data$node.comp < 2,], mapping = aes(x = long, y = lat, r = 2)) +
  geom_point(data = spring.data[spring.data$node.comp == 2,], mapping = aes(x = long, y = lat, fill = single.reg), shape= 21, cex = 5, colour = "black",  show.legend = F)+
  scale_fill_manual(values = c("West" = "#D55E00", "Central" = "#009E73", "East" = "#0072B2"), name = "Region") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", text = element_text(size = 12), legend.key = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## Panel ----
plot_grid(fall.gplot.comp, spring.gplot.comp)

# Figure 2: Fall and spring migratory network communities ----

## Fall network communities ----
fall.com.plot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = , fill = as.factor(label.prop.comm)), shape=21, size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "none", text = element_text(size = 12), legend.key = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank())+
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## Spring network communities -----
spring.com.plot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = , fill = as.factor(infomap.comm)), shape=21, size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "none", text = element_text(size = 12), legend.key = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank())+
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## Panel ----
plot_grid(fall.com.plot, spring.com.plot)

# Figure 3: Fall and springcentrality metrics ---- 

## Betweenness centrality in fall network ----
fall.gplot.betw <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = betweenness), shape=21, size  = 4)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Betweenness", 
                       guide = guide_colorbar(frame.colour = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.4), text = element_text(size = 12), legend.key = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank())

## Betweenness centrality in spring network ----
spring.gplot.betw <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, lwd = weight, colour = edge.type),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = betweenness), shape=21, size  = 4)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Betweenness", 
                       guide = guide_colorbar(frame.colour = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.4), text = element_text(size = 12), legend.key = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank())

## Bridge strength in fall network ----
fall.gplot.bridge.betw <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = bridge.betweenness), shape=21, size  = 4)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Bridge strength", 
                       guide = guide_colorbar(frame.colour = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.4), text = element_text(size = 12), legend.key = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank())

## Bridge strength in spring network ----
spring.gplot.bridge.betw <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = bridge.betweenness), shape=21, size  = 4)+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Bridge Strength", 
                       guide = guide_colorbar(frame.colour = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.4), text = element_text(size = 12), legend.key = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank())

# create panel
plot_grid(fall.gplot.betw , spring.gplot.betw,
          fall.gplot.bridge.betw , spring.gplot.bridge.betw)
