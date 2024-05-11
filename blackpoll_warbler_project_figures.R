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
#library(maptools)
library(ebirdst)
library(scatterpie)
library(ggnewscale)
library(cowplot)
library(ggpubr)
library(patchwork)
library(purrr)
library(stringr)
library(gt)
library(shadowtext)
library(PieGlyph)
library(ggarchery)
library(viridis)

# network specific libraries ----
library(igraph)
library(ggnetwork)
library(intergraph)
library(networktools)
library(clustAnalytics)
library(pls)
library(tnet)

# ebird
library(ebirdst)

# set ebirddist access key
# set_ebirdst_access_key("bmedjn18aoku")

# Will need to run the network analysis and construction scripts ----
#source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/Network_Analysis.R")
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Load required data for the fall
fall.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.graph.edge.list.txt", directed = TRUE)
meta.fall.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.csv")
fall.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.edge.weights.csv")
#equinox_region <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Mapping_components/Data/Fall_Equinox_affected_regionV6.shp")
bpw_range <-  read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Birdlife_international_species_distribution/SppDataRequest.shp")
fall.stat <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stationary.data.csv")
fall.move <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.all.locations.csv")

# Load required data for the spring 
spring.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.graph.edge.list.txt", directed = TRUE)
meta.spring.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.node.metadata.csv")
spring.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.edge.weights.csv")
spring.stat <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.stationary.data.csv")
spring.move <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.all.locations.csv")

# add cluster numbers tp fall and spring graph node attributes--- 
V(fall.graph)$cluster.num <- meta.fall.ab$X
V(spring.graph)$cluster.num <- meta.spring.ab$X

# add abundance weights to fall and spring graph node attributes--- 
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

# Number of individuals using each node 
V(fall.graph)$n.individuals <- meta.fall.ab$n.individuals.at.cluster
V(spring.graph)$n.individuals <- meta.spring.ab$n.individuals.at.cluster

# Node type (breeding, stopover or nonbreeding)
V(fall.graph)$node.type <- meta.fall.ab$node.type
V(spring.graph)$node.type <- meta.spring.ab$node.type

# Node population composition by number of individuals and relative abundance 
V(fall.graph)$West <- meta.fall.ab$prop.ab.western
V(fall.graph)$Central <- meta.fall.ab$prop.ab.central
V(fall.graph)$East <- meta.fall.ab$prop.ab.eastern
V(fall.graph)$Northwest <- meta.fall.ab$prop.ab.northwest

V(spring.graph)$West <- meta.spring.ab$prop.ab.western
V(spring.graph)$Central <- meta.spring.ab$prop.ab.central
V(spring.graph)$East <- meta.spring.ab$prop.ab.eastern
V(spring.graph)$Northwest  <- meta.spring.ab$prop.ab.northwest

# Node population composition by number of individuals 
V(fall.graph)$n.West <- meta.fall.ab$n.western
V(fall.graph)$n.Central <- meta.fall.ab$n.central
V(fall.graph)$n.East <- meta.fall.ab$n.eastern
V(fall.graph)$n.Northwest <- meta.fall.ab$n.northwestern

V(spring.graph)$n.West <- meta.spring.ab$n.western
V(spring.graph)$n.Central <- meta.spring.ab$n.central
V(spring.graph)$n.East <- meta.spring.ab$n.eastern
V(spring.graph)$n.Northwest  <- meta.spring.ab$n.northwestern

# Node positions 
V(fall.graph)$long <- meta.fall.ab$Lon.50.
V(fall.graph)$lat <- meta.fall.ab$Lat.50.
 
V(spring.graph)$long <- meta.spring.ab$Lon.50.
V(spring.graph)$lat <- meta.spring.ab$Lat.50.

#Node composition by region 
fall.comp <- meta.fall.ab %>% rowwise() %>%
  mutate(comp = length(which(c(prop.ab.central, 
                               prop.ab.eastern, 
                               prop.ab.western,
                               prop.ab.northwestern)==0)),
         single.reg = case_when(comp == 3 & prop.ab.western != 0 ~ "West",
                                comp == 3 & prop.ab.central != 0 ~ "Central",
                                comp == 3 & prop.ab.eastern != 0 ~ "East",
                                comp == 3 & prop.ab.northwestern != 0 ~ "Northwest",
                                .default = NA)) 

spring.comp <- meta.spring.ab %>% rowwise() %>%
  mutate(comp = length(which(c(prop.ab.central, 
                               prop.ab.eastern, 
                               prop.ab.western,
                               prop.ab.northwestern)==0)),
         single.reg = case_when(comp == 3 & prop.ab.western != 0 ~ "West",
                                comp == 3 & prop.ab.central != 0 ~ "Central",
                                comp == 3 & prop.ab.eastern != 0 ~ "East",
                                comp == 3 & prop.ab.northwestern  != 0 ~ "Northwest",
                                .default = NA)) 

V(fall.graph)$node.comp <- fall.comp$comp
V(spring.graph)$node.comp <- spring.comp$comp

V(fall.graph)$single.reg <- fall.comp$single.reg
V(spring.graph)$single.reg <- spring.comp$single.reg
  
# Fall network structure
undirected.fall.graph <- as.undirected(fall.graph, mode = "collapse",
                                       edge.attr.comb = list(weight = "sum", edge.type = "ignore"))

fall.label.prop <- concensusCluster(graph = undirected.fall.graph, thresh = 0.5, algiter = 1000)
fall.infomap <- cluster_infomap(fall.graph)
fall.walktrap <- cluster_walktrap(fall.graph, steps =6)

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
spring.walktrap <- cluster_walktrap(spring.graph, steps = 4)

modularity(spring.graph, spring.label.prop$`community structure`$membership)
modularity(spring.graph, spring.infomap$membership)
modularity(spring.graph, spring.walktrap$membership, directed = T)

V(spring.graph)$label.prop.comm <- spring.label.prop$`community structure`$membership
V(spring.graph)$infomap.comm <- spring.infomap$membership
V(spring.graph)$walktrap.comm <- spring.walktrap$membership

# Fall and spring betweenness centrality 

# Fall migratory network without spring edges 
fall.e <- which(E(fall.graph)$edge.type == "spring")
fall.graph.disc <- fall.graph - edge(fall.e)

#betweenness calculation 
V(fall.graph)$betweenness <- betweenness(fall.graph.disc, directed = T, weights = 1/E(fall.graph.disc)$weight) 
V(fall.graph)$betweenness.unweighted <- 1/betweenness(fall.graph.disc, directed = T, weights = NULL) 

#betweenness calculation (Opshal et al.)
fall.edge.list <- cbind(get.edgelist(fall.graph.disc), E(fall.graph.disc)$weight)
fall.net <- as.tnet(fall.edge.list, type = "weighted one-mode tnet")

V(fall.graph)$betweenness.TO <- betweenness_w(fall.edge.list, directed = T, alpha = 0.5)[,2]

# spring migratory network without fall edges  
spring.e <- which(E(spring.graph)$edge.type == "fall")
spring.graph.disc <- spring.graph - edge(spring.e)

#betweenness calculation 
V(spring.graph)$betweenness <- betweenness(spring.graph.disc, directed =T, weights = 1/E(spring.graph.disc)$weight) 
V(spring.graph)$betweenness.unweighted <- 1/betweenness(spring.graph.disc, directed = T, weights = NULL) 

#betweenness calculation (Opshal et al.)
spring.edge.list <- cbind(get.edgelist(spring.graph.disc), E(spring.graph.disc)$weight)
spring.net <- as.tnet(spring.edge.list, type = "weighted one-mode tnet")

V(spring.graph)$betweenness.TO <- betweenness_w(spring.edge.list, directed = T, alpha = 0.5)[,2]

# fall and spring eigenvector centrality coefficient
V(fall.graph)$eigen <- eigen_centrality(as.undirected(fall.graph.disc))$vector
V(spring.graph)$eigen <- eigen_centrality(as.undirected(spring.graph.disc))$vector

# Fall and spring degree centrality 

V(fall.graph)$degree <- degree(fall.graph.disc)
V(spring.graph)$degree <- degree(spring.graph.disc)

V(fall.graph)$in.degree <- degree(fall.graph.disc, mode = "in")
V(spring.graph)$in.degree <- degree(spring.graph.disc, mode = "in")

V(fall.graph)$out.degree <- degree(fall.graph.disc, mode = "out")
V(spring.graph)$out.degree <- degree(spring.graph.disc, mode = "out")

# fall and spring use by time 

# Extract abundance and time spent data for the fall
fall.breed.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.abundance.per.bird.csv")
fall.stat.ab <- merge(fall.stat, fall.breed.ab[,c("ab.unit", "geo_id")], by = "geo_id")

# Summed time spent in each node for the fall
fall.use.time <- fall.stat.ab %>% 
  mutate(duration = ifelse(site_type != "Stopover", 0, duration)) %>%
  # filter(!(cluster == 5 & study.site %in% c("Quebec","Mount Mansfield, Vermont, USA")),
  #        !(cluster == 7 & study.site %in% c("Nova Scotia, Canada"))) %>%
  group_by(cluster, geo_id) %>%
  summarize(time.per.node.ind =sum(duration)) %>%
  summarize(time.per.node = mean(time.per.node.ind))

# Summed time spent in each node weighed by relative abundance for the spring 
fall.use.timeab <- fall.stat.ab %>%
  mutate(duration = ifelse(site_type != "Stopover", 0, duration),
         ab.unit = ifelse(site_type != "Stopover", 0, ab.unit )) %>%
  # filter(!(cluster == 5 & study.site %in% c("Quebec","Mount Mansfield, Vermont, USA")),
  #        !(cluster == 7 & study.site %in% c("Nova Scotia, Canada"))) %>%
  group_by(cluster, geo_id) %>%
  summarize(time.per.node = sum(duration), ab.units = unique(ab.unit)) %>%
  mutate(time.per.ab.ind = time.per.node * ab.units) %>%
  group_by(cluster) %>%
  summarize(time.per.ab = sum(time.per.ab.ind))

V(fall.graph)$time.spent <- fall.use.time$time.per.node
V(fall.graph)$time.spent.ab <- fall.use.timeab$time.per.ab/(max(fall.use.timeab$time.per.ab))
  
# Extract abundance and time spent data for the spring
spring.breed.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.abundance.per.bird.csv")
spring.stat.ab <- merge(spring.stat, spring.breed.ab[,c("ab.unit", "geo_id")], by = "geo_id")

# Summed time spent in each node for the spring
spring.use.time <- spring.stat.ab %>% 
  mutate(duration = ifelse(site_type != "Stopover", 0, duration)) %>%
  group_by(cluster, geo_id) %>%
  summarize(time.per.node.ind =sum(duration)) %>%
  summarize(time.per.node = mean(time.per.node.ind))

# Summed time spent in each node weighed by relative abundance for the spring 
spring.use.timeab <- spring.stat.ab %>%
   mutate(duration = ifelse(site_type != "Stopover", 0, duration),
          ab.unit = ifelse(site_type != "Stopover", 0, ab.unit )) %>%
  group_by(geo_id, cluster) %>%
  summarize(time.per.node = sum(duration), ab.units = unique(ab.unit)) %>%
  mutate(time.per.ab.ind = time.per.node * ab.units) %>%
  group_by(cluster) %>%
  summarize(time.per.ab = sum(time.per.ab.ind))

V(spring.graph)$time.spent <- spring.use.time$time.per.node
V(spring.graph)$time.spent.ab <- spring.use.timeab$time.per.ab/max(spring.use.timeab$time.per.ab)

# fall and spring bridge centrality
fall.graph.brd <- fall.graph.disc 
spring.graph.brd <- spring.graph.disc 

#E(fall.graph.brd)$weight  <- 1/E(fall.graph.disc)$weight
#E(spring.graph.brd)$weight <- 1/E(spring.graph.disc)$weight

V(fall.graph)$bridge.strength <- bridge(fall.graph.brd,  nodes =as.character(V(fall.graph.brd)), communities = V(fall.graph)$walktrap.comm , directed = T)$`Bridge Strength`
V(spring.graph)$bridge.strength <- bridge(spring.graph.brd, nodes =as.character(V(spring.graph.brd)), communities = V(spring.graph)$wakltrap.comm, directed = T)$`Bridge Strength`
V(fall.graph)$bridge.betweenness <- bridge(fall.graph.brd,  nodes =as.character(V(fall.graph.brd)), communities = V(fall.graph)$walktrap.comm , directed = T)$`Bridge Betweenness`
V(spring.graph)$bridge.betweenness <- bridge(spring.graph.brd, nodes =as.character(V(spring.graph.brd)), communities = V(spring.graph)$wakltrap.comm, directed = T)$`Bridge Betweenness`
V(fall.graph)$bridge.indegree <- bridge(fall.graph.brd,  nodes =as.character(V(fall.graph.brd)), communities = V(fall.graph)$walktrap.comm , directed = T)$`Bridge Indegree`
V(spring.graph)$bridge.indegree <- bridge(spring.graph.brd, nodes =as.character(V(spring.graph.brd)), communities = V(spring.graph)$wakltrap.comm, directed = T)$`Bridge Indegree`

# fall and spring participation coefficient

# Code from Dai Shizuka: https://dshizuka.github.io/networkanalysis/example_usairports.html

fall.edges <- get.data.frame(fall.graph.disc)

p.fall <- vector(length=vcount(fall.graph.disc ))
for(i in 1:vcount(fall.graph.disc )){
  nei <- neighbors(fall.graph.disc , i, mode  = "all") #get all neighbors of node i 
  nei.w <- fall.edges %>% filter((from == i & to %in% nei) | (to == i & from %in% nei)) %>%
  #nei.w <- fall.edges %>% filter(to == i & from %in% nei) %>%
    mutate(nei = ifelse(from != i, from, to)) %>% group_by(nei) %>% summarize(t.weight =sum(weight))
  nei.w$comm <- fall.walktrap$membership[names=nei.w$nei]
  weight.vec <- nei.w %>% group_by(comm) %>% summarize(comm.weight = sum(t.weight))
  p.fall [i] <- 1-sum((weight.vec$comm.weight/sum(weight.vec$comm.weight))^2)
  #if (nrow(nei.w) == 0){p.fall [i] <- 0} 
}

spring.edges <- get.data.frame(spring.graph.disc )

p.spring <- vector(length=vcount(spring.graph.disc ))
for(i in 1:vcount(spring.graph.disc )){
  nei <- neighbors(spring.graph.disc , i, mode = "all") #get all neighbors of node i 
  nei.w <- spring.edges %>% filter((from == i & to %in% nei) | (to == i & from %in% nei)) %>%
  #nei.w <- spring.edges %>% filter(to == i & from %in% nei) %>%
     mutate(nei = ifelse(from != i, from, to)) %>% group_by(nei) %>% summarize(t.weight =sum(weight))
  nei.w$comm <- spring.walktrap$membership[names=nei.w$nei]
  weight.vec <- nei.w %>% group_by(comm) %>% summarize(comm.weight = sum(t.weight))
  p.spring [i] <- 1-sum((weight.vec$comm.weight/sum(weight.vec$comm.weight))^2)
  #if (nrow(nei.w) == 0){p.spring [i] <- 0} 
}

V(fall.graph)$participation.coef <-  p.fall
V(spring.graph)$participation.coef <-  p.spring

#Save graphs with their attribute data ----
fall.gdata <- igraph::as_data_frame(fall.graph, what = "vertices") %>% mutate(cluster = seq(1:vcount(fall.graph)))
spring.gdata <- igraph::as_data_frame(spring.graph, what = "vertices") %>% mutate(cluster = seq(1:vcount(spring.graph)))

# write.csv(fall.gdata, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/fall.graph.data.csv")
# write.csv(spring.gdata, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/spring.graph.data.csv")

# Figures 1 and 2: Fall and spring migratory network node types and stationary location clusters ----
America <- wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),]

# Create the equinox region for the fall network 

# load the equinox polygon
equipol <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Manual_stat_site_clustering/Layers/equipol.shp")

# load blackpoll warbler range polygon
bpw.range <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Birdlife_international_species_distribution/SppDataRequest.shp") %>%
  dplyr::filter(seasonal %in% c(2,3,4)) %>% st_union() %>% st_transform(crs(wrld_simpl))

# Find its intersection with the Blackpoll warbler's breeding range 
sf_use_s2(FALSE)
equi.region <- st_intersection(st_as_sf(America), bpw.range) %>% st_intersection(equipol)
bpw.range.full <- st_intersection(st_as_sf(America), st_union(bpw.range))

## fall node types ---- 
fall.ggnet <- ggnetwork(fall.graph, layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), scale = F)
fall.gplot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = equi.region, fill = "#D9D5B2", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(6, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.3), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = node.type), shape=21)+
  scale_size(range = c(0.8, 4), guide = "none")+
  scale_fill_manual(values=c("Breeding"  = "#440154FF", "Stopover" = "#FDE725FF", "Nonbreeding" = "#21908CFF"), name = "Node type")+
  ggtitle("(a) Fall migration network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.2),
        legend.position = c(0.18, 0.4), legend.key = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = NA),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        title = element_text(size = 8),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.line=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin= unit(c(0,0,1,0), "pt"))+
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## fall stationary location clusters 
fall.clustplot<- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_errorbar(data = fall.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.1, alpha = 0.3, color = "black") +
  geom_errorbar(data = fall.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.1, alpha = 0.3, color = "black") +
  #geom_path(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = fall.stat[(fall.stat$site_type!= "Breeding"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = as.factor(cluster)), cex = 1, shape = 21, col = "white", stroke = 0.1) +
  #geom_text(data = meta.fall.ab[meta.fall.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex), cex = 3, fontface = "bold")+
  geom_shadowtext(data = meta.fall.ab[meta.fall.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex), cex = 3, fontface = "bold", col = "black", bg.colour = "white")+
  labs(colour = "Cluster") +
  theme_bw() +
  ggtitle("(c) Fall stationary location clusters") + 
  theme(text = element_text(size = 6), legend.position = "None",
        axis.line=element_blank(),
        axis.text =element_blank(),
        axis.ticks=element_blank(),
        axis.title =element_blank(),
        title = element_text(size = 8),
        axis.ticks.length = unit(0, "pt"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.margin= unit(c(2,0,0,0), "pt"))

## Spring node types ----
spring.ggnet <- ggnetwork(spring.graph, layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), scale = F)
spring.gplot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(6, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.3)), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = node.type), shape=21)+
  scale_size(range = c(0.8, 4), breaks = c(0.1, 0.2, 0.3), name = "Node weight")+
  scale_fill_manual(values=c("Breeding"  = "#440154FF", "Stopover" = "#FDE725FF", "Nonbreeding" = "#21908CFF"), guide = "none")+
  ggtitle("(b) Spring migration network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.2),
        legend.position = c(0.15, 0.4), legend.key = element_rect(fill = "white", colour = NA),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        title = element_text(size = 8),
        legend.background = element_rect(fill = NA),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.line=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin= unit(c(0,0,1,0), "pt"))

## spring stationary location clusters 
spring.clustplot<- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_errorbar(data = spring.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.1, alpha = 0.3, color = "black") +
  geom_errorbar(data = spring.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.1, alpha = 0.3, color = "black") +
  #geom_path(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = spring.stat[(spring.stat$site_type!= "Breeding"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = as.factor(cluster)), cex = 1, shape = 21, col = "white", stroke = 0.1) +
  #geom_text(data = meta.spring.ab[meta.spring.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex), cex = 3, fontface = "bold")+
  geom_shadowtext(data = meta.spring.ab[meta.spring.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex), cex = 3, fontface = "bold", col = "black", bg.colour = "white")+
  labs(colour = "Cluster") +
  theme_bw() +
  ggtitle("(d) Spring stationary location clusters") + 
  theme(text = element_text(size = 8), legend.position = "None",
        axis.line=element_blank(),
        axis.text =element_blank(),
        axis.ticks=element_blank(),
        axis.title =element_blank(),
        title = element_text(size = 8),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.margin = unit(c(2,0,0,0), "pt"))

## Panel ----
nodes.fig <- (fall.gplot | spring.gplot)/ (fall.clustplot |spring.clustplot) #+
 # plot_annotation(tag_levels = 'a') &
 # theme(plot.tag.position = c(0.05, 0.95),
       # plot.tag = element_text(face = 'bold', size = 10))

ggsave(plot = nodes.fig, filename = "nodes.figure.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures",
       units = "cm", width = 25*1.2, height = 12*1.2, dpi = "print", bg = "white")

# Figure 2: Node population composition ---- 

## fall population composition accounting for abundance during the fall migration  ----
fall.data <- igraph::as_data_frame(fall.graph, "vertices")
fall.gplot.comp <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_pie_glyph(slices= c( "Northwest", "West", "Central", "East"), colour = "black", data = fall.data[fall.data$node.comp < 3,], mapping = aes(x = long, y = lat, radius = node.weight)) +
  scale_radius(range = c(0.6, 3), unit = "mm", guide = "none")+
  geom_point(data = fall.data[fall.data$node.comp == 3,], mapping = aes(x = long, y = lat, fill = single.reg, size = node.weight), shape= 21, colour = "black",  show.legend = F)+
  scale_size(range = c(1.2, 6), guide = "none")+
  scale_fill_manual(values = c("West" = "#D55E00", "Central" = "#009E73", "East" = "#0072B2", "Northwest" = "#F0E442"), name = "Breeding origin") +
  ggtitle("(a) Fall node use")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.15, 0.58), text = element_text(size = 10), legend.key = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))+ 
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## Spring population composition accounting for abundance during the spring migration ----
spring.data <- igraph::as_data_frame(spring.graph, "vertices")
spring.gplot.comp <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  geom_pie_glyph(slices= c("West", "Central", "East", "Northwest"), colour = "black", data = spring.data, mapping = aes(x = long, y = lat, radius = node.weight)) +
  scale_radius(range = c(0.6, 3), unit = "mm", guide = "none")+
  geom_point(data = spring.data[spring.data$node.comp == 3,], mapping = aes(x = long, y = lat, fill = single.reg, size = node.weight), shape= 21, colour = "black")+
  geom_point(data = spring.data, mapping = aes(x = long, y = lat, fill = single.reg, size = node.weight), colour = NA, shape= 21)+
  scale_size(range = c(1.2, 7), breaks = c(0.1, 0.2, 0.3), name = "Node weight")+
  scale_fill_manual(values = c("West" = "#D55E00", "Central" = "#009E73", "East" = "#0072B2", "Northwest" = "#F0E442"), guide = "none") +
  ggtitle("(b) Spring node use")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.15, 0.55), text = element_text(size = 10), legend.key = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.line=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin= unit(c(0,0,0,0), "pt"))

# # Fall population composition of the nonbreeding nodes
# fall.nbr.node.comp <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.nbr.node.composition.csv") 
# fl.nbr.node.comp <- ggplot(st_as_sf(America))+
#   geom_sf(colour = "black", fill = "#F7F7F7") +
#   coord_sf(xlim = c(-80, -50),ylim = c(-2, 15)) +
#   geom_pie_glyph(slices = c("Eastern.Region","Northwestern.Region", "Western.Region", "Central.Region"), colour = "black", data = fall.nbr.node.comp, mapping = aes(x = Lon.50., y = Lat.50., radius = tot.abundance))+
#   scale_radius(range = c(6, 12), unit = "mm", guide = "none")+
#   geom_point(data = fall.nbr.node.comp[fall.nbr.node.comp$reg.no == "single reg",], mapping = aes(x = Lon.50., y = Lat.50., size = tot.abundance, fill = single_reg), shape= 21,  show.legend = F)+
#   scale_size(range = c(2, 40), guide = "none")+
#   scale_fill_manual(values = c("Northwestern.Region" = "#F0E442", "Western.Region" = "#D55E00", "Central.Region" = "#009E73", "Eastern.Region" = "#0072B2"), name = "Breeding Region") +
#   ggtitle("Fall")+
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
#           legend.position = "None", text = element_text(size = 12), legend.key = element_blank(),
#           axis.title =element_blank(),
#           axis.text =element_blank(),
#           axis.ticks =element_blank(),
#           axis.line=element_blank(),
#           axis.ticks.length = unit(0, "pt"),
#           plot.margin= unit(c(0,0,0,0), "pt"))+ 
#     guides(fill = guide_legend(override.aes = list(size = 5)), )
# 
# # Spring population composition of the nonbreeding nodes
# spring.nbr.node.comp <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/spring.nbr.node.composition.csv") 
# spr.nbr.node.comp.plot <- ggplot(st_as_sf(America))+
#   geom_sf(colour = "black", fill = "#F7F7F7") +
#   coord_sf(xlim = c(-80, -50),ylim = c(-2, 15)) +
#   geom_pie_glyph(slices = c("Eastern.Region","Northwestern.Region", "Western.Region", "Central.Region"), colour = "black", data = spring.nbr.node.comp, mapping = aes(x = Lon.50., y = Lat.50., radius = tot.abundance))+
#   scale_radius(range = c(6, 12), unit = "mm", guide = "none")+
#   geom_point(data = spring.nbr.node.comp[spring.nbr.node.comp$reg.no == "single reg",], mapping = aes(x = Lon.50., y = Lat.50., size = tot.abundance, fill = single_reg), shape= 21,  show.legend = F)+
#   scale_size(range = c(1.5, 44), guide = "none")+
#   scale_fill_manual(values = c("Northwestern.Region" = "#F0E442", "Western.Region" = "#D55E00", "Central.Region" = "#009E73", "Eastern.Region" = "#0072B2"), name = "Breeding Region") +
#   ggtitle("Spring")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
#         legend.position = "None", text = element_text(size = 12), legend.key = element_blank(),
#         axis.title =element_blank(),
#         axis.text =element_blank(),
#         axis.ticks =element_blank(),
#         axis.ticks.length = unit(0, "pt"),
#         plot.margin = unit(c(0,0,0,0), "pt"))+ 
#   guides(fill = guide_legend(override.aes = list(size = 5)), )
# 
# # Fall population composition of the nonbreeding nodes (stopover use)
# fall.stp.node.comp <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stp.node.composition.csv") 
# fl.stp.node.comp <- ggplot(st_as_sf(America))+
#   geom_sf(colour = "black", fill = "#F7F7F7") +
#   coord_sf(xlim = c(-80, -50),ylim = c(-2, 15)) +
#   geom_pie_glyph(slices = c("Eastern.Region","Northwestern.Region", "Western.Region", "Central.Region"), colour = "black", data = fall.stp.node.comp, mapping = aes(x = Lon.50., y = Lat.50., radius = tot.abundance))+
#   scale_radius(range = c(0.6, 3), unit = "mm", guide = "none")+
#   geom_point(data = fall.stp.node.comp[fall.stp.node.comp$reg.no == "single reg",], mapping = aes(x = Lon.50., y = Lat.50., size = tot.abundance, fill = single_reg), shape= 21,  show.legend = F)+
#   scale_size(range = c(1.2, 6), guide = "none")+
#   scale_fill_manual(values = c("Northwestern.Region" = "#F0E442", "Western.Region" = "#D55E00", "Central.Region" = "#009E73", "Eastern.Region" = "#0072B2"), name = "Breeding Region") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
#         legend.position = "None", text = element_text(size = 12), legend.key = element_blank(),
#         axis.title =element_blank(),
#         axis.text =element_blank(),
#         axis.ticks =element_blank(),
#         axis.line=element_blank(),
#         axis.ticks.length = unit(0, "pt"),
#         plot.margin= unit(c(0,0,0,0), "pt"))+ 
#   guides(fill = guide_legend(override.aes = list(size = 5)), )
# 
# # Spring population composition of the nonbreeding nodes (stopover use)
# spring.stp.node.comp <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.stp.node.composition.csv") 
# spr.stp.node.comp.plot <- ggplot(st_as_sf(America))+
#   geom_sf(colour = "black", fill = "#F7F7F7") +
#   coord_sf(xlim = c(-80, -50),ylim = c(-2, 15)) +
#   geom_pie_glyph(slices = c("Eastern.Region","Northwestern.Region", "Western.Region", "Central.Region"), colour = "black", data = spring.stp.node.comp, mapping = aes(x = Lon.50., y = Lat.50., radius = tot.abundance))+
#   scale_radius(range = c(0.6, 3), unit = "mm", guide = "none")+
#   geom_point(data = spring.stp.node.comp[spring.stp.node.comp$reg.no == "single reg",], mapping = aes(x = Lon.50., y = Lat.50., size = tot.abundance, fill = single_reg), shape= 21,  show.legend = F)+
#   scale_size(range = c(1.2, 6), guide = "none")+
#   scale_fill_manual(values = c("Northwestern.Region" = "#F0E442", "Western.Region" = "#D55E00", "Central.Region" = "#009E73", "Eastern.Region" = "#0072B2"), name = "Breeding Region") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
#         legend.position = "None", text = element_text(size = 12), legend.key = element_blank(),
#         axis.title =element_blank(),
#         axis.text =element_blank(),
#         axis.ticks =element_blank(),
#         axis.ticks.length = unit(0, "pt"),
#         plot.margin = unit(c(0,0,0,0), "pt"))+ 
#   guides(fill = guide_legend(override.aes = list(size = 5)), )

## Zoom-in onto movements in the nonbreeding range 
load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.sub.graph.R")
fall.nbr.node.comp <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.nbr.node.composition.csv") 
meta.fall.sub <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.sub.csv")
fall.sub.ggnet <- ggnetwork(fall.graph.sub.weighed, layout = as.matrix(meta.fall.sub[, c("Lon.50.", "Lat.50.")]),  scale = F)
fall.gplot.comp.nbr <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-80, -54),ylim = c(-0, 13)) +
  geom_edges(data = fall.sub.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, lwd = ab.weight),
             colour = adjustcolor("black", alpha = 0.5),
             arrow = arrow(length = unit(5, "mm"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 1), guide = "none") +
  geom_pie_glyph(slices = c("Eastern.Region","Northwestern.Region", "Western.Region", "Central.Region"), colour = "black", data = fall.nbr.node.comp, mapping = aes(x = Lon.50., y = Lat.50., radius = tot.abundance))+
  scale_radius(range = c(1.2, 3), unit = "mm", guide = "none")+
  geom_point(data = fall.nbr.node.comp[fall.nbr.node.comp$reg.no == "single reg",], mapping = aes(x = Lon.50., y = Lat.50., size = tot.abundance, fill = single_reg), shape= 21,  show.legend = F)+
  scale_size(range = c(1.2, 10), guide = "none")+
  scale_fill_manual(values = c("Northwestern.Region" = "#F0E442", "Western.Region" = "#D55E00", "Central.Region" = "#009E73", "Eastern.Region" = "#0072B2"), name = "Breeding Region") +
  ggtitle("(c) Edges in South America and \n final nonbreeding node use")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", text = element_text(size = 10), legend.key = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))+ 
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## Zoom-in onto movements in the nonbreeding range 
load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.sub.graph.R")
spring.nbr.node.comp <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/spring.nbr.node.composition.csv") 
spring.ggnet.nbr <- spring.ggnet%>% filter(yend < 13, y < 13)
spring.gplot.comp.nbr <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-80, -54),ylim = c(-0, 13)) +
  geom_edges(data = spring.ggnet.nbr , mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(5, "mm"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 1), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_pie_glyph(slices = c("Eastern.Region","Northwestern.Region", "Western.Region", "Central.Region"), colour = "black", data = spring.nbr.node.comp, mapping = aes(x = Lon.50., y = Lat.50., radius = tot.abundance))+
  scale_radius(range = c(1.2, 3), unit = "mm", guide = "none")+
  geom_point(data = spring.nbr.node.comp[spring.nbr.node.comp$reg.no == "single reg",], mapping = aes(x = Lon.50., y = Lat.50., size = tot.abundance, fill = single_reg), shape= 21,  show.legend = F)+
  scale_size(range = c(1.2, 10), guide = "none")+
  scale_fill_manual(values = c("Northwestern.Region" = "#F0E442", "Western.Region" = "#D55E00", "Central.Region" = "#009E73", "Eastern.Region" = "#0072B2"), name = "Breeding Region") +
  ggtitle("(d) Edges in South America and \n initial nonbreeding node use")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", text = element_text(size = 10), legend.key = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.line=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin= unit(c(0,0,0,0), "pt"))

## Panel ----
node.comp.fig <- (fall.gplot.comp |spring.gplot.comp)

p1 <- fall.gplot.comp + inset_element(fall.gplot.comp.nbr, 0.01,  -0.05, 0.45, 0.45)
p2 <- spring.gplot.comp + inset_element(spring.gplot.comp.nbr,  0.01,  -0.05, 0.45, 0.45)

node.comp.insert <- (p1 | p2)

#ggsave(plot = node.comp.fig, filename = "Node.comp.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
#       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")

png("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures/Node.comp.png", units = "cm", width = 25*1.2, height = 8*1.2, res = 500)
node.comp.fig
dev.off()

png("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures/Node.comp.insert.png", units = "cm", width = 25*1.2, height = 12*1.2, res = 400)
node.comp.insert
dev.off()

# Figure 3: Fall and spring migratory network communities ----

## Fall network communities ----
fall.com.plot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = , fill = as.factor(walktrap.comm)), shape=21, size = 4)+
  scale_fill_viridis(discrete = T, begin= 0.2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.2),
        legend.position = c(0.15, 0.4), legend.key = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = NA),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))+ 
  guides(fill=guide_legend(title="Fall communities"))

## Spring network communities -----
spring.com.plot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = , fill = as.factor(walktrap.comm)), shape=21, size = 4)+
  scale_fill_viridis(discrete = T, begin= 0.2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.2),
        legend.position = c(0.17, 0.4), legend.key = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = NA),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))+ 
  guides(fill=guide_legend(title="Spring communities"))

## Panel ----
communities.fig <- (fall.com.plot | spring.com.plot)

ggsave(plot = communities.fig, filename = "communities.figure.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")

# # Table 1: Significance of the spring and fall network communities ----
# 
# ## Significance of the Fall network communities ----
# 
# # create multiple rewired versions of the network, apply the clustering method, then calculate scoring functions
# iter = 100
# 
# rand.data = NULL
# 
# for (i in seq(1,iter)){
# 
#   # Create randomized graph
#   rewire.fall <- rewireCpp(g = fall.graph, weight_sel="max_weight", Q = 1)
# 
#   # cluster analysis
#   rewire.cluster <- cluster_walktrap(rewire.fall)
#   rewire.comms <-  rewire.cluster$membership
# 
#   # calculate scoring function
#   rewire.score <- scoring_functions(rewire.fall, com = rewire.comms, weighted = T, type = "global")
# 
#   # Build data.frame with results
#   if (is.null(rand.data)){
#     rand.data <- as.data.frame(rewire.score)
#   }else{
#     rand.data <-rbind(rand.data,  as.data.frame(rewire.score))
#   }
# 
#   print(paste0("progress: ", as.character(i), " to ", as.character(iter)))
# }
# 
# # test whether scoring functions differ between the observed and randomized networks
# score.test.fall <- as.data.frame(NULL)
# 
# # scoring functions for the observed graph
# ob.score <- as.data.frame(scoring_functions(fall.graph, V(fall.graph)$walktrap.comm, weighted = T, type = "global"))
# 
# for (i in seq(1, length(colnames(rand.data)))){
# 
#   if (length(unique(rand.data[,i])) > 1){
#     test.results <- t.test(rand.data[,i], mu = ob.score[,i])
# 
#     score.test.fall[i,"score.function"] <- colnames(rand.data)[i]
#     score.test.fall[i,"observed.score"] <- ob.score[,i]
#     score.test.fall[i,"random.function.mean"] <- mean(rand.data[,i], na.rm = T)
#     score.test.fall[i, "random.function.se"] <- sd(rand.data[,i], na.rm = T)/sqrt(length(rand.data[,i][!is.na(rand.data[,i])]))
#     score.test.fall[i,"statistic"] <- test.results$statistic
#     score.test.fall[i,"parameter"] <- test.results$parameter
#     score.test.fall[i,"p-value"] <- test.results$p.value
#   }else{
#     score.test.fall[i,"score.function"] <- colnames(rand.data)[i]
#     score.test.fall[i,"observed.score"] <- ob.score[,i]
#     score.test.fall[i,"random.function.mean"] <- mean (rand.data[,i])
#     score.test.fall[i, "random.function.se"] <- sd(rand.data[,i], na.rm = T)/sqrt(length(rand.data[,i][!is.na(rand.data[,i])]))
#     score.test.fall[i,"statistic"] <- NA
#     score.test.fall[i,"parameter"] <- NA
#     score.test.fall[i,"p-value"] <- NA
#   }
# }
# 
# # Save the fall score test results
# write.csv(score.test.fall, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Table_data/fall.network.significance.test.csv")
# 
# ## Significance of the Spring network communities ----
# 
# # create multiple rewired versions of the network, apply the clustering method, then calculate scoring functions
# iter = 100
# 
# rand.data = NULL
# 
# for (i in seq(1,iter)){
# 
#   # Create randomized graph
#   rewire.spring <- rewireCpp(g = spring.graph, weight_sel="max_weight", Q = 1)
# 
#   # cluster analysis
#   rewire.cluster <- cluster_walktrap(rewire.spring)
#   rewire.comms <-  rewire.cluster$membership
# 
#   # calculate scoring function
#   rewire.score <- scoring_functions(rewire.spring, com = rewire.comms, weighted = T, type = "global")
# 
#   # Build data.frame with results
#   if (is.null(rand.data)){
#     rand.data <- as.data.frame(rewire.score)
#   }else{
#     rand.data <-rbind(rand.data,  as.data.frame(rewire.score))
#   }
# 
#   print(paste0("progress: ", as.character(i), " to ", as.character(iter)))
# }
# 
# # test whether scoring functions differ between the observed and randomized networks
# score.test.spring <- as.data.frame(NULL)
# 
# # scoring functions for the observed graph
# ob.score <- as.data.frame(scoring_functions(spring.graph, V(spring.graph)$walktrap.comm, weighted = T, type = "global"))
# 
# for (i in seq(1, length(colnames(rand.data)))){
# 
#   if (length(unique(rand.data[,i])) > 1){
#     test.results <- t.test(rand.data[,i], mu = ob.score[,i])
# 
#     score.test.spring [i,"score.function"] <- colnames(rand.data)[i]
#     score.test.spring [i,"observed.score"] <- ob.score[,i]
#     score.test.spring [i,"random.function.mean"] <- mean(rand.data[,i], na.rm = T)
#     score.test.spring [i, "random.function.se"] <- sd(rand.data[,i], na.rm = T)/sqrt(length(rand.data[,i][!is.na(rand.data[,i])]))
#     score.test.spring [i,"statistic"] <- test.results$statistic
#     score.test.spring [i,"parameter"] <- test.results$parameter
#     score.test.spring [i,"p-value"] <- test.results$p.value
#   }else{
#     score.test.spring[i,"score.function"] <- colnames(rand.data)[i]
#     score.test.spring[i,"observed.score"] <- ob.score[,i]
#     score.test.spring[i,"random.function.mean"] <- mean (rand.data[,i])
#     score.test.spring[i, "random.function.se"] <- sd(rand.data[,i], na.rm = T)/sqrt(length(rand.data[,i][!is.na(rand.data[,i])]))
#     score.test.spring[i,"statistic"] <- NA
#     score.test.spring[i,"parameter"] <- NA
#     score.test.spring[i,"p-value"] <- NA
#   }
# }
# 
# # Save the spring score test results
# write.csv(score.test.spring, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Table_data/spring.network.significance.test.csv")
# 
# 
# # Table 2: Stability of the spring and fall network communities ----
# 
# ## Assessment of community stability in the fall ------
# 
# # Stability metrics for observed network 
# obs.stab <- boot_alg_list(g = undirected.spring.graph, alg_list = list(walktrap
#               = cluster_walktrap))
# 
# # create multiple rewired versions of the network, apply the bootstrapping procedure to each one 
# iter = 100
# 
# rand.data.stab.fall = NULL
# 
# for (i in seq(1,iter)){
# 
#   # Create randomized graph
#   rewire.fall <- rewireCpp(g = fall.graph, weight_sel="max_weight", Q = 1)
# 
#   # boostrap stability
#   rewire.stab.score <- boot_alg_list(g = rewire.fall, alg_list = list(walktrap = cluster_walktrap))
# 
#   # Build data.frame with results
#   if (is.null(rand.data.stab.fall)){
#     rand.data.stab.fall <- as.data.frame(rewire.stab.score )
#   }else{
#     rand.data.stab.fall <-rbind(rand.data.stab.fall,  as.data.frame(rewire.stab.score ))
#   }
# 
#   print(paste0("progress: ", as.character(i), " to ", as.character(iter)))
# }

# Figure 4: Metrics of node importance ---- 

## Betweenness centrality in fall network ----
fall.gplot.betw <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  #geom_text(data = fall.ggnet, mapping = aes(x = x, y = y, label = weight), nudge_x = 4)+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = betweenness), shape=21, size  = 3)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Betweenness \ncentrality", begin  = 0.3,
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0),  max(fall.ggnet$betweenness, spring.ggnet$betweenness)))+
  ggtitle("Fall network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(-0.18, 0.5), legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## Betweenness centrality in spring network ----
spring.gplot.betw <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, lwd = weight, colour = edge.type),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10),
             curvature = 0.3)+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = betweenness), shape=21, size  = 3)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Betweenness \ncentrality", begin   = 0.3,
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0),  max(fall.ggnet$betweenness, spring.ggnet$betweenness)))+
  ggtitle("Spring network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## Time-adjusted node weight ----
fall.gplot.metric2 <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent.ab), shape=21, size  = 3)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Time-adjusted \nnode weight", begin  = 0.3,
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0), max(fall.ggnet$time.spent.ab, spring.ggnet$time.spent.ab)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(-0.18, 0.5), legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## Time-adjusted node weight ----
spring.gplot.metric2 <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent.ab), shape=21, size  = 3)+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  #geom_text(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, label = participation.coef))+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Time-adjusted \nnode weight", begin  = 0.3,
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0), max(fall.gdata$time.spent.ab, fall.gdata$time.spent.ab)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## create panel ----
metrics.fig <- (fall.gplot.betw | spring.gplot.betw)/ (fall.gplot.metric2| spring.gplot.metric2 )

ggsave(plot = metrics.fig, filename = "nodes.metrics.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")


# Figure 4.2: average time spent at each node 

# Time spent in the fall network 
fall.gplot.time.spent <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  #geom_text(data = fall.ggnet, mapping = aes(x = x, y = y, label = weight), nudge_x = 4)+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent), shape=21, size  = 3)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Average time spent \n(days)", begin  = 0.3, 
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0),  max(fall.ggnet$time.spent, spring.ggnet$timespent)))+
  ggtitle("(a) Fall network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        legend.position = c(0.18, 0.4), legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## Time spent spring network ----
spring.gplot.time.spent<- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, lwd = weight, colour = edge.type),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent), shape=21, size  = 3)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Time spent", begin  = 0.3,
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0),  max(fall.ggnet$time.spent, spring.ggnet$time.spent)))+
  ggtitle("(b) Spring network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## create panel ----
metrics.time.fig <- (fall.gplot.time.spent | spring.gplot.time.spent)

ggsave(plot = metrics.time.fig, filename = "nodes.metrics.time.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
       units = "cm", width = 20*1.2, height = 8*1.2, dpi = "print", bg = "white")

# Table 2: Table of node characteristics ----

## Fall node characteristics
fall.char <- fall.gdata %>% dplyr::select("cluster", "node.type", "node.weight", "n.individuals",
                                   "betweenness", "bridge.indegree", "time.spent", "time.spent.ab") %>% arrange(factor(node.type, levels = c("Breeding","Stopover","Nonbreeding"))) 

write_csv(fall.char, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Table_data/Fall_node_table_data.csv")

## Spring node characteristics
spring.char <- spring.gdata %>% dplyr::select("cluster", "node.type", "node.weight", "n.individuals",
                  "betweenness", "bridge.indegree", "time.spent", "time.spent.ab") %>% arrange(factor(node.type, levels = c("Breeding","Stopover","Nonbreeding")))

write_csv(spring.char, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Table_data/Spring_node_table_data.csv")

# # Figure 5 Node population composition by network community ----
# 
# ## Fall data and plot ----  
# 
# ## Create a dataframe with the proportion of individuals from each section fall community
# fall.comms.breed <- data.frame(cluster = V(fall.graph), community = V(fall.graph)$label.prop.comm) %>%
#   filter(cluster %in% fall.breed$cluster) %>% merge(fall.breed.ab, by  = "cluster") %>%
#   mutate(community = paste0("community", community))
# 
# fall.stat.ab <- merge(fall.stat, fall.comms.breed[,c("geo_id", "ab.unit", "community")], by = "geo_id")
# 
# fall.ab.by.comm <- fall.stat.ab %>% group_by(cluster, community) %>%
#   summarize(comm.ab.units = sum(ab.unit)) %>% ungroup() %>%
#   complete(cluster, community, fill = list(comm.ab.units = 0))
# 
# 
# # Convert data from wide to long
# fall.ab.by.comm <- fall.ab.by.comm %>% pivot_wider(names_from = community, values_from = comm.ab.units)
# 
# # Merge abundance data with the node metadata
# meta.fall.ab <- merge(meta.fall, fall.ab.by.comm, by.x = "vertex", by.y = "cluster") 
# 
# # Also add a column with the overall abundance at each site 
# fall.stat.ab.per.site <- fall.stat.ab %>% group_by(cluster) %>% 
#   summarize(r.abundance.at.cluster = sum (ab.unit))
# 
# #create a column that can be converted to a numeric vector
# meta.fall.ab <- transform(meta.fall.ab, num.reg.ab.vector = asplit(cbind(meta.fall.ab[,unique(fall.comms.breed$community)]), 1))
# 
# # Plot of proportional node use during the fall migration 
# 
# # Create vector of vertex shapes 
# meta.fall.ab <- meta.fall.ab %>% rowwise() %>%
#   mutate(shape_single = length(which(c(community1, community2)==0))) %>%
#   ungroup() %>%
#   mutate(shape_single = ifelse(shape_single == 1 & node.type != "Breeding", "circle", "none")) %>%
#   mutate(shape_single_breeding = ifelse(node.type == "Breeding", "square", "none"))%>%
#   mutate(shape_multiple = ifelse(shape_single == "none" & shape_single_breeding == "none", "pie", "none")) %>%
#   mutate(shape_colour_single = case_when(shape_single != "none" & community1 != 0 ~ "#D55E00",
#                                          shape_single != "none" & community2 != 0 ~ "#009E73",
#                                          .default = NA)) %>%
#   mutate(shape_colour_single_breeding = case_when(shape_single_breeding != "none" & community1 != 0 ~ "#D55E00",
#                                                   shape_single_breeding != "none" & community2 != 0 ~ "#009E73",
#                                                   .default = NA))
# 
# # Create a palette for site use by community 
# reg.ab.palette <- list(c("#D55E00", "#009E73"))
# 
# # Prepare edge colours for spring and fall edges (spring edges should not appear in this plot)
# edge.cols.fall <- fall.con.ab %>% mutate(col = case_when(
#   edge.type == "fall" ~ adjustcolor("darkgray", alpha.f = 0.9),
#   edge.type == "spring" ~ NA))
# 
# plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
#      xlim = c(-170, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)
# 
# plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
#      xlim = c(-170, -35), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)
# 
# plot(fall.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
#      vertex.shape = meta.fall.ab$shape_single_breeding, vertex.color = meta.fall.ab$shape_colour_single_breeding,
#      edge.arrow.width = 0,edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,
#      layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
#      ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)), 
#      edge.color = edge.cols.fall$col, add = T)
# 
# plot(fall.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
#      vertex.shape = meta.fall.ab$shape_single, vertex.color = meta.fall.ab$shape_colour_single,
#      edge.arrow.width = 0,edge.width = 0, edge.arrow.size = 0, edge.arrow.width = 0,
#      layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
#      ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)), add = T)
# 
# plot(fall.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
#      vertex.shape = meta.fall.ab$shape_multiple, vertex.pie = meta.fall.ab$num.reg.ab.vector,
#      vertex.pie.color = reg.ab.palette,edge.width = 0,edge.arrow.size = 0, edge.arrow.width = 0,
#      layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
#      ylim = c(-15, 70), vertex.label = NA, vertex.label.dist = 30, add = T)
# 
# ## Spring data and plot ----  
# 
# ## Create a dataframe with the proportion of individuals from each section spring community
# spring.comms.breed <- data.frame(cluster = V(spring.graph), community = V(spring.graph)$infomap.comm) %>%
#   filter(cluster %in% spring.breed$cluster) %>% merge(spring.breed.ab, by  = "cluster") %>%
#   mutate(community = paste0("community", community))
# 
# spring.comms.breed$community
# 
# spring.stat.ab <- merge(spring.stat, spring.comms.breed[,c("geo_id", "ab.unit", "community")], by = "geo_id")
# 
# spring.ab.by.comm <- spring.stat.ab %>% group_by(cluster, community) %>%
#   summarize(comm.ab.units = sum(ab.unit)) %>% ungroup() %>%
#   complete(cluster, community, fill = list(comm.ab.units = 0))
# 
# 
# # Convert data from wide to long
# spring.ab.by.comm <- spring.ab.by.comm %>% pivot_wider(names_from = community, values_from = comm.ab.units)
# 
# # Merge abundance data with the node metadata
# meta.spring.ab <- merge(meta.spring, spring.ab.by.comm, by.x = "vertex", by.y = "cluster") 
# 
# # Also add a column with the overall abundance at each site 
# spring.stat.ab.per.site <- spring.stat.ab %>% group_by(cluster) %>% 
#   summarize(r.abundance.at.cluster = sum (ab.unit))
# 
# #create a column that can be converted to a numeric vector
# meta.spring.ab <- transform(meta.spring.ab, num.reg.ab.vector = asplit(cbind(meta.spring.ab[,unique(spring.comms.breed$community)]), 1))
# 
# # Plot of proportional node use during the spring migration 
# 
# # Create vector of vertex shapes 
# meta.spring.ab <- meta.spring.ab %>% rowwise() %>%
#   mutate(shape_single = length(which(c(community1, community2, community3)==0))) %>%
#   ungroup() %>%
#   mutate(shape_single = ifelse(shape_single == 2 & node.type != "Breeding", "circle", "none")) %>%
#   mutate(shape_single_breeding = ifelse(node.type == "Breeding", "square", "none"))%>%
#   mutate(shape_multiple = ifelse(shape_single == "none" & shape_single_breeding == "none", "pie", "none")) %>%
#   mutate(shape_colour_single = case_when(shape_single != "none" & community1 != 0 ~ "#D55E00",
#                                          shape_single != "none" & community2 != 0 ~ "#009E73",
#                                          shape_single != "none" & community3 != 0 ~ "#0072B2",
#                                          .default = NA)) %>%
#   mutate(shape_colour_single_breeding = case_when(shape_single_breeding != "none" & community1 != 0 ~ "#D55E00",
#                                                   shape_single_breeding != "none" & community2 != 0 ~ "#009E73",
#                                                   shape_single_breeding != "none" & community3 != 0 ~ "#0072B2",
#                                                   .default = NA))
# 
# # Create a palette for site use by community 
# reg.ab.palette <- list(c("#D55E00","#0072B2","#009E73"))
# 
# # Prepare edge colours for spring and spring edges (spring edges should not appear in this plot)
# edge.cols.spring <- spring.con.ab %>% mutate(col = case_when(
#   edge.type == "spring" ~ adjustcolor("darkgray", alpha.f = 0.9),
#   edge.type == "spring" ~ NA))
# 
# plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
#      xlim = c(-170, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)
# 
# plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
#      xlim = c(-170, -35), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)
# 
# plot(spring.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
#      vertex.shape = meta.spring.ab$shape_single_breeding, vertex.color = meta.spring.ab$shape_colour_single_breeding,
#      edge.arrow.width = 0,edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,
#      layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
#      ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)), 
#      edge.color = edge.cols.spring$col, add = T)
# 
# plot(spring.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
#      vertex.shape = meta.spring.ab$shape_single, vertex.color = meta.spring.ab$shape_colour_single,
#      edge.arrow.width = 0,edge.width = 0, edge.arrow.size = 0, edge.arrow.width = 0,
#      layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
#      ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)), add = T)
# 
# plot(spring.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
#      vertex.shape = meta.spring.ab$shape_multiple, vertex.pie = meta.spring.ab$num.reg.ab.vector,
#      vertex.pie.color = reg.ab.palette,edge.width = 0,edge.arrow.size = 0, edge.arrow.width = 0,
#      layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
#      ylim = c(-15, 70), vertex.label = NA, vertex.label.dist = 30, add = T)


#  Figure 6: nonbreeding regions for the MC metric (supplementary) ----

# Load the data for the first nonbreeding regions and generate the first nonbreeding sites 
proj <- '+proj=aeqd +lat_0=0 +lon_0=-74'

fall.nbr <- fall.stat %>% group_by(geo_id) %>% 
  filter(period == "Non-breeding period", sitenum == max(sitenum), !is.na(StartTime), !is.na(EndTime)) %>%
  arrange(geo_id) 

fall.nbr.sf <- st_as_sf(fall.nbr, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld_simpl), remove = F)
fall.nbr.sf <- st_transform(fall.nbr.sf, st_crs(wrld_simpl)) 

fall.nbr.regions <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data/fall.nbr.regions.shp") 

## Fall nonbreeding regions plot ----
fall.nbr.sf$Breeding_region_MC <- factor(fall.nbr.sf$Breeding_region_MC , levels = c("Eastern Region", "Central Region", "Western Region", "Northwestern Region"))
First.nbr.regions <- ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7")+
  geom_sf(data = fall.nbr.regions, col = "black",fill = "lightgray", alpha = 0.9)+
  geom_errorbar(data = fall.nbr.sf, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), col = "black", linewidth = 0.5)+
  geom_errorbar(data = fall.nbr.sf, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), col = "black", linewidth = 0.5)+
  geom_sf(data = fall.nbr.sf, aes(fill = Breeding_region_MC), shape = 21, cex = 3)+
  scale_fill_manual(values = c("Northwestern Region" = "#F0E442",
                               "Central Region" = "#009E73",
                               "Eastern Region" = "#0072B2",
                               "Western Region"  = "#D55E00"), name = "Breeding origin") +
  coord_sf(xlim = c(-90, -35), ylim = c(-15, 20))+
  theme_bw()+
  ggtitle("(a) Fall nonbreeding sites")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        legend.position = c(0.8, 0.8),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

# Load the data for the second nonbreeding regions and generate the second nonbreeding sites 
spring.nbr <- spring.stat %>% group_by(geo_id) %>% 
  filter(period == "Non-breeding period") %>% filter(sitenum == min(sitenum))

spring.nbr.sf <- st_as_sf(spring.nbr, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld_simpl), remove = F)
spring.nbr.sf <- st_transform(spring.nbr.sf, st_crs(proj)) 

spring.nbr.sf <- st_as_sf(spring.nbr, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld_simpl), remove = F)
spring.nbr.sf <- st_transform(spring.nbr.sf, st_crs(proj)) 

spring.nbr.regions <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data/spring.nbr.regions.shp") 

## spring nonbreeding regions plot ----
spring.nbr.sf$Breeding_region_MC <- factor(spring.nbr.sf$Breeding_region_MC , levels = c("Eastern Region", "Central Region", "Western Region", "Northwestern Region"))
second.nbr.regions <- ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7")+
  geom_sf(data = spring.nbr.regions, col = "black", fill = "lightgray", alpha = 0.9)+
  geom_errorbar(data = spring.nbr.sf, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), col = "black", linewidth = 0.5)+
  geom_errorbar(data = spring.nbr.sf, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), col = "black", linewidth = 0.5)+
    geom_sf(data = spring.nbr.sf, aes(fill= Breeding_region_MC), shape = 21, cex = 3)+
  scale_fill_manual(values = c("Northwestern Region" = "#F0E442",
                               "Central Region" = "#009E73",
                               "Eastern Region" = "#0072B2",
                               "Western Region"  = "#D55E00"), name = "Breeding origin") +
  coord_sf(xlim = c(-90, -35), ylim = c(-15, 20))+
  theme_bw()+
  ggtitle("(b) Spring nonbreeding sites")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        legend.position = "None",
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## create panel ----
MC.nbr.regions.fig <- (First.nbr.regions  | second.nbr.regions) 
ggsave(plot = MC.nbr.regions.fig, filename = "MC.nbr.regions.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")

# Figure 7: Abundance propagation regions ---- 

#Load blackpoll warbler reference dataset to get breeding site
bpw.ref <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv") %>%
  st_as_sf(coords = c("deploy.longitude", "deploy.latitude"), crs = crs(wrld_simpl))

# Load abundance propagation region polygons 
br.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/bpw_abundance_regions_adjusted.shp") %>%
  st_transform(crs(wrld_simpl)) %>%
  st_cast("MULTIPOLYGON")

# load blackpoll warbler range polygons
bpw.range <- load_ranges(path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports",
                         species = "bkpwar",
                           resolution = "27km",
                           smoothed = T) %>%
  dplyr::filter(season %in% c("breeding", "nonbreeding")) %>% st_transform(crs(wrld_simpl))

# polygon for america 
sf_use_s2(FALSE)
bpw.range <- st_intersection(st_as_sf(America) , bpw.range)

# Plot 
ab.prop.regions.fig <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = bpw.range, aes(fill = season),col = NA, alpha = 0.7) +
  scale_fill_discrete(labels = c("Breeding range", "Nonbreeding range"), name = "", guide = guide_legend(order = 3)) +
  geom_sf(data = br.regions, aes(col = "Abundance propagation regions"), fill = NA, linewidth = 0.5) +
  scale_colour_manual(values = c("Abundance propagation regions" = "black"), name = "")+
  geom_sf_label(data = br.regions, aes(label = region, stroke = region), nudge_y = c(-10, -12,-13 ,-13), nudge_x = c(10, 10, -4,0), cex =4)+
  new_scale_fill() +
  geom_sf(data =  bpw.ref, aes(fill = "black"), col = "white", shape = 21, cex = 3)+
  scale_fill_manual(values = c("black"), labels = c("Geolocator deployment sites"), name = "", guide = guide_legend(order = 1))+
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70))+
  geom_point(data = spring.stat[spring.stat$geo_id == "WRMA04173" & spring.stat$sitenum == 5,],
             aes(x = Lon.50., y = Lat.50.), shape = 4, cex = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        legend.title=element_text(size=11),
        legend.text=element_text(size=11),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black"),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-18, "pt"),
        legend.position = c(0.22, 0.3),
        legend.margin=margin(c(-5,5,9,5)),
        legend.key = element_rect(colour = "transparent", fill = "white"))

## Save the plot ----
ggsave(plot = ab.prop.regions.fig, filename = "abundance.prop.regions.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")

# # Figure 8 individual migratory tracks for the fall migration  ----
# 
# # We create a list of plots
# fall.ind <- list()
# 
# fall.stat <- fall.stat %>% arrange(geo_id, StartTime)
# 
# # We will use fall.stat
# for (i in unique(fall.stat$geo_id)){
#   
#   ind.stat <- fall.stat[fall.stat$geo_id == i,] %>% filter(!is.na(StartTime))
#   #ind.move <- fall.move[fall.move$geo_id == i,]
#   ind.data.stat <- fall.stat[fall.stat$geo_id == i & fall.stat$sitenum !=0,]
#   
#   fall.ind[[i]] <- ggplot(st_as_sf(America))+
#     geom_sf(colour = "black", fill = "white") +
#     coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
#     #coord_sf(xlim = c(min(ind.stat$Lon.50.)-15, max(ind.stat$Lon.50.)+15),ylim = c(min(ind.stat$Lat.50.)-15, max(ind.stat$Lat.50.)+15)) +
#     geom_errorbar(data = ind.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
#     geom_errorbar(data = ind.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
#     geom_path(data = ind.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5, col = "blue") +
#     geom_point(data =  ind.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = sitenum), col = "black", pch = 21, cex = 1.5)+
#     scale_fill_continuous(low = "yellow", high = "purple")+
#     ggtitle(i) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
#           legend.position = "none",
#           plot.title=element_text(size=8, vjust=-1),
#           axis.title =element_blank(),
#           axis.text =element_blank(),
#           axis.ticks =element_blank(),
#           axis.ticks.length = unit(0, "pt"),
#           legend.spacing = unit(-5, "pt"),
#           plot.margin = unit(c(0,0,0,0), "pt"),
#           legend.key = element_rect(colour = "transparent", fill = "white"))
#   
# }
# 
# # Create a panel of plots 
# # function from : https://stackoverflow.com/questions/66688668/automatically-assemble-plots-for-patchwork-from-a-list-of-ggplots 
# plot_a_list <- function(plots, nrows, ncols) {
#   
#   patchwork::wrap_plots(plots, 
#                         nrow = nrows, ncol = ncols)
# }
# 
# fall.ind.panel1 <- plot_a_list(fall.ind[1:24], 6, 4)
# fall.ind.panel2 <- plot_a_list(fall.ind[25:44], 6, 4)
# 
# ## Save the plots ----
# ggsave(plot = fall.ind.panel1, filename = "individual.fall.movements1.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
#        units = "cm", width = 24*1.2, height = 30*1.2, dpi = "print", bg = "white")
# 
# ggsave(plot = fall.ind.panel2, filename = "individual.fall.movements2.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
#        units = "cm", width = 24*1.2, height = 30*1.2, dpi = "print", bg = "white")
# 
# # Figure 9 individual migratory tracks for the spring migration  ----
# 
# # We create a list of plots
# spring.ind <- list()
# 
# spring.stat <- spring.stat %>% arrange(geo_id, StartTime)
# 
# # We will use spring.stat
# for (i in unique(spring.stat$geo_id)){
#   
#   ind.stat <- spring.stat[spring.stat$geo_id == i,] %>% filter(!is.na(StartTime))
#   ind.move <- spring.move[spring.move$geo_id == i,]
#   ind.data.stat <- spring.stat[spring.stat$geo_id == i & spring.stat$sitenum !=0,]
#   
#   spring.ind[[i]] <- ggplot(st_as_sf(America))+
#     geom_sf(colour = "black", fill = "white") +
#     coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
#     geom_errorbar(data = ind.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
#     geom_errorbar(data = ind.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
#     geom_path(data = ind.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5, col = "blue") +
#     geom_point(data =  ind.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = sitenum), col = "black", pch = 21, cex = 1.5)+
#     scale_fill_continuous(low = "yellow", high = "purple")+
#     ggtitle(i) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
#           plot.title=element_text(size=8, vjust=-1),
#           axis.title =element_blank(),
#           axis.text =element_blank(),
#           axis.ticks =element_blank(),
#           axis.ticks.length = unit(0, "pt"),
#           legend.spacing = unit(-5, "pt"),
#           legend.position = "none",
#           plot.margin = unit(c(0,0,0,0), "pt"),
#           legend.key = element_rect(colour = "transparent", fill = "white"))
#   
# }
# 
# # Create a panel of plots 
# spring.ind.panel1 <- plot_a_list(spring.ind[1:24], 6, 4)
# spring.ind.panel2 <- plot_a_list(spring.ind[25:35], 6, 4)
# 
# ## Save the plots ----
# ggsave(plot = spring.ind.panel1, filename = "individual.spring.movements1.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
#        units = "cm", width = 24*1.2, height = 30*1.2, dpi = "print", bg = "white")
# 
# ggsave(plot = spring.ind.panel2, filename = "individual.spring.movements2.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
#        units = "cm", width = 24*1.2, height = 30*1.2, dpi = "print", bg = "white")


# Figure 10 individual migratory tracks for full year  ----

# We create a list of plots
loc.ind <- list()

#Load locations processed during the netwrok construction
geo.all <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/All.locations.csv") %>% arrange(geo_id, StartTime) %>%
  group_by(geo_id) %>%
  mutate(period = ifelse(sitenum == 1 & geo_id != "WRMA04173", "Breeding", period),
         period = ifelse(sitenum == max(sitenum) & geo_id != "WRMA04173" & Recorded_North_South_mig == "South and partial North", "Failure", period),
         period = ifelse(sitenum == max(sitenum) & geo_id == "WRMA04173", "Breeding", period))

# loop through the locations an create the plots
for (i in unique(geo.all$geo_id)){

  ind.data <- geo.all[geo.all$geo_id == i,] %>% filter(!is.na(StartTime))
  ind.data.stat <- geo.all[geo.all$geo_id == i & geo.all$sitenum !=0,]

  loc.ind[[i]] <- ggplot(st_as_sf(America))+
    geom_sf(colour = "black", fill = "white") +
    coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
    #coord_sf(xlim = c(min(ind.stat$Lon.50.)-15, max(ind.stat$Lon.50.)+15),ylim = c(min(ind.stat$Lat.50.)-15, max(ind.stat$Lat.50.)+15)) +
    geom_errorbar(data = ind.data.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
    geom_errorbar(data = ind.data.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
    geom_path(data = ind.data, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, col = period), alpha = 0.9, col = "firebrick") +
    geom_point(data =  ind.data.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period, pch = period, col = period), cex = 2.5)+
    scale_shape_manual(values=c("Post-breeding migration" = 21 , "Non-breeding period"  = 22, "Pre-breeding migration" = 21, "Breeding"  = 24, "Failure" = 4)) +
    scale_colour_manual(values=c("Post-breeding migration" = "black" , "Non-breeding period"  = "white", "Pre-breeding migration" = "black", "Breeding"  = "white", "Failure" = "black")) +
    scale_fill_manual(values=c("Post-breeding migration" = "#FDE725FF" , "Non-breeding period"  = "black", "Pre-breeding migration" = "#21908CFF", "Breeding"  = "black"))+
    #scale_fill_continuous(low = "yellow", high = "purple")+
    ggtitle(i) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
          plot.title=element_text(size=8, vjust=-1),
          legend.position = "None",
          axis.title =element_blank(),
          axis.text =element_blank(),
          axis.ticks =element_blank(),
          axis.ticks.length = unit(0, "pt"),
          legend.spacing = unit(-5, "pt"),
          plot.margin = unit(c(0,0,0,0), "pt"),
          legend.key = element_rect(colour = "transparent", fill = "white"))
   #if (i == first(unique(geo.all$geo_id))){theme(legend.position = c(0.5, 0.2))} else {theme(legend.position = "None")}
}

# Create a panel of plots
# function from : https://stackoverflow.com/questions/66688668/automatically-assemble-plots-for-patchwork-from-a-list-of-ggplots
plot_a_list <- function(plots, nrows, ncols) {

  patchwork::wrap_plots(plots,
                        nrow = nrows, ncol = ncols)
}

loc.ind.panel1 <- plot_a_list(loc.ind[1:24], 6, 4)
loc.ind.panel2 <- plot_a_list(loc.ind[25:45], 6, 4)

## Save the plots ----
ggsave(plot = loc.ind.panel1, filename = "individual.movements1.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures",
       units = "cm", width = 24*1.2, height = 30*1.2, dpi = "print", bg = "white")

ggsave(plot = loc.ind.panel2, filename = "individual.movements2.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures",
       units = "cm", width = 24*1.2, height = 30*1.2, dpi = "print", bg = "white")


# Figure 11 group migratory tracks for full year  ----


## Migratory track for eastern breeders in the fall ----
east.fall.data <- geo.all %>% filter(Breeding_region_MC == "Eastern Region") %>%
  group_by(geo_id) %>%
  filter(StartTime <= StartTime[which(NB_count == 1)])
  
east.fall.mig.routes <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "white") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = east.fall.data[east.fall.data$sitenum > 0,], aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5., group = geo_id), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_errorbar(data = east.fall.data[east.fall.data$sitenum > 0,], aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5., group = geo_id), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_path(data = east.fall.data[east.fall.data$sitenum > 0,], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, col = period), alpha = 0.9, col = "firebrick") +
  geom_point(data =  east.fall.data[east.fall.data$sitenum > 0,], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period, pch = period, col = period), cex = 2.5)+
  scale_shape_manual(values=c("Post-breeding migration" = 21 , "Non-breeding period"  = 22, "Pre-breeding migration" = 21, "Breeding"  = 24, "Failure" = 4)) +
  scale_colour_manual(values=c("Post-breeding migration" = "black" , "Non-breeding period"  = "white", "Pre-breeding migration" = "black", "Breeding"  = "white", "Failure" = "black")) +
  scale_fill_manual(values=c("Post-breeding migration" = "#FDE725FF" , "Non-breeding period"  = "black", "Pre-breeding migration" = "#21908CFF", "Breeding"  = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        plot.title=element_text(size=8, vjust=-1),
        legend.position = "None",
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-5, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.key = element_rect(colour = "transparent", fill = "white"))

## Migratory track for eastern breeders in the spring ----
east.spring.data <- geo.all %>% filter(Breeding_region_MC == "Eastern Region") %>%
  group_by(geo_id) %>%
  filter(StartTime >= StartTime[which(NB_count == max(NB_count, na.rm = T))])

east.spring.mig.routes <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "white") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = east.spring.data [east.spring.data$sitenum > 0,], aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5., group = geo_id), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_errorbar(data = east.spring.data [east.spring.data$sitenum > 0,], aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5., group = geo_id), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_path(data = east.spring.data [east.spring.data$sitenum > 0,], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, col = period), alpha = 0.9, col = "firebrick") +
  geom_point(data =  east.spring.data [east.spring.data $sitenum > 0,], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period, pch = period, col = period), cex = 2.5)+
  scale_shape_manual(values=c("Post-breeding migration" = 21 , "Non-breeding period"  = 22, "Pre-breeding migration" = 21, "Breeding"  = 24, "Failure" = 4)) +
  scale_colour_manual(values=c("Post-breeding migration" = "black" , "Non-breeding period"  = "white", "Pre-breeding migration" = "black", "Breeding"  = "white", "Failure" = "black")) +
  scale_fill_manual(values=c("Post-breeding migration" = "#FDE725FF" , "Non-breeding period"  = "black", "Pre-breeding migration" = "#21908CFF", "Breeding"  = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        plot.title=element_text(size=8, vjust=-1),
        legend.position = "None",
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-5, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.key = element_rect(colour = "transparent", fill = "white"))

## Migratory track for western breeders in the fall ----
west.fall.data <- geo.all %>% filter(Breeding_region_MC %in% c("Western Region", "Northwestern Region", "Central Region")) %>%
  group_by(geo_id) %>%
  filter(StartTime <= StartTime[which(NB_count == 1)])

west.fall.mig.routes <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "white") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = west.fall.data [west.fall.data $sitenum > 0,], aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5., group = geo_id), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_errorbar(data = west.fall.data [west.fall.data $sitenum > 0,], aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5., group = geo_id), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_path(data = west.fall.data [west.fall.data $sitenum > 0,], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, col = period), alpha = 0.9, col = "firebrick") +
  geom_point(data =  west.fall.data [west.fall.data $sitenum > 0,], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = Breeding_region_MC, pch = period, col = period), cex = 2.5)+
  scale_shape_manual(values=c("Post-breeding migration" = 21 , "Non-breeding period"  = 22, "Pre-breeding migration" = 21, "Breeding"  = 24, "Failure" = 4)) +
  scale_colour_manual(values=c("Post-breeding migration" = "black" , "Non-breeding period"  = "white", "Pre-breeding migration" = "black", "Breeding"  = "white", "Failure" = "black")) +
  #scale_fill_manual(values=c("Post-breeding migration" = "#FDE725FF" , "Non-breeding period"  = "black", "Pre-breeding migration" = "#21908CFF", "Breeding"  = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        plot.title=element_text(size=8, vjust=-1),
        legend.position = "None",
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-5, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.key = element_rect(colour = "transparent", fill = "white"))

## Migratory track for western breeders in the spring ----
west.spring.data <- geo.all %>% filter(Breeding_region_MC %in% c("Western Region", "Northwestern Region")) %>%
  group_by(geo_id) %>%
  filter(StartTime >= StartTime[which(NB_count == max(NB_count, na.rm = T))])

west.spring.mig.routes <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "white") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = west.spring.data [west.spring.data $sitenum > 0,], aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5., group = geo_id), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_errorbar(data = west.spring.data [west.spring.data $sitenum > 0,], aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5., group = geo_id), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_path(data = west.spring.data, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, col = period), alpha = 0.9, col = "darkgray") +
  geom_point(data =  west.spring.data [west.spring.data $sitenum > 0,], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = Breeding_region_MC, pch = period, col = period), cex = 2.5, alpha = 0.8)+
  scale_shape_manual(values=c("Post-breeding migration" = 21 , "Non-breeding period"  = 22, "Pre-breeding migration" = 21, "Breeding"  = 24, "Failure" = 4)) +
  scale_colour_manual(values=c("Post-breeding migration" = "black" , "Non-breeding period"  = "white", "Pre-breeding migration" = "black", "Breeding"  = "white", "Failure" = "black")) +
  #scale_fill_manual(values=c("Post-breeding migration" = "#FDE725FF" , "Non-breeding period"  = "black", "Pre-breeding migration" = "#21908CFF", "Breeding"  = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        plot.title=element_text(size=8, vjust=-1),
        legend.position = "None",
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-5, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.key = element_rect(colour = "transparent", fill = "white"))

# Figure 12 threshold data and phenology dates ----
geo.all <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/All.locations.csv") %>% arrange(geo_id, StartTime)
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"

# load threshold paths 
tpaths <- findThresLocData()
timings <- geo.all %>% group_by(geo_id) %>% summarize(fall.br.departure = as.Date(first(br.departure)), 
                                                        spring.br.arrival = as.Date(first(br.arrival)),
                                                        nbr.departure = as.Date(first(nbr.departure)),
                                                        nbr.arrival = as.Date(first(nbr.arrival)))

# We create a list of plots
dates.ind.lon <- list()
dates.ind.lat <- list()

# loop through the locations an create the plots
for (i in unique(tpaths$geo_id)){
  
  #individual data
  ind.data <- tpaths[tpaths$geo_id == i,]
  
  #individual timing
  i.timing <- data.frame(event = c("br.departure", "br.arrival", "nbr.departure", "nbr.arrival"),
                            time = as.numeric(as.Date(c(timings[timings$geo_id == i,]$fall.br.departure,
                                     timings[timings$geo_id == i,]$spring.br.arrival,
                                     timings[timings$geo_id == i,]$nbr.departure,
                                     timings[timings$geo_id == i,]$nbr.arrival))))
  
  dates.ind.lon[[i]] <- ggplot(data = ind.data, aes(y = lon, x = as.numeric(as.Date(Twilight)), group = 1))+
    geom_line()+
    geom_vline(data =  i.timing, aes(xintercept = time, group = event, colour = event))+
    scale_colour_manual(values=c("br.departure" = "red" , "br.arrival"  = "red", "nbr.departure" = "blue", "nbr.arrival"  = "blue"))+
    ggtitle(i) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
          plot.title=element_text(size=8, vjust=-1),
          legend.position = "None",
          axis.title =element_blank(),
          axis.text =element_blank(),
          axis.ticks =element_blank(),
          axis.ticks.length = unit(0, "pt"),
          legend.spacing = unit(-5, "pt"),
          plot.margin = unit(c(0,0,0,0), "pt"),
          legend.key = element_rect(colour = "transparent", fill = "white"))
  
  dates.ind.lat[[i]] <- ggplot(data = ind.data, aes(y = lat, x = as.numeric(as.Date(Twilight)), group = 1))+
    geom_line()+
    geom_vline(data =  i.timing, aes(xintercept = time, group = event, colour = event))+
    scale_colour_manual(values=c("br.departure" = "red" , "br.arrival"  = "red", "nbr.departure" = "blue", "nbr.arrival"  = "blue"))+
    ggtitle(i) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
          plot.title=element_text(size=8, vjust=-1),
          legend.position = "None",
          axis.title =element_blank(),
          axis.text =element_blank(),
          axis.ticks =element_blank(),
          axis.ticks.length = unit(0, "pt"),
          legend.spacing = unit(-5, "pt"),
          plot.margin = unit(c(0,0,0,0), "pt"),
          legend.key = element_rect(colour = "transparent", fill = "white"))
}

# Create a panel of plots
# function from : https://stackoverflow.com/questions/66688668/automatically-assemble-plots-for-patchwork-from-a-list-of-ggplots
plot_a_list <- function(plots, nrows, ncols) {
  
  patchwork::wrap_plots(plots,
                        nrow = nrows, ncol = ncols)
}

loc.ind.panel1 <- plot_a_list(dates.ind.lon[1:24], 6, 4)
loc.ind.panel2 <- plot_a_list(dates.ind.lon[25:47], 6, 4)
loc.ind.panel3 <- plot_a_list(dates.ind.lat[1:24], 6, 4)
loc.ind.panel4 <- plot_a_list(dates.ind.lat[25:47], 6, 4)

## Save the plots ----
ggsave(plot = loc.ind.panel1, filename = "Movement.timings.longitude1.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures",
       units = "cm", width = 24*1.2, height = 30*1.2, dpi = "print", bg = "white")

ggsave(plot = loc.ind.panel2, filename = "Movement.timings.longitude2.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures",
       units = "cm", width = 24*1.2, height = 30*1.2, dpi = "print", bg = "white")

#Figure 14 stopovers in the nonbreeding range ----

# Fall stopovers in the nonbreeding range
fall.nbr.stp <- fall.stat %>% filter(Lat.50. < 13) %>% group_by(geo_id) %>%
  filter(length(geo_id) > 1)

fall.nbr.stp.plot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "white") +
  coord_sf(xlim = c(-80, -45),ylim = c(-5, 15)) +
  #geom_errorbar(data = fall.nbr.stp , aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, width = 0, alpha = 0.2, color = "black") +
  #geom_errorbar(data = fall.nbr.stp , aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, width = 0, alpha = 0.2, color = "black") +
  geom_path(data = fall.nbr.stp , mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, col = geo_id), col = adjustcolor("black", 0.5)) +
  geom_point(data =  fall.nbr.stp , mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period), cex = 2.5, pch= 21)+
  #geom_text(data =  fall.nbr.stp , mapping = aes(x = Lon.50., y = Lat.50., label = round(duration)), cex = 2.5)+
  scale_fill_manual(values = c("purple", "yellow"), labels = c("wintering site", "stopover"), name = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        plot.title=element_text(size=8, vjust=-1),
        legend.position = c(0.75, 0.75),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-5, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
  
# spring stopovers in the nonbreeding range
spring.nbr.stp <- spring.stat %>% filter(Lat.50. < 13) %>% group_by(geo_id) %>%
  filter(length(geo_id) > 1)

spring.nbr.stp.plot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "white") +
  coord_sf(xlim = c(-80, -45),ylim = c(-5, 15)) +
  geom_errorbar(data = spring.nbr.stp , aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, width = 0, alpha = 0.2, color = "black") +
  geom_errorbar(data = spring.nbr.stp , aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, width = 0, alpha = 0.2, color = "black") +
  geom_path(data = spring.nbr.stp , mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, col = geo_id), col = adjustcolor("black", 0.5)) +
  geom_point(data =  spring.nbr.stp , mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period), cex = 2.5, pch= 21)+
  geom_text(data =  spring.nbr.stp , mapping = aes(x = Lon.50., y = Lat.50., label = round(duration)), cex = 2.5)+
  scale_fill_manual(values = c("Non-breeding period" = "purple", "Pre-breeding migration" = "yellow"))+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
      plot.title=element_text(size=8, vjust=-1),
      legend.position = "None",
      axis.title =element_blank(),
      axis.text =element_blank(),
      axis.ticks =element_blank(),
      axis.ticks.length = unit(0, "pt"),
      legend.spacing = unit(-5, "pt"),
      plot.margin = unit(c(0,0,0,0), "pt"),
      legend.key = element_rect(colour = "transparent", fill = "white"))


# Figure 15 time spent at nonbreeding sites  ----

# Run the script for nonbreeding movements 
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Blackpoll_warbler_mapping_scripts/Blackpoll_nonbreeding_movements.R")

NB.move <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/nonbreeding.movements.csv")
NB.stat.mean <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/nonbreeding.mean.csv")

## Create new column measuring the stage of the nonbreeding season for individual birds ----
NB.move <- NB.move%>% mutate(nbr.stage = as.numeric(difftime(anytime(move.start), anytime(nbr.arrival), units = "days"))/
                                           as.numeric(difftime(anytime(nbr.departure), anytime(nbr.arrival), units = "days"))) 

## Create new column measuring the stage of the nonbreeding movement base on average start of the nonbreeding season 
avg.nbr.arrival <- as.Date(mean(yday(NB.move$nbr.arrival), na.rm = T))
NB.move <- NB.move %>% group_by(geo_id) %>% mutate(avg.nbr.arrival = anytime(avg.nbr.arrival))
year(NB.move$avg.nbr.arrival) <- ifelse(NB.move$geo_id != "WRMA04173" , year(NB.move$deploy.on.date), year(NB.move$deploy.on.date)-1)

NB.move <- NB.move%>% mutate(nbr.stage.common = as.numeric(difftime(anytime(move.start), anytime(avg.nbr.arrival), units = "days"))/
                               as.numeric(difftime(anytime(nbr.departure), anytime(avg.nbr.arrival), units = "days"))) 

## Nonbreeding movement directions ----
nbr.move.plot <- ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7", lwd = 0.3) +
  coord_sf(xlim = c(-95, -48),ylim = c(-8, 15)) +
  geom_point(data = NB.stat.mean, mapping = aes(x =  mean.lon, y =  mean.lat,fill = "darkgray"), colour = "black", cex = 3, shape = 21, stroke = 0.5) +
  scale_fill_manual(values = c("darkgray"),label = c("Stationary individuals"), name = "") +
  new_scale_fill()+
  geom_arrowsegment(data = NB.move, mapping = aes(x = start.lon, y = start.lat, xend = end.lon, yend = end.lat,
                    col = nbr.stage.common , 
                    fill = nbr.stage.common,
                    linetype = equinox.nbr.move),
            arrows = arrow(end = "last", type = "closed", length = unit(0.1, "inches")), arrow_positions = 1, lwd = 0.6)+
  scale_fill_viridis( begin = 0, end = 0.9, breaks = c("Early" = 0.25, "Middle" = 0.50, "Late" = 0.75),
                      name = "Timing")+
  scale_color_viridis( begin = 0, end = 0.9, guide = "none")+
  scale_linetype_manual(values = c("dashed", "solid"), guide = "none")+
  #geom_text(data = NB.move, mapping = aes(x = Lon.50., y = Lat.50., label = geo_id), cex = 2.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        text = element_text(size = 14),
        plot.title=element_text(size=8, vjust=-1),
        legend.position = c(0.14, 0.4),
        axis.title =element_blank(),
        #axis.text =element_blank(),
        #axis.ticks =element_blank(),
        #axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-5, "pt"),
        #plot.margin = unit(c(0,0,0,0), "pt"),
        legend.key = element_rect(colour = "transparent", fill = "white"))+
          guides(fill = guide_colourbar(order=1))

## Save the plot ----
ggsave(plot = nbr.move.plot, filename = "nbr.movements.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")


# Flight by individuals with ID V8757-096 -----

geo.id <- "V8757_096"

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# import lig data 
lig <- readLig(paste0(dir,"/Raw_light_data_", geo.id, ".lig"), skip = 1)

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/", geo.id, "_twl_times.csv"))

#load the adjusted threshold path path x0_ad
load(file = paste0(dir,"/", geo.id, "adjusted_initial_path_raw.csv"))

#Fall transoceanic flight
start <- "2012-10-25"
end <- "2012-11-12"

#first flight
f1.start <- "2012-11-01"
f1.end <- "2012-11-05"


# # plot  of light transitions 
# ggplot()+
#   geom_line(aes(x = lig$Date[lig$Date > start & lig$Date < end],
#                 y = lig$Light[lig$Date > start & lig$Date < end]))+
#   annotate("rect", xmin = anytime(f1.start), xmax = anytime(f1.end), ymin = min(lig$Light)-2, ymax = max(lig$Light)+2,
#            alpha = .1,fill = "yellow")

# Plot lat, lon and light transitions  
jpeg("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures/light_transitions_V8757_096.png",
            width = 2124 , height = 1090, quality = 100, res = 300)

par(cex.lab=1, cex.axis=1, mfrow=c(3,1), 
    mar = c(1,5,0.1,5), ps = 12)

plot(lig$Date[lig$Date > start & lig$Date < end], lig$Light[lig$Date > start & lig$Date < end], type = "o",
     ylab = "Light level", xaxt='n')
rect(anytime(f1.start), min(lig$Light)-2, anytime(f1.end), max(lig$Light)+2, col = alpha("yellow", 0.3), lty=0)

plot(anytime(twl$Twilight[twl$Twilight> start & twl$Twilight < end], tz = "UTC"), x0_ad[,1][twl$Twilight > start & twl$Twilight < end],
     ylab = "Longitude", xaxt='n')
rect(anytime(f1.start), min(x0_ad[,1])-2, anytime(f1.end), max(x0_ad[,1])+2, col = alpha("yellow", 0.3), lty=0)

par(mar = c(5,5,0.1,5))

plot(anytime(twl$Twilight[twl$Twilight> start & twl$Twilight < end], tz = "UTC"), x0_ad[,2][twl$Twilight > start & twl$Twilight < end],
     ylab = "Latitude", xlab = "Time")
rect(anytime(f1.start), min(x0_ad[,2])-2, anytime(f1.end), max(x0_ad[,2])+2, col = alpha("yellow", 0.3), lty=0)

par(cex.lab= 1, cex.axis= 1)

dev.off()
  