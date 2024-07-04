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
fall.walktrap <- cluster_walktrap(fall.graph, steps = 5) #cluster_walktrap(fall.graph.disc, steps =6)

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
V(fall.graph)$betweenness <- betweenness(fall.graph.disc, weights = 1/E(fall.graph.disc)$weight) 
V(fall.graph)$betweenness.unweighted <- betweenness(fall.graph.disc, directed = T, weights = NULL) 

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
Lakes <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/World_map/ne_110m_lakes/ne_110m_lakes.shp")%>%
  st_transform(crs = crs(wrld_simpl))

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

## create ggnet objects for plotting
fall.ggnet <- ggnetwork(fall.graph, layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), scale = F)
spring.ggnet <- ggnetwork(spring.graph, layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), scale = F)

# Create raster of north America ----
wrld_rast <- terra::rast("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/World_raster/NE1_50M_SR_W/NE1_50M_SR_W/NE1_50M_SR_W.tif")
plot(wrld_rast, xlim = c(-175, -20), ylim = c(-60, 80))

# Create polygons for America and American lakes
America <- wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),]
Lakes <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/World_map/ne_110m_lakes/ne_110m_lakes.shp")%>%
  st_transform(crs = crs(wrld_simpl))

# Plot of the breeding and nonbreeding range ----

#load reference data
ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")

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
bpw.range <- st_intersection(st_as_sf(America) , bpw.range) %>% mutate(colour = case_when(season == "breeding" ~ "red",
                                                                                          season == "nonbreeding" ~ "purple"))

# Raster imagery of America with overlayed polygons showing the breeding and non-breeding ranges
jpeg("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures/range.png")

plot(wrld_rast, xlim = c(-175, -20), ylim = c(-60, 80))
plot(America, col = NA, border = "black",  add = T)
plot(bpw.range, col = adjustcolor(bpw.range$colour, alpha = 0.5), add  = T)

dev.off()

# map of the geolocator deployment sites ----

#load reference data
ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")

# load blackpoll warbler range polygons
bpw.range <- load_ranges(path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports",
                         species = "bkpwar",
                         resolution = "27km",
                         smoothed = T) %>%
  dplyr::filter(season %in% c("breeding", "nonbreeding")) %>% st_transform(crs(wrld_simpl))

dep.sites <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = bpw.range, aes(fill = season),col = NA, alpha = 0.7) +
  scale_fill_discrete(labels = c("Breeding range", "Nonbreeding range"), name = "", guide = guide_legend(order = 3)) +
  new_scale_fill() +
  # all breeding sites included in the fall network
  geom_point(data =  fall.stat[fall.stat$sitenum == 1,], aes(fill = "black", x = Lon.50., y = Lat.50.), col = "white", shape = 21, cex = 3)+
  scale_fill_manual(values = c("black"), labels = c("Geolocator deployment sites"), name = "", guide = guide_legend(order = 1))+
  # Deployment location for WRMA04173
  geom_point(data =  ref.data[ref.data$geo.id == "WRMA04173",], aes(fill = "black", x = mod.deploy.lon, y = mod.deploy.lat), col = "white", shape = 21, cex = 3)+
  scale_fill_manual(values = c("black"), labels = c("Geolocator deployment sites"), name = "", guide = guide_legend(order = 1))+
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70))+
  # estimated breeding location for WRMA04173
  geom_point(data = spring.stat[spring.stat$geo_id == "WRMA04173" & spring.stat$sitenum == 5,],
             aes(x = Lon.50., y = Lat.50.), shape = 4, cex = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = NA, fill = NA),
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
ggsave(plot = dep.sites, filename = "deployment.sites.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")


# Location clusters 

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

# Fall location clusters
fall.clustplot<- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  geom_sf(data = equi.region, fill = "#D9D5B2", lwd = 0.2, alpha = 0.5) +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_errorbar(data = fall.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.2, alpha = 0.3, color = "black") +
  geom_errorbar(data = fall.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.2, alpha = 0.3, color = "black") +
  #geom_path(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = fall.stat[(fall.stat$site_type!= "Breeding"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = as.factor(cluster)), cex = 3, shape = 21, col = "white", stroke = 0.1) +
  #geom_text(data = meta.fall.ab[meta.fall.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex), cex = 3, fontface = "bold")+
  geom_shadowtext(data = meta.fall.ab[meta.fall.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex), cex = 5, fontface = "bold", col = "black", bg.colour = "white")+
  #labs(colour = "Cluster") +
  theme_bw() +
  #ggtitle("(c) Fall stationary location clusters") + 
  theme(text = element_text(size = 12), legend.position = "None",
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

# Spring location clusters 
spring.clustplot<- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_errorbar(data = spring.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.2, alpha = 0.3, color = "black") +
  geom_errorbar(data = spring.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.2, alpha = 0.3, color = "black") +
  #geom_path(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = spring.stat[(spring.stat$site_type!= "Breeding"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = as.factor(cluster)), cex = 3, shape = 21, col = "white", stroke = 0.1) +
  #geom_text(data = meta.spring.ab[meta.spring.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex), cex = 3, fontface = "bold")+
  geom_shadowtext(data = meta.spring.ab[meta.spring.ab$node.type != "Breeding",], mapping = aes(x = Lon.50., y = Lat.50., label = vertex), cex = 5, fontface = "bold", col = "black", bg.colour = "white")+
  #labs(colour = "Cluster") +
  theme_bw() +
  #ggtitle("(d) Spring stationary location clusters") + 
  theme(text = element_text(size = 12), legend.position = "None",
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
clust.fig <- (fall.clustplot |spring.clustplot)

ggsave(plot = clust.fig, filename = "cluster.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures",
       units = "cm", width = 24*1.2, height = 12*1.2, dpi = "print", bg = "white")

# Plot of blackpoll warbler abundance 

# Load blackpoll warbler abundance raster
bpw.ab <- load_raster(path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports",
                           species = "bkpwar",
                           product = "abundance",
                           period = "seasonal",
                           resolution = "9km")

bpw.ab <- project(bpw.ab, crs(wrld_simpl))
breed.ab.raster <- bpw.ab %>% as.data.frame(bpw.ab$breeding, xy = T)  %>% filter()#mutate(breeding = ifelse(breeding == 0, NA, breeding))

# Plot 
ab.fig <- ggplot(st_as_sf(America))+
   geom_sf(colour = "black", fill = "#F7F7F7") +
  scale_fill_continuous(low="thistle2", high="darkred", 
                        guide="colorbar",na.value="white")+
  geom_tile(breed.ab.raster, mapping = aes(x = x, y = y, fill = breeding), na.rm = T)+
  geom_sf(colour = "black", fill = "NA") +
  scale_colour_manual(values = c("Abundance propagation regions" = "black"), name = "")+
  new_scale_fill()+
  geom_sf_label(data = br.regions, aes(label = region), nudge_y = c(-10, -8,-13 ,-13), nudge_x = c(10, 10, -4,0), cex =4)+
  geom_sf(data = br.regions, aes(col = "Abundance propagation regions"), fill = NA, linewidth = 0.5) +
   # all breeding sites included in the fall network
  new_scale_fill()+
  geom_point(data =  fall.stat[fall.stat$sitenum == 1,], aes(fill = "black", x = Lon.50., y = Lat.50.), col = "white", shape = 21, cex = 3)+
  scale_fill_manual(values = c("black"), labels = c("Geolocator deployment sites"), name = "", guide = guide_legend(order = 1))+
  # Deployment location for WRMA04173
  geom_point(data =  ref.data[ref.data$geo.id == "WRMA04173",], aes(fill = "black", x = mod.deploy.lon, y = mod.deploy.lat), col = "white", shape = 21, cex = 3)+
  scale_fill_manual(values = c("black"), labels = c("Geolocator deployment sites"), name = "", guide = guide_legend(order = 1))+
  coord_sf(xlim = c(-170, -40),ylim = c(-5, 70))+
  # estimated breeding location for WRMA04173
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
        legend.position = "None", #c(0.22, 0.3),
        legend.margin=margin(c(-5,5,9,5)),
        legend.key = element_rect(colour = "transparent", fill = "white"))

ggsave(plot = ab.fig, filename = "abundance.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures",
       units = "cm", width = 12*1.2, height = 12*1.2, dpi = "print", bg = "white")


# Bertweenness centrality on 
g <- make_graph(c(1,2,2,3,3,4,1,4,1,5,1,6,5,7,1,7,1,8,7,8,3,8,8,9,9,10,10,11,9,12), directed = T)
set.seed(0)

plot(g, edge.arrow.size = 0)
E(g)$weight <- runif(min = 1, max = 20, n = ecount(g))

sg <- betweenness(g, directed = T, weights = 1/E(g)$weights)
sg.palette  <- hcl.colors(n = length(seq(0, max(sg)+1,1)), palette = "Viridis", rev = T) 
names(sg.palette) <- seq(0, max(sg)+1,1)

plot(g, edge.width = E(g)$weight/2, edge.arrow.size = 0,
     vertex.color = sg.palette[as.character(round(sg, digits = 0))],
     vertex.size = 20, vertex.label = NA)

legend_image <- as.raster(matrix(hcl.colors(n = length(seq(0, max(sg), 1)), palette = "Viridis", rev = F) , ncol=1))
plot(c(0,2),c(-0.01,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x= 1.25, y = seq(0,1,l=5), labels = round(seq(0,max(sg),l=5), digits = 1), cex = 1)
rasterImage(legend_image, 0.75, 0, 1,1)

