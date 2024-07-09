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

# Figure 3: Plot of the breeding and nonbreeding range ----

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

# Figure 4: map of the geolocator deployment sites ----

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


# Betweenness centrality in a sample network  
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

# Figure 5: ample of light level data during a flight over the carribean ----
geo.id <- "V8296_004"

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# import lig data 
lig <- readLig(paste0(dir,"/Raw_light_data_", geo.id, ".lig"), skip = 1)

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/", geo.id, "_twl_times.csv"))

#load the adjusted threshold path path x0_ad
load(file = paste0(dir,"/", geo.id, "adjusted_initial_path_raw.csv"))

#Fall transoceanic flight
start <- "2019-09-23"
end <- "2019-10-21"

#first flight
f1.start <- "2019-09-30"
f1.end <- "2019-10-02"

#second flight
f2.start <- "2019-10-13 12:00"
f2.end <- "2019-10-14"

# Plot lat, lon and light transitions  
# jpeg("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures/light_transitions_V8757_096.png",
#      width = 2124 , height = 1090, quality = 100, res = 300)

par(cex.lab=1.4)
par(cex.axis=1.4)
par(mfrow=c(2,1), mar = c(5,5,0.1,5))
plot(jitter(as.numeric(lig$Date[lig$Date > start & lig$Date < end])), lig$Light[lig$Date > start & lig$Date < end], type = "o",
     ylab = "Light level", xlab = "Time")
rect(anytime(f1.start), min(lig$Light)-2, anytime(f1.end), max(lig$Light)+2, col = alpha("yellow", 0.2), lty=0)
rect(anytime(f2.start), min(lig$Light)-2, anytime(f2.end), max(lig$Light)+2, col = alpha("yellow", 0.2), lty=0)

plot(anytime(twl$Twilight[twl$Twilight> start & twl$Twilight < end]), x0_ad[,1][twl$Twilight > start & twl$Twilight < end],
     ylab = "Longitude", xlab = "Time")
rect(anytime(f1.start), min(x0_ad[,1])-2, anytime(f1.end), max(x0_ad[,1])+2, col = alpha("yellow", 0.2), lty=0)
rect(anytime(f2.start), min(x0_ad[,1])-2, anytime(f2.end), max(x0_ad[,1])+2, col = alpha("yellow", 0.2), lty=0)

# plot(anytime(twl$Twilight[twl$Twilight > start & twl$Twilight < end]), x0_ad[,2][twl$Twilight > start & twl$Twilight < end],
#      ylab = "Latitude", xlab = "Time")
# rect(anytime(f1.start), min(x0_ad[,2])-2, anytime(f1.end), max(x0_ad[,2])+2, col = alpha("yellow", 0.2), lty=0)
# rect(anytime(f2.start), min(x0_ad[,2])-2, anytime(f2.end), max(x0_ad[,2])+2, col = alpha("yellow", 0.2), lty=0)
par(cex.lab= 1)
par(cex.axis= 1)

#dev.off()

# Figure 6: the fall and spring networks  ----

## fall node types ---- 
fall.ggnet <- ggnetwork(fall.graph, layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), scale = F)
fall.gplot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  geom_sf(data = equi.region, fill = "#D9D5B2", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(6, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.3), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = node.type), shape=21)+
  scale_size(range = c(1.6, 8), guide = "none")+
  scale_fill_manual(values=c("Breeding"  = "#440154FF", "Stopover" = "#FDE725FF", "Nonbreeding" = "#21908CFF"), name = "Node type")+
  ggtitle("(a) Fall migration network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.2),
        legend.position = c(0.18, 0.4), legend.key = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = NA),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        title = element_text(size = 14),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.line=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin= unit(c(6,6,6,6), "pt"))+
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## Spring node types ----
spring.ggnet <- ggnetwork(spring.graph, layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), scale = F)
spring.gplot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(6, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.3)), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = node.type), shape=21)+
  scale_size(range = c(1.6, 8), breaks = c(0.1, 0.2, 0.3), name = "Node weight")+
  scale_fill_manual(values=c("Breeding"  = "#440154FF", "Stopover" = "#FDE725FF", "Nonbreeding" = "#21908CFF"), guide = "none")+
  ggtitle("(b) Spring migration network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.2),
        legend.position = c(0.15, 0.4), legend.key = element_rect(fill = "white", colour = NA),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        title = element_text(size = 14),
        legend.background = element_rect(fill = NA),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.line=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin= unit(c(6,6,6,6), "pt"))

# Panel ----
nodes.fig <- (fall.gplot | spring.gplot) 

ggsave(plot = nodes.fig, filename = "nodes.figure.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures",
       units = "cm", width = 25*1.2, height = 12*1.2, dpi = "print", bg = "white")

## fall population composition accounting for abundance during the fall migration  ----
fall.data <- igraph::as_data_frame(fall.graph, "vertices")
fall.gplot.comp <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
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
        legend.position = c(0.15, 0.58), text = element_text(size = 14), legend.key = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(6,6,6,6), "pt"))+ 
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## Spring population composition accounting for abundance during the spring migration ----
spring.data <- igraph::as_data_frame(spring.graph, "vertices")
spring.gplot.comp <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
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
        legend.position = c(0.15, 0.55), text = element_text(size = 14), legend.key = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.line=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin= unit(c(6,6,6,6), "pt"))

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
        legend.position = "None", text = element_text(size = 12), legend.key = element_blank(),
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
        legend.position = "None", text = element_text(size = 12), legend.key = element_blank(),
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

png("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures/Node_comp.png", units = "cm", width = 28*1.2, height = 15*1.2, res = 400)
node.comp.insert
dev.off()

# Figure 4: Fall and spring migratory network communities ----

## Fall network communities ----
fall.com.plot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
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
        plot.margin = unit(c(6,6,6,6), "pt"))+ 
  guides(fill=guide_legend(title="Fall communities"))

## Spring network communities -----
spring.com.plot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
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
        plot.margin = unit(c(6,6,6,6), "pt"))+ 
  guides(fill=guide_legend(title="Spring communities"))

## Panel ----
communities.fig <- (fall.com.plot | spring.com.plot)

ggsave(plot = communities.fig, filename = "communities.figure.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")

#  Figure 5: nonbreeding regions for the MC metric (supplementary) ----

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
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        title = element_text(size = 16),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(6,6,6,6), "pt"))

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
        title = element_text(size = 16),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(6,6,6,6), "pt"))

## create panel ----
MC.nbr.regions.fig <- (First.nbr.regions  | second.nbr.regions) 
ggsave(plot = MC.nbr.regions.fig, filename = "MC.nbr.regions.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")


# Figure 6: Metrics of node importance ---- 

## Betweenness centrality in fall network ----
fall.gplot.betw <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  #geom_text(data = fall.ggnet, mapping = aes(x = x, y = y, label = weight), nudge_x = 4)+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = betweenness.TO), shape=21, size  = 5)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Betweenness \ncentrality", begin  = 0.3,
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0),  max(fall.ggnet$betweenness.TO, spring.ggnet$betweenness.TO)))+
  ggtitle("Fall network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.5), legend.key = element_blank(),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        plot.title = element_text(size=16),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(6,6,6,6), "pt"))

## Betweenness centrality in spring network ----
spring.gplot.betw <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, lwd = weight, colour = edge.type),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10),
             curvature = 0)+
  #geom_nodelabel(data = spring.ggnet, aes(label = cluster.num, x = x, y = y), nudge_x =  2, nudge_y =  2)+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = betweenness.TO), shape=21, size  = 5)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Betweenness \ncentrality", begin   = 0.3,
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0),  max(fall.ggnet$betweenness.TO, spring.ggnet$betweenness.TO)))+
  ggtitle("Spring network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", legend.key = element_blank(),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        plot.title = element_text(size=16),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(6,6,6,6), "pt"))

betw.fig <- (fall.gplot.betw | spring.gplot.betw)

ggsave(plot = betw.fig, filename = "nodes.betweenness.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")

## Fall time-adjusted node weight ----
fall.gplot.metric2 <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent.ab), shape=21, size  = 5)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Time-adjusted \nnode weight", begin  = 0.3,
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0), max(fall.ggnet$time.spent.ab, spring.ggnet$time.spent.ab)))+
  ggtitle("Fall network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.5), legend.key = element_blank(),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        plot.title = element_text(size=16),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(6,6,6,6), "pt"))

## Spring time-adjusted node weight ----
spring.gplot.metric2 <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent.ab), shape=21, size  = 5)+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Time-adjusted \nnode weight", begin  = 0.3,
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0), max(fall.gdata$time.spent.ab, fall.gdata$time.spent.ab)))+
  ggtitle("Spring network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", legend.key = element_blank(),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        plot.title = element_text(size=16),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(6,6,6,6), "pt"))

## create panel ----
ats.fig <- (fall.gplot.metric2| spring.gplot.metric2)

ggsave(plot = ats.fig, filename = "nodes.adjtimsp.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")


# Sample plot of longitude and latitude for a bird with the equinox highlighted 

# Get time during the equinox 
equi.t <- twl[twl$Twilight > "2019-09-09" & twl$Twilight < "2019-10-06" | twl$Twilight > "2020-03-09" & twl$Twilight < "2020-04-06",]
lat.equi <- x0_ad[,2][which(twl$Twilight > "2019-09-09" & twl$Twilight < "2019-10-06" | twl$Twilight > "2020-03-09" & twl$Twilight < "2020-04-06")]

# plot of longitude and latitude 
par(mfrow = c(2,1), mar = c(1, 4, 3, 3))
plot(twl$Twilight, x0_ad[,1], ylab = "longitude", xlab = NA, type = "o")
# abline(v = anytime(arr.nbr))
# abline(v = anytime(dep.nbr))
abline(v = fall.equi, col = "orange")
abline(v = spring.equi, col = "orange")
par(mar = c(3, 4, 3, 3))
plot(twl$Twilight, x0_ad[,2], ylab = "latitude", xlab = NA, type = "o")
points(equi.t$Twilight, lat.equi, col = "red")
# abline(v = anytime(arr.nbr))
# abline(v = anytime(dep.nbr))
abline(v = fall.equi, col = "orange")
abline(v = spring.equi, col = "orange")


# Figure 7 stopovers in the nonbreeding range ----

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
  scale_fill_manual(values = c("#009E73", "yellow"), labels = c("wintering site", "stopover"), name = "")+
  ggtitle("South American stopovers in the fall") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        plot.title=element_text(size=14, vjust=-1),
        legend.text=element_text(size=12),
        legend.position = c(0.75, 0.75),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-5, "pt"),
        plot.margin = unit(c(6,6,6,6), "pt"),
        legend.key = element_rect(colour = "transparent", fill = "white"))


# spring stopovers in the nonbreeding range
spring.nbr.stp <- spring.stat %>% filter(Lat.50. < 13) %>% group_by(geo_id) %>%
  filter(length(geo_id) > 1)


# addition of non-breeding movements with a predominent northward direction 
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Blackpoll_warbler_mapping_scripts/Blackpoll_nonbreeding_movements.R")

NB.move <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/nonbreeding.movements.csv")
NB.move.stp <- NB.move %>% filter(month(move.start) %in% c(3, 4))

spring.nbr.stp.plot <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "white") +
  coord_sf(xlim = c(-80, -45),ylim = c(-5, 15)) +
  #geom_errorbar(data = spring.nbr.stp , aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, width = 0, alpha = 0.2, color = "black") +
  #geom_errorbar(data = spring.nbr.stp , aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, width = 0, alpha = 0.2, color = "black") +
  geom_path(data = spring.nbr.stp , mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, col = geo_id), col = adjustcolor("black", 0.5)) +
  geom_point(data =  spring.nbr.stp , mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period), cex = 2.5, pch= 21)+
  #geom_text(data =  spring.nbr.stp , mapping = aes(x = Lon.50., y = Lat.50., label = round(duration)), cex = 2.5)+
  scale_fill_manual(values = c("#009E73", "yellow"), labels = c("Wintering site", "Stopover"), name = "")+
  geom_segment(data = NB.move.stp, mapping = aes(x = start.lon, y = start.lat, xend = end.lon, yend = end.lat), col = adjustcolor("black", 0.5)) +
  new_scale_fill()+
  geom_point(data = NB.move.stp, mapping = aes(x = start.lon, y = start.lat), col = "black", fill = "#009E73", cex = 2.5, pch= 21) +
  geom_point(data = NB.move.stp, mapping = aes(x = end.lon, y= end.lat, fill = "orange"),  col = "black", cex = 2.5, pch= 21) +
  scale_fill_manual(values = c("orange"), labels = c("Extended stopover"), name = "")+
  ggtitle("South American stopovers in the spring") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        plot.title=element_text(size=14, vjust=-1),
        legend.text=element_text(size=12),
        legend.position = c(0.8, 0.8),
        legend.background = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-30, "pt"),
        plot.margin = unit(c(6,6,6,6), "pt"),
        legend.key = element_rect(colour = "transparent", fill = "white"))

## create panel ----
nbr.stops <- (fall.nbr.stp.plot | spring.nbr.stp.plot) 
ggsave(plot = nbr.stops, filename = "Nonbreeding.stopovers.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_defense/Presentation_figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")

# Figure 8 individual migratory tracks for full year  ----

#Load locations processed during the netwrok construction
geo.all <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/All.locations.csv") %>% arrange(geo_id, StartTime) %>%
  group_by(geo_id) %>%
  mutate(period = ifelse(sitenum == 1 & geo_id != "WRMA04173", "Breeding", period),
         period = ifelse(sitenum == max(sitenum) & geo_id != "WRMA04173" & Recorded_North_South_mig == "South and partial North", "Failure", period),
         period = ifelse(sitenum == max(sitenum) & geo_id == "WRMA04173", "Breeding", period))

i <- "V8296_004"

# get location data for inddividual i
ind.data <- geo.all[geo.all$geo_id == i,] %>% filter(!is.na(StartTime))
ind.data.stat <- geo.all[geo.all$geo_id == i & geo.all$sitenum !=0,]

# get probability raster for individual i
load(paste0("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", i ,"/", i,"_SGAT_GroupedThreshold_fit.R"))

xlim <- c(-170, -40) 
ylim <- c(-40, 75)

r <- rast(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmin = xlim[1],
          xmax = xlim[2], ymin = ylim[1], ymax = ylim[2], crs = proj4string(wrld_simpl))
s <- slices(type = "intermediate", breaks = NULL, mcmc = fit, grid = r)
sk <- slice(s, sliceIndices(s))
sk.df <- as.data.frame(sk, xy = T, na.rm = T) %>% filter(lyr.1 > quantile(sk$lyr.1, probs = 0.95))

individual.track <- ggplot(st_as_sf(America))+
    scale_fill_gradient(low = adjustcolor("lightblue", alpha = 0.4), high = "#04364B", name = "Probability")+ 
    geom_tile(data = sk.df, aes(x = x, y = y, fill = lyr.1))+
    new_scale_fill()+
    geom_sf(colour = "black", fill = NA) +
    coord_sf(xlim = c(-100, -45),ylim = c(-5, 50)) +
    #coord_sf(xlim = c(min(ind.stat$Lon.50.)-15, max(ind.stat$Lon.50.)+15),ylim = c(min(ind.stat$Lat.50.)-15, max(ind.stat$Lat.50.)+15)) +
    geom_errorbar(data = ind.data.stat[ind.data.stat$duration >=2,], aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
    geom_errorbar(data = ind.data.stat[ind.data.stat$duration >=2,], aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
    geom_path(data = ind.data, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, col = period), alpha = 0.9, col = "firebrick", linewidth = 1.5) +
    geom_point(data =  ind.data.stat[ind.data.stat$duration >=2,], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period, pch = period, col = period), cex = 5)+
    scale_shape_manual(values=c("Post-breeding migration" = 21 , "Non-breeding period"  = 22, "Pre-breeding migration" = 21, "Breeding"  = 24, "Failure" = 4)) +
    scale_colour_manual(values=c("Post-breeding migration" = "black" , "Non-breeding period"  = "white", "Pre-breeding migration" = "black", "Breeding"  = "white", "Failure" = "black")) +
    scale_fill_manual(values=c("Post-breeding migration" = "#FDE725FF" , "Non-breeding period"  = "black", "Pre-breeding migration" = "#21908CFF", "Breeding"  = "black"))+
    #scale_fill_continuous(low = "yellow", high = "purple")+
    #ggtitle(first(geo.all[geo.all$geo_id == i,]$Standard.geo.id)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
          plot.title=element_text(size=8, vjust=-1),
         legend.position = "None",
          axis.title =element_blank(),
          axis.text =element_blank(),
          axis.ticks =element_blank(),
          axis.ticks.length = unit(0, "pt"),
          legend.spacing = unit(-5, "pt"),
          plot.margin = unit(c(6,6,6,6), "pt"),
          legend.key = element_rect(colour = "transparent", fill = "white")) +
    if (i == "V8757_096"){theme(panel.border = element_rect(colour = "firebrick", fill=NA, size=1))}
  #if (i == first(unique(geo.all$geo_id))){theme(legend.position = c(0.5, 0.2))} else {theme(legend.position = "None")}


