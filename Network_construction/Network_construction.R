# clear objects from workspace
rm(list=ls())

#Load libraries
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

source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis/Geolocator_analysis_helper_functions.R")

# Import data ##################################################################

# path to reference data file 
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

# location data 
geo.all <- findLocData(geo.ids = c("V8757_010",
                                   "V8296_004",
                                   "V8296_005",
                                   "V8296_006",
                                   "V8757_055",
                                   "V8757_018",
                                   "V8757_021",
                                   "V8296_015",
                                   "V8296_017",
                                   "V8296_026",
                                   "V8296_025",
                                   "V8296_007",
                                   "V8296_008",
                                   "V8757_019",
                                   "V8757_096",
                                   "V8757_134",
                                   "V8757_029",
                                   "V8757_078",
                                   "blpw09",
                                   "blpw12",
                                   "3254_001",
                                   "4068_014",
                                   "blpw14",
                                   "3254_003",
                                   "3254_008",
                                   "3254_011",
                                   "3254_057",
                                   "blpw15",
                                   "blpw25",
                                   "4105_008",
                                   "4105_009",
                                   "4105_016",
                                   "4105_017",
                                   "4210_002",
                                   "4210_004",
                                   "4210_006",
                                   "4210_010",
                                   "WRMA04173",
                                   "A",
                                   "B",
                                   "C",
                                   #"E",
                                   "D"), check_col_length = F, ref_path = ref_path)

# Define nodes as breeding, stopover, or non-breeding ##########################

geo.all <- geo.all %>% group_by(geo_id) %>% mutate(site_type = case_when(
  (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ "Breeding",
  sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ "Breeding",
  (sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ "Breeding",
  period == "Non-breeding period" & (duration >= 14 | sitenum == 1 | sitenum == max(sitenum)) ~ "Nonbreeding",
  .default = "Stopover"))

# For geolocators that were deployed in close proximity, assign one location to the geolocator deployment site
# Here this location will be the median of all geolocator deployment locations for that site. 

# we first add new columns with the median deployment locations in our reference dataset (metadata for each geolcator)
ref_data <- ref_data %>% group_by(study.site) %>% 
  mutate(mod.deploy.lon = median(deploy.longitude))%>%
  mutate(mod.deploy.lat = median(deploy.latitude))

#View(ref_data %>% group_by(study.site, geo.id) %>% summarise(disp = mean(deploy.longitude)))

# Now we can modify our location dataset
#geo.all <- merge(geo.all, ref_data[,c("geo.id", "mod.deploy.lon", "mod.deploy.lat")], by.x = "geo_id", by.y = "geo.id")

geo.all <- geo.all %>% group_by(geo_id) %>% mutate(Lon.50. = case_when(
  (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ mod.deploy.lon,
  sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ mod.deploy.lon,
  #(sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ mod.deploy.lon,
  .default = Lon.50.)) %>%
  mutate(Lat.50. = case_when(
    (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ mod.deploy.lat,
    sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ mod.deploy.lat,
   # (sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ mod.deploy.lat,
    .default = Lat.50.))

# Check the location of breeding and nonbreeding sites   
#View(geo.all[,c("geo_id", "Lon.50.", "Lat.50.", "mod.deploy.lon", "mod.deploy.lat")])

# plot locations of deployment sites (there should be 10) 
dep.sites <- geo.all[(geo.all$sitenum == 1),]
plot(wrld_simpl, col = "lightgray", xlim = c(-170, -30), ylim = c(-20, 70))
points(dep.sites$Lon.50., dep.sites$Lat.50., col = "darkred")

# plot locations of all stationary sites 
stat.sites <- geo.all[(geo.all$sitenum > 0),]
plot(wrld_simpl, col = "lightgray", xlim = c(-170, -30), ylim = c(-20, 70))
points(stat.sites$Lon.50., stat.sites$Lat.50., col = "darkred")

# Count stationary sites by season for each geolocator #########################
geo.all <- geo.all %>% arrange(geo_id, StartTime)
num <- 1
geo.all$NB_count <- NA
geo.id.ind <- geo.all$geo_id[1]

for (i in seq(1, nrow(geo.all))){
  
  if (geo.all$site_type[i] == "Nonbreeding" & geo.all$geo_id[i] == geo.id.ind){
    geo.all$NB_count[i] <- num
    num <- num + 1 
   }
  
  if (geo.all$site_type[i] == "Nonbreeding" & geo.all$geo_id[i] != geo.id.ind){
    geo.id.ind <- geo.all$geo_id[i]
    num <-1 
    geo.all$NB_count[i] <- num
    num <- num + 1 
  }
}

#View(geo.all[,c("geo_id","StartTime", "sitenum", "site_type", "NB_count")])

#Save the new columns in the reference data
write.csv(ref_data,"C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv", row.names=FALSE)

################################################################################
# Fall migration network #######################################################
################################################################################

# Clustering breeding and nonbreeding location #################################

# Either import a file with manual clusters, or create cluster in R


# Create clusters for fall stopover locations and first nonbreeding location in R 
# first we must filter our location data to retain only those relevant to the fall migration 

fall.stat <- geo.all %>% filter(sitenum > 0, site_type %in% c("Stopover","Nonbreeding"),
                                period %in% c("Post-breeding migration","Non-breeding period"),
                                Recorded_North_South_mig %in% c("Both", "South and partial North", "South"))

#get the timing of the first nonbreeding area
fall.timings.nb <- geo.all %>% filter(NB_count == 1) %>% dplyr::select(NB.first.site.arrival = StartTime)
fall.stat <- merge(fall.stat, fall.timings.nb, by = "geo_id")

#only retain the stopovers and the first nonbreeding sites occupied
fall.stat <- fall.stat %>% group_by(geo_id) %>% filter(StartTime <= NB.first.site.arrival)

# code by stack overflow user jlhoward (https://stackoverflow.com/questions/21095643/approaches-for-spatial-geodesic-latitude-longitude-clustering-in-r-with-geodesic)
# Calculates the geodesic distance between points and creates a distance matrix
geo.dist = function(df) {
  require(geosphere)
  d <- function(i,z){         # z[1:2] contain long, lat
    dist <- rep(0,nrow(z))
    dist[i:nrow(z)] <- distHaversine(z[i:nrow(z),1:2],z[i,1:2])
    return(dist)
  }
  dm <- do.call(cbind,lapply(1:nrow(df),d,df))
  return(as.dist(dm))
}

dist.matrix <- geo.dist(fall.stat[,c("Lon.50.","Lat.50.")])
fall.clust <- hclust(dist.matrix)
plot(fall.clust)
fall.stat$cluster  <- cutree(fall.clust, k = 16)

# Add breeding sites with separate clusters
fall.breed.n <- geo.all  %>% filter(site_type == "Breeding",
                                      period %in% c("Post-breeding migration","Non-breeding period")) %>% 
  group_by(study.site) %>% summarize(geo_per_site = n())

fall.breed <- geo.all  %>% filter(site_type == "Breeding",
                                    period %in% c("Post-breeding migration","Non-breeding period")) %>%
  arrange(study.site, geo_id, StartTime)

fall.breed$cluster <- rep(seq(max(fall.stat$cluster) + 1, max(fall.stat$cluster) + nrow(fall.breed.n)), fall.breed.n$geo_per_site)

fall.stat <- rbind(fall.stat, fall.breed) %>% arrange(geo_id, StartTime)

# plot stopover and nonbreeding nodes

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_path(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = site_type)) +
  theme_bw() +
  theme(text = element_text(size = 16))

# plot the clusters 

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_path(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster))) +
  labs(colour = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 16))

# Generate the network from our location data and clusters #####################

# Add a column with inter-cluster connections to our location dataset
fall.stat <- fall.stat %>% mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# Create a dataframe with the edge info 
fall.edge.df <- fall.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                            Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                            Lat.2.5., Lat.97.5., EndTime,
                                            sitenum, duration, period,study.site,
                                            Range_region, NB_count, period,
                                            site_type) %>%
filter(!is.na(next.cluster))

# calculate node locations as the median locations of the points included within 
# Note, we use the 50% quantile of the posterior distribution 
fall.node.lon <- fall.stat %>% group_by(cluster) %>% 
  summarize(node.lon = median(Lon.50.))

fall.node.lat <- fall.stat %>% group_by(cluster) %>% 
  summarize(node.lat = median(Lat.50.))

# We also assess which use was the more common for a node (breeding, nonbreeding, or stopover)
# The most common use will be the one assigned to the node
fall.node.type <- fall.stat %>% group_by(cluster, site_type) %>%
  summarize(use = n()) %>% filter(use == max(use)) %>% mutate(site_type_num = case_when(
    site_type == "Breeding" ~ 1,
    site_type == "Nonbreeding" ~ 2,
    site_type == "Stopover" ~ 3
  ))

# Create a layout with node locations 
meta <- data.frame("vertex" = seq(min(fall.edge.df$cluster), max(fall.edge.df$cluster)), 
                   "Lon.50." = fall.node.lon$node.lon,
                   "Lat.50." = fall.node.lat$node.lat,
                   "node.type" = fall.node.type$site_type,
                   "node.type.num" = fall.node.type$site_type_num)

location <- as.matrix(meta[, c("Lon.50.", "Lat.50.")])

# Create the fall network
fall.graph <- graph_from_data_frame(fall.edge.df, directed = T, vertices = meta)

# Colour palette for site type
type.palette <- rainbow(3)  

# plot the fall network over North and South America
plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = 1, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.color = type.palette[meta$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta$node.type.num)],
       pch = 16)

# Add edge weights to the fall network to show the flow of individuals #########

# We start with edges based on the number of individuals moving between nodes

# list of connections and the number of times they occur
fall.con <- fall.edge.df %>% group_by(cluster, next.cluster) %>% 
  summarize(weight = n()) 

# Create a fall network with weighed edges
fall.graph.weighed <- graph_from_data_frame(fall.con, directed = T, vertices = meta)
is_weighted(fall.graph.weighed)

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con$weight/1.5, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.3, 0.3), nrow(fall.con)),
     vertex.color = type.palette[meta$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta$node.type.num)],
       pch = 16)

################################################################################
# Spring migration network #####################################################
################################################################################

# Clustering breeding and nonbreeding location #################################

# Either import a file with manual clusters, or create cluster in R

# Create clusters for spring stopover locations and last nonbreeding locations in R 
spring.stat <- geo.all %>% filter(sitenum > 0, site_type %in% c("Stopover","Nonbreeding"),
                                period %in% c("Pre-breeding migration","Non-breeding period"),
                                Recorded_North_South_mig %in% c("Both","North", "South and partial North"))

#get the timing of the lastnonbreeding area
spring.timings.nb <- geo.all %>% group_by(geo_id) %>% filter(NB_count == max(NB_count, na.rm = T)) %>% dplyr::select(NB.last.site.arrival = StartTime)
spring.stat <- merge(spring.stat, spring.timings.nb, by = "geo_id")

#only retain the stopovers and the first nonbreeding sites occupied
spring.stat <- spring.stat %>% group_by(geo_id) %>% filter(StartTime >= NB.last.site.arrival)

# code by stack overflow user jlhoward (https://stackoverflow.com/questions/21095643/approaches-for-spatial-geodesic-latitude-longitude-clustering-in-r-with-geodesic)
# Calculates the geodesic distance between points and creates a distance matrix
geo.dist = function(df) {
  require(geosphere)
  d <- function(i,z){         # z[1:2] contain long, lat
    dist <- rep(0,nrow(z))
    dist[i:nrow(z)] <- distHaversine(z[i:nrow(z),1:2],z[i,1:2])
    return(dist)
  }
  dm <- do.call(cbind,lapply(1:nrow(df),d,df))
  return(as.dist(dm))
}

dist.matrix <- geo.dist(spring.stat[,c("Lon.50.","Lat.50.")])
spring.clust <- hclust(dist.matrix)
plot(spring.clust)
spring.stat$cluster  <- cutree(spring.clust, k = 16)

# Add breeding sites with separate clusters
spring.breed.n <- geo.all  %>% filter(site_type == "Breeding",
                                 period %in% c("Pre-breeding migration","Non-breeding period")) %>% 
  group_by(study.site) %>% summarize(geo_per_site = n())

spring.breed <- geo.all  %>% filter(site_type == "Breeding",
                                      period %in% c("Pre-breeding migration","Non-breeding period"))%>%
  arrange(study.site, geo_id, StartTime) 

spring.breed$cluster <- rep(seq(max(spring.stat$cluster) + 1, max(spring.stat$cluster) + nrow(spring.breed.n)), spring.breed.n$geo_per_site)

spring.stat <-rbind(spring.stat, spring.breed) %>% arrange(geo_id, StartTime)

# plot stopover and nonbreeding nodes 

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_path(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = site_type))+ 
  theme_bw() +
  theme(text = element_text(size = 16))

# plot the clusters 
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_path(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster))) +
  labs(colour = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 16))

# Generate the network from our location data and clusters #####################

# Add a column with inter-cluster connections to our location dataset
spring.stat  <- spring.stat %>% mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# Create a dataframe with the edge info 
spring.edge.df <- spring.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                            Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                            Lat.2.5., Lat.97.5., EndTime,
                                            sitenum, duration, period,study.site,
                                            Range_region, NB_count, period,
                                            site_type) %>%
filter(!is.na(next.cluster))

# calculate node locations as the median locations of the points included within 
# Note, we use the 50% quantile of the posterior distribution 
spring.node.lon <- spring.stat %>% group_by(cluster) %>% 
  summarize(node.lon = median(Lon.50.))

spring.node.lat <- spring.stat %>% group_by(cluster) %>% 
  summarize(node.lat = median(Lat.50.))

# We also assess which use was the more common for a node (breeding, nonbreeding, or stopover)
# The most common use will be the one assigned to the node
spring.node.type <- spring.stat %>% group_by(cluster, site_type) %>%
  summarize(use = n()) %>% filter(use == max(use)) %>% mutate(site_type_num = case_when(
    site_type == "Breeding" ~ 1,
    site_type == "Nonbreeding" ~ 2,
    site_type == "Stopover" ~ 3
  ))

# Create a layout with node locations 
meta <- data.frame("vertex" = seq(min(spring.stat$cluster), max(spring.stat$cluster)), 
                   "Lon.50." = spring.node.lon$node.lon,
                   "Lat.50." = spring.node.lat$node.lat,
                   "node.type.num" = spring.node.type$site_type_num)

location <- as.matrix(meta[, c("Lon.50.", "Lat.50.")])

# Create the spring network
spring.graph <- graph_from_data_frame(spring.edge.df, directed = T, vertices = meta)

# Colour palette for site type
type.palette <- rainbow(3)  

# plot the fall network over North and South America
plot(wrld_simpl, xlim = c(-170, -30),ylim = c(-15, 70))
plot(spring.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = 1, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.color = type.palette[meta$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta$node.type.num)],
       pch = 16)

# Add edge weights to the spring network to show the flow of individuals #########

# We start with edges based on the number of individuals moving between nodes

# list of connections and the number of times they occur
spring.con <- spring.edge.df %>% group_by(cluster, next.cluster) %>% 
  summarize(weight = n()) 

# Create a new spring network
spring.graph.weighed <- graph_from_data_frame(spring.con, directed = T, vertices = meta)
is_weighted(spring.graph.weighed)

plot(wrld_simpl, xlim = c(-170, -30),ylim = c(-15, 70))
plot(spring.graph.weighed, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con$weight/1.5, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.3, 0.3), nrow(spring.con)),
     vertex.color = type.palette[meta$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta$node.type.num)],
       pch = 16)

################################################################################
# Nonbreeding movements network ################################################
################################################################################

# Clustering breeding and nonbreeding location #################################

# Either import a file with manual clusters, or create cluster in R

# Create clusters for nonbreeding sites in R 
NB.stat <- geo.all %>% filter(sitenum > 0, site_type %in% c("Nonbreeding"),
                                  period %in% c("Non-breeding period"),
                                  Recorded_North_South_mig %in% c("Both", "South and partial North"))

# #get the timing of the first nonbreeding area
# NB.start.time <- geo.all %>% filter(NB_count == 1) %>% dplyr::select(NB.first.site.arrival = StartTime)
# NB.stat <- merge(NB.stat, NB.start.time, by = "geo_id")
# 
# #get the timing of the lastnonbreeding area
# NB.end.time <- geo.all %>% group_by(geo_id) %>% filter(NB_count == max(NB_count, na.rm = T)) %>% dplyr::select(NB.last.site.arrival = StartTime)
# NB.stat <- merge(NB.stat.stat, NB.end.time.nb, by = "geo_id")

# code by stack overflow user jlhoward (https://stackoverflow.com/questions/21095643/approaches-for-spatial-geodesic-latitude-longitude-clustering-in-r-with-geodesic)
# Calculates the geodesic distance between points and creates a distance matrix
geo.dist = function(df) {
  require(geosphere)
  d <- function(i,z){         # z[1:2] contain long, lat
    dist <- rep(0,nrow(z))
    dist[i:nrow(z)] <- distHaversine(z[i:nrow(z),1:2],z[i,1:2])
    return(dist)
  }
  dm <- do.call(cbind,lapply(1:nrow(df),d,df))
  return(as.dist(dm))
}

dist.matrix <- geo.dist(NB.stat[,c("Lon.50.","Lat.50.")])
NB.clust <- hclust(dist.matrix)
plot(NB.clust)
NB.stat$cluster  <- cutree(NB.clust, k = 5)

# plot stopover and nonbreeding nodes 

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-90, -30),ylim = c(-15, 20))+
  geom_errorbar(data = NB.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), color = "red", width=1, alpha = 0.2) + 
  geom_point(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = site_type)) +
  geom_path(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5,
            arrow = arrow(end = "last", type = "open", length = unit(0.10, "inches"))) 

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-90, -30),ylim = c(-15, 20)) +
  geom_point(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = cluster)) +
  geom_path(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5,
            arrow = arrow(end = "last", type = "open", length = unit(0.10, "inches"))) +
   scale_colour_gradientn(colours=rainbow(5))

# Generate the network from our location data and clusters #####################

# Add a column with inter-cluster connections to our location dataset
NB.stat  <- NB.stat %>% mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# Create a dataframe with the edge info 
NB.edge.df <- NB.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                                Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                                Lat.2.5., Lat.97.5., EndTime,
                                                sitenum, duration, period,study.site,
                                                Range_region, NB_count, period,
                                                site_type) %>%
  filter(!is.na(next.cluster))

# calculate node locations as the median locations of the points included within 
# Note, we use the 50% quantile of the posterior distribution 
NB.node.lon <- NB.stat %>% group_by(cluster) %>% 
  summarize(node.lon = median(Lon.50.))

NB.node.lat <- NB.stat %>% group_by(cluster) %>% 
  summarize(node.lat = median(Lat.50.))

# We also assess which use was the more common for a node (breeding, nonbreeding, or stopover)
# The most common use will be the one assigned to the node
NB.node.type <- NB.stat %>% group_by(cluster, site_type) %>%
  summarize(use = n()) %>% filter(use == max(use)) %>% mutate(site_type_num = case_when(
    site_type == "Nonbreeding" ~ 1,
    site_type == "Stopover" ~ 2
  ))

# Create a layout with node locations 
meta <- data.frame("vertex" = seq(min(NB.stat$cluster), max(NB.stat$cluster)), 
                   "Lon.50." = NB.node.lon$node.lon,
                   "Lat.50." = NB.node.lat$node.lat,
                   "node.type.num" = NB.node.type$site_type_num)

location <- as.matrix(meta[, c("Lon.50.", "Lat.50.")])

# Create the spring network
NB.graph <- graph_from_data_frame(NB.edge.df, directed = T, vertices = meta)

# Colour palette for site type
type.palette <- rainbow(3)  

# plot the fall network over North and South America
plot(NB.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = 1, edge.arrow.size = 0.5, edge.arrow.width = 0.5,  
     layout = location, rescale = F, asp = 0, xlim = c(-90, -30),
     ylim = c(-15, 20), vertex.color = type.palette[meta$node.type.num])
plot(wrld_simpl, add = T)

# Add edge weights to the spring network to show the flow of individuals #########

# We start with edges based on the number of individuals moving between nodes

# list of connections and the number of times they occur
NB.con <- NB.edge.df %>% group_by(cluster, next.cluster) %>% 
  summarize(weight = n()) 

# Create a new spring network
NB.graph.weighed <- graph_from_data_frame(NB.con, directed = T, vertices = meta)
is_weighted(NB.graph.weighed)

plot(wrld_simpl, xlim = c(-90, -30), ylim = c(-15, 20))
plot(NB.graph.weighed, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = NB.con$weight/1.5, edge.arrow.size = 0.5, edge.arrow.width = 1.2,  
     layout = location, rescale = F, asp = 0, edge.curved = rep(c(-0.3, 0.3), nrow(NB.con)),
     vertex.color = type.palette[meta$node.type.num], add = T)

################################################################################
# Assess and propagate the relative abundance of blackpoll warblers from eBird
################################################################################

# create and export geolocator deployment sites
geo.breed <- geo.all %>% filter(site_type == "Breeding" &
                                  period == "Post-breeding migration") %>%
  dplyr::select(geo_id, Lon.50., Lat.50., sitenum, site_type, period, study.site)
breed.sites <- st_as_sf(geo.breed, coords = c("Lon.50.", "Lat.50."))

# set crs
st_crs(breed.sites) <- st_crs(wrld_simpl)

# export breeding sites
#st_write(breed.sites, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/bpw_breeding_sites.shp")

# import breeding sites and join with breeding site data 
abundance.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/bpw_abundance_regions.shp")
abundance.regions <- st_join(abundance.regions, breed.sites)

# import breeding season abundance data
bpw.fall.ab <- load_raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports/2021/bkpwar",
                           product = "abundance",
                           period = "seasonal",
                           resolution = "lr")

bpw.fall.ab <- terra::project(bpw.fall.ab, crs(abundance.regions))

# extract the abundance for each region,
ab.extract <- terra::extract(bpw.fall.ab$breeding, abundance.regions, fun = sum, na.rm=TRUE)
ab.extract$ID <- abundance.regions$geo_id
ab.extract$breedregionname <- abundance.regions$breedregio

# Create a dataframe with the abundance per region
ab.per.region <- merge(as.data.frame(abundance.regions), ab.extract, by.x = "id", by.y = "ID") %>%
  dplyr::select(-geometry)


