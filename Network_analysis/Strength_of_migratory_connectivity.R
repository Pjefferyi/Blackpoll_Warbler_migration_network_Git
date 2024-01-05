# Calculate the strength of migratory connectivity for the blackpoll warbler

# migconnectivity package 
#install.packages("devtools")
# usethis::create_github_token()
# usethis::edit_r_environ()
#devtools::install_github("SMBC-NZP/MigConnectivity")
library(MigConnectivity)

# load required packages 
library(tidyverse)
library(lubridate)
library(anytime)
library(ggplot2)
library(geosphere)
library(terra)
library(sf)
library(maptools)
library(ebirdst)
library(geosphere)


# Data preparation for MC between Breedign site and first nonbreeding site ######

# run the network preparation script
#source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Network_construction.R")

# projection for the analysis: Azimuthal Equidistant projection 
proj <- '+proj=aeqd +lat_0=0 +lon_0=-74'
#map <- spTransform(wrld_simpl, CRS(proj))
#plot(map)

# fall network breeding points
fall.br <- fall.breed %>% arrange(geo_id)
fall.br.sf <- st_as_sf(fall.br, coords = c("Lon.50.", "Lat.50."), crs = crs(wrld_simpl))
fall.br.sf <- st_transform(fall.br.sf, CRS(proj)) 

# fall network nonbreeding points
fall.stat <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stationary.locations.csv")
fall.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.graph.edge.list.txt", directed = TRUE)
E(fall.graph)$weight <- fall.con.ab$weight

undirected.fall.graph <- as.undirected(fall.graph, mode = "collapse",
                                       edge.attr.comb = list(weight = "sum", edge.type = "ignore"))

fall.label.prop <- concensusCluster(graph = undirected.fall.graph, thresh = 0.5, algiter = 1000)
V(fall.graph)$label.prop.comm <- fall.label.prop$`community structure`$membership

# merge fall community detection data with position data 
comm.data <- data.frame(cluster = V(fall.graph), community =  V(fall.graph)$label.prop.comm)
fall.stat <- merge(fall.stat, comm.data, by = "cluster")

fall.nbr <- fall.stat %>% group_by(geo_id) %>% 
  filter(period == "Non-breeding period", sitenum == max(sitenum), !is.na(StartTime), !is.na(EndTime)) %>%
  arrange(geo_id) 

fall.nbr.sf <- st_as_sf(fall.nbr, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld_simpl))
fall.nbr.sf <- st_transform(fall.nbr.sf, st_crs(proj)) 
#st_write(fall.nbr.sf, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data/bpw_fall_nbr_sitesV1.shp", append = F)

# Origin sites  (polygons)
br.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/bpw_abundance_regions_adjusted.shp") %>%
  st_transform(CRS(proj)) %>%
  st_cast("MULTIPOLYGON")

# Target sites(polygons)

# Uncomment to import polygons prepared in QGIS
# fall.nbr.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data/NonbreedingregionsV3.shp") %>%
#   st_transform(CRS(proj))%>%
#   st_cast("MULTIPOLYGON")

# creating target sites in R based on the nonbreeding nodes
geometry <- st_sfc(lapply(seq(1, length(unique(fall.nbr.sf$cluster))), function(x) st_geometrycollection()), crs = st_crs(proj))#st_sf(crs = st_crs(proj))
fall.nbr.regions <- st_sf(id = seq(1, length(unique(fall.nbr.sf$cluster))), geometry = geometry)

iter = 1

for (i in unique(fall.nbr.sf$cluster)){
  pt.sb <- fall.nbr.sf[fall.nbr.sf$cluster== i,]
  
  poly <- st_convex_hull(st_union(pt.sb))
  poly <- st_buffer(poly, dist = 80000)
  poly <- st_transform(poly, st_crs(proj))
  poly <- st_difference(poly, fall.nbr.regions)[1]
  #st_combine(c(fall.nbr.regions, poly))
  fall.nbr.regions$geometry[iter] <- poly
  fall.nbr.regions$cluster[iter] <- i
  iter = iter +1
}

# geoBIas (in meters)
# We will use the uncertainty for the threshold location estimates because the ouput of the 
# GroupThresholdModel provides only one location for stationary periods
Thresh.loc.data <- findThresLocData()
Thresh.loc.data$Twilight <- anytime(Thresh.loc.data$Twilight)

# We calculate spatial bias as the mean distance from the geolocator deployment site while the bird was known to be at that location
geo.bias.dists <- Thresh.loc.data %>% group_by(geo_id) %>%
  filter(Twilight > min(Twilight) + days(1) & Twilight < min(Twilight) + days(15)) %>%
  mutate(lat.geo.bias = distHaversine(cbind(deploy.longitude,lat), cbind(deploy.longitude, deploy.latitude))) %>%
  mutate(lon.geo.bias = distHaversine(cbind(lon,deploy.latitude), cbind(deploy.longitude, deploy.latitude)))

# bias in longitude and latitude estimates, in meters 
geo.bias <- c(lon.bias = mean(geo.bias.dists$lon.geo.bias),
              lat.bias = mean(geo.bias.dists$lat.geo.bias))

# Geolocator variance covariance at origin sites 
#geo.vcov <- cov(fall.nbr[, c("Lon.50.","Lat.50.")])
geo.vcov <- cov(st_coordinates(fall.nbr.sf))

#number of bootstrap samples 
nSamples <- 100

## Verify that the breeding and first nonbreeding sites and points makes sense

# Nonbreeding
plot(spTransform(wrld_simpl, CRS(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/6, 19659061/6),
     main = "Blackpoll warbler nonbreeding regions")
plot(as_Spatial(fall.nbr.regions), col = "yellow", add = T)
plot(as_Spatial(fall.nbr.sf), col = "blue", add  = T)

# Nonbreeding regions ggplot
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "lightgray")+
  geom_sf(data = fall.nbr.regions, aes(fill = as.factor(community)), col = "black", alpha = 0.5)+
  scale_fill_discrete(name = "Nonbreeding regions") +
  geom_sf(data = fall.nbr.sf, aes(col = "Individual nonbreeding locations (1 bird each)"), shape = 4, cex = 3)+
  scale_colour_manual(values = c("Individual nonbreeding locations (1 bird each)" = "blue"), name = "") +
  coord_sf(xlim = c(-90, -30),ylim = c(-15, 20))+
  theme_bw()+
  theme()+
  ggtitle("Blackpoll warbler nonbreeding regions")
  
# Breeding
plot(spTransform(wrld_simpl, CRS(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/8, 19659061/2))
plot(as_Spatial(br.regions), add = T)
plot(as_Spatial(fall.br.sf), add  = T)

# Breeding regions ggplot
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "lightgray")+
  geom_sf(data = st_transform(br.regions), aes(fill = region), col = "black", alpha = 0.5)+
  scale_fill_discrete(name = "Breeding regions") +
  geom_sf(data = st_transform(fall.br.sf), aes(col = "Geolocator deployment sites"), shape = 4, cex = 3)+
  scale_colour_manual(values = c("Geolocator deployment sites" = "blue"), name = "") +
  coord_sf(xlim = c(-170, -30),ylim = c(40, 90))+
  theme_bw()+
  theme()+
  ggtitle("Blackpoll warbler breeding regions")

# Estimation and resampling of uncertainty for transition probabilities (psi) ################################################################################
GL_psi_nbr1 <- estTransition(isGL=TRUE,
                        geoBias = geo.bias,
                        geoVCov = geo.vcov,
                        targetSites = fall.nbr.regions,
                        originSites = br.regions,
                        originPoints = fall.br.sf,
                        targetPoints = fall.nbr.sf,
                        verbose = 2,
                        nSamples = nSamples,
                        resampleProjection = CRS(proj),
                        maxTries = 1000)

# Save the output of estTransition
save(GL_psi_nbr1, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/estTransition_ouput_nbr1.R")
# load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/estTransition_ouput.R")

# Retrieve the relative abundance in each breeding site/region for the fall network ################################################################################

# import breeding season abundance data
bpw.breed.ab <- load_raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports/2021/bkpwar",
                           product = "abundance",
                           period = "seasonal",
                           resolution = "lr") %>%
 terra::project(crs(proj))

# # Plot of breeding regions over relative abundance surface
# plot(bpw.breed.ab$breeding, , xlim = c(-19312101/3, 19889532/6), ylim = c(-19793221/8, 19659061/2))
# plot(spTransform(wrld_simpl, CRS(proj)), add = T)
# plot(as_Spatial(br.regions), add = T)
# plot(as_Spatial(fall.br.sf), add  = T)

# extract the relative abundance for each site 
br.site.ab <- terra::extract(bpw.breed.ab$breeding, br.regions, fun = sum,
                             na.rm=TRUE, ID  = T) %>%
  mutate(breeding = breeding/sum(breeding),
         region = br.regions$region)
  
# Compute distance matrixes ####################################################
fall.br.centroids <- centroid(as_Spatial(br.regions))
fall.br.dists <- distFromPos(fall.br.centroids, "plane")

fall.nbr.centroids <- centroid(as_Spatial(fall.nbr.regions))
fall.nbr.dists <- distFromPos(fall.nbr.centroids, "plane")

# Estimate the strength of migratory connectivity metric #######################

MC_metric_nbr1 <- estStrength(targetDist = fall.nbr.dists , # targetSites distance matrix
                         originDist = fall.br.dists , # originSites distance matrix
                         targetSites = fall.nbr.sf, # Non-breeding target sites
                         originSites = fall.br.sf, # Breeding origin sites
                         psi = GL_psi_nbr1,
                         originRelAbund = br.site.ab$breeding,
                         nSamples = 1000,
                         sampleSize = nrow(fall.br.sf))

# Save output
 save(MC_metric_nbr1, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/MC_metric_nbr1.R")
# load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/MC_metric_nbr1.R")

# Estimate the Mantel correlation ##############################################

rM_metric_nbr1 <- estMantel(isGL= T,#Logical vector: light-level GL(T)/GPS(F)
                                 geoBias = geo.bias, # Geolocator location bias
                                 geoVCov = geo.vcov, # Location covariance matrix
                                 targetSites = fall.nbr.regions, # Non-breeding target sites
                                 originPoints = fall.br.sf, # Capture Locations
                                 targetPoints = fall.nbr.sf, # Target locations
                                 verbose = 1,   # output options
                                 nBoot = 100, # This is set low for example
                                 maxTries = 1000,
                                 resampleProjection = CRS(proj))

#Save output 
save(rM_metric_nbr1, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/rM_metric_nbr1.R")

# Data preparation for MC between Breeding site and second nonbreeding site ######

# fall network breeding points
spring.br <- spring.breed %>% arrange(geo_id)
spring.br.sf <- st_as_sf(spring.br, coords = c("Lon.50.", "Lat.50."), crs = crs(wrld_simpl))
spring.br.sf <- st_transform(spring.br.sf, CRS(proj)) 

# spring network nonbreeding points
spring.stat <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/spring.stationary.locations.csv")
spring.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/spring.graph.edge.list.txt", directed = TRUE)
E(spring.graph)$weight <- spring.con.ab$weight

spring.infomap <- cluster_infomap(graph = spring.graph)
V(spring.graph)$infomap.comm <- spring.infomap$membership

# merge spring community detection data with position data 
comm.data <- data.frame(cluster = V(spring.graph), community =  V(spring.graph)$infomap.comm)
spring.stat <- merge(spring.stat, comm.data, by = "cluster")

spring.nbr <- spring.stat %>% group_by(geo_id) %>% 
  filter(period == "Non-breeding period") %>% filter(sitenum == min(sitenum))

spring.nbr.sf <- st_as_sf(spring.nbr, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld_simpl))
spring.nbr.sf <- st_transform(spring.nbr.sf, st_crs(proj)) 
#st_write(spring.nbr.sf, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data/bpw_spring_nbr_sitesV1.shp", append = F)

# Origin sites  (polygons)
spring.br.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/bpw_abundance_regions_adjusted.shp") %>%
  st_transform(CRS(proj)) %>%
  st_cast("MULTIPOLYGON")

# Target sites(polygons)

# Uncomment to import polygons prepared in QGIS
# spring.nbr.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data/NonbreedingregionsV3.shp") %>%
#   st_transform(CRS(proj))%>%
#   st_cast("MULTIPOLYGON")

# creating target sites in R based on community clustering analysis 
geometry <- st_sfc(lapply(seq(1:max(spring.nbr.sf$community)), function(x) st_geometrycollection()), crs = st_crs(proj))#st_sf(crs = st_crs(proj))
spring.nbr.regions <- st_sf(id = seq(1:max(spring.nbr.sf$community)), geometry = geometry)

for (i in unique(spring.nbr.sf$community)){
  pt.sb <- spring.nbr.sf[spring.nbr.sf$community== i,]
  
  poly <- st_convex_hull(st_union(pt.sb))
  poly <- st_buffer(poly, dist = 80000)
  poly <- st_transform(poly, st_crs(proj))
  poly <- st_difference(poly, spring.nbr.regions)[1]
  #st_combine(c(spring.nbr.regions, poly))
  spring.nbr.regions$geometry[i] <- poly
  spring.nbr.regions$community[i] <- i
}

## Verify that the second nonbreeding sites and points makes sense

# Nonbreeding
plot(spTransform(wrld_simpl, CRS(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/6, 19659061/6),
     main = "Blackpoll warbler nonbreeding regions")
plot(as_Spatial(spring.nbr.regions), col = "yellow", add = T)
plot(as_Spatial(spring.nbr.sf), col = "blue", add  = T)

# Nonbreeding regions ggplot
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "lightgray")+
  geom_sf(data = spring.nbr.regions, aes(fill = as.factor(community)), col = "black", alpha = 0.5)+
  scale_fill_discrete(name = "Nonbreeding regions") +
  geom_sf(data = spring.nbr.sf, aes(col = "Individual nonbreeding locations (1 bird each)"), shape = 4, cex = 3)+
  scale_colour_manual(values = c("Individual nonbreeding locations (1 bird each)" = "blue"), name = "") +
  coord_sf(xlim = c(-90, -30),ylim = c(-15, 20))+
  theme_bw()+
  theme()+
  ggtitle("Blackpoll warbler nonbreeding regions")

# Breeding
plot(spTransform(wrld_simpl, CRS(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/8, 19659061/2))
plot(as_Spatial(br.regions), add = T)
plot(as_Spatial(spring.br.sf), add  = T)

# Breeding regions ggplot
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "lightgray")+
  geom_sf(data = st_transform(br.regions), aes(fill = region), col = "black", alpha = 0.5)+
  scale_fill_discrete(name = "Breeding regions") +
  geom_sf(data = st_transform(spring.br.sf), aes(col = "Geolocator deployment sites"), shape = 4, cex = 3)+
  scale_colour_manual(values = c("Geolocator deployment sites" = "blue"), name = "") +
  coord_sf(xlim = c(-170, -30),ylim = c(40, 90))+
  theme_bw()+
  theme()+
  ggtitle("Blackpoll warbler breeding regions")

# Estimation and resampling of uncertainty for transition probabilities (psi) ################################################################################
GL_psi_nbr2 <- estTransition(isGL=TRUE,
                        geoBias = geo.bias,
                        geoVCov = geo.vcov,
                        targetSites = spring.nbr.regions,
                        originSites = spring.br.regions,
                        originPoints = spring.br.sf,
                        targetPoints = spring.nbr.sf,
                        verbose = 2,
                        nSamples = nSamples,
                        resampleProjection = CRS(proj),
                        maxTries = 1000)

# Save the output of estTransition
save(GL_psi_nbr2, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/estTransition_ouput_nbr2.R")
#load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/estTransition_ouput_nbr2.R")

# Retrieve the relative abundance in each breeding site/region for the fall network ################################################################################

# # Plot of breeding regions over relative abundance surface
# plot(bpw.breed.ab$breeding, , xlim = c(-19312101/3, 19889532/6), ylim = c(-19793221/8, 19659061/2))
# plot(spTransform(wrld_simpl, CRS(proj)), add = T)
# plot(as_Spatial(br.regions), add = T)
# plot(as_Spatial(fall.br.sf), add  = T)

# extract the relative abundance for each site 
br.site.ab <- terra::extract(bpw.breed.ab$breeding, br.regions, fun = sum,
                                  na.rm=TRUE, ID  = T) %>%
  mutate(breeding = breeding/sum(breeding),
         region = br.regions$region)

# Compute distance matrixes ####################################################
spring.br.centroids <- centroid(as_Spatial(br.regions))
spring.br.dists <- distFromPos(spring.br.centroids, "plane")

spring.nbr.centroids <- centroid(as_Spatial(spring.nbr.regions))
spring.nbr.dists <- distFromPos(spring.nbr.centroids, "plane")

# Estimate the strength of migratory connectivity metric #######################
MC_metric_nbr2 <- estStrength(targetDist = spring.nbr.dists , # targetSites distance matrix
                              originDist = spring.br.dists , # originSites distance matrix
                              targetSites = spring.nbr.sf, # Non-breeding target sites
                              originSites = spring.br.sf, # Breeding origin sites
                              psi = GL_psi_nbr2,
                              originRelAbund = br.site.ab$breeding,
                              nSamples = 1000,
                              sampleSize = nrow(spring.br.sf))

#Save output
save(MC_metric_nbr2, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/MC_metric_nbr2.R")

# Estimate the Mantel correlation ##############################################

rM_metric_nbr2 <- estMantel(isGL= T,#Logical vector: light-level GL(T)/GPS(F)
                            geoBias = geo.bias, # Geolocator location bias
                            geoVCov = geo.vcov, # Location covariance matrix
                            targetSites = spring.nbr.regions, # Non-breeding target sites
                            originPoints = spring.br.sf, # Capture Locations
                            targetPoints = spring.nbr.sf, # Target locations
                            verbose = 1,   # output options
                            nBoot = 100, # This is set low for example
                            maxTries = 1000,
                            resampleProjection = CRS(proj))

#Save output 
save(rM_metric_nbr2, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/rM_metric_nbr2.R")

