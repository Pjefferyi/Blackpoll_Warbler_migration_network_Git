# Calculate the strength of migratory connectivity for the blackpoll warbler

# install the migconnectivity package 
# install.packages("devtools")
# usethis::create_github_token()
# usethis::edit_r_environ()
# devtools::install_github("SMBC-NZP/MigConnectivity")
library(MigConnectivity)

# load other required packages 
library(tidyverse)
library(lubridate)
library(anytime)
library(ggplot2)
library(geosphere)
library(terra)
library(sf)
library(ebirdst)
library(geosphere)
library(rnaturalearth)
library(igraph)

# Load helper functions 
source("Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Load country boundary polygons 
wrld <- as_Spatial(ne_countries(type = "countries", scale = "large"))

# First assessment of migratory connectivity ----

# Includes all individuals available irrespective of whether they were tracked during both the spring and fall

## Data preparation to estimate MC between the breeding site and first nonbreeding site ----

## Load required data generated when building the network
fall.breed <- read.csv("Network_construction/Fall.abundance.per.bird.csv")

# projection for the analysis: Azimuthal Equidistant projection 
proj <- '+proj=aeqd +lat_0=0 +lon_0=-74'
#map <- spTransform(wrld, CRS(proj))
#plot(map)

# fall network breeding points
fall.br <- fall.breed %>% arrange(geo_id)
fall.br.sf <- st_as_sf(fall.br, coords = c("Lon.50.", "Lat.50."), crs = crs(wrld))
fall.br.sf <- st_transform(fall.br.sf, crs(proj)) 

# fall network nonbreeding points
fall.stat <- read.csv("Network_construction/Fall.stationary.locations.csv")
fall.con.ab <- read.csv("Network_construction/Fall.edge.weights.csv")

fall.graph <- read_graph("Network_construction/Fall.graph.edge.list.txt", directed = TRUE)
E(fall.graph)$weight <- fall.con.ab$weight

# filter to get nonbreeding sites 
fall.nbr <- fall.stat %>% group_by(geo_id) %>% 
  filter(period == "Non-breeding period", sitenum == max(sitenum), !is.na(StartTime), !is.na(EndTime)) %>%
  arrange(geo_id)

fall.nbr.sf <- st_as_sf(fall.nbr, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld), remove = F)
fall.nbr.sf <- st_transform(fall.nbr.sf, st_crs(proj)) 

# Origin sites(abundance region polygons)
br.regions <- read_sf("Analysis_input_data/Breeding_regions/bpw_abundance_regions.shp") %>%
  st_transform(crs(proj)) %>%
  st_cast("MULTIPOLYGON")

# creating target sites in R based on the nonbreeding nodes

# First check whether there will be any isolated nonbreeding sites and add them to the nearest group
grp.sizes <- fall.nbr.sf %>% group_by(cluster) %>% 
  summarize(count = n(), med.lon = median(Lon.50.), med.lat = median(Lat.50.)) %>%
  mutate(cluster.change = NA)

fall.nbr.sf.temp <- fall.nbr.sf

# If there are any small groups of isolated sites, we merge them with the nearest polygons
for (i in seq(1, nrow(grp.sizes))){
  
  dist.dtfi <- grp.sizes
  
  dist.dtfi$lon.sites <- grp.sizes$med.lon[i]
  dist.dtfi$lat.sites <- grp.sizes$med.lat[i]
  
  dist.dtfi <- dist.dtfi %>% rowwise() %>% mutate(dist = distHaversine(c(lon.sites,lat.sites), c(med.lon, med.lat))) %>%
    filter(dist != 0)
  
  if (grp.sizes$count[i] < 4){
    
    min.index <- which(dist.dtfi$dist == min(dist.dtfi$dist))
    
    fall.nbr.sf.temp[fall.nbr.sf$cluster == grp.sizes$cluster[i],]$cluster <- dist.dtfi$cluster[min.index]
  }
  
}

fall.nbr.sf$cluster <- fall.nbr.sf.temp$cluster

# build the destination polyongs using buffered convex shells  
geometry <- st_sfc(lapply(seq(1, length(unique(fall.nbr.sf$cluster))), function(x) st_geometrycollection()), crs = st_crs(proj))#st_sf(crs = st_crs(proj))
fall.nbr.regions <- st_sf(id = seq(1, length(unique(fall.nbr.sf$cluster))), geometry = geometry)

iter = 1

for (i in sort(unique(fall.nbr.sf$cluster))){
  pt.sb <- fall.nbr.sf[fall.nbr.sf$cluster== i,]
  
  poly <- st_convex_hull(st_union(pt.sb))
  poly <- st_buffer(poly, dist = 80000)
  poly <- st_transform(poly, st_crs(proj))
  #st_combine(c(fall.nbr.regions, poly))
  fall.nbr.regions$geometry[iter] <- poly
  fall.nbr.regions$cluster[iter] <- i
  iter = iter +1
}

# we remove overlaps between the polygons
fall.nbr.regions <- st_difference(fall.nbr.regions)

# Save the first nonbreeding sites shapefile and the associated nonbreeding regions 
write_sf(fall.nbr.regions, "Analysis_output_data/Migratory_connectivity_regions/fall.nbr.regions.shp")
write.csv(data.frame(fall.nbr.sf), "Analysis_output_data/Migratory_connectivity_regions/fall.nbr.sf.csv")

# Uncertainty and error estimation 
# We will use daily location estimates from the grouped threshold model ran without stationary periods 
locs <- findThresModData() 
locs <- locs %>% group_by(geo_id) %>% filter(Time1 > IH.calib.start & Time1  < IH.calib.end)

# Modelled locations, projected
Thresh.mod.data.loc <- locs %>% st_as_sf(coords = c("Lon.mean",
                                                              "Lat.mean"),
                                                   crs = 4326) %>%
  st_transform(proj)

# origin locations, projected
Thresh.mod.data.or <- locs %>% st_as_sf(coords = c("deploy.longitude",
                                                                 "deploy.latitude"),
                                                      crs = 4326) %>%
  st_transform(proj)

# We calculate spatial bias as the mean distance from the geolocator deployment site while the bird was known to be at that location
lon.errors <- rep(NA, nrow(fall.nbr.sf))
lat.errors <- rep(NA, nrow(fall.nbr.sf))

for (i in seq(1, length(fall.nbr.sf$geo_id))){
  
  mod.locs <- Thresh.mod.data.loc %>% filter(geo_id == fall.nbr.sf$geo_id[i])
  or.locs <-  Thresh.mod.data.or %>% filter(geo_id == fall.nbr.sf$geo_id[i])
  
  lon.error <- mean(st_coordinates(mod.locs)[,1] - st_coordinates(or.locs)[,1])
  lat.error <- mean(st_coordinates(mod.locs)[,2] - st_coordinates(or.locs)[,2])
  
  lon.errors[i] <- lon.error
  lat.errors[i] <- lat.error
}

# Geolocator error and variance/covariance at origin sites
mod <- lm(cbind(lon.errors, lat.errors) ~ 1)
geo.bias <- coef(mod)
geo.vcov <- vcov(mod)

# Number of bootstrap samples 
nSamples <- 100

## Verify that the breeding and first nonbreeding sites and points makes sense

# Nonbreeding
plot(spTransform(wrld, crs(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/6, 19659061/6),
     main = "Blackpoll warbler nonbreeding regions")
plot(as_Spatial(fall.nbr.regions), col = "yellow", add = T)
plot(as_Spatial(fall.nbr.sf), col = "blue", add  = T)

# Nonbreeding regions ggplot
ggplot(st_as_sf(wrld))+
  geom_sf(colour = "black", fill = "lightgray")+
  geom_sf(data = fall.nbr.regions, aes(fill = as.factor(cluster)), col = "black", alpha = 0.5)+
  scale_fill_discrete(name = "Nonbreeding regions") +
  #geom_sf(data = fall.nbr.sf, aes(col = "Individual nonbreeding locations (1 bird each)"), shape = 4, cex = 3)+
  geom_sf(data = fall.nbr.sf, shape = 4, cex = 3)+
  #scale_colour_manual(values = c("Individual nonbreeding locations (1 bird each)" = "blue"), name = "") +
  coord_sf(xlim = c(-100, -30),ylim = c(-15, 20))+
  theme_bw()+
  theme()+
  ggtitle("Blackpoll warbler nonbreeding regions")
  
# Breeding
plot(spTransform(wrld, crs(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/8, 19659061/2))
plot(as_Spatial(br.regions), add = T)
plot(as_Spatial(fall.br.sf), add  = T)

# Breeding regions ggplot
ggplot(st_as_sf(wrld))+
  geom_sf(colour = "black", fill = "lightgray")+
  geom_sf(data = st_transform(br.regions), aes(fill = region), col = "black", alpha = 0.5)+
  scale_fill_discrete(name = "Breeding regions") +
  geom_sf(data = st_transform(fall.br.sf), aes(col = "Geolocator deployment sites"), shape = 4, cex = 3)+
  scale_colour_manual(values = c("Geolocator deployment sites" = "blue"), name = "") +
  coord_sf(xlim = c(-170, -30),ylim = c(40, 90))+
  theme_bw()+
  theme()+
  ggtitle("Blackpoll warbler breeding regions")

### Estimation and resampling of uncertainty for transition probabilities (psi) ################################################################################
GL_psi_nbr1 <- estTransition(isGL=TRUE,
                        geoBias = geo.bias,
                        geoVCov = geo.vcov,
                        targetSites = fall.nbr.regions,
                        originSites = br.regions,
                        originPoints = fall.br.sf,
                        targetPoints = fall.nbr.sf,
                        verbose = 2,
                        nSamples = nSamples,
                        resampleProjection = crs(proj),
                        maxTries = 1000)

# Save the output of estTransition
save(GL_psi_nbr1, file = "Analysis_output_data/MC_estimation/estTransition_ouput_nbr1.R")
# load("Analysis_output_data/MC_estimation/estTransition_ouput_nbr1.R")

### Retrieve the relative abundance in each breeding site/region for the fall network ################################################################################

# import breeding season abundance data
bpw.breed.ab <- load_raster(path = "Analysis_input_data/eBird_imports",
                            species = "bkpwar",
                            product = "abundance",
                            period = "seasonal",
                            resolution = "9km") %>%
 terra::project(crs(proj))

# # Plot of breeding regions over relative abundance surface
# plot(bpw.breed.ab$breeding, , xlim = c(-19312101/3, 19889532/6), ylim = c(-19793221/8, 19659061/2))
# plot(spTransform(wrld, CRS(proj)), add = T)
# plot(as_Spatial(br.regions), add = T)
# plot(as_Spatial(fall.br.sf), add  = T)

# extract the relative abundance for each site 
br.site.ab <- terra::extract(bpw.breed.ab$breeding, br.regions, fun = sum,
                             na.rm=TRUE, ID  = T) %>%
  mutate(breeding = breeding/sum(breeding),
         region = br.regions$region)
  
### Compute distance matrixes ####################################################
fall.br.centroids <- centroid(as_Spatial(br.regions))
fall.br.dists <- distFromPos(fall.br.centroids, "plane")

fall.nbr.centroids <- centroid(as_Spatial(fall.nbr.regions))
fall.nbr.dists <- distFromPos(fall.nbr.centroids, "plane")

### Estimate the strength of migratory connectivity metric #######################
MC_metric_nbr1 <- estStrength(targetDist = fall.nbr.dists , # targetSites distance matrix
                         originDist = fall.br.dists , # originSites distance matrix
                         targetSites = fall.nbr.sf, # Non-breeding target sites
                         originSites = fall.br.sf, # Breeding origin sites
                         psi = GL_psi_nbr1,
                         originRelAbund = br.site.ab$breeding,
                         nSamples = 1000,
                         sampleSize = nrow(fall.br.sf))

# Save output
save(MC_metric_nbr1, file = "Analysis_output_data/MC_estimation/MC_metric_nbr1.R")

### Estimate the Mantel correlation ##############################################

rM_metric_nbr1 <- estMantel(isGL= T,#Logical vector: light-level GL(T)/GPS(F)
                                 geoBias = geo.bias, # Geolocator location bias
                                 geoVCov = geo.vcov, # Location covariance matrix
                                 targetSites = fall.nbr.regions, # Non-breeding target sites
                                 originPoints = fall.br.sf, # Capture Locations
                                 targetPoints = fall.nbr.sf, # Target locations
                                 verbose = 1,   # output options
                                 nBoot = 100, # This is set low for example
                                 maxTries = 1000,
                                 resampleProjection = crs(proj))

#Save output 
save(rM_metric_nbr1, file = "Analysis_output_data/MC_estimation/rm_metric_nbr1.R")


## Data preparation for MC between Breeding site and second nonbreeding site ----

## Load required data generated when building the network
spring.breed <- read.csv("Network_construction/spring.abundance.per.bird.csv")

# spring network breeding points
spring.br <- spring.breed %>% arrange(geo_id)
spring.br.sf <- st_as_sf(spring.br, coords = c("Lon.50.", "Lat.50."), crs = crs(wrld))
spring.br.sf <- st_transform(spring.br.sf, crs(proj)) 

# spring network nonbreeding points
spring.stat <- read.csv("Network_construction/spring.stationary.locations.csv")
spring.graph <- read_graph("Network_construction/spring.graph.edge.list.txt", directed = TRUE)
E(spring.graph)$weight <- spring.con.ab$weight

# Filter to get the last nonbreeding site used before the spring migration 
spring.nbr <- spring.stat %>% group_by(geo_id) %>% 
  filter(period == "Non-breeding period") %>% filter(sitenum == min(sitenum))

spring.nbr.sf <- st_as_sf(spring.nbr, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld), remove = F)
spring.nbr.sf <- st_transform(spring.nbr.sf, st_crs(proj)) 

# Uncertainty and error estimation for the spring 
# We again use daily location estimates from the grouped threshold model ran without stationary periods 

# We calculate spatial bias as the mean distance from the geolocator deployment site while the bird was known to be at that location
lon.errors <- rep(NA, nrow(spring.nbr.sf))
lat.errors <- rep(NA, nrow(spring.nbr.sf))

for (i in seq(1, length(spring.nbr.sf$geo_id))){
  
  mod.locs <- Thresh.mod.data.loc %>% filter(geo_id == spring.nbr.sf$geo_id[i])
  or.locs <-  Thresh.mod.data.or %>% filter(geo_id == spring.nbr.sf$geo_id[i])
  
  lon.error <- mean(st_coordinates(mod.locs)[,1] - st_coordinates(or.locs)[,1])
  lat.error <- mean(st_coordinates(mod.locs)[,2] - st_coordinates(or.locs)[,2])
  
  lon.errors[i] <- lon.error
  lat.errors[i] <- lat.error
}

# Geolocator error and variance/covariance at origin sites
mod <- lm(cbind(lon.errors, lat.errors) ~ 1)
geo.bias <- coef(mod)
geo.vcov <- vcov(mod)

# Origin sites  (polygons)
spring.br.regions <- read_sf("Analysis_input_data/Breeding_regions/bpw_abundance_regions.shp") %>%
  st_transform(crs(proj)) %>%
  st_cast("MULTIPOLYGON")

# creating target sites in R based on the nonbreeding nodes

# First check whether there will be any isolated nonbreeding site and add them to the nearest group
grp.sizes <- spring.nbr.sf %>% group_by(cluster) %>% 
  summarize(count = n(), med.lon = median(Lon.50.), med.lat = median(Lat.50.))

spring.nbr.sf.temp <- spring.nbr.sf

# If there are any small groups of isolated sites, we merge them with the nearest polygons
for (i in seq(1, nrow(grp.sizes))){
  
  dist.dtfi <- grp.sizes
  
  dist.dtfi$lon.sites <- grp.sizes$med.lon[i]
  dist.dtfi$lat.sites <- grp.sizes$med.lat[i]
  
  dist.dtfi <- dist.dtfi %>% rowwise() %>% mutate(dist = distHaversine(c(lon.sites,lat.sites), c(med.lon, med.lat))) %>%
    filter(dist != 0)
  
  if (grp.sizes$count[i] < 4){
    
    min.index <- which(dist.dtfi$dist == min(dist.dtfi$dist))
    
    spring.nbr.sf.temp[spring.nbr.sf$cluster == grp.sizes$cluster[i],]$cluster <- dist.dtfi$cluster[min.index]
  }
  
}

spring.nbr.sf$cluster <- spring.nbr.sf.temp$cluster

geometry <- st_sfc(lapply(seq(1, length(unique(spring.nbr.sf$cluster))), function(x) st_geometrycollection()), crs = st_crs(proj))#st_sf(crs = st_crs(proj))
spring.nbr.regions <- st_sf(id = seq(1, length(unique(spring.nbr.sf$cluster))), geometry = geometry)

iter = 1

for (i in unique(spring.nbr.sf$cluster)){
  pt.sb <- spring.nbr.sf[spring.nbr.sf$cluster== i,]
  
  poly <- st_convex_hull(st_union(pt.sb))
  poly <- st_buffer(poly, dist = 80000)
  poly <- st_transform(poly, st_crs(proj))
  #st_combine(c(spring.nbr.regions, poly))
  spring.nbr.regions$geometry[iter] <- poly
  spring.nbr.regions$cluster[iter] <- i
  iter = iter +1
}

# we remove overlaps between the polygons
spring.nbr.regions <- st_difference(spring.nbr.regions)

# Save the second nonbreeding sites shapefile and the associated nonbreeding regions 
write_sf(spring.nbr.regions, "Analysis_output_data/Migratory_connectivity_regions/spring.nbr.regions.shp")
write.csv(data.frame(spring.nbr.sf), "Analysis_output_data/Migratory_connectivity_regions/spring.nbr.sf.csv")

## Verify that the second nonbreeding sites and points makes sense

# Nonbreeding
plot(spTransform(wrld, crs(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/6, 19659061/6),
     main = "Blackpoll warbler nonbreeding regions")
plot(as_Spatial(spring.nbr.regions), col = "yellow", add = T)
plot(as_Spatial(spring.nbr.sf), col = "blue", add  = T)

# Nonbreeding regions ggplot
ggplot(st_as_sf(wrld))+
  geom_sf(colour = "black", fill = "lightgray")+
  geom_sf(data = spring.nbr.regions, aes(fill = as.factor(cluster)), col = "black", alpha = 0.5)+
  scale_fill_discrete(name = "Nonbreeding regions") +
  geom_sf(data = spring.nbr.sf, aes(col = "Individual nonbreeding locations (1 bird each)"), shape = 4, cex = 3)+
  scale_colour_manual(values = c("Individual nonbreeding locations (1 bird each)" = "blue"), name = "") +
  coord_sf(xlim = c(-90, -30),ylim = c(-15, 20))+
  theme_bw()+
  theme()+
  ggtitle("Blackpoll warbler nonbreeding regions")

# Breeding
plot(spTransform(wrld, crs(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/8, 19659061/2))
plot(as_Spatial(br.regions), add = T)
plot(as_Spatial(spring.br.sf), add  = T)

# Breeding regions ggplot
ggplot(st_as_sf(wrld))+
  geom_sf(colour = "black", fill = "lightgray")+
  geom_sf(data = st_transform(br.regions), aes(fill = region), col = "black", alpha = 0.5)+
  scale_fill_discrete(name = "Breeding regions") +
  geom_sf(data = st_transform(spring.br.sf), aes(col = "Geolocator deployment sites"), shape = 4, cex = 3)+
  scale_colour_manual(values = c("Geolocator deployment sites" = "blue"), name = "") +
  coord_sf(xlim = c(-170, -30),ylim = c(40, 90))+
  theme_bw()+
  theme()+
  ggtitle("Blackpoll warbler breeding regions")

# Generate an indicator to signal that WRMA04173 was captured in the nonbreeding range
capture.loc <- ifelse(spring.nbr.sf$geo_id != "WRMA04173", "origin", "target")

### Estimation and resampling of uncertainty for transition probabilities (psi) ################################################################################
GL_psi_nbr2 <- estTransition(isGL=TRUE,
                        geoBias = geo.bias,
                        geoVCov = geo.vcov,
                        targetSites = spring.nbr.regions,
                        originSites = spring.br.regions,
                        originPoints = spring.br.sf,
                        targetPoints = spring.nbr.sf,
                        captured = capture.loc,
                        verbose = 2,
                        nSamples = nSamples,
                        resampleProjection = crs(proj),
                        maxTries = 1000)

# Save the output of estTransition
save(GL_psi_nbr2, file = "Analysis_output_data/MC_estimation/estTransition_ouput_nbr2.R")
#load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/estTransition_ouput_nbr2.R")

### Retrieve the relative abundance in each breeding site/region for the fall network ################################################################################

# # Plot of breeding regions over relative abundance surface
# plot(bpw.breed.ab$breeding, , xlim = c(-19312101/3, 19889532/6), ylim = c(-19793221/8, 19659061/2))
# plot(spTransform(wrld, CRS(proj)), add = T)
# plot(as_Spatial(br.regions), add = T)
# plot(as_Spatial(fall.br.sf), add  = T)

# extract the relative abundance for each site 
br.site.ab <- terra::extract(bpw.breed.ab$breeding, br.regions, fun = sum,
                                  na.rm=TRUE, ID  = T) %>%
  mutate(breeding = breeding/sum(breeding),
         region = br.regions$region)

### Compute distance matrixes ####################################################
spring.br.centroids <- centroid(as_Spatial(br.regions))
spring.br.dists <- distFromPos(spring.br.centroids, "plane")

spring.nbr.centroids <- centroid(as_Spatial(spring.nbr.regions))
spring.nbr.dists <- distFromPos(spring.nbr.centroids, "plane")

### Estimate the strength of migratory connectivity metric #######################
MC_metric_nbr2 <- estStrength(targetDist = spring.nbr.dists , # targetSites distance matrix
                              originDist = spring.br.dists , # originSites distance matrix
                              targetSites = spring.nbr.sf, # Non-breeding target sites
                              originSites = spring.br.sf, # Breeding origin sites
                              psi = GL_psi_nbr2,
                              originRelAbund = br.site.ab$breeding,
                              nSamples = 1000,
                              sampleSize = nrow(spring.br.sf))

#Save output
save(MC_metric_nbr2, file = "Analysis_output_data/Migratory_connectivity_regions/MC_metric_nbr2.R")
#load("Analysis_output_data/Migratory_connectivity_regions/MC_metric_nbr2.R")

### Estimate the Mantel correlation ##############################################

rM_metric_nbr2 <- estMantel(isGL= T,#Logical vector: light-level GL(T)/GPS(F)
                            geoBias = geo.bias, # Geolocator location bias
                            geoVCov = geo.vcov, # Location covariance matrix
                            targetSites = spring.nbr.regions, # Non-breeding target sites
                            originPoints = spring.br.sf, # Capture Locations
                            targetPoints = spring.nbr.sf, # Target locations
                            verbose = 1,   # output options
                            nBoot = 100, # This is set low for example
                            maxTries = 1000,
                            resampleProjection = crs(proj))

# Save output 
save(rM_metric_nbr2, file = "Analysis_output_data/Migratory_connectivity_regions/rM_metric_nbr2.R")

# Second assessment of migratory connectivity ----

# Includes only individuals tracked in during both the spring and fall migration

## Estimate MC between the breeding site and first nonbreeding site ----

# This is to verify that the change of migratory connectivity is not due to a difference in our sample 
fall.nbr.sf.test <- fall.nbr.sf %>% filter(geo_id %in% spring.nbr.sf$geo_id)
fall.br.sf.test <- fall.br.sf %>% filter(geo_id %in% spring.br.sf$geo_id)

# creating target sites in R based on the nonbreeding nodes

# First check whether there will be any isolated nonbreeding site and add them to the nearest group
grp.sizes <- fall.nbr.sf.test %>% group_by(cluster) %>% 
  summarize(count = n(), med.lon = median(Lon.50.), med.lat = median(Lat.50.)) %>%
  mutate(cluster.change = NA)

fall.nbr.sf.temp <- fall.nbr.sf.test

for (i in seq(1, nrow(grp.sizes))){
  
  dist.dtfi <- grp.sizes
  
  dist.dtfi$lon.sites <- grp.sizes$med.lon[i]
  dist.dtfi$lat.sites <- grp.sizes$med.lat[i]
  
  dist.dtfi <- dist.dtfi %>% rowwise() %>% mutate(dist = distHaversine(c(lon.sites,lat.sites), c(med.lon, med.lat))) %>%
    filter(dist != 0)
  
  if (grp.sizes$count[i] < 4){
    
    min.index <- which(dist.dtfi$dist == min(dist.dtfi$dist))
    
    fall.nbr.sf.temp[fall.nbr.sf$cluster == grp.sizes$cluster[i],]$cluster <- dist.dtfi$cluster[min.index]
  }
  
}

fall.nbr.sf.test$cluster <- fall.nbr.sf.temp$cluster

geometry <- st_sfc(lapply(seq(1, length(unique(fall.nbr.sf$cluster))), function(x) st_geometrycollection()), crs = st_crs(proj))#st_sf(crs = st_crs(proj))
fall.nbr.regions <- st_sf(id = seq(1, length(unique(fall.nbr.sf$cluster))), geometry = geometry)

iter = 1

for (i in sort(unique(fall.nbr.sf.test$cluster))){
  pt.sb <- fall.nbr.sf.test[fall.nbr.sf.test$cluster== i,]
  
  poly <- st_convex_hull(st_union(pt.sb))
  poly <- st_buffer(poly, dist = 80000)
  poly <- st_transform(poly, st_crs(proj))
  #st_combine(c(fall.nbr.regions, poly))
  fall.nbr.regions$geometry[iter] <- poly
  fall.nbr.regions$cluster[iter] <- i
  iter = iter +1
}

# we remove overlaps between the polygons
fall.nbr.regions.test <- st_difference(fall.nbr.regions)

# Uncertainty and error estimation 
# We will use daily location estimates from the grouped threshold model ran without stationary periods 
locs <- findThresModData() 
locs <- locs %>% group_by(geo_id) %>% filter(Time1 > IH.calib.start & Time1  < IH.calib.end)

# modelled locations, projected
Thresh.mod.data.loc <- locs %>% st_as_sf(coords = c("Lon.mean",
                                                    "Lat.mean"),
                                         crs = 4326) %>%
  st_transform(proj)

# origin locations, projected
Thresh.mod.data.or <- locs %>% st_as_sf(coords = c("deploy.longitude",
                                                   "deploy.latitude"),
                                        crs = 4326) %>%
  st_transform(proj)

# We calculate spatial bias as the mean distance from the geolocator deployment site while the bird was known to be at that location
lon.errors <- rep(NA, nrow(fall.nbr.sf.test))
lat.errors <- rep(NA, nrow(fall.nbr.sf.test))

for (i in seq(1, length(fall.nbr.sf.test$geo_id))){
  
  mod.locs <- Thresh.mod.data.loc %>% filter(geo_id == fall.nbr.sf.test$geo_id[i])
  or.locs <-  Thresh.mod.data.or %>% filter(geo_id == fall.nbr.sf.test$geo_id[i])
  
  lon.error <- mean(st_coordinates(mod.locs)[,1] - st_coordinates(or.locs)[,1])
  lat.error <- mean(st_coordinates(mod.locs)[,2] - st_coordinates(or.locs)[,2])
  
  lon.errors[i] <- lon.error
  lat.errors[i] <- lat.error
}

# Geolocator error and variance/covariance at origin sites
mod <- lm(cbind(lon.errors, lat.errors) ~ 1)
geo.bias <- coef(mod)
geo.vcov <- vcov(mod)

#number of bootstrap samples 
nSamples <- 100

## Verify that the breeding and first nonbreeding sites and points makes sense

# Nonbreeding
plot(spTransform(wrld, crs(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/6, 19659061/6),
     main = "Blackpoll warbler nonbreeding regions")
plot(as_Spatial(fall.nbr.regions), col = "yellow", add = T)
plot(as_Spatial(fall.nbr.sf.test), col = "blue", add  = T)

### Estimation and resampling of uncertainty for transition probabilities (psi) ################################################################################
GL_psi_nbr3 <- estTransition(isGL=TRUE,
                             geoBias = geo.bias,
                             geoVCov = geo.vcov,
                             targetSites = fall.nbr.regions.test,
                             originSites = br.regions,
                             originPoints = fall.br.sf.test,
                             targetPoints = fall.nbr.sf.test,
                             verbose = 2,
                             nSamples = nSamples,
                             resampleProjection = crs(proj),
                             maxTries = 1000)

x <- as.data.frame(GL_psi_nbr1$psi$mean)
colnames(x) <- fall.nbr.regions$cluster
rownames(x) <- br.regions$region

# Save the output of estTransition
save(GL_psi_nbr3, file = "Analysis_output_data/MC_estimation/estTransition_ouput_nbr3.R")
# load("Analysis_output_data/MC_estimation/estTransition_ouput_nbr3.R")

### Retrieve the relative abundance in each breeding site/region for the fall network ################################################################################

# import breeding season abundance data
bpw.breed.ab <- load_raster(path = "Analysis_input_data/eBird_imports",
                            species = "bkpwar",
                            product = "abundance",
                            period = "seasonal",
                            resolution = "9km") %>%
  terra::project(crs(proj))

# # Plot of breeding regions over relative abundance surface
# plot(bpw.breed.ab$breeding, , xlim = c(-19312101/3, 19889532/6), ylim = c(-19793221/8, 19659061/2))
# plot(spTransform(wrld, CRS(proj)), add = T)
# plot(as_Spatial(br.regions), add = T)
# plot(as_Spatial(fall.br.sf), add  = T)

# extract the relative abundance for each site 
br.site.ab <- terra::extract(bpw.breed.ab$breeding, br.regions, fun = sum,
                             na.rm=TRUE, ID  = T) %>%
  mutate(breeding = breeding/sum(breeding),
         region = br.regions$region)

### Compute distance matrixes ####################################################
fall.br.centroids <- centroid(as_Spatial(br.regions))
fall.br.dists <- distFromPos(fall.br.centroids, "plane")

fall.nbr.centroids <- centroid(as_Spatial(fall.nbr.regions))
fall.nbr.dists <- distFromPos(fall.nbr.centroids, "plane")

### Estimate the strength of migratory connectivity metric #######################
MC_metric_nbr3 <- estStrength(targetDist = fall.nbr.dists , # targetSites distance matrix
                              originDist = fall.br.dists , # originSites distance matrix
                              targetSites = fall.nbr.sf, # Non-breeding target sites
                              originSites = fall.br.sf, # Breeding origin sites
                              psi = GL_psi_nbr1,
                              originRelAbund = br.site.ab$breeding,
                              nSamples = 1000,
                              sampleSize = nrow(fall.br.sf))

# Save output
save(MC_metric_nbr3, file = "Analysis_output_data/MC_estimation/MC_metric_nbr3.R")
# load("Analysis_output_data/MC_estimation/MC_metric_nbr3.R")

### Estimate the Mantel correlation ##############################################
rM_metric_nbr3 <- estMantel(isGL= T,#Logical vector: light-level GL(T)/GPS(F)
                            geoBias = geo.bias, # Geolocator location bias
                            geoVCov = geo.vcov, # Location covariance matrix
                            targetSites = fall.nbr.regions, # Non-breeding target sites
                            originPoints = fall.br.sf, # Capture Locations
                            targetPoints = fall.nbr.sf, # Target locations
                            verbose = 1,   # output options
                            nBoot = 100, # This is set low for example
                            maxTries = 1000,
                            resampleProjection = crs(proj))


## Estimate MC between the breeding site and second nonbreeding site ----

# spring network breeding points
spring.br.sf.test <- spring.br.sf %>% filter(geo_id %in% fall.nbr.sf$geo_id)

# Spring network nonbreeding points 
spring.nbr.sf.test <- spring.nbr.sf %>% filter(geo_id %in% fall.br.sf$geo_id)

# Uncertainty and error estimation for the spring 
# We again use daily location estimates from the grouped threshold model ran without stationary periods 

# We calculate spatial bias as the mean distance from the geolocator deployment site while the bird was known to be at that location
lon.errors <- rep(NA, nrow(spring.nbr.sf.test))
lat.errors <- rep(NA, nrow(spring.nbr.sf.test))

for (i in seq(1, length(spring.nbr.sf.test$geo_id))){
  
  mod.locs <- Thresh.mod.data.loc %>% filter(geo_id == spring.nbr.sf.test$geo_id[i])
  or.locs <-  Thresh.mod.data.or %>% filter(geo_id == spring.nbr.sf.test$geo_id[i])
  
  lon.error <- mean(st_coordinates(mod.locs)[,1] - st_coordinates(or.locs)[,1])
  lat.error <- mean(st_coordinates(mod.locs)[,2] - st_coordinates(or.locs)[,2])
  
  lon.errors[i] <- lon.error
  lat.errors[i] <- lat.error
}

# Geolocator error and variance/covariance at origin sites
mod <- lm(cbind(lon.errors, lat.errors) ~ 1)
geo.bias <- coef(mod)
geo.vcov <- vcov(mod)

# Origin sites  (polygons)
spring.br.regions <- read_sf("Analysis_input_data/Breeding_regions/bpw_abundance_regions.shp") %>%
  st_transform(crs(proj)) %>%
  st_cast("MULTIPOLYGON")

# creating target sites in R based on the nonbreeding nodes

# First check whether there will be any isolated nonbreeding site and add them to the nearest group
grp.sizes <- spring.nbr.sf.test %>% group_by(cluster) %>% 
  summarize(count = n(), med.lon = median(Lon.50.), med.lat = median(Lat.50.))

spring.nbr.sf.temp <- spring.nbr.sf.test

for (i in seq(1, nrow(grp.sizes))){
  
  dist.dtfi <- grp.sizes
  
  dist.dtfi$lon.sites <- grp.sizes$med.lon[i]
  dist.dtfi$lat.sites <- grp.sizes$med.lat[i]
  
  dist.dtfi <- dist.dtfi %>% rowwise() %>% mutate(dist = distHaversine(c(lon.sites,lat.sites), c(med.lon, med.lat))) %>%
    filter(dist != 0)
  
  if (grp.sizes$count[i] < 4){
    
    min.index <- which(dist.dtfi$dist == min(dist.dtfi$dist))
    
    spring.nbr.sf.temp[spring.nbr.sf.test$cluster == grp.sizes$cluster[i],]$cluster <- dist.dtfi$cluster[min.index]
  }
  
}

spring.nbr.sf.test$cluster <- spring.nbr.sf.temp$cluster

geometry <- st_sfc(lapply(seq(1, length(unique(spring.nbr.sf.test$cluster))), function(x) st_geometrycollection()), crs = st_crs(proj))#st_sf(crs = st_crs(proj))
spring.nbr.regions <- st_sf(id = seq(1, length(unique(spring.nbr.sf.test$cluster))), geometry = geometry)

iter = 1

for (i in unique(spring.nbr.sf$cluster)){
  pt.sb <- spring.nbr.sf.test[spring.nbr.sf.test$cluster== i,]
  
  poly <- st_convex_hull(st_union(pt.sb))
  poly <- st_buffer(poly, dist = 80000)
  poly <- st_transform(poly, st_crs(proj))
  #st_combine(c(spring.nbr.regions, poly))
  spring.nbr.regions$geometry[iter] <- poly
  spring.nbr.regions$cluster[iter] <- i
  iter = iter +1
}

# we remove overlaps between the polygons
spring.nbr.regions.test <- st_difference(spring.nbr.regions)

## Verify that the second nonbreeding sites and points makes sense

# Nonbreeding
plot(spTransform(wrld, crs(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/6, 19659061/6),
     main = "Blackpoll warbler nonbreeding regions")
plot(as_Spatial(spring.nbr.regions.test), col = "yellow", add = T)
plot(as_Spatial(spring.nbr.sf.test), col = "blue", add  = T)

### Estimation and resampling of uncertainty for transition probabilities (psi) ################################################################################
GL_psi_nbr4 <- estTransition(isGL=TRUE,
                             geoBias = geo.bias,
                             geoVCov = geo.vcov,
                             targetSites = spring.nbr.regions.test,
                             originSites = spring.br.regions,
                             originPoints = spring.br.sf.test,
                             targetPoints = spring.nbr.sf.test,
                             verbose = 2,
                             nSamples = nSamples,
                             resampleProjection = crs(proj),
                             maxTries = 1000)

x <- as.data.frame(GL_psi_nbr4$psi$mean)
colnames(x) <- spring.nbr.regions$cluster
rownames(x) <- spring.br.regions$region

# Save the output of estTransition
save(GL_psi_nbr4, file = "Analysis_output_data/MC_estimation/estTransition_ouput_nbr4.R")
#load("Analysis_output_data/MC_estimation/estTransition_ouput_nbr4.R")

### Retrieve the relative abundance in each breeding site/region for the fall network ################################################################################

# # Plot of breeding regions over relative abundance surface
# plot(bpw.breed.ab$breeding, , xlim = c(-19312101/3, 19889532/6), ylim = c(-19793221/8, 19659061/2))
# plot(spTransform(wrld, CRS(proj)), add = T)
# plot(as_Spatial(br.regions), add = T)
# plot(as_Spatial(fall.br.sf), add  = T)

# extract the relative abundance for each site 
br.site.ab <- terra::extract(bpw.breed.ab$breeding, br.regions, fun = sum,
                             na.rm=TRUE, ID  = T) %>%
  mutate(breeding = breeding/sum(breeding),
         region = br.regions$region)

### Compute distance matrixes ####################################################
spring.br.centroids <- centroid(as_Spatial(br.regions))
spring.br.dists <- distFromPos(spring.br.centroids, "plane")

spring.nbr.centroids <- centroid(as_Spatial(spring.nbr.regions))
spring.nbr.dists <- distFromPos(spring.nbr.centroids, "plane")

### Estimate the strength of migratory connectivity metric #######################
MC_metric_nbr4 <- estStrength(targetDist = spring.nbr.dists , # targetSites distance matrix
                              originDist = spring.br.dists , # originSites distance matrix
                              targetSites = spring.nbr.sf, # Non-breeding target sites
                              originSites = spring.br.sf, # Breeding origin sites
                              psi = GL_psi_nbr4,
                              originRelAbund = br.site.ab$breeding,
                              nSamples = 1000,
                              sampleSize = nrow(spring.br.sf))

#Save output
save(MC_metric_nbr4, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/MC_metric_nbr4.R")
#load("Analysis_output_data/MC_estimation/MC_metric_nbr4.R")

### Estimate the Mantel correlation ##############################################
rM_metric_nbr4 <- estMantel(isGL= T,#Logical vector: light-level GL(T)/GPS(F)
                            geoBias = geo.bias, # Geolocator location bias
                            geoVCov = geo.vcov, # Location covariance matrix
                            targetSites = spring.nbr.regions, # Non-breeding target sites
                            originPoints = spring.br.sf, # Capture Locations
                            targetPoints = spring.nbr.sf, # Target locations
                            verbose = 1,   # output options
                            nBoot = 100, # This is set low for example
                            maxTries = 1000,
                            resampleProjection = crs(proj))

#Save output 
save(rM_metric_nbr4, file = "Analysis_output_data/MC_estimation/rM_metric_nbr4.R")
#load("Analysis_output_data/MC_estimation/rM_metric_nbr4.R")
