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

################################################################################
# Data preparation 
################################################################################

# run the network preparation script
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Network_construction.R")

# projection for the analysis: Azimuthal Equidistant projection 
proj <- '+proj=aeqd +lat_0=0 +lon_0=-74'
#map <- spTransform(wrld_simpl, CRS(proj))
#plot(map)

# fall network breeding points
fall.br <- fall.breed %>% arrange(geo_id)
fall.br.sf <- st_as_sf(fall.br, coords = c("Lon.50.", "Lat.50."), crs = crs(wrld_simpl))
fall.br.sf <- st_transform(fall.br.sf, CRS(proj)) 

# fall network nonbreeding points
fall.nbr<- fall.stat %>% group_by(geo_id) %>% 
  filter(sitenum == max(sitenum)) %>%
  arrange(geo_id) 

fall.nbr.sf <- st_as_sf(fall.nbr, coords = c("Lon.50.", "Lat.50."), crs = crs(wrld_simpl))
fall.nbr.sf <- st_transform(fall.nbr.sf, CRS(proj)) 
#st_write(fall.nbr.sf, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data/bpw_nbr_sites.shp")

# Origin sites 
fall.br.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/bpw_abundance_regions.shp") %>%
  st_transform(CRS(proj)) %>%
  st_cast("MULTIPOLYGON")

# Target sites 
fall.nbr.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data/Nonbreedingregions.shp") %>%
  st_transform(CRS(proj))%>%
  st_cast("MULTIPOLYGON")

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

geo.bias <- c(lon.bias = mean(geo.bias.dists$lon.geo.bias),
              lat.bias = mean(geo.bias.dists$lat.geo.bias))

# Geolocator variance covariance at origin sites 
#geo.vcov <- cov(fall.nbr[, c("Lon.50.","Lat.50.")])
geo.vcov <- cov(st_coordinates(fall.nbr.sf))

#number of bootstrap samples 
nSamples <- 100

## Verify that the breeding and nonbreeding sites and points for the fall makes sense

# # Nonbreeding 
# plot(spTransform(wrld_simpl, CRS(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/6, 19659061/6))
# plot(as_Spatial(fall.nbr.regions), add = T)
# plot(as_Spatial(fall.nbr.sf), add  = T)
# 
# # Breeding
# plot(spTransform(wrld_simpl, CRS(proj)), xlim = c(-19312101/6, 19889532/6), ylim = c(-19793221/8, 19659061/2))
# plot(as_Spatial(fall.br.regions), add = T)
# plot(as_Spatial(fall.br.sf), add  = T)

################################################################################
### Estimation and resampling of uncertainty for transition probabilities (psi)
################################################################################

# GL_psi <- estTransition(isGL=TRUE,
#                         geoBias = geo.bias,
#                         geoVCov = geo.vcov,
#                         targetSites = fall.nbr.regions,
#                         originSites = fall.br.regions,
#                         originPoints = fall.br.sf,
#                         targetPoints = fall.nbr.sf,
#                         verbose = 2,
#                         nSamples = nSamples,
#                         resampleProjection = CRS(proj),
#                         maxTries = 10000)

# Save the output of estTransition
#save(GL_psi, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/estTransition_ouput.R")
load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_analysis/estTransition_ouput.R")

################################################################################
# Retrieve the relative abundance in each breeding site/region for the fall network
################################################################################

# import breeding season abundance data
bpw.fall.ab <- load_raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports/2021/bkpwar",
                           product = "abundance",
                           period = "seasonal",
                           resolution = "lr") %>%
 terra::project(crs(proj))

# # Plot of breeding regions over relative abundance surface 
# plot(bpw.fall.ab$breeding, , xlim = c(-19312101/3, 19889532/6), ylim = c(-19793221/8, 19659061/2))
# plot(spTransform(wrld_simpl, CRS(proj)), add = T)
# plot(as_Spatial(fall.br.regions), add = T)
# plot(as_Spatial(fall.br.sf), add  = T)

# extract the relative abundance for each site 
fall.br.site.ab <- terra::extract(bpw.fall.ab$breeding, fall.br.regions, fun = sum,
                             na.rm=TRUE, ID  = T)

################################################################################
# Compute distance matrixes
################################################################################


