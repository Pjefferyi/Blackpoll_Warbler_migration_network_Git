loadpackages=function(packages){  for(p in packages){
  if(!require(p,character.only=T)){install.packages(p)}
  # IF require returns FALSE, the package is missing and will be installed
  library(p,character.only=T,quietly=T,verbose=F)}} # next, it calls the package with library

loadpackages(c(
  "ggplot2", # Makes ggplots!
  "sf",   # "Simple Features", handles vector data
  "raster", # For working with Raster Data
  "rgdal", # 'Geospatial' Data Abstraction Library
  # note: rgdal will be retired by the end of 2023
  "tmap", # Thematic Map Visualization
  "ggsn", # Several scale bars for maps
  "ggrepel", # Labels on ggplots
  "maptools", # Manipulating geographic data
  # note: maptools will be retired by the end of 2023
  "plyr", # The split-apply-combine paradigm
  "data.table", # Works with data.frames
  "dplyr", # Data manipulation
  "purrr", # Functional Programming Tools
  "devtools", # Download custom R packages
  "spData", # Spatial Datasets
  "terra",
  "tmap",
  "scales",
  "zoo"
))

if(!require("spDataLarge",character.only=T)){devtools::install_github("Nowosad/spDataLarge")}
# devtools::install_github installs a package from its github directory


################################################################################
# MAPPING MOVEBANK MIGRATION TRACKS ############################################
################################################################################

par(mfrow = c(1,1))

# import Movebank data  
deluca_2015 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Movebank_data/Blackpoll Warbler eastern North America (data from DeLuca et al. 2015).csv")
deluca_2019 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Movebank_data/A boreal songbird's migration across North America and the Atlantic Ocean, Setophaga striata.csv")

# Plot tracks for bird from Deluca et al. 2019 based on their Geolocator ID
locs <- deluca_2019 %>%
  filter(tag.local.identifier == "3254-011")

# some of the latitudes for this individuals are NA
locs <- locs %>% 
 # group_by(event.id) %>% 
  mutate(location.lat = ifelse(is.na(location.lat), na.locf(location.lat), location.lat))
  #ungroup 

#create track of points 
track <- st_as_sf(locs, coords = c('location.long', 'location.lat'), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#creare track as lines 
trackline <- track %>%
  st_coordinates() %>%
  st_linestring()

tm_shape(world, bbox = c(-180, -25, -30, 70)) +
  tm_borders(lwd = 1) +
  tm_shape(track) +
  tm_dots(col = "black", size = 0.05)

plot(wrld_simpl,  xlim=xlim, ylim=ylim)
points(locs$location.long, locs$location.lat, pch = 16, cex = 0.5, col = "firebrick")
lines(locs$location.long, locs$location.lat)

# Plot tracks for bird from Deluca et al. 2019 based on their Geolocator ID
locs <- deluca_2015 %>%
  filter(is.na(gls.light.level) & individual.local.identifier == "D")

# some of the latitudes for this individuals are NA
locs <- locs %>% 
  mutate(location.lat = ifelse(is.na(location.lat), na.locf(location.lat), location.lat)) %>%
  mutate(location.long = ifelse(is.na(location.long), na.locf(location.long), location.long))

#create track of points 
track <- st_as_sf(locs, coords = c('location.long', 'location.lat'), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#create track as lines 
trackline <- track %>%
  st_coordinates() %>%
  st_linestring()

#plot points with lines 
plot(wrld_simpl,  xlim=xlim, ylim=ylim)
points(locs$location.long, locs$location.lat, pch = 16, cex = 0.5, col = "firebrick")
lines(locs$location.long, locs$location.lat)

