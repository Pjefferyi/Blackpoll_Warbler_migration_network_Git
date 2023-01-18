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
  "scales"
))

if(!require("spDataLarge",character.only=T)){devtools::install_github("Nowosad/spDataLarge")}
# devtools::install_github installs a package from its github directory


################################################################################
# MAPPING MOVEBANK MIGRATION TRACKS ############################################
################################################################################

# import Movebank data  
deluca_2015 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/raw_data/Deluca_et_al_2015/Blackpoll Warbler eastern North America (data from DeLuca et al. 2015).csv")
deluca_2019 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/raw_data/Deluca_et_al_2019/Movebank/A boreal songbird's migration across North America and the Atlantic Ocean, Setophaga striata.csv")

# Plot tracks for bird with geolocator 4105-008 (Deluca et al. 2019)
locs <- deluca_2019 %>%
  filter(tag.local.identifier == "3254-003")

# some of the latitudes for this individuals are NA
locs <- locs %>% 
  group_by(event.id) %>% 
  mutate(location.lat = ifelse(is.na(location.lat), na.locf(location.lat), location.lat)) %>%
  ungroup 

track <- st_as_sf(locs, coords = c('location.long', 'location.lat'), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

trackline <- track %>%
  st_coordinates() %>%
  st_linestring()


tm_shape(world, bbox = c(-180, -25, -30, 70)) +
  tm_borders(lwd = 1) +
  tm_shape(track) +
  tm_dots(col = "black", size = 0.05)

plot(wrld_simpl)
lines(locs$location.long, locs$location.lat)
