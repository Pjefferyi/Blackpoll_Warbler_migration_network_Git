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

#Blank map for the whole of the blackpoll warbler's range
wholerange <- tm_shape(shp = world, bbox = c(-180, -25, -30, 70)) + tm_borders()

#Blank map for the blackpoll warbler's breeding range
Breedingrange <- tm_shape(shp = world, bbox = c(-180, 30, -30, 70)) + tm_borders()

#Blank map for the blackpoll warbler's overwintering range
overwinteringrange <- tm_shape(shp = world, bbox = c(-180, -25, -30, 20)) + tm_borders()

#map of the blackpoll warbler's relative abundance in the breeding range
bpbbox = tmaptools::bb(matrix(c(
 -6000000, #xmin
 -2000000, #ymin  
 4000000,#xmax
 5000000  #ymax
),2,2))

breedingraster <- raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geo_spatial_data/bkpwar_abundance_seasonal_breeding_mean_2021.tif")
breedingraster[breedingraster == 0] <- NA

#"+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
#"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
projcrs <- crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

crs(breedingraster) <- projcrs

#convert dataframe of sampled sites to layer of points 
sites_df <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Sampling_sites.csv")
sites_points <- st_as_sf(sites_df, coords = c('Longitude', 'Latitude'), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

tm_shape(breedingraster, bbox = bpbbox,  raster.warp = FALSE) + 
    tm_raster(palette = "Reds", alpha = 1, title = "Relative abundance", style = "cont", legend.show = FALSE) +
  tm_shape(shp = world[(world$name_long %in% c("Canada", "United States")),]) +
  tm_borders(lwd = 1) +
  tm_shape(sites_points) +
  tm_dots(col = "black", size = 0.6)
#tm_grid()# +
#tm_layout(legend.position = c("left","BOTTOM"),
#legend.text.size = 0.8,
#legend.hist.size = 0.5,
#legend.title.size = 1)

#map of the blackpoll warbler's abundance in the nonbreeding range 
bpbbox_rangemap = tmaptools::bb(matrix(c(
  -10000000, #xmin
  -2500000, #ymin  
  -2000000,#xmax
  2000000  #ymax
),2,2))

rangeraster <- raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geo_spatial_data/bkpwar_abundance_seasonal_nonbreeding_mean_2021.tif")

crs(rangeraster) <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

tm_shape(rangeraster, bbox = bpbbox_rangemap,  raster.warp = FALSE) +
  tm_raster(palette = "Blues", alpha = 1, title = "Relative abundance", style = "cont") +
  tm_shape(shp = world, bbox = c(-180, -25, -30, 70)) +
  tm_borders(lwd = 1) +
  tm_layout(legend.position = c("right","TOP"))
          
#map of the blackpoll warbler's nonbreeding range (gray)
bpbbox_rangemap = tmaptools::bb(matrix(c(
  -10000000, #xmin
  -2500000, #ymin  
  -2000000,#xmax
  2000000  #ymax
),2,2))

rangeraster <- raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geo_spatial_data/bkpwar_abundance_seasonal_nonbreeding_mean_2021.tif")
rangeraster[rangeraster == 0] <- NA
rangeraster[rangeraster > 0] <- 1

crs(rangeraster) <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

tm_shape(rangeraster, bbox = bpbbox_rangemap,  raster.warp = FALSE) +
  tm_raster(col = 'layer', palette = "gray", alpha = 0.5, legend.show = FALSE) +
  tm_shape(shp = world, bbox = c(-180, -25, -30, 70)) +
  tm_borders(lwd = 1)

#map of the blackpoll warbler's breeding range (gray)
bpbbox_rangemap_breeding = tmaptools::bb(matrix(c(
  -20000000, #xmin
  5000000, #ymin  
  -2000000,#xmax
  12000000  #ymax
),2,2))

rangeraster_breeding <- raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geo_spatial_data/bkpwar_abundance_seasonal_breeding_mean_2021.tif")
rangeraster_breeding[rangeraster_breeding == 0] <- NA
rangeraster_breeding[rangeraster_breeding > 0] <- 1

crs(rangeraster_breeding) <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

tm_shape(rangeraster_breeding, bbox = bpbbox_rangemap_breeding,  raster.warp = FALSE) +
  tm_raster(col = 'layer', palette = "gray", alpha = 0.5, legend.show = FALSE) +
  tm_shape(shp = world, bbox = c(-180, -25, -30, 70)) +
  tm_borders(lwd = 1)


#map of the blackpoll warbler's breeding and nonbreeding range (red and blue)
bpbbox_fullrangemap = tmaptools::bb(matrix(c(
  -20000000, #xmin
  -2500000, #ymin  
  -2000000,#xmax
  12000000  #ymax
),2,2))

tm_shape(rangeraster_breeding, bbox = bpbbox_fullrangemap,  raster.warp = FALSE) +
  tm_raster(col = 'layer', palette = "blue", alpha = 0.1, legend.show = FALSE) +
  tm_shape(rangeraster, bbox = bpbbox_rangemap,  raster.warp = FALSE) +
  tm_raster(col = 'layer', palette = "red", alpha = 0.1, legend.show = FALSE) +
  tm_shape(shp = world, bbox = c(-180, -25, -30, 70)) +
  tm_borders(lwd = 1)


# map of the blackpoll warbler's breeding, nonbreeding and stopover range
bpbbox_fullrangemap = tmaptools::bb(matrix(c(
  -20000000, #xmin
  -2500000, #ymin  
  -2000000,#xmax
  12000000  #ymax
),2,2))

crs(rangeraster_breeding) <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

full_rangeraster <- raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geo_spatial_data/bkpwar_abundance_seasonal_full-year_mean_2021.tif")
full_rangeraster[full_rangeraster == 0] <- NA
full_rangeraster[full_rangeraster > 0] <- 1

crs(full_rangeraster) <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

tm_shape(full_rangeraster, bbox = bpbbox_fullrangemap,  raster.warp = FALSE) +
  tm_raster(col = 'layer', palette = "red", alpha = 0.1, legend.show = FALSE) +
  tm_shape(world, bbox = c(-180, -25, -30, 70)) +
  tm_borders(lwd = 1)

#Range geopackage
data <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Maps_&_Figures/Range_Map/bkpwar_range_2021/bkpwar_range_2021.gpkg")
