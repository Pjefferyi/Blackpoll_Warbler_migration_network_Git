# Packages

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
  "ebirdst" #Extract and process eBird status and trends data
))

# Set eBird access key (key already set)
#set_ebirdst_access_key("k187tqe2nh49")

# Directory for eBird data 
ebird.dir <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports"

#Download eBird data for blackpoll warbler data (already downloaded)
#ebirdst_download("bkpwar", path = ebird.dir, tifs_only = TRUE)

################################################################################
# Create an empty raster with the desired dimensions and resolutions ###########
################################################################################

n = 20 # there will be n*n raster cells per grid squares of 1 degree of latitude and longitude 

# Boundaries of the raster
xlim = c(-120, -60) 
ylim = c(-22, 65)

# create empty raster with desired resolution
r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
           xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))

# create a raster for the stationary period, in this case by giving land a value of 1 and sea NA
land.mask = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
              rasterize(wrld_simpl, r, 1, silent = TRUE), 
              rasterize(elide(wrld_simpl, shift = c(360, 0)), r, 1, silent = TRUE))


################################################################################
# Create a range mask using abundance data available in eBird Status and Trends
###############################################################################

# Load data (this data is in a .tif format)
abundance <- raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geo_spatial_data/bkpwar_abundance_seasonal_full-year_mean_2021.tif")

# Project abundance data to the same parameters as the land mask 
abundance_resamp <- projectRaster(abundance, land.mask, method = "ngb")

# Set NA in cells where abundance is 0 or NaN
abundance_resamp[abundance_resamp == 0 | is.nan(abundance_resamp)] <- NA

# areas where abundance is greater than 1 get a value of 1 
abundance_resamp[abundance_resamp > 0 ] <- 1

#Create range mask 
range.mask1 <- abundance_resamp * land.mask

################################################################################
# Create a range mask using the range polygons available in eBird status and 
# trends 
################################################################################
#load range data
EB.data <- load_ranges(get_species_path("bkpwar", ebird.dir), resolution = "mr", smoothed = T)

#Dissolve the polygons 
EB.range <- EB.data %>%
  group_by("species_code") %>% 
  dplyr::summarise()

# this needs work 

################################################################################
# Create a range mask using the data provided by Birdlife international
################################################################################

#load data 
BLI.data <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Birdlife_international_species_distribution/SppDataRequest.shp")

#Dissolve the polygons 
BLI.range <- BLI.data %>%
  group_by("sci_name") %>% 
  dplyr::summarise()

plot(BLI.range)


