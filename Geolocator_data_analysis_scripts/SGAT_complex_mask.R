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

# Set eBird access key
set_ebirdst_access_key("k187tqe2nh49")

# Directory for eBird data 
ebird.dir <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports"

#Download eBird data for blackpoll warbler data
ebirdst_download("bkpwar", path = ebird.dir, tifs_only = TRUE)

################################################################################
# Create an empty raster with the desired dimensions and resolutions ###########
################################################################################

n = 10 # there will be n*n raster cells per grid squares of 1 degree of latitude and longitude 

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
range.poly <- load_ranges(get_species_path("bkpwar", ebird.dir), smoothed = T)
range.poly$prediction_year <- 1 
range.poly$start_date <- 1 
range.raster <- rasterize(range.poly, land.mask)
 
xbin = seq(xmin(maski),xmax(maski),length=ncol(maski)+1)
ybin = seq(ymin(maski),ymax(maski),length=nrow(maski)+1)

maskx <- as.matrix(maski)

## Define the log prior for x and z
log.prior <- function(p) {
  f <- mask(p)
  ifelse(is.na(f), log(1), f) 
}

#modified matrix subset approach 
maski[cbind(.bincode(lat.calib,ybin), .bincode(lon.calib,xbin))]

#original matrix subset approach
maski[cbind(length(ybin) -.bincode(lat.calib,ybin), .bincode(lon.calib,xbin))]

#breeding site location 
maski[cbind(length(ybin) -.bincode(x0[,2],ybin), .bincode(x0[,1],xbin))] <-2
maski[cbind(.bincode(x0[,2],ybin), .bincode(x0[,1],xbin))] <- 1

plot(maski)

getValues(mask)
