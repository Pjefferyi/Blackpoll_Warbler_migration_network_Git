# clear objects from workspace
#rm(list=ls())

#Load libraries
library(tidyverse)
library(tidyterra)
library(igraph)
library(lubridate)
library(anytime)
library(ggplot2)
library(geosphere)
library(terra)
library(sf)
library(maptools)
library(ebirdst)

source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

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
                                   "D"), check_col_length = F, edits = F)

# Define nodes as breeding, stopover, or non-breeding ##########################

# Define the different site types 
geo.all <- geo.all %>% group_by(geo_id) %>% mutate(site_type = case_when(
  (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ "Breeding",
  sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ "Breeding",
  (sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ "Breeding",
  (sitenum == max(sitenum)) & Recorded_North_South_mig == "South and partial North" & path.elongated == T ~ "Breeding",
  period == "Non-breeding period" & (duration >= 14 | sitenum == 1 | sitenum == max(sitenum)) ~ "Nonbreeding",
  .default = "Stopover"))

# Define stopover types ########################################################

geo.all <- geo.all %>% mutate(site_type = case_when(
  site_type == "Stopover" & period == "Post-breeding migration" ~ "Fall stopover",
  site_type == "Stopover" & period == "Pre-breeding migration" ~ "Spring stopover",
  site_type == "Stopover" & period == "Non-breeding period" ~ "Nonbreeding stopover",
  .default = site_type
))

# plot geolocator paths ########################################################

ggplot(spData::world[(spData::world$continent %in% c("North America", "South America")),]) +
  geom_sf() +
  coord_sf() +
  geom_path(data = geo.all[(geo.all$sitenum) > 0,], aes(x = Lon.50., y = Lat.50., group = geo_id, colour = geo_id), linewidth = 0.2)+
  geom_point(data = geo.all[(geo.all$sitenum) > 0,], aes(x = Lon.50., y = Lat.50., group = geo_id, colour = geo_id), size = 0.3)+
  theme(legend.position = "None")

# Find paths to probability density estimates for all birds with geolocator data #####################

# Create a list of path to all files with fit data 
paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data",
                    pattern = "SGAT_GroupedThreshold_fit.R", recursive = T, full.names = T )

#remove duplicates
duplicated(paths)

# Create a list of geolocator names for the slices
spl.paths <- strsplit(paths, "/")
geonames <- seq(1:length(paths))

for (i in seq(1:length(paths))){
  
  geonames[i] <- spl.paths[[i]][11]
}

# Set the name of each file path (vector entries) to be a geo.id
names(paths) <- geonames

# Parameters to set prior to the data extraction #####################################################

# site for which we will extrcat probability density rasters 
site.type <- "Spring stopover"

# A grid with the dimensions of our rasters 
xlim = c(-180, 180)
ylim = c(-90, 90)

r <- raster(ncol = sum(diff(xlim))*2, nrow = sum(diff(ylim))*2, xmn = xlim[1], xmx = xlim[2],
            
            ymn = ylim[1], ymx = ylim[2])

# For each bird extract distribution ##################################################################

dist.list <- vector(mode='list', length= length(geonames))
names(dist.list) <- geonames

for (bird in (geonames)){
  
  if (bird %in% geo.all$geo_id){
  
    # load the model fit file 
    load(paths[bird])
    
    # Create a slices object 
    s <- slices(type = "intermediate", breaks = NULL, mcmc = fit, grid = r)
    
    ##Create a raster with the indexes and their time intervals 
    #index <- sliceIndices(s)
    #s.df <- as.data.frame(index)
    
    #int.func <- function(ind.ob) sliceInterval(s, k = c(ind.ob))
    #s.df$start <- anytime(unlist(lapply(lapply(sliceIndices(s), int.func), '[[',1)), asUTC  = T, tz = "UTC")
    #s.df$end <- anytime(unlist(lapply(lapply(sliceIndices(s), int.func), '[[',2)),  asUTC  = T,tz = "UTC")
    
    geo.info <- geo.all[(geo.all$geo_id == bird),]
    time.index <- which(geo.info$site_type == site.type)
    
    sk <-  SGAT::slice(s, k = time.index-1) 
    
    #sliceInterval(s, k = time.index-1, mcmc = s$mcmc)
    #sliceIndices(s, s$mcmc)
    
    dist.list[bird] <- sk
    
  }
}

# remove null elements 
dist.list <-dist.list[!sapply(dist.list ,is.null)]

# plot rasters
# for (i in seq(1, length(dist.list))){
#  plot(dist.list[[i]], main = names(dist.list[i]))
# }

# Convert the rasters to binary surface (1/0), where pixels with a value of 1 fall in the 95% credibility interval

# We create a function

binSurface <- function(surface, prob){
  
  surface <- rast(surface)
  
  # 95% probability 
  problim <- quantile(values(surface), probs = prob, na.rm = T)
  
  #edit raster surface
  values(surface)[values(surface) < problim] <- 0
  values(surface)[values(surface) >= problim] <- 1
  values(surface)[is.na(values(surface))] <- 0

  
  #return raster
  return(surface)
}

#function test case 
#y <- binSurface(dist.list[[1]], prob = c(0.95))
#plot(y)

# apply the function to all the rasters 
for (i in seq(1, length(dist.list))){
  dist.list[[i]] <- binSurface(dist.list[[i]], prob = c(0.95)) 
}

# convert the raster list of a terra raster dataset 
surface.sds <- sds(dist.list)

# sum the rasters 
surface.sum <- app(surface.sds, fun = sum)

#apply a water mask and turn zos into NA
surface.sum <- terra::mask(x = surface.sum, mask = vect(wrld_simpl))

# # plot the summarized raster surface, in this case it shows the estimated bumber of bird in each pixel
# ggplot(st_as_sf(wrld_simpl))+
#   geom_sf(fill = "white", colour = "darkgray") +
#   coord_sf(xlim = c(-88, -30),ylim = c(30, -15)) +
#   geom_spatraster(data = surface.sum) +
#   #scale_fill_continuous(low = adjustcolor("white", alpha = 1), high = "firebrick", na.value = "white")+
#   scale_fill_gradientn(colours = rev(terrain.colors(255)), na.value = "white")+
#   geom_sf(fill = NA, colour = "black") +
#   coord_sf(xlim = c(-88, -30),ylim = c(30, -15)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #axis.title.x=element_blank(),
#         #axis.text.x=element_blank(),
#         #axis.ticks.x=element_blank(),
#         #axis.title.y=element_blank(),
#         #axis.text.y=element_blank(),
#         #axis.ticks.y=element_blank(),
#         legend.position = c(0.75, 0.75)) +
#   theme(text = element_text(size = 14))+
#   labs(fill = "#Birds")

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(fill = "white", colour = "darkgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(80, -15)) +
  geom_spatraster(data = surface.sum) +
  #scale_fill_gradient2(low = adjustcolor("white", alpha = 0), high = "firebrick", na.value = "white")+
  scale_fill_gradientn(colours = rev(terrain.colors(255)), na.value = "white")+
  geom_sf(fill = NA, colour = "black") +
  coord_sf(xlim = c(-170, -30),ylim = c(80, -15)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.1, 0.25)) +
  theme(text = element_text(size = 14))+
  labs(fill = "#Birds")

# Notes 
load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/3254_001/3254_001_SGAT_GroupedThreshold_summary.csv")
nrow(sm)

nrow(geo.all[(geo.all$geo_id == bird),])


s <- slices(type = c("intermediate"), breaks = "week", mcmc = fit, grid = r)
sk <- slice(s, sliceIndices(s))

plot(sk, useRaster = F, legend = F, breaks = c(seq(0, quantile(sk, prob = 0.95), length = 100), maxValue(sk)),
     
     col = c(rep("transparent", 1), rev(topo.colors(99))))
