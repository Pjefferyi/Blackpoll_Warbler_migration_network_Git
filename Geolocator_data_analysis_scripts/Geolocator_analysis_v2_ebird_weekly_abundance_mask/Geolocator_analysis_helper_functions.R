#Helper function for the geolocator primary analysis

#load packages
require(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(remotes)
library(anytime)
library(lubridate)
library(ebirdst)

#load spatial packages 
library(ggmap)
library(terra)
library(spData)
library(sf)

#load and install and load geolocator packages 
# install_github("eldarrak/FLightR")
library(FLightR)
# install_github("SLisovski/TwGeos")
library(TwGeos)
# install_github("SWotherspoon/SGAT")
library(SGAT)
# install_github("MTHallworth/LLmig")
library(LLmig)
# install_github("SLisovski/GeoLocTools")
library(GeoLocTools)
setupGeolocation()

# Set eBirdst key 
#set_ebirdst_access_key("mk5l4atjq2bg", overwrite = TRUE)

# thresholdOverLight############################################################

#Function to visualize the threshold in a plot of light level over time 

thresholdOverLight <- function(data, threshold, span = c()){
  
  col = colorRampPalette(c('black',"purple",'orange'))(50)[as.numeric(cut(data$Light[2000:5000],breaks = 50))]
  
  par(mfrow = c(1, 1), mar = c(2, 2, 2, 2) )
  with(data[span[1]:span[2],], plot(Date, Light, type = "o", pch=16,  col = col, cex = 0.5)) 
  abline(h=threshold, col="orange", lty = 2, lwd = 2)
}

# Shift time recordings  #######################################################

# Function to identify the span of a time shift in a geolocator relative to the 
# expect Greenwich Mean time

# This function works by finding the difference between noon (or midnight) 
# in the times recorded by the geolocator in the breeding grounds, and the expected times
# based on the output of TwGeos' twilight() function

# Noon and midnight occur when the sun crosses the meridian, so at the halfway point 
# of every day and night, respectively. 

shiftSpan <- function(twl, lig, period, est.zenith, dep.lon, dep.lat){
  
  # Get a subset of the observed twilights 
  ob_twl_sub <- subset(twl, twl$Twilight > period[1] & twl$Twilight < period[2])
  
  check.dates <- rep(seq(from = date(min(ob_twl_sub$Twilight)), to = date(max(ob_twl_sub$Twilight)), by = "day"), each = 2)
  
  # Check that there are no missing days 
  if (F %in% (date(ob_twl_sub$Twilight) == check.dates)) {
    
    print("One or more days twilights are missing in twl")
    
  }

  # # Sort the subset of data so that sunsrise is always reported before sunrise
  # # otherwise they may not match with the expected data 
  # ob_twl_sub <- ob_twl_sub[order(date(ob_twl_sub$Twilight), -ob_twl_sub$Rise),]
  
  # Convert the lig file import to a series of days and rises  
  dates <- seq(from = min(lig$Date), to = max(lig$Date), by = "day")
  rise <- rep(c(TRUE, FALSE), length(dates))
  
  # Generate the expected twilight times for the location (with time in GMT)
  exp_twl <-  data.frame(twilight =
                           twilight(rep(dates, each = 2), 
                                    lon = dep.lon,
                                    lat = dep.lat,
                                    zenith = est.zenith, # adjust zenith to match observed and known twilights
                                    rise = rise),
                         rise = rise)
  
  exp_twl_sub <- subset(exp_twl, exp_twl$twilight > period[1] & exp_twl$twilight < period[2])
  
  #There are some cases where we must change the subset of observed times
  if ((ob_twl_sub$Rise[1] == exp_twl_sub$rise[1] & exp_twl_sub$twilight[1] > exp_twl_sub$twilight[2])|
    (ob_twl_sub$Rise[1] != exp_twl_sub$rise[1] & exp_twl_sub$twilight[1] < exp_twl_sub$twilight[2])){

  # expected twilights must be adjusted
  exp_twl_sub <- exp_twl_sub[-1,]

  exp_twl_sub[nrow(exp_twl_sub) +1,] <- exp_twl[(exp_twl$twilight > period[2]),][1,]

  # Else, if sunrise and sunset occur within the span of a single calendar day (based on GMT),
  # we will calculate the shift using the timing of midnight
  } 
  
  # Check that the two list of twilights are the same length 
  if (nrow(ob_twl_sub) != nrow(exp_twl_sub)) {
    warning("The observed and predicted lists of twilights have different lengths")
  }
  
  t1 <- mean(ob_twl_sub$Twilight)
  t2 <- mean(exp_twl_sub$twilight)
  
  #measure the time shift: the time between the observed and expected noons or midnights 
  shift <- t1 - t2
  
  #return a list with the time shift, expected time, and observed time
  return(list(shift = shift, observed_times = ob_twl_sub,
              expected_times = exp_twl_sub,
              mean_observed_meridian_time = mean(ob_twl_sub$Twilight),
              mean_expected_meridian_time = mean(exp_twl_sub$twilight)))
}


# findLocData ##################################################################

# Function to recover location data from the folders allocated to each geolocator 

findLocData <- function(geo.ids = c(), check_col_length = F, ref_path = NA){

  # Create a list of path to all files with location data 
  paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data",
                     pattern = "SGAT_GroupedThreshold_summary.csv", recursive = T, full.names = T )
  
  location_set <- data.frame()
  
  # Check the number of columns in the dataset for each geolocator 
  if (check_col_length == T){
    for (f in paths){
      load(file = f)
      print(f)
      print(ncol(sm))
    }
    return()
  }
  
  # For loop to load the data from each file and add it to dataset 
  for (f in paths){
    
    #add the geolocator to the larger dataset if it is in geo_ids
    load(file = f)
    if (length(geo.ids) == 0 | sm$geo_id[1] %in% geo.ids) {
      location_set <- rbind(location_set, sm) 
    }
  }
  

  # Add reference data 
  if (!is.na(ref_path)){
    
    ref.data <- read.csv(ref_path)
    
    # Only retain relevant rows
    # ref.data <- ref.data[,c("geo.id",
    #                         "deploy.latitude",
    #                         "deploy.longitude",
    #                         "study.site",
    #                         "Range_region")]
    
    #Join the region and location data (inner join)
    location_set <- merge(location_set, ref.data, by.x = "geo_id", by.y = "geo.id")
  }
  
  return(location_set)
}

# Test calls  for findLocData ##################################################

# check length of the dataframes with the data for each geolocator 
# r1 <- findLocData(geo.ids = c(), check_col_length = T)

# Extract location data for specific geolocators 
# r2 <- findLocData(geo.ids = c("V8757_010", "V8296_004"), check_col_length = F)

# Extract location data for specific geolocators with reference information 

# path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"   
# r3 <- findLocData(geo.ids = c("V8757_010", "V8296_004"), check_col_length = F, ref_path = path)

# MapLocData ###################################################################

# Function to map vector data from multiple geolocator analyses

plotLocVec <- function(data, stati_only = F, timing = c("Post-breeding migration",
                                                        "Pre-breeding migration",
                                                        "Non-breeding migration")
                       , er_bars = F, legend = T){
  
  # whether to use only stationary locations
  if (stati_only == T){
  data <- data[(data$sitenum > 0),]
  }
  
  data <- data[(data$period %in% timing),]
  
  #Create the map 
  ggplot(spData::world[(spData::world$continent %in% c("North America", "South America")),]) +
    geom_sf() +
    coord_sf() +
    {if(er_bars ==  T)geom_errorbar(data = data, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, alpha = 0.3, color = "black")} + 
    {if(er_bars ==  T)geom_errorbar(data = data, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, alpha = 0.3, color = "black")} + 
    geom_point(data = data, mapping = aes(x = Lon.50., y = Lat.50., color = geo_id), size = 1.1) +
    #geom_path(data = data, mapping = aes(x = Lon.50., y = Lat.50., color = geo_id), linewidth = 0.3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    {if(legend ==  F)theme(legend.position = "none")}
    
}
    
# Test calls  for mapLocData ###################################################

#call for all geolocators

# ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"   
# 
# geo.all <- findLocData(geo.ids = c("V8757_010",
#                               "V8296_004",
#                               "V8296_005",
#                               "V8296_006",
#                               "V8757_055",
#                               "V8757_018",
#                               "V8757_021",
#                               "V8296_015",
#                               "V8296_017",
#                               "V8296_026",
#                               "V8296_025",
#                               "V8296_007",
#                               "V8296_008",
#                               "V8757_019",
#                               "V8757_096",
#                               "V8757_134",
#                               "V8757_029",
#                               "V8757_078",
#                               "blw09",
#                               "blpw12",
#                               "3254_001",
#                               "4068_014",
#                               "blpw14",
#                               "3254_003",
#                               "3254_008",
#                               "3254_011",
#                               "3254_057",
#                               "blpw15",
#                               "blpw25",
#                               "4105_008",
#                               "4105_009",
#                               "4105_016",
#                               "4105_017",
#                               "4210_002",
#                               "4210_004",
#                               "4210_006",
#                               "4210_010",
#                               "A",
#                               "B",
#                               "C",
#                               "D",
#                               "WRMA04173"), check_col_length = F, ref_path = ref_path)
# 
# 
# geo.church <- geo.all[(geo.all$study.site == "Churchill, Manitoba"),]
# geo.nome <- geo.all[(geo.all$study.site == "Nome, Alaska"),]
# geo.denali <- geo.all[(geo.all$study.site == "Denali, Alaska"),]
# geo.whitehorse <- geo.all[(geo.all$study.site == "Whitehorse, Yukon"),]
# geo.west <- geo.all[(geo.all$study.site %in% c("Whitehorse, Yukon", "Nome, Alaska", "Denali, Alaska")),]
# geo.quebec <- geo.all[(geo.all$study.site == "Quebec"),]
# geo.2015 <- geo.all[(geo.all$Study == "Deluca et al. 2015"),]
# 
# geo.sample <- findLocData(geo.ids = c("V8757_134",  "blpw14", "4210_004", "A", "3254_057", "V8757_029"), check_col_length = F, ref_path = ref_path)
# 
# geos <- geo.all

# plotLocVec(data = geos, er_bars =  T, stati_only = T, legend = T,timing = c("Post-breeding migration", "Non-breeding period"))
# plotLocVec(data = geos, er_bars =  T, stati_only = T, legend = T,timing = c("Pre-breeding migration", "Non-breeding period"))
# plotLocVec(data = geos, er_bars =  T,stati_only = T, legend = T,timing = c("Non-breeding period"))
# plotLocVec(data = geos, er_bars =  T,stati_only = T, legend = T, timing = c("Post-breeding migration", "Pre-breeding migration", "Non-breeding period"))


# findSlicesData ##################################################################

# Function to recover thedata required to obtain the density estimates for each geolocator 

findSlicesData <- function(periods, xlim = c(-170, -40), ylim = c(-40, 75)  ){
  
  # Create a list of path to all files with fit data 
  paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data",
                      pattern = "SGAT_GroupedThreshold_fit.R", recursive = T, full.names = T )
  
  # Create a list of geolocator names for the slices
  spl.paths <- strsplit(paths, "/")
  geonames <- seq(1:length(paths))
  
  for (i in seq(1:length(paths))){
    
    geonames[i] <- spl.paths[[i]][11]
  }
  
  #list to store the density rasters 
  dens.list <- list() 
  
  #Create a grid to project density rasters 
  r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1],
              xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # extract the density data 
  for (f in seq(1:length(paths))) {
    
    # load the model fit file 
    load(paths[f])
    
    # Create a slices object 
    s <- slices(type = "intermediate", breaks = NULL, mcmc = fit, grid = r)
    
    # Retrieve the index for the period of interest
    load(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geonames[f] ,"/", geonames[f] ,"_SGAT_GroupedThreshold_summary.csv"))
    period.index <- which(sm$period %in% periods)
    
    # generate the density estimate  
    sk <- slice(s, sliceIndices(s)[period.index-1]) # must had -1 here so that the index coresponds with the bins 
    
    # save the density raster
    dens.list[[geonames[f]]] <-sk 
  }
  
  return(dens.list)
}

# Arguments for findSlicesData #################################################

# periods: a vecator of periods for which the density estimates will be extracted
# e.g. c("Non-breeding period", "Pre-breeding migration")

# Test calls  for findSlicesData ###############################################

# periods = c("Non-breeding period")
# nbr.dens <- findSlicesData(periods)
# plot(nbr.dens[["A"]])

# plotBreedSites ##################################################################

# Function to plot the breeding sites for all blackpoll warblers in the analysis

plotBreedSites <- function(){
  
  #load location data
  geo.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")
  
  #Create and crop base sf object
  range_crop <- st_crop(spData::world[(spData::world$continent %in% c("North America", "South America")),],
                        xmin = -190, xmax = -10, ymin = -10, ymax = 60)
  
  #Create the map of breeding sites 
  brd.sites <- ggplot(range_crop) +
    geom_sf() +
    coord_sf() +
    geom_point(data = geo.data, mapping = aes(x = deploy.longitude, y = deploy.latitude), size = 1.4, col = "firebrick") +
    xlab("Longitude") +
    ylab("Latitude")+
    theme_bw() 
  
  return(brd.sites)
}  


# earthseaMask ####################################################################

# Function to create a spatial mask
earthseaMask <- function(xlim, ylim, n = 2, pacific=FALSE, index) {
  
  if (pacific) { wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1
  rs = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
             rasterize(wrld_simpl, r, 1, silent = TRUE), 
             rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  
  #load weekly rasters of blackpoll warbler abundance
  ab.ras <- load_raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports/2021/bkpwar",
                        product = "occurrence",
                        period = "weekly",
                        resolution = "lr")
  
  names(ab.ras) <- as.numeric(strftime(names(ab.ras), format = "%W"))
  
  #project the abundance rasters
  ab.ras.pr <- terra::project(ab.ras, crs(proj4string(wrld_simpl)), method = "near") 
  
  values(ab.ras.pr)[is.nan(values(ab.ras.pr))] <- NA
  values(ab.ras.pr)[!is.na(values(ab.ras.pr))] <- 1
  
  # get bincodes linking geolocator twilight measurement times to the weeks of  
  # each abundance layer 
  t.datetime <- (twl %>% filter(group != lag(group, default = -1)))$Twilight
  t.weeks <- week(t.datetime)
  t.code <- .bincode(t.weeks, as.numeric(names(ab.ras)))
  
  xbin = seq(xmin(ab.ras.pr),xmax(ab.ras.pr),length=ncol(ab.ras.pr)+1)
  ybin = seq(ymin(ab.ras.pr),ymax(ab.ras.pr),length=nrow(ab.ras.pr)+1)
  ab.arr <- as.array(ab.ras.pr)
  
  # If the bird becomes stationary, the prior for a location is selected from an eBird 
  # relative abundance raster for the week during which the bird arrived at the location 
  # While the bird is moving, the prior always has a value of 0
  function(p){
    
    ifelse(stationary, 
           ab.arr[cbind(length(ybin)-.bincode(p[,2],ybin), .bincode(p[,1],xbin), t.code)],
           0)
  }
}

# Test calls  for earthseaMask  ###############################################

# arguments must be defined using geolocator data
# mask <- earthseaMask(xlim, ylim, n = 10, index=index)
