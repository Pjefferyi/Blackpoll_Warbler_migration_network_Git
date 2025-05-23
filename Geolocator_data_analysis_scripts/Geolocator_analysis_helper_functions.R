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
library(MASS)

#load spatial packages 
library(ggmap)
library(terra)
library(spData)
library(sf)
library(cluster)
library(raster)

#load and install and load geolocator packages 
library(GeoLight)
# install_github("eldarrak/FLightR")
library(FLightR)
# install_github("SLisovski/TwGeos")
library(TwGeos)
# install_github("SWotherspoon/SGAT")
library(SGAT)
# install_github("SLisovski/GeoLocTools")
library(GeoLocTools)
#install_github("slisovski/invMovement")


# Set eBirdst key 
# set_ebirdst_access_key("mk5l4atjq2bg", overwrite = TRUE)

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
# expected Greenwich Mean time

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

# Function to retrieve the location data obtained using SGAT's GroupedThresholdModel for each geolocator 
# If the output of the model were edited to account the results of the light data review
# then the function will automatically extract the edited location data

findLocData <- function(geo.ids = NULL, check_col_length = F, edits = T){
  
  # Create a list of path to all files with location data
  folder_paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", full.names = T)
  geo_names <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data")
  
  location_set <- data.frame()
  
  # We must detect any geolocators where the geolocator data was edited
  # based on the result of the light data analysis to detect carribean stopovers. 
  ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")
  with_edits = ref.data[(ref.data$Fall_carrib_edits == T),]$geo.id
  
  # If no vector of geo.ids was provided, extract data for all geolocators in dataset
  if (is.null(geo.ids)){
    geo.ids = unique(ref.data$geo.id)
  }
  
  # Check the number of columns in the dataset for each geolocator 
  if (check_col_length == T){
    for (i in seq(1:length(folder_paths))){
      if (geo_names[i] %in% with_edits & geo_names[i] %in% geo.ids){
        file_path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_GroupedThreshold_summary_fall_edit.csv")
        load(file = file_path)
        print(ncol(sm.fall.edit))
      }
      if (geo_names[i] %in% geo.ids){
        file_path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_GroupedThreshold_summary.csv")
        load(file = file_path)
        print(folder_paths[i])
        print(ncol(sm))
      }
    }
    return()
  }
  
  for (i in seq(1:length(folder_paths))){
    # load the data from each file and add it to dataset if it is in geo_ids 
    if (geo_names[i] %in% geo.ids | is.null(geo.ids)){
      # some of the data has edits during during the fall transoceanic flight
      if (geo_names[i] %in% with_edits & edits == T){
        file_path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_GroupedThreshold_summary_fall_edit.csv")
        load(file = file_path)
        location_set <- rbind(location_set, sm.fall.edit)
      } else {
        file_path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_GroupedThreshold_summary.csv")
        load(file = file_path)
        location_set <- rbind(location_set, sm)
      } 
    }
  }
  
  # Add reference data to location data 
  location_set <- merge(location_set, ref.data, by.x = "geo_id", by.y = "geo.id")
  
  return(location_set)
}

# Test calls  for findLocData ##################################################

# check length of the location dataframe of a specific geolocator
# r1 <- findLocData(geo.ids = c("V8757_055"), check_col_length = T)

# check length of the dataframes with the data for each geolocator 
# r2 <- findLocData(check_col_length = T)

# Extract location data for specific geolocators 
# r3 <- findLocData(geo.ids = c("V8757_010", "V8296_004"), check_col_length = F)

# Extract location data for all geolocators 
# r4 <- findLocData()

# findThresLocData ##################################################################

# Function to retrieve threshold location data for each geolocator 

findThresLocData <- function(geo.ids = NULL){
  
  # Create a list of path to all files with location data 
  folder_paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", full.names = T)
  geo_names <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data")
  
  location_set <- data.frame()
  
  #load the geolocator reference data/metadata
  ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")
  
  # If no vector of geo.ids was provided, extract data for all geolocators in dataset
  if (is.null(geo.ids)){
    geo.ids = unique(ref.data$geo.id)
  }
  
  for (i in seq(1:length(folder_paths))){
    # load the data from each file and add it to dataset if it is in geo_ids 
    if (geo_names[i] %in% geo.ids | is.null(geo.ids)){
      
      print(geo_names[i])
      
      #extract twilight and threshold location data 
      thresloc.path <- paste0(folder_paths[i], "/",geo_names[i],"_initial_path_raw.csv")
      twl.path <- paste0(folder_paths[i], "/",geo_names[i],"_twl_times.csv")
      
      load(file = thresloc.path)
      twl <- read.csv(file = twl.path)
      
      # merge twilight data (for timing) with threshold location data 
      time.loc.data <- cbind(twl[c("X", "Twilight")], as.data.frame(x0_r)) %>% 
        mutate(geo_id = geo_names[i])
      
      location_set <- rbind(location_set, time.loc.data)
    }
  }
  
  # Add reference data 
  location_set <- merge(location_set, ref.data, by.x = "geo_id", by.y = "geo.id")
  
  return(location_set)
}

# findThresLocData ##################################################################

# Function to retrieve threshold location data for each geolocator 

findThresModData <- function(geo.ids = NULL){
  
  # Create a list of path to all files with location data 
  folder_paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", full.names = T)
  geo_names <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data")
  
  location_set <- data.frame()
  
  #load the geolocator reference data/metadata
  ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")
  
  # If no vector of geo.ids was provided, extract data for all geolocators in dataset
  if (is.null(geo.ids)){
    geo.ids = unique(ref.data$geo.id)
  }
  
  for (i in seq(1:length(folder_paths))){
    # load the data from each file and add it to dataset if it is in geo_ids 
    if (geo_names[i] %in% geo.ids | is.null(geo.ids)){
      
      print(geo_names[i])
      
      #extract threshold location data 
      thresmod.path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_Threshold_summary.csv")
      
      load(file = thresmod.path)
      
      # Add acolumn with the geo ID
      sm.s$geo_id <- geo_names[i]
      
      location_set <- rbind(location_set, sm.s)
    }
  }
  
  # Add reference data 
  location_set <- merge(location_set, ref.data, by.x = "geo_id", by.y = "geo.id")
  
  return(location_set)
}

# Test calls  for findThresLocData #############################################

# Extract threshold location data for specific geolocators
r1 <- findThresLocData(geo.ids = c("V8757_010"))

# Extract threshold location data for all geolocators
#r2 <- findThresLocData()

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

# #call for all geolocators
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
#                               "WRMA04173"), check_col_length = F)
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
# geo.sample <- findLocData(geo.ids = c("V8757_134",  "blpw14", "4210_004", "A", "3254_057", "V8757_029"), check_col_length = F)
# 
# geos <- geo.all
# 
# plotLocVec(data = geos, er_bars =  T, stati_only = T, legend = F,timing = c("Post-breeding migration", "Non-breeding period"))
# plotLocVec(data = geos, er_bars =  F, stati_only = F, legend = T,timing = c("Pre-breeding migration", "Non-breeding period"))
# plotLocVec(data = geos, er_bars =  T,stati_only = T, legend = T,timing = c("Non-breeding period"))
# plotLocVec(data = geos, er_bars =  T,stati_only = T, legend = T, timing = c("Post-breeding migration", "Pre-breeding migration", "Non-breeding period"))

# findSlicesData ##################################################################

# Function to recover the data required to obtain the density estimates for each geolocator 

findSlicesData <- function(periods, xlim = c(-170, -40), ylim = c(-40, 75)){
  
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
  r <- rast(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmin = xlim[1],
              xmax = xlim[2], ymin = ylim[1], ymax = ylim[2], crs = proj4string(wrld_simpl))
  
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
    dens.list[[geonames[f]]] <- sk 
  }
  
  return(dens.list)
}

# Arguments for findSlicesData #################################################

# periods: a vector of periods for which the density estimates will be extracted
# e.g. c("Non-breeding period", "Pre-breeding migration")

# Test calls  for findSlicesData ###############################################

# periods = c("Non-breeding period", "Pre-breeding migration", "Post-breeding migration")
# nbr.dens <- findSlicesData(periods)
# plot(nbr.dens[["V8757_029"]])
# plot(wrld_simpl, add = T)

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
  
  if (pacific) {wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1
  rs = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
             rasterize(wrld_simpl, r, 1, silent = TRUE), 
             rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  
  #load polygon of blackpoll's range
  load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Birdlife_int_Full_blackpoll_range_polygon.R")
  
  #rasterize the polygon 
  range.raster <- rasterize(BLI.range.poly, rs)
  
  #Update the stationary mask 
  rs <- range.raster * rs
  
  # make the movement raster the same resolution as the stationary raster, but allow the bird to go anywhere by giving all cells a value of 1
  rm = rs; rm[] = 1
  
  # stack the movement and stationary rasters on top of each other
  mask = stack(rs, rm)
  
  xbin = seq(xmin(mask),xmax(mask),length=ncol(mask)+1)
  ybin = seq(ymin(mask),ymax(mask),length=nrow(mask)+1)
  mask = as.array(mask)[,,sort(unique(index)),drop=FALSE]
  
  function(p) mask[cbind(length(ybin)-.bincode(p[,2],ybin), .bincode(p[,1],xbin), index)]
}

# Test calls  for earthseaMask  ###############################################

# arguments must be defined using geolocator data
#mask <- earthseaMask(xlim, ylim, n = 10, index=index)

# runGeoScripts ################################################################

#Function to run all geolocator analysis scripts

runGeoScripts <- function(scripts = c()){
  
  # Create a list of path to all files with location data 
  paths <- list.files("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis",
                      pattern = "LLG_analysis", recursive = T, full.names = T )
  
  # For loop to load the data from each file and add it to dataset 
  for (i in paths){
    
    script.name <- i
    
    if(script.name %in% scripts){
      print(i)
      source(i, local = F)
    }
  }
}

# Test calls to runGeoScripts ##################################################

paths <- list.files("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis",
                    pattern = "LLG_analysis", recursive = T, full.names = T)
 
# runGeoScripts(scripts = paths)

# insertLoc ####################################################################

# Function to add stopovers over the carribean that were undetected due to the equinox and that were identified by inspecting plots of longitude estimates and light-level measurements  

insertLoc <- function(data, lat.at.loc, start.date, end.date, period, thresh.locs, twl, geo_id, sep1, sep2){
  
  # threshold location estimates 
  thresh.loc.es <- data.frame(x = thresh.locs[,1], datetime = twl$Twilight) 
  
  # find the longitude for the stopover 
  lon.at.loc <- filter(thresh.loc.es, thresh.loc.es$datetime >= start.date & thresh.loc.es$datetime <= end.date)$x
  
  # Add values to the row that will be inserted
  new.loc <- data.frame(StartTime = start.date,
                        EndTime = end.date,
                        Lon.50. = mean(lon.at.loc),
                        Lon.2.5. = quantile(lon.at.loc, 0.025)[[1]],
                        Lon.97.5. = quantile(lon.at.loc, 0.975)[[1]],
                        Lat.50. = lat.at.loc,
                        Lat.2.5. = NA,
                        Lat.97.5. = NA,
                        sitenum = -1,
                        duration = as.numeric(difftime(end.date, start.date), unit = "days"),
                        geo_id = geo_id,
                        period = period)
  
  # remove overlapping locations recorded while bird was moving
  data.mod <- data  %>% filter(!(StartTime > start.date & EndTime <  anytime(end.date) + days(1) & duration == 0))#remove overlapping locations when the bird was moving
  data.mod <- data  %>% filter(!(StartTime > start.date & EndTime < anytime(end.date) + days(1))) # remove  stationary locations estimated during the Caribbean stopover 
  
  # Adjust stationary location times  
  data.mod <- data.mod %>% mutate(StartTime = if_else((StartTime < start.date & EndTime >  anytime(end.date) + days(1) & (EndTime - anytime(end.date)) > (anytime(start.date) - StartTime)) |
                                                        (StartTime > start.date & StartTime <  anytime(end.date) + days(1) & EndTime >  anytime(end.date) + days(1)),
                                                      anytime(end.date) + sep2,
                                                      StartTime)) %>%
    mutate(EndTime = if_else((StartTime < start.date & EndTime > start.date & EndTime <  anytime(end.date) + days(1) & (EndTime - anytime(end.date)) < (anytime(start.date) - StartTime)) |
                               (EndTime > start.date & EndTime >  anytime(end.date) + days(1) & StartTime < start.date),
                             anytime(start.date) - sep1,
                             EndTime))
  
  #Add the new row (location) and order data by date
  data.mod <- rbind(data.mod, new.loc) %>% ungroup() %>% arrange(StartTime)
  
  # recalculate durations
  data.mod$duration <- as.numeric(difftime(data.mod$EndTime, data.mod$StartTime), unit = "days") 
  
  #remove any locations with a negative duration (this is any short stop that was caugth in the buffer around the carribean stopover)
  data.mod <- data.mod %>% filter(duration >= 0)
  
  # edit the site numbers 
  data.mod[(data.mod$sitenum != 0), ]$sitenum <- seq(1, nrow(data.mod[(data.mod$sitenum != 0), ]), 1)
  
  print(paste0("Estimated stopover longitude: ", as.character(mean(lon.at.loc))))
  return(data.mod)
}

# clusterLocs ####################################################################

# Function to cluster stationary location estimates
# Clusters are made smaller until a certain minimum value is reached
# This mimunum value should be selected on the basis of geolocator uncertainty

# This function is based on the methods described by Lagasse et al. 2022: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0270957#sec015
# It returns a number of clusters: k 

clusterLocs <- function(locs, maxdiam = 300, lon.only = F){
  
  # locs: Input locations should be in a dafatframe with reference information 
  # maxdiam: the maximum diameter of clusters (measured as the median distance between the two furthest location in each cluster) in km
  
  if (lon.only == F){
  
    # function to Calculate the geodesic distance between points and creates a distance matrix
    geo.dist = function(df) {
      require(geosphere)
      d <- function(i,z){         # z[1:2] contain long, lat
        dist <- rep(0,nrow(z))
        dist[i:nrow(z)] <- distHaversine(z[i:nrow(z),1:2],z[i,1:2])
        return(dist/1000)
      }
      dm <- do.call(cbind,lapply(1:nrow(df),d,df))
      return(as.dist(dm))
    }
    
    # calculate distance matrix for all stationary locations
    dist.matrix <- geo.dist(locs[,c("Lon.50.","Lat.50.")])
    
    # set intial values for the number of clusters and their maximum diameter
    runs = 1
    k = 1
    diam = Inf 
    
    # increase the number of clusters until the mean diameter of the cluster is less than the maximum value
    while (diam > maxdiam){
      
      # cluster stationary location and add locations to dataframe 
      clust <- pam(dist.matrix, k, diss  = T, cluster.only =  T)
      locs$cluster <- clust
      
      # for each cluster, calculate the diameter
      diam.list <- c()
      
      for (i in unique(clust)){
        # Get the stationary locations in each subset 
        clust.subset <- locs %>% dplyr::filter(cluster == i)
        
        if (nrow(clust.subset) == 1){
          subset.dist.matrix <- 0
        } else {
          
          # Generate distance matrix
          subset.dist.matrix <- geo.dist(clust.subset [,c("Lon.50.","Lat.50.")])
          
          # Calculate the cluster diameter (here defined as the distance between the two furthest locations in the cluster)
          diam.list <- append(diam.list, max(subset.dist.matrix))
        }
      }
      
      diam <- median(diam.list)
      print(paste("median diameter =", as.character(diam), "km"))
      print(paste("Run #", as.character(runs)))
      print(paste("k =", as.character(k)))
      runs <-  runs + 1 
      k = k + 1 
    }
    
    return(list(clusters = clust, k = k-1))
    
  }else{
    
    # function to Calculate the geodesic distance between points and creates a distance matrix
    geo.dist = function(df) {
      require(geosphere)
      d <- function(i,z){         # z[1:2] contain long, lat
        dist <- rep(0,nrow(z))
        dist[i:nrow(z)] <- distHaversine(z[i:nrow(z),1:2],z[i,1:2])
        return(dist/1000)
      }
      dm <- do.call(cbind,lapply(1:nrow(df),d,df))
      return(as.dist(dm))
    }
    
    # calculate distance matrix for all stationary locations
    locs$zerocol <- 0 
    dist.matrix <- geo.dist(locs[,c("Lon.50.","zerocol")])
    
    # set intial values for the number of clusters and their maximum diameter
    runs = 1
    k = 1
    diam = Inf 
    
    # increase the number of clusters until the mean diameter of the cluster is less than the maximum value
    while (diam > maxdiam){
      
      # cluster stationary location and add locations to dataframe 
      clust <- pam(dist.matrix, k, diss  = T, cluster.only =  T)
      locs$cluster <- clust
      
      # for each cluster, calculate the diameter
      diam.list <- c()
      
      for (i in unique(clust)){
        # Get the stationary locations in each subset 
        clust.subset <- locs %>% dplyr::filter(cluster == i)
        
        if (nrow(clust.subset) == 1){
          subset.dist.matrix <- 0
        } else {
          
          # Generate distance matrix
          subset.dist.matrix <- geo.dist(clust.subset [,c("Lon.50.","zerocol")])
          
          # Calculate the cluster diameter (here defined as the distance between the two furtherst locations in the cluster)
          diam.list <- append(diam.list, max(subset.dist.matrix))
        }
      }
      
      diam <- median(diam.list)
      print(paste("median diameter =", as.character(diam), "km"))
      print(paste("Run #", as.character(runs)))
      print(paste("k =", as.character(k)))
      runs <-  runs + 1 
      k = k + 1 
    }
    
    return(list(clusters = clust, k = k-1))
  }
}

# test cases for clusterlocs ###################################################

# # path to reference data file 
# ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
# ref_data <- read.csv(ref_path)
# 
# # location data 
# geo.all <- findLocData(geo.ids = c("V8757_010",
#                                    "V8296_004",
#                                    "V8296_005",
#                                    "V8296_006",
#                                    "V8757_055",
#                                    "V8757_018",
#                                    "V8757_021",
#                                    "V8296_015",
#                                    "V8296_017",
#                                    "V8296_026",
#                                    "V8296_025",
#                                    "V8296_007",
#                                    "V8296_008",
#                                    "V8757_019",
#                                    "V8757_096",
#                                    "V8757_134",
#                                    "V8757_029",
#                                    "V8757_078",
#                                    "blpw09",
#                                    "blpw12",
#                                    "3254_001",
#                                    "4068_014",
#                                    "blpw14",
#                                    "3254_003",
#                                    "3254_008",
#                                    "3254_011",
#                                    "3254_057",
#                                    "blpw15",
#                                    "blpw25",
#                                    "4105_008",
#                                    "4105_009",
#                                    "4105_016",
#                                    "4105_017",
#                                    "4210_002",
#                                    "4210_004",
#                                    "4210_006",
#                                    "4210_010",
#                                    "WRMA04173",
#                                    "A",
#                                    "B",
#                                    "C",
#                                    #"E",
#                                    "D"), check_col_length = F)
# 
# # First extract stationary locations for the fall 
# geo.all <- geo.all %>% group_by(geo_id) %>% mutate(site_type = case_when(
#   (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ "Breeding",
#   sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ "Breeding",
#   (sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ "Breeding",
#   (sitenum == max(sitenum)) & Recorded_North_South_mig == "South and partial North" & path.elongated == T ~ "Breeding",
#   period == "Non-breeding period" & (duration >= 14 | sitenum == 1 | sitenum == max(sitenum)) ~ "Nonbreeding",
#   .default = "Stopover"))
# 
# fall.stat <- geo.all %>% filter(sitenum > 0, site_type %in% c("Stopover","Nonbreeding"),
#                                 period %in% c("Post-breeding migration","Non-breeding period"),
#                                 Recorded_North_South_mig %in% c("Both", "South and partial North", "South"))
# 
# # Run clustering function
# k <- clusterLocs(locs = fall.stat, maxdiam = 700)
