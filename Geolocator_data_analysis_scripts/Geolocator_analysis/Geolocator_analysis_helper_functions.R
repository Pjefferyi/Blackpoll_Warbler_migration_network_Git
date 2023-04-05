#Helper function for the geolocator primary analysis

#load packages
require(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(remotes)
library(anytime)
library(lubridate)

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

findLocData <- function(geo.ids = c(), check_col_length = F){

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
  
  return(location_set)
}

# test calls  for findLocData ##################################################

# check length of the dataframes with the data for each geolocator 
# r1 <- findLocData(geo.ids = c(), check_col_length = T)

# Extract location data for specific geolocators 
# r2 <- findLocData(geo.ids = c("V8757_010", "V8296_004"), check_col_length = F)

# MapLocData ###################################################################

# Function to map vector data from multiple geolocator analyses

plotLocVec <- function(data, stati_only = F, timing = c("Post-breeding migration",
                                                        "Pre-breeding migration",
                                                        "Non-breeding migration")){
  # whether to use only stationary locations
  if (stati_only == T){
  data <- data[(data$sitenum > 0),]
  }
  
  data <- data[(data$period %in% timing),]
  
  #Create the map 
  ggplot(spData::world[(spData::world$continent %in% c("North America", "South America")),]) +
    geom_sf() +
    coord_sf() +
    geom_point(data = data, mapping = aes(x = Lon.50., y = Lat.50., color = geo_id), size = 1.2) +
    geom_path(data = data, mapping = aes(x = Lon.50., y = Lat.50., color = geo_id), linewidth = 0.4) +
    theme_bw() 
}
    
# test calls  for mapLocData ###################################################

# call for all geolocators 
r2 <- findLocData(geo.ids = c("V8757_010",
                                 "V8296_004", 
                                 "V8757_055",
                                 "V8757_018",
                                  "V8296_015",
                                  "V8296_017",
                                 "V8296_007",
                                 "V8296_008",
                                "V8757_019",
                                "V8757_096",
                              "V8757_134",
                              "V8757_029"), check_col_length = F)

# Call for geolocators from Quebec
#r2 <- findLocData(geo.ids = c("V8296_004", "V8757_018", "V8296_007", "V8296_015", "V8296_017"), check_col_length = F)

# Call for geolocators from Newfoundland 
#r2 <- findLocData(geo.ids = c("V8757_096", "V8757_134"), check_col_length = F)

# Call for geolocators from British Columbia
#r2 <- findLocData(geo.ids = c("V8757_019", "V8757_010", "V8757_029"), check_col_length = F)

# Custom call
# r2 <- findLocData(geo.ids = c("V8757_055"), check_col_length = F)

plotLocVec(data = r2, stati_only = T, timing = c("Post-breeding migration", "Non-breeding period"))
plotLocVec(data = r2, stati_only = T, timing = c("Pre-breeding migration", "Non-breeding period"))  
plotLocVec(data = r2, stati_only = T, timing = c("Non-breeding period"))  
plotLocVec(data = r2, stati_only = T, timing = c("Post-breeding migration", "Pre-breeding migration", "Non-breeding period"))
 
