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


# Function: visualize threshold in plot of light level   
thresholdOverLight <- function(data, threshold, span = c()){
  
  col = colorRampPalette(c('black',"purple",'orange'))(50)[as.numeric(cut(data$Light[2000:5000],breaks = 50))]
  
  par(mfrow = c(1, 1), mar = c(2, 2, 2, 2) )
  with(data[span[1]:span[2],], plot(Date, Light, type = "o", pch=16,  col = col, cex = 0.5)) 
  abline(h=threshold, col="orange", lty = 2, lwd = 2)
}


# Function Identify the span of a time shift in a geolocator relative to the expect Greenwich Meam time

shiftSpan <- function(twl, lig, period, est.zenith, dep.lon, dep.lat){

  # Get a subset of the observed twilights 
  ob_twl_sub <- subset(twl, twl$Twilight > period[1] & twl$Twilight < period[2])
  
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
  
  #Get a subset of the expected twilight times 
  exp_twl_sub <- subset(exp_twl, exp_twl$twilight > period[1] & exp_twl$twilight < period[2])
  
  #measure the time shift 
  shift <- mean(ob_twl_sub$Twilight) - mean(exp_twl_sub$twilight)
  
  #return a list with the time shift, expected time, and observed time
  return(list(shift = shift, observed_times = ob_twl_sub,
              expected_times = exp_twl_sub,
              mean_observed_noon = mean(ob_twl_sub$Twilight),
              mean_expected_noon = mean(exp_twl_sub$twilight)))
}

# This function returns a difference of 54.03504 mins for V8296-005
# This function returns a difference of 


shiftSpan(twl = twl2,
          lig = lig,
          period = as.POSIXct(c("2019-06-19", "2019-06-22"), tz = "UTC"),
          est.zenith = 92,
          dep.lon = lon.calib,
          dep.lat = lat.calib)


shiftSpan(twl = twl2,
          lig = lig,
          period = as.POSIXct(c("2019-08-01", "2019-08-03"), tz = "UTC"),
          est.zenith = 92,
          dep.lon = lon.calib,
          dep.lat = lat.calib)


