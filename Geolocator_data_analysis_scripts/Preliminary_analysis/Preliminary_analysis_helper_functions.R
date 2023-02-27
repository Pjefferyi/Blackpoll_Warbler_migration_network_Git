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


# thresholdOverLight############################################################

#Function to visualize the threshold in a plot of light level over time 

thresholdOverLight <- function(data, threshold, span = c()){
  
  col = colorRampPalette(c('black',"purple",'orange'))(50)[as.numeric(cut(data$Light[2000:5000],breaks = 50))]
  
  par(mfrow = c(1, 1), mar = c(2, 2, 2, 2) )
  with(data[span[1]:span[2],], plot(Date, Light, type = "o", pch=16,  col = col, cex = 0.5)) 
  abline(h=threshold, col="orange", lty = 2, lwd = 2)
}


# shiftSpan ####################################################################

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
  
  dates <- rep(seq(from = min(twl$Twilight), to = max(twl$Twilight), by = "day"), each = 2)
  
  # Check that there are no missing days 
  if days(ob_twl_sub$Twilight) == rep()

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

  #There are some cases where we must change the subset of observed times
  if ((ob_twl_sub$Rise[1] == exp_twl$rise[1] & exp_twl$twilight[1] > exp_twl$twilight[2])|
    (ob_twl_sub$Rise[1] != exp_twl$rise[1] & exp_twl$twilight[1] < exp_twl$twilight[2])){

  # expected twilights must be adjusted
  exp_twl_sub <- subset(exp_twl, exp_twl$twilight > period[1] & exp_twl$twilight < period[2])

  exp_twl_sub <- exp_twl_sub[-1,]

  exp_twl_sub[nrow(exp_twl_sub) +1,] <- exp_twl[(exp_twl$twilight > period[2]),][1,]

  # Else, if sunrise and sunset occur within the span of a single calendar day (based on GMT),
  # we will calculate the shift using the timing of midnight
  } else {
  
  # expected twilights do not have to be adjusted 
  exp_twl_sub <- subset(exp_twl, exp_twl$twilight > period[1] & exp_twl$twilight < period[2])
  
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



