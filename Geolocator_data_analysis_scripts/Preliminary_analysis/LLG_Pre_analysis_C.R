# source: Deluca et al. 2015 
# tag number: C
# site: Bon portage island Nova Scotia

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

# load and install and load geolocator packages 
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

# geo deployment location 
lat.calib <- 43.464762
lon.calib <- -65.7499


# time of deployment
deploy <- anytime("2013-06-16  12:04:45", tz = "GMT")

# data directory
dir <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data"

###############################################################################
#DATA EXTRACTION ##############################################################
###############################################################################

# import lig data 
lig <- read.csv(paste0(dir,"/raw_data/Deluca_et_al_2015/Blackpoll Warbler eastern North America (data from DeLuca et al. 2015).csv"))

#remove rows with processed data
lig <- lig[(is.na(lig$comments) == TRUE),]  %>%
  #get data for individual C
    filter(individual.local.identifier == "C") %>%
  # rename columns  
    rename(Date = timestamp, Light = gls.light.level) %>%
  #remove rows before and after deployment time 
    filter(Date > deploy)

#convert Dates to as.POXIct format
lig$Date <- anytime(lig$Date, tz = "UTC")

###############################################################################
#TWILIGHT ANNOTATION ##########################################################
###############################################################################

threshold <- 1.5 

# visualize threshold over light levels  
col = colorRampPalette(c('black',"purple",'orange'))(50)[as.numeric(cut(lig[2000:5000,2],breaks = 50))]

thresholdOverLight(lig, threshold, span =c(30000, 35000))

# plot light levels over the deployment period 
offset <- 12 # adjusts the y-axis to put night (dark shades) in the middle

lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

#Detect twilight times, for now do not edit twilight times  
twl <- preprocessLight(lig, 
                       threshold = threshold,
                       offset = offset, 
                       lmax = 64,         # max. light value
                       gr.Device = "x11", # MacOS version (and windows)
                       dark.min = 60)

# Adjust sunset times by 120 second sampling interval
twl <- twilightAdjust(twilights = twl, interval = 120)

# Automatically adjust or mark false twilights 
twl <- twilightEdit(twilights = twl, 
                    window = 4,           
                    outlier.mins = 90,    
                    stationary.mins = 45, 
                    plot = TRUE)

# Visualize light and twilight time-series
lightImage(lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))