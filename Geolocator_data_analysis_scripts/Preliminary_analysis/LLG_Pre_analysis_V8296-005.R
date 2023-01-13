# source: 
# tag number: V8296 005
# site: Quebec

#load packages
require(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(remotes)
library(anytime)

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
lat.calib <- 47.38323	
lon.calib <- -71.09479

# time of deployment
deploy <- anytime("2019-06-18	18:45:00", tz = "GMT")
  
# data directory
dir <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data"

###############################################################################
#DATA EXTRACTION ##############################################################
###############################################################################

# import lig data 
lig <- readLig(paste0(dir,"/raw_data/BLPW_geo_2020/breed/ML6740 V8296 005/ML6740 V8296 005 reconstructed_000.lig"))

#remove the first row 
lig <- lig[2:nrow(ligdata),]

#remove rows before and after deployment time 
lig <- lig[(lig$Date > deploy),]

lig$Date <- lig$Date - 1*60*60

###############################################################################
#TWILIGHT ANNOTATION ##########################################################
###############################################################################

threshold <- 1.5 

# visualize threshold over light levels  
thresholdOverLight(lig, threshold, span =c(30000, 35000))

# plot light levels over the deployment period 
offset <- 12 # adjusts the y-axis to put night (dark shades) in the middle

lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

#Detect twilight times 
twl <- preprocessLight(lig, 
                       threshold = threshold,
                       offset = offset, 
                       lmax = 64,         # max. light value
                       gr.Device = "x11") # MacOS version (and windows)

head(twl)
