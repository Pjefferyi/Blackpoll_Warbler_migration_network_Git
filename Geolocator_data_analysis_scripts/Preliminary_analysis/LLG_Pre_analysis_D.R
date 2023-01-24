# source: Deluca et al. 2015 
# tag number: D
# site: Seal Island, Nova Scotia, Canada

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

geo.id <- "D"

# geo deployment location 
lat.calib <- 43.408825
lon.calib <- -66.014189

# time of deployment
deploy.start <- anytime("2013-08-25", tz = "GMT")
deploy.end <- anytime("2014-04-15", tz = "GMT")

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

###############################################################################
#DATA EXTRACTION ##############################################################
###############################################################################

# import light data 
lig <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Movebank_data/Blackpoll Warbler eastern North America (data from DeLuca et al. 2015).csv")

#remove rows with processed data
lig <- lig[(is.na(lig$comments) == TRUE),]  %>%
  #get data for individual D
    filter(individual.local.identifier == "D") %>%
  # rename columns  
    rename(c( "timestamp" = "Date", "gls.light.level" = "Light")) %>%
  #remove rows before and after deployment time 
    filter(Date > deploy.start) %>%
  #Convert light levels to log
    mutate(Light = log(Light + 0.0001) + abs(min(log(Light+0.0001)))) 

#convert Dates to as.POXIct format
lig$Date <- anytime(lig$Date, tz = "UTC")

#adjust time 
lig$Date <- lig$Date - 4*60*60

###############################################################################
#TWILIGHT ANNOTATION ##########################################################
###############################################################################

threshold <- 1

# visualize threshold over light levels  
thresholdOverLight(lig, threshold, span =c(2000, 2500))

# plot light levels over the deployment period 
offset <- 18 # adjusts the y-axis to put night (dark shades) in the middle

# open jpeg
jpeg(paste0(dir, "/D_light_plot.png"), width = 1024, height = 990)

lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 2))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

dev.off()

#Detect twilight times, for now do not edit twilight times  
twl <- preprocessLight(lig, 
                       threshold = threshold,
                       offset = offset, 
                       lmax = 2,         # max. light value
                       gr.Device = "x11", # MacOS version (and windows)
                       dark.min = 60)


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

# Save the twilight times 
#write.csv(twl, paste0(dir,"/Pre_analysis_D_twl_times.csv"))

###############################################################################
# SGAT ANALYSIS ###############################################################
###############################################################################

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/Pre_analysis_D_twl_times.csv"))
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage(tagdata = lig,
            offset = offset,     
            zlim = c(0, 2))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

#calibration period before the migration 
tm.calib <- as.POSIXct(c(deploy.start, deploy.start + days(14)), tz = "UTC")

abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

# subset of twilight times that are within the calibration period
d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])

# perform the calibration and verify the fit with the gamma or log distribution
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

#parameters of the error distribution 
zenith  <- calib[1] 
zenith0 <- calib[2]

alpha <- calib[3:4]

# Alternative calibration ######################################################

startDate <- "2013-11-20"
endDate   <- "2014-02-10"

start = min(which(as.Date(twl$Twilight) == startDate))
end = max(which(as.Date(twl$Twilight) == endDate))

(zenith_sd <- findHEZenith(twl, tol=0.01, range=c(start,end)))

# Movement model ###############################################################

#this movement model should be based on the estimated migration speed of the blackpoll warbler 
beta  <- c(2.2, 0.08)
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

# Initial Path #################################################################
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith_sd, tol=0.01)

x0 <- path$x
z0 <- trackMidpts(x0)

# open jpeg
jpeg(paste0(dir, "/D_Threshold_path.png"), width = 1024, height = 990)

data(wrld_simpl)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(path$x, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()

dev.off()
