# source: Deluca et al. 2019
# tag number: 3254-011
# site: Nome Alaska

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

geo.id <- "3254-011"

# geo deployment location 
lat.calib <- 64.86398
lon.calib <- -163.69373

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# time of deployment
deploy.start <- anytime("2016-06-18", tz = "GMT")

###############################################################################
#DATA EXTRACTION ##############################################################
###############################################################################

# import lig data 
lig <- readLig(paste0(dir,"/ML6440 V3254 011 reconstructed_000.lig"), skip = 1)

#remove rows before and after deployment time 
lig <- lig[(lig$Date > deploy.start),]

lig <- lig[(is.na(lig$Light) == FALSE),]

###############################################################################
#TWILIGHT ANNOTATION ##########################################################
###############################################################################

threshold <- 1.5 

# visualize threshold over light levels  
thresholdOverLight(lig, threshold, span =c(30000, 35000))

# plot light levels 
offset <- 20 # adjusts the y-axis to put night (dark shades) in the middle

# open jpeg
jpeg(paste0(dir, "/3254-011_light_plot.png"), width = 1024, height = 990)

lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 64))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

dev.off()

#Detect twilight times, for now do not edit twilight times  
#twl <- preprocessLight(lig, 
#                       threshold = threshold,
#                       offset = offset, 
#                       lmax = 64,         # max. light value
#                       gr.Device = "x11", # MacOS version (and windows)
#                       dark.min = 60)

# Adjust sunset times by 120 second sampling interval
#twl <- twilightAdjust(twilights = twl, interval = 120)

# Automatically adjust or mark false twilights 
#twl <- twilightEdit(twilights = twl, 
#                    window = 6,           
#                    outlier.mins = 90,    
#                    stationary.mins = 45, 
#                    plot = TRUE)

# Visualize light and twilight time-series
#lightImage(lig, offset = 19)
#tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
#              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

# Save the twilight times 
#write.csv(twl, paste0(dir,"/Pre_analysisis_ML6440_V3254_011_twl_times.csv"))

###############################################################################
# SGAT ANALYSIS ###############################################################
###############################################################################

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/Pre_analysis_V3254_011_twl_times.csv"))
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

#calibration period before the migration 
tm.calib <- as.POSIXct(c("2016-08-01", "2016-08-10"), tz = "UTC")

abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

# subset of twilight times that are within the calibration period
d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])

# perform the calibration and verify the fit with the gamma or log distribution
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

# In-habitat calibration not possible during the breeding period due to shading/ incomplete darkness (high latitude)

# Alternative calibration approach =============================================

startDate <- "2016-11-01"
endDate   <- "2017-01-12"

start = min(which(as.Date(twl$Twilight) == startDate))
end = max(which(as.Date(twl$Twilight) == endDate))

(zenith_sd <- findHEZenith(twl, tol=0.01, range=c(start,end)))

#convert to geolight format
geo_twl <- export2GeoLight(twl)

# this is just to find places where birds have been for a long time, would not use these parameters for stopover identification, detailed can be found in grouped model section
cL <- changeLight(twl=geo_twl, quantile=0.8, summary = F, days = 10, plot = T)
# merge site helps to put sites together that are separated by single outliers.
mS <- mergeSites(twl = geo_twl, site = cL$site, degElevation = 90-zenith0, distThreshold = 500)

#specify which site is the stationary one
site           <- mS$site[mS$site>0] # get rid of movement periods
stationarySite <- which(table(site) == max(table(site))) # find the site where bird is the longest

#find the dates that the bird arrives and leaves this stationary site
start <- min(which(mS$site == stationarySite))
end   <- max(which(mS$site == stationarySite))

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
jpeg(paste0(dir, "/3254_011_Threshold_path.png"), width = 1024, height = 990)

# plot of threshold path  
data(wrld_simpl)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(path$x, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()

dev.off()

