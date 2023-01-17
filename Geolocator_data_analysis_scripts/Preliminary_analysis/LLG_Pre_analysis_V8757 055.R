# source: Unpublished data 
# tag number: V8757 055
# site: Northwestern territories

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
lat.calib <- 62.475799
lon.calib <- -114.699932

# data directory
dir <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data"

# time of deployment
deploy <- anytime("2019-06-27", tz = "GMT")

###############################################################################
#DATA EXTRACTION ##############################################################
###############################################################################

# import lig data 
lig <- readLig(paste0(dir,"/raw_data/BLPW_geo_2020/breed/ML6740 V8757 055/ML6740 V8757 055 reconstructed_000.lig"), skip = 1)

#adjust year
lig$Date <- lig$Date %m+% years(7)

#remove rows before and after deployment time 
lig <- lig[(lig$Date > deploy),]

#adjust time 
lig$Date <- lig$Date + 5.5 *60*60

###############################################################################
#TWILIGHT ANNOTATION ##########################################################
###############################################################################

threshold <- 1.5 

# visualize threshold over light levels  
thresholdOverLight(lig, threshold, span =c(30000, 35000))

# plot light levels 
offset <- 20 # adjusts the y-axis to put night (dark shades) in the middle

lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 64))

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
                    window = 6,           
                    outlier.mins = 90,    
                    stationary.mins = 45, 
                    plot = TRUE)

# Visualize light and twilight time-series
lightImage(lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

# Save the twilight times 
write.csv(twl, paste0(dir,"/intermediate_data/Twilight_times/Pre_analysis_ML6740_V8757_055_twl_times.csv"))

###############################################################################
# SGAT ANALYSIS ###############################################################
###############################################################################

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/intermediate_data/Twilight_times/Pre_analysis_ML6740_V8757_055_twl_times.csv"))
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

#calibration period before the migration 
tm.calib <- as.POSIXct(c("2019-08-10", "2019-09-01"), tz = "UTC")

abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

# subset of twilight times that are within the calibration period
d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])

# perform the calibration and verify the fit with the gamma or log distribution
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

#parameters of the error distribution 
zenith  <- calib[1] 
zenith0 <- calib[2]

alpha <- calib[3:4]

# Movement model ###############################################################

#this movement model should be based on the estimated migration speed of the blackpoll warbler 
beta  <- c(2.2, 0.08)
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

# Initial Path #################################################################
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol=0.01)

x0 <- path$x
z0 <- trackMidpts(x0)

data(wrld_simpl)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(path$x, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()
