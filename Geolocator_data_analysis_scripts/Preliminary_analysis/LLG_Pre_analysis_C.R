# source: Deluca et al. 2015 
# tag number: C
# site: Bon Portage Island, Nova Scotia, Canada

#load packages
require(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(remotes)
library(anytime)
library(lubridate)
library(parallel)

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

geo.id <- "C"

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# geo deployment location 
lat.calib <- 49.45833
lon.calib <- -56.225

# time of deployment
deploy.start <- anytime("2013-06-16", tz = "GMT")
deploy.end <- anytime("2014-04-15", tz = "GMT")

#Equinox times
fall.equi <- anytime("2013-09-22", asUTC = T, tz = "GMT")
spring.equi <- anytime("2014-03-20", asUTC = T, tz = "GMT")

#Find number of cores available for analysis
Threads= detectCores()-1

###############################################################################
#DATA EXTRACTION ##############################################################
###############################################################################

# import light data 
lig <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Movebank_data/Blackpoll Warbler eastern North America (data from DeLuca et al. 2015).csv")

#remove rows with processed data
lig <- lig[(is.na(lig$comments) == TRUE),]  %>%
  #get data for individual C
  filter(individual.local.identifier == "C") %>%
  # rename columns  
  rename(c( Date = timestamp, Light = gls.light.level)) %>%
  #remove rows before and after deployment time 
  filter(Date > deploy.start) %>%
  #Convert light levels to log
  mutate(Light = log(Light + 0.0001) + abs(min(log(Light+0.0001)))) 

#convert Dates to as.POXIct format
lig$Date <- anytime(lig$Date, tz = "UTC")

#parameter to visualize the data 
offset <- 18 # adjusts the y-axis to put night (dark shades) in the middle

#Threshold light level 
threshold <- 1.5 

# visualize threshold over light levels  
thresholdOverLight(lig, threshold, span =c(0, 25000))

# FIND TIME SHIFT ##############################################################

#This geolocator has a a time shift visible on this plot
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 2))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# we will do an initial twilight annotation to find identify the time interval
# by which we need to shift time
# There should be not need to edit, delete or insert twilights for this
twl_in <- preprocessLight(lig,
                          threshold = threshold,
                          offset = offset,
                          lmax = 2,         # max. light value
                          gr.Device = "x11", # MacOS version (and windows)
                          dark.min = 60)

#write.csv(twl_in, paste0(dir,"/Pre_analysis_C_twl_times_initial.csv"))
twl_in <- read.csv(paste0(dir,"/Pre_analysis_C_twl_times_initial.csv"))
twl_in$Twilight <- as.POSIXct(twl_in$Twilight, tz = "UTC")

# Period over which to calculate the time shift. It should be while the bird is 
# still in the breeding grounds 
period <- as.POSIXct(c("2013-06-28", "2013-07-19"), tz = "UTC")

#plot the period over the light image 
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 2))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

abline(v = period, lwd = 2, lty = 2, col = "orange")

# calculate the time shift
shift <- shiftSpan(twl = twl_in, lig = lig, period = period, est.zenith = 92,
                   dep.lon = lon.calib,
                   dep.lat = lat.calib)

# verify the output 
# shift 

#adjust time 
lig$Date <- lig$Date - (shift$shift)

#TWILIGHT ANNOTATION ##########################################################

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
lightImage(lig, offset = offset)
tsimagePoints(twl$Twilight, offset = offset, pch = 16, cex = 0.5,
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
tm.calib <- as.POSIXct(c(deploy.start + days(1), deploy.start + days(14)), tz = "UTC")

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

# adjust the zenith angles calculated from the breeding sites 
zenith0_ad <- zenith0 + abs(zenith - zenith_sd)
zenith_ad  <- zenith_sd

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



twl_in2 <- twl_in[order(date(twl_in$Twilight), -twl_in$Rise),]
