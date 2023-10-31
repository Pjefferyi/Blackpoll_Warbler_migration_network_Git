# source: unpublished data 
# tag number: V7638-010
# site: Denali, Alaska

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

# install and load geolocator packages 
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

# clear object from workspace
rm(list=ls())

# Load helper functions 
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis/Geolocator_analysis_helper_functions.R")

geo.id <- "V7638_010"

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# read file with consolidated geolocator data
ref_data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")

#Assign geolocator deployment site 
lat.calib <- ref_data$deploy.latitude[which(ref_data$geo.id == geo.id)]
lon.calib <- ref_data$deploy.longitude[which(ref_data$geo.id == geo.id)]

# time of deployment (from reference file)
deploy.start <- anytime(ref_data$deploy.on.date[which(ref_data$geo.id == geo.id)], asUTC = T, tz = "GMT")

# time of recovery (estimate from light data)
deploy.end <- anytime(ref_data$deploy.off.date[which(ref_data$geo.id == geo.id)], asUTC = T, tz = "GMT")

#Equinox times
fall.equi <- anytime("2018-09-22", asUTC = T, tz = "GMT")
spring.equi <- anytime("2018-03-20", asUTC = T, tz = "GMT")

#DATA EXTRACTION ##############################################################

# import lig data 
lig <- readLig(paste0(dir,"/Raw_light_data_", geo.id, ".lig"), skip = 1)

#remove rows before and after deployment time 
lig <- lig[(lig$Date > deploy.start),]
lig <- lig[(lig$Date < deploy.end),]

#parameter to visualize the data 
offset <- -4# adjusts the y-axis to put night (dark shades) in the middle

#Threshold light level 
threshold <- 1.5 

# visualize threshold over light levels  
thresholdOverLight(lig, threshold, span =c(45000, 55000))

# FIX CLOCK DRIFT ##############################################################

# This geolocator has no visible time shift, but it has visisble clock drift 
lightImage( tagdata = lig,
            offset = offset,
            zlim = c(0, 64))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# we will do an initial twilight annotation to find identify the time interval
# by which we need to shift time
# There should be not need to edit, delete or insert twilights for this
twl_in <- preprocessLight(lig,
                          threshold = threshold,
                          offset = offset,
                          lmax = 64,         # max. light value
                          gr.Device = "x11", # MacOS version (and windows)
                          dark.min = 60)

#write.csv(twl_in, paste0(dir,"/", geo.id,"_twl_times_initial.csv"))
twl_in <- read.csv(paste0(dir,"/", geo.id,"_twl_times_initial.csv"))
twl_in$Twilight <- as.POSIXct(twl_in$Twilight, tz = "UTC")

# calculate the difference between noon before the departure frot he nonbreeding grounds and after the return in the spring 

# period prior to the departure  
premig.period <- as.POSIXct(c("2018-08-01", "2018-08-30"), tz = "UTC")
postmig.period <- as.POSIXct(c("2019-06-05", "2019-06-25"), tz = "UTC")

# plot the period over the light image
lightImage( tagdata = lig,
            offset = offset,
            zlim = c(0, 64))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

abline(v = premig.period, lwd = 2, lty = 2, col = "orange")
abline(v = postmig.period, lwd = 2, lty = 2, col = "orange")

# calculate the time shift
shift <- shiftSpan(twl = twl_in, lig = lig, period = period, est.zenith = 92,
                   dep.lon = lon.calib,
                   dep.lat = lat.calib)


# verify the that the time shift measured makes sense
shift

#adjust  the based on the time shift detected
#lig$Date <- lig$Date - (shift$shift)
#lig$Date <- lig$Date + 5.5*60*60


#TWILIGHT ANNOTATION ##########################################################

#open jpeg
jpeg(paste0(dir, "/", geo.id, "_light_plot.png"), width = 1024, height = 990)

lightImage( tagdata = lig,
            offset = offset,
            zlim = c(0, 20))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

abline( h = 10, col = "red")

dev.off()

# Detect twilight times, for now do not edit twilight times
twl <- preprocessLight(lig,
                       threshold = threshold,
                       offset = offset,
                       lmax = 64,         # max. light value
                       gr.Device = "x11", # MacOS version (and windows)
                       dark.min = 120)

# Adjust sunset times by 120 second sampling interval
twl <- twilightAdjust(twilights = twl, interval = 120)

# Save the twilight times
# write.csv(twl, paste0(dir,"/",geo.id , "_twl_times.csv"))

###############################################################################
# SGAT ANALYSIS ###############################################################
###############################################################################

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/", geo.id, "_twl_times.csv"))
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

# Automatically adjust or mark false twilights
twl <- twilightEdit(twilights = twl,
                    window = 4,
                    outlier.mins = 35,
                    stationary.mins = 25,
                    plot = TRUE)

# Visualize light and twilight time-series
lightImage(lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

# calibration is only possible before the migration due to the high latitude of this breeding site  
tm.calib1 <- as.POSIXct(c("2018-07-25", "2018-08-22"), tz = "UTC")

abline(v = tm.calib1, lwd = 2, lty = 2, col = "orange")

# subset of twilight times that are within the calibration period
d_calib <- subset(twl, Twilight>=tm.calib1[1] & Twilight<=tm.calib1[2])

# perform the calibration and verify the fit with the gamma or log distribution
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

#parameters of the error distribution 
zenith  <- calib[1] 
zenith0 <- calib[2]

alpha <- calib[3:4]

# Alternative calibration #######################################################

#convert to geolight format
geo_twl <- export2GeoLight(twl)

# this is just to find places where birds have been for a long time, would not use these parameters for stopover identification, detailed can be found in grouped model section
cL <- changeLight(twl=geo_twl, quantile=0.9, summary = F, days = 10, plot = T)
# merge site helps to put sites together that are separated by single outliers.
mS <- mergeSites(twl = geo_twl, site = cL$site, degElevation = 90-zenith, distThreshold = 1000)

#specify which site is the stationary one
site           <- mS$site[mS$site>0] # get rid of movement periods
stationarySite <- which(table(site) == max(table(site))) # find the site where bird is the longest

#find the dates that the bird arrives and leaves this stationary site
start <- min(which(mS$site == stationarySite))
end   <- max(which(mS$site == stationarySite))

(zenith_sd <- findHEZenith(twl, tol=0.01, range=c(start,end)))

# the Hill-ekstrom zenith and in-habitat zenith angles differ by more than 0.5 degrees

# adjust the zenith angles calculated from the breeding sites for the non-breeding sites
zenith0_ad <- zenith0 + abs(zenith - zenith_sd)
zenith_ad  <- zenith_sd

# Find approximate  timing of arrival and departure from the nonbreeding grounds 
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol=0)

x0_r<- path$x
z0 <- trackMidpts(x0_r)

#Save raw path (no linear interpolation around the equinox)
save(x0_r, file = paste0(dir,"/", geo.id, "_initial_path_raw.csv"))

# Check the following times of arrival and departure using a plot 
arr.nbr <- "2018-10-16" 
dep.nbr <- "2019-04-29" 

# open jpeg
jpeg(paste0(dir, "/", geo.id, "_LatLon_scatterplot.png"), width = 1024, height = 990)

par(mfrow = c(2,1))
plot(twl$Twilight, x0_r[,1], ylab = "longitude")
abline(v = anytime(arr.nbr))
abline(v = anytime(dep.nbr))
abline(v = fall.equi, col = "orange")
abline(v = spring.equi, col = "orange")
plot(twl$Twilight, x0_r[,2], ylab = "latitude")
abline(v = anytime(arr.nbr))
abline(v = anytime(dep.nbr))
abline(v = fall.equi, col = "orange")
abline(v = spring.equi, col = "orange")

dev.off()

# Using approximate timings of arrival and departure from the breeding grounds
zenith_twl_zero <- data.frame(Date = twl$Twilight) %>%
  mutate(zenith = case_when(Date < anytime(arr.nbr) ~ zenith0,
                            Date > anytime(arr.nbr) & Date < anytime(dep.nbr) ~ zenith0_ad,
                            Date > anytime(dep.nbr) ~ zenith0_ad))

zeniths0 <- zenith_twl_zero$zenith

zenith_twl_med <- data.frame(Date = twl$Twilight) %>%
  mutate(zenith = case_when(Date < anytime(arr.nbr) ~ zenith,
                            Date > anytime(arr.nbr) & Date < anytime(dep.nbr) ~ zenith_sd,
                            Date > anytime(dep.nbr) ~ zenith_sd))

zeniths_med <- zenith_twl_med$zenith

# plot longitudes and latitudes with the new zenith angles 
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zeniths_med, tol= 0)

x0_ad <- path$x
z0 <- trackMidpts(x0_ad)

#Save raw path (no linear interpolation around the equinox)
save(x0_ad, file = paste0(dir,"/", geo.id, "adjusted_initial_path_raw.csv"))

# open jpeg
jpeg(paste0(dir, "/", geo.id, "_LatLon_scatterplot_adjusted.png"), width = 1024, height = 990)

par(mfrow = c(2,1))
plot(twl$Twilight, x0_ad[,1], ylab = "longitude")
abline(v = anytime(arr.nbr))
abline(v = anytime(dep.nbr))
abline(v = fall.equi, col = "orange")
abline(v = spring.equi, col = "orange")
plot(twl$Twilight, x0_ad[,2], ylab = "latitude")
abline(v = anytime(arr.nbr))
abline(v = anytime(dep.nbr))
abline(v = fall.equi, col = "orange")
abline(v = spring.equi, col = "orange")

dev.off()

# Initial Path #################################################################
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zeniths_med, tol=0.11)

x0 <- path$x
z0 <- trackMidpts(x0)

# open jpeg
jpeg(paste0(dir, "/", geo.id, "_Threshold_path.png"), width = 1024, height = 990)

par(mfrow=c(1,1))
data(wrld_simpl)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(path$x[], pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()

dev.off()

#Save initial path 
save(x0, file = paste0(dir,"/", geo.id, "_initial_path.csv"))