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
lon.calib <- -65.7499


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
lig$Date <- anytime(lig$Date, tz = "UTC", asUTC = T)

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
period <- as.POSIXct(c("2013-06-16", "2013-06-30"), tz = "UTC")

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
lightImage(lig, offset = offset, zlim = c(0, 2))
tsimagePoints(twl$Twilight, offset = offset, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

# Save the twilight times 
# write.csv(twl, paste0(dir,"/Pre_analysis_C_twl_times.csv"))

###############################################################################
# SGAT ANALYSIS ###############################################################
###############################################################################

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/Pre_analysis_C_twl_times.csv"))
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage(tagdata = lig,
           offset = offset,     
           zlim = c(0, 2))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

#calibration period before the migration 
tm.calib <- as.POSIXct(c(deploy.start + days(1), deploy.start + days(40)), tz = "UTC")

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
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith_sd, tol=0.18)

x0 <- path$x
z0 <- trackMidpts(x0)

# open jpeg
jpeg(paste0(dir, "/C_Threshold_path.png"), width = 1024, height = 990)

data(wrld_simpl)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(path$x, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()

dev.off()

# Define known locations #######################################################

#we set the location of geolocator deploymentlocation for the MCMC sampler (there is no recovery location) 

fixedx <- rep(F, nrow(x0))
fixedx[1:2] <- T # first two location estimates

x0[fixedx, 1] <- lon.calib
x0[fixedx, 2] <- lat.calib

z0 <- trackMidpts(x0) # we need to update the z0 locations

# Land mask ####################################################################

earthseaMask <- function(xlim, ylim, n = 2, pacific=FALSE) {
  
  if (pacific) { wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1 and sea NA
  mask = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
               rasterize(wrld_simpl, r, 1, silent = TRUE), 
               rasterize(elide(wrld_simpl, shift = c(360, 0)), r, 1, silent = TRUE))
  
  #load polygon of blackpoll's range
  #load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_Full_blackpoll_range_polygon.R.R")
  load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Birdlife_int_Full_blackpoll_range_polygon.R")
  
  #rasterize the polygon 
  #range.raster <- rasterize(range.poly, mask)
  range.raster <- rasterize(BLI.range.poly, mask)
  
  #Update the land mask 
  mask <- range.raster * mask
  
  xbin = seq(xmin(mask),xmax(mask),length=ncol(mask)+1)
  ybin = seq(ymin(mask),ymax(mask),length=nrow(mask)+1)
  
  function(p) mask[cbind(length(ybin) -.bincode(p[,2],ybin),.bincode(p[,1],xbin))]
  
}

xlim <- range(x0[,1]+c(-5,5))
ylim <- range(x0[,2]+c(-5,5))

mask <- earthseaMask(xlim, ylim, n = 4)

## Define the log prior for x and z
log.prior <- function(p) {
  f <- mask(p)
  #ifelse(is.na(f), log(1), f)  # if f is the relative abundance within a grid square 
  ifelse(is.na(f), log(1), log(2)) # if f indicates the distribution of the blackpoll warbler 
  
}

# Run the Estelle model ########################################################

#Define the model
model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "ModifiedGamma",
                        alpha = alpha,
                        beta = beta,
                        logp.x = log.prior, logp.z = log.prior, 
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = fixedx)

#Define the error distribution around each location 
proposal.x <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0))

fit <- estelleMetropolis(model, proposal.x, proposal.z, iters = 1000, thin = 20)

# We tune the proposals 
x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "ModifiedGamma",
                        alpha = alpha,
                        beta = beta,
                        logp.x = log.prior, logp.z = log.prior, 
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = fixedx)

x.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl))
z.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl) - 1)

# Fit multiple runs to tune the proposals
for (k in 1:3) {
  fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x), 
                           z0 = chainLast(fit$z), iters = 300, thin = 20)
  
  x.proposal <- mvnorm(chainCov(fit$x), s = 0.2)
  z.proposal <- mvnorm(chainCov(fit$z), s = 0.2)
}

# Check that the chain is well mixed 
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!fixedx, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!fixedx, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)

#Final Run 
x.proposal <- mvnorm(chainCov(fit$x), s = 0.25)
z.proposal <- mvnorm(chainCov(fit$z), s = 0.25)

fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x), 
                         z0 = chainLast(fit$z), iters = 1000, thin = 20)

#Summarize the results
sm <- locationSummary(fit$z, time=fit$model$time)
head(sm)

#Save the output of the estelle model 
#save(sm, file = paste0(dir,"/Pre_analysis_C_SGAT_estelle_summary.csv"))
#save(fit, file = paste0(dir,"/Pre_analysis_C_SGAT_estelle_fit.R"))

#load the output of the estelle model 
#load(sm, file = paste0(dir,"/Pre_analysis_C_SGAT_estelle_summary.csv"))
#load(fit, file = paste0(dir,"/Pre_analysis_C_SGAT_estelle_fit.R"))

# open jpeg
jpeg(paste0(dir, "/C__Estelle_path.png"), width = 1024 , height = 990)

#Plot the results
par(mfrow=c(1,1))
# empty raster of the extent
r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = proj4string(wrld_simpl))

s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- slice(s, sliceIndices(s))

plot(sk, useRaster = F,col = rev(viridis::viridis(50)))
plot(wrld_simpl, xlim=xlim, ylim=ylim,add = T, bg = adjustcolor("black",alpha=0.1))

#plot location track. Locations in blue occured during the fall equinox 
lines(sm[,"Lon.50%"], sm[,"Lat.50%"], 
      col = ifelse(sm$Time1 > spring.equi - days(15) & sm$Time1 < spring.equi + days(15), adjustcolor("blue", alpha.f = 0.6), adjustcolor("firebrick", alpha.f = 0.6)),
      type = "o", pch = 16)

#close jpeg
dev.off()

# Plot of mean longitude and latitude

# open jpeg
jpeg(paste0(dir, "3254_011__mean_lon_lat.png"), width = 1024 , height = 990)

par(mfrow=c(2,1),mar=c(4,4,1,1))

plot(sm$Time1, sm$"Lon.50%", ylab = "Longitude", xlab = "", yaxt = "n", type = "n", ylim = c(min(sm$Lon.mean) - 10, max(sm$Lon.mean) + 10))
axis(2, las = 2)
polygon(x=c(sm$Time1,rev(sm$Time1)), y=c(sm$`Lon.2.5%`,rev(sm$`Lon.97.5%`)), border="gray", col="gray")
lines(sm$Time1,sm$"Lon.50%", lwd = 2)
abline(v = fall.equi, lwd = 2, lty = 2, col = "orange")
abline(v = spring.equi, lwd = 2, lty = 2, col = "orange")

plot(sm$Time1,sm$"Lat.50%", type="n", ylab = "Latitude", xlab = "", yaxt = "n", ylim = c(min(sm$Lat.mean) - 10, max(sm$Lat.mean) + 10))
axis(2, las = 2)
polygon(x=c(sm$Time1,rev(sm$Time1)), y=c(sm$`Lat.2.5%`,rev(sm$`Lat.97.5%`)), border="gray", col="gray")
lines(sm$Time1, sm$"Lat.50%", lwd = 2)
abline(v = fall.equi, lwd = 2, lty = 2, col = "orange")
abline(v = spring.equi, lwd = 2, lty = 2, col = "orange")

#close jpeg
dev.off()

# Identify stopover areas using median longitude and latitude
sm <- sm %>% mutate(stationary = ifelse(abs(lead(Lon.mean) - Lon.mean) < 2 & abs(lead(Lat.mean) - Lat.mean) < 2, 1, 0)) 

par(mfrow=c(2,1))

plot(sm$Time1, sm$"Lon.50%", ylab = "Longitude", xlab = "", yaxt = "n", type = "n", ylim = c(min(sm$Lon.mean) - 10, max(sm$Lon.mean) + 10))
axis(2, las = 2)
polygon(x=c(sm$Time1,rev(sm$Time1)), y=c(sm$`Lon.2.5%`,rev(sm$`Lon.97.5%`)), border="gray", col="gray")
lines(sm$Time1,sm$"Lon.50%", lwd = 2)
abline(v = fall.equi, lwd = 2, lty = 2, col = "orange")
abline(v = spring.equi, lwd = 2, lty = 2, col = "orange")
points(sm$Time1, sm$"Lon.50%", col = ifelse(sm$stationary == 1, "blue", "red"), cex = 1.2)
grid()

plot(sm$Time1,sm$"Lat.50%", type="n", ylab = "Latitude", xlab = "", yaxt = "n", ylim = c(min(sm$Lat.mean) - 10, max(sm$Lat.mean) + 10))
axis(2, las = 2)
polygon(x=c(sm$Time1,rev(sm$Time1)), y=c(sm$`Lat.2.5%`,rev(sm$`Lat.97.5%`)), border="gray", col="gray")
lines(sm$Time1, sm$"Lat.50%", lwd = 2)
abline(v = fall.equi, lwd = 2, lty = 2, col = "orange")
abline(v = spring.equi, lwd = 2, lty = 2, col = "orange")
points(sm$Time1, sm$"Lat.50%", col = ifelse(sm$stationary == 1, "blue", "red"), cex = 1.2)
grid()

par(mfrow=c(1,1))
# empty raster of the extent
r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = proj4string(wrld_simpl))

s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- slice(s, sliceIndices(s))

plot(sk, useRaster = F,col = rev(viridis::viridis(50)))
plot(wrld_simpl, xlim=xlim, ylim=ylim,add = T, bg = adjustcolor("black",alpha=0.1))

lines(sm[,"Lon.50%"], sm[,"Lat.50%"], 
      col = ifelse(sm$stationary == 1, "blue", "red"),
      type = "o", pch = 16)




