# source: Deluca et al. 2019
# tag number: WRMA04173
# site: Colombia

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

geo.id <- "WRMA04173"

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# read file with consolidated geolocator data
ref_data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")

#Assign geolocator deployment site 
lat.calib <- ref_data$deploy.latitude[which(ref_data$geo.id == geo.id)]
lon.calib <- ref_data$deploy.longitude[which(ref_data$geo.id == geo.id)]

# time of deployment (from reference file)
#deploy.start <- anytime("", asUTC = T, tz = "GMT")

# time of recovery (estimate from light data)
#deploy.end <- anytime("", asUTC = T, tz = "GMT")

#Equinox times
spring.equi <- anytime("2020-03-19", asUTC = T, tz = "GMT")

#DATA EXTRACTION ##############################################################

# import lig data 
lig <- readLig(paste0(dir,"/Raw_light_data_", geo.id, ".lig"), skip = 1)

#remove rows before and after deployment time 
#lig <- lig[(lig$Date > deploy.start),]

#parameter to visualize the data 
offset <- 20 # adjusts the y-axis to put night (dark shades) in the middle

#Threshold light level 
threshold <- 1.5 

# visualize threshold over light levels  
thresholdOverLight(lig, threshold, span =c(0, 25000))

# FIND TIME SHIFT ##############################################################

#This geolocator has no time shift 

#TWILIGHT ANNOTATION ##########################################################

#open jpeg
jpeg(paste0(dir, "/", geo.id, "_light_plot.png"), width = 1024, height = 990)

lightImage( tagdata = lig,
            offset = offset,
            zlim = c(0, 20))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))
dev.off()

# #Detect twilight times, for now do not edit twilight times
# twl <- preprocessLight(lig,
#                        threshold = threshold,
#                        offset = offset,
#                        lmax = 64,         # max. light value
#                        gr.Device = "x11", # MacOS version (and windows)
#                        dark.min = 60)
# 
# # Adjust sunset times by 120 second sampling interval
# twl <- twilightAdjust(twilights = twl, interval = 120)
# 
# # Visualize light and twilight time-series
# lightImage(lig, offset = 19)
# tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
#               col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

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

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

#calibration period before the migration 
tm.calib <- as.POSIXct(c("2020-03-15", "2020-04-15"), tz = "UTC")

abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

# subset of twilight times that are within the calibration period
d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])

# perform the calibration and verify the fit with the gamma or log distribution
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

#parameters of the error distribution 
zenith  <- calib[1] 
zenith0 <- calib[2]

alpha <- calib[3:4]

# Alternative calibration #######################################################

# Alternative Calibration is not necessary here 

# Find approximate  timing of arrival and departure from the nonbreeding grounds 
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol= 0)

x0_r<- path$x
z0 <- trackMidpts(x0_r)

#Save raw path (no linear interpolation around the equinox)
save(x0_r, file = paste0(dir,"/", geo.id, "_initial_path_raw.csv"))

# Check the following times of arrival and departure using a plot 
dep.nbr <- "2020-06-01" 

# open jpeg
jpeg(paste0(dir, "/", geo.id, "_LatLon_scatterplot.png"), width = 1024, height = 990)

par(mfrow = c(2,1))
plot(twl$Twilight, x0_r[,1], ylab = "longitude")
abline(v = anytime(dep.nbr))
plot(twl$Twilight, x0_r[,2], ylab = "latitude")
abline(v = anytime(dep.nbr))
abline(v = spring.equi, col = "orange")

dev.off()

# Movement model ###############################################################

#this movement model should be based on the estimated migration speed of the blackpoll warbler 
beta  <- c(0.7, 0.05)
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

# Initial Path #################################################################
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol=0.05)

#Adjusted tol until second stopover was located over North Carolina rather than further South. 
x0 <- path$x
z0 <- trackMidpts(x0)

# open jpeg
jpeg(paste0(dir, "/", geo.id, "_Threshold_path.png"), width = 1024, height = 990)

par(mfrow=c(1,1))
data(wrld_simpl)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(path$x, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()

dev.off()

#Save initial path 
save(x0, file = paste0(dir,"/", geo.id, "_initial_path.csv"))

# Define known locations #######################################################

#we set the location of geolocator deployment and recovery as fixed locations for the MCMC sampler 

fixedx <- rep(F, nrow(x0))
fixedx[1:2] <- T # first two location estimates

# We only fix the fist location because this track is incomplete

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
  load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Birdlife_int_Full_blackpoll_range_polygon.R")
  #load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_Full_blackpoll_range_polygon.R")
  
  #rasterize the polygon 
  range.raster <- rasterize(BLI.range.poly, mask)
  #range.raster <- rasterize(range.poly, mask)
  
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
  #ifelse(f | is.na(f), log(2), log(1)) #original function from Lisovski et al. 2020
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
#save(sm, file = paste0(dir,"/", geo.id, "_SGAT_estelle_summary.csv"))
#save(fit, file = paste0(dir,"/", geo.id, "_SGAT_estelle_fit.R"))

# Plot Results #################################################################

# open jpeg
jpeg(paste0(dir, "/", geo.id,"_Estelle_path.png"), width = 1024 , height = 990)

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
      col = ifelse(sm$Time1 > fall.equi - days(10) & sm$Time1 < fall.equi + days(10), adjustcolor("blue", alpha.f = 0.6), adjustcolor("firebrick", alpha.f = 0.6)),
      type = "o", pch = 16)

#close jpeg
dev.off()

# Plot of mean longitude and latitude

# open jpeg
jpeg(paste0(dir, "/", geo.id,"_mean_lon_lat.png"), width = 1024 , height = 990)

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

#SGAT Group model analysis ####################################################

# group twilight times were birds were stationary 
geo_twl <- export2GeoLight(twl)

# Often it is necessary to play around with quantile and days
# quantile defines how many stopovers there are. the higher, the fewer there are
# days indicates the duration of the stopovers 
cL <- changeLight(twl=geo_twl, quantile=0.9, summary = F, days = 2, plot = T)

# merge site helps to put sites together that are separated by single outliers.
mS <- mergeSites(twl = geo_twl, site = cL$site, degElevation = 90-zenith, distThreshold = 500)

##back transfer the twilight table and create a group vector with TRUE or FALSE according to which twilights to merge 
twl.rev <- data.frame(Twilight = as.POSIXct(geo_twl[,1], geo_twl[,2]), 
                      Rise     = c(ifelse(geo_twl[,3]==1, TRUE, FALSE), ifelse(geo_twl[,3]==1, FALSE, TRUE)),
                      Site     = rep(mS$site,2))
twl.rev <- subset(twl.rev, !duplicated(Twilight), sort = Twilight)

grouped <- rep(FALSE, nrow(twl.rev))
grouped[twl.rev$Site>0] <- TRUE 
grouped[c(1:3, (length(grouped)-2):length(grouped))] <- TRUE

# Create a vector which indicates which numbers sites
g <- makeGroups(grouped)

# Add data to twl file
twl$group <- c(g, g[length(g)])

# Add behavior vector
behaviour <- c()
for (i in 1:max(g)){
  behaviour<- c(behaviour, which(g==i)[1])
}
stationary <- grouped[behaviour]
sitenum <- cumsum(stationary==T)
sitenum[stationary==F] <- 0

# Initiate the model ###########################################################

#set initial path
x0 <- cbind(tapply(path$x[,1],twl$group,median), 
            tapply(path$x[,2],twl$group,median))


#set fixed locations 
fixedx <- rep_len(FALSE, length.out = nrow(x0))
fixedx[1] <- TRUE # We only fix the first location because this is an incomplete track

x0[fixedx,1] <- lon.calib
x0[fixedx,2] <- lat.calib

z0 <- trackMidpts(x0)

# plot stationary locations ####################################################
dtx0 <- as.data.frame(x0)
names(dtx0) <- c("x", "y")

par(mfrow=c(1,1))
data(wrld_simpl)
plot(dtx0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(dtx0, pch=19, col=  "cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
points(dtx0[sitenum > 0,], pch = 16, cex = 1, col = "green")
box()

# Movement model ###############################################################

# Here the model only reflects speed during active flight 
beta  <- c(2.2, 0.08)
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

# Create a Land mask for the group model #######################################
earthseaMask <- function(xlim, ylim, n = 2, pacific=FALSE, index) {
  
  if (pacific) { wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1
  rs = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
             rasterize(wrld_simpl, r, 1, silent = TRUE), 
             rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  
  #load polygon of blackpoll's range
  load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Birdlife_int_Full_blackpoll_range_polygon.R")
  #load("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_Full_blackpoll_range_polygon.R")
  
  
  #rasterize the polygon 
  range.raster <- rasterize(BLI.range.poly, rs)
  #range.raster <- rasterize(range.poly, rs)
  
  #Update the stationary mask 
  rs <- range.raster * rs
  
  # make the movement raster the same resolution as the stationary raster, but allow the bird to go anywhere by giving all cells a value of 1
  rm = rs; rm[] = 1
  
  # stack the movement and stationary rasters on top of each other
  mask = stack(rs, rm)
  
  xbin = seq(xmin(mask),xmax(mask),length=ncol(mask)+1)
  ybin = seq(ymin(mask),ymax(mask),length=nrow(mask)+1)
  mask = as.array(mask)[,,sort(unique(index)),drop=FALSE]
  
  function(p) mask[cbind(length(ybin)-.bincode(p[,2],ybin), .bincode(p[,1],xbin), index)]
}

#create the mask using the function 

xlim <- range(x0[,1])+c(-5,5)
ylim <- range(x0[,2])+c(-5,5)

index <- ifelse(stationary, 1, 2)

# testing #################
#  dtsm <- sm[,c("Lon.50.","Lat.50.")]
# # dtsm$index <- index
# # dtx0$index <- index
# # 
#  i <- dtsm[,1:2]
#  logp(i) 
############################

mask <- earthseaMask(xlim, ylim, n = 10, index=index)

# We will give locations on land a higher prior 
## Define the log prior for x and z
logp <- function(p) {
  f <- mask(p)
  ifelse(is.na(f), -1000, log(2))
}

# Define the Estelle model ####################################################
model <- groupedThresholdModel(twl$Twilight,
                               twl$Rise,
                               group = twl$group, #This is the group vector for each time the bird was at a point
                               twilight.model = "ModifiedGamma",
                               alpha = alpha,
                               beta =  beta,
                               x0 = x0, # median point for each group (defined by twl$group)
                               z0 = z0, # middle points between the x0 points
                               zenith = zenith0,
                               logp.x = logp,# land sea mask
                               fixedx = fixedx)


# define the error shape
x.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(x0))
z.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(z0))

# Fit the model
fit <- estelleMetropolis(model, x.proposal, z.proposal, iters = 1000, thin = 20)

#Tuning ########################################################################

# use output from last run
x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

model <- groupedThresholdModel(twl$Twilight, 
                               twl$Rise, 
                               group = twl$group,
                               twilight.model = "ModifiedGamma",
                               alpha = alpha, 
                               beta =  beta,
                               x0 = x0, z0 = z0,
                               logp.x = logp,
                               missing=twl$Missing,
                               zenith = zenith0,
                               fixedx = fixedx)

for (k in 1:3) {
  x.proposal <- mvnorm(chainCov(fit$x), s = 0.3)
  z.proposal <- mvnorm(chainCov(fit$z), s = 0.3)
  fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x),
                           z0 = chainLast(fit$z), iters = 300, thin = 20)
}

## Check if chains mix
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!fixedx, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!fixedx, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)

#Final model run ###############################################################
x.proposal <- mvnorm(chainCov(fit$x), s = 0.3)
z.proposal <- mvnorm(chainCov(fit$z), s = 0.3)

fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x),
                         z0 = chainLast(fit$z), iters = 2000, thin = 20, chain = 1)

#Summarize results #############################################################

# sm <- locationSummary(fit$x, time=fit$model$time)
sm <- SGAT2Movebank(fit$x, time = twl$Twilight, group = twl$group)

#create a plot of the stationary locations #####################################

# open jpeg
jpeg(paste0(dir, "/", geo.id,"_grouped_threshold_model_map.png"), width = 1024 , height = 990)

par(mfrow=c(1,1))
colours <- c("black",colorRampPalette(c("blue","yellow","red"))(max(twl.rev$Site)))
data(wrld_simpl)

# empty raster of the extent
r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
            xmx = xlim[2]+10, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = proj4string(wrld_simpl))

s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- slice(s, sliceIndices(s))

plot(sk, useRaster = F,col = c("transparent", rev(viridis::viridis(50))))
plot(wrld_simpl, xlim=xlim, ylim=ylim,add = T, bg = adjustcolor("black",alpha=0.1))

with(sm[sitenum>0,], arrows(`Lon.50.`, `Lat.2.5.`, `Lon.50.`, `Lat.97.5.`, length = 0, lwd = 2.5, col = "firebrick"))
with(sm[sitenum>0,], arrows(`Lon.2.5.`, `Lat.50.`, `Lon.97.5.`, `Lat.50.`, length = 0, lwd = 2.5, col = "firebrick"))
lines(sm[,"Lon.50."], sm[,"Lat.50."], col = adjustcolor("black", alpha = 0.6), lwd = 2)
points(sm[,"Lon.50."], sm[,"Lat.50."], col = ifelse(sm$StartTime > spring.equi - days(20) & sm$StartTime < spring.equi + days(20), "blue", "darkorchid4"), lwd = 2)

points(sm[,"Lon.50."], sm[,"Lat.50."], pch=21, bg=colours[sitenum+1], 
       cex = ifelse(sitenum>0, 3, 0), col = "firebrick", lwd = 2.5)

# Use this to number the stationary locations in the order they were use by the bird 
# points(sm[,"Lon.50."], sm[,"Lat.50."], pch=as.character(sitenum),
#        cex = ifelse(sitenum>0, 1, 0))

#The text in the symbols indicates the estimated number of days spent at each stopover location 
text(sm[,"Lon.50."], sm[,"Lat.50."], ifelse(sitenum>0, as.integer(((sm$EndTime - sm$StartTime)/86400)), ""), col="black") 

#Show dates
#text(sm[,"Lon.50."], sm[,"Lat.50."], ifelse(sitenum>0, as.character(sm$StartTime), ""), col="red", pos = 1) 

#Close jpeg
dev.off()

#plot of longitude and latitude ################################################

# open jpeg
jpeg(paste0(dir, "/", geo.id,"_grouped_threshold_model_lon_lat.png"), width = 1024 , height = 990)

par(mfrow=c(2,1))

plot(sm$StartTime, sm$"Lon.50.", ylab = "Longitude", xlab = "", yaxt = "n", type = "n", ylim = c(min(sm$Lon.50.) - 10, max(sm$Lon.50.) + 10))
axis(2, las = 2)
polygon(x=c(sm$StartTime,rev(sm$StartTime)), y=c(sm$`Lon.2.5.`,rev(sm$`Lon.97.5.`)), border="gray", col="gray")
lines(sm$StartTim,sm$"Lon.50.", lwd = 2)
abline(v = spring.equi, lwd = 2, lty = 2, col = "orange")

#Add points for stopovers 
points(sm$StartTime, sm$"Lon.50.", pch=21, bg=colours[sitenum+1], 
       cex = ifelse(sitenum>0, 3, 0), col = "firebrick", lwd = 2.5)

#The text in the symbols indicates the estimated number of days spent at each stopover location 
text(sm$StartTime, sm$"Lon.50.", ifelse(sitenum>0, as.integer(((sm$EndTime - sm$StartTime)/86400)), ""), col="black") 

plot(sm$StartTime, sm$"Lat.50.", ylab = "Latitude", xlab = "", yaxt = "n", type = "n", ylim = c(min(sm$Lat.50.) - 10, max(sm$Lat.50.) + 10))
axis(2, las = 2)
polygon(x=c(sm$StartTime,rev(sm$StartTime)), y=c(sm$`Lat.2.5.`,rev(sm$`Lat.97.5.`)), border="gray", col="gray")
lines(sm$StartTim,sm$"Lat.50.", lwd = 2)
abline(v = spring.equi, lwd = 2, lty = 2, col = "orange")

#Add points for stopovers 
points(sm$StartTime, sm$"Lat.50.", pch=21, bg=colours[sitenum+1], 
       cex = ifelse(sitenum>0, 3, 0), col = "firebrick", lwd = 2.5)

#The text in the symbols indicates the estimated number of days spent at each stopover location 
text(sm$StartTime, sm$"Lat.50.", ifelse(sitenum>0, as.integer(((sm$EndTime - sm$StartTime)/86400)), ""), col="black")

#Close jpeg
dev.off()

#Extract Stationary locations ##################################################
sm$sitenum <- sitenum
sm$duration <- as.numeric(difftime(sm$EndTime, sm$StartTime), unit = "days")
stat.loc <- sm[sitenum > 0, ]

#plot only stationary locations
par(mfrow=c(1,1))

data(wrld_simpl)
plot(wrld_simpl, xlim=xlim, ylim=ylim, col = "grey95")
points(sm$Lon.50., sm$Lat.50., pch = 16, cex = 0, col = "firebrick", type = "o")
points(stat.loc$Lon.50., stat.loc$Lat.50., pch = 16, cex = 1.5, col = "firebrick")

# add column with geolocator ID
sm$geo_id <- geo.id

#add a column that categorizes the locations (based on the groupthreshold model output)
sm <- sm %>% mutate(period= case_when(StartTime <= anytime("2020-04-30 10:29:15", asUTC = T, tz = "GMT")  ~ "Non-breeding period",
                                      StartTime > anytime("2020-04-30 10:29:15", asUTC = T, tz = "GMT") ~ "Pre-breeding migration"))

#Save the output of the model 
#save(sm, file = paste0(dir,"/", geo.id,"_SGAT_GroupedThreshold_summary.csv"))
#save(fit, file = paste0(dir,"/", geo.id,"_SGAT_GroupedThreshold_fit.R"))

#load the output of the model 
#load(file = paste0(dir,"/", geo.id,"_SGAT_GroupedThreshold_summary.csv"))
#load(file = paste0(dir,"/", geo.id,"_SGAT_GroupedThreshold_fit.R"))

# Record details for the geolocator analysis 
writeLines(c(paste("Median zenith angle in Breeding grounds =", zenith),
             paste("Zero deviation angle in Breeding grounds =", zenith0),
             paste("Hill Ekstrom zenith angle =", zenith_sd)),
           paste0(dir, "/", geo.id,"_analysis_details.txt", sep = ""))

# Examine twilights ############################################################

#load the raw threshold path path x0_r
load(file = paste0(dir,"/", geo.id, "_initial_path_raw.csv"))

#plot of light transitions throughout the fall migration 
par(mfrow=c(2,1))
plot(twl$Twilight, x0_r[,1], type = "o", ylab = "longitude", xlab = "time")
rect(anytime("2013-10-13"), min(x0_r[,1])-2, anytime("2013-10-15"), max(x0_r[,1])+2, col = alpha("orange", 0.2), lty=0)
rect(anytime("2013-10-17"), min(x0_r[,1])-2, anytime("2013-10-18"), max(x0_r[,1])+2, col = alpha("orange", 0.2), lty=0)

plot(twl$Twilight, x0_r[,2], type = "o", ylab = "latitude", xlab = "time")
rect(anytime("2013-10-13"), min(x0_r[,2])-2, anytime("2013-10-15"), max(x0_r[,2])+2, col = alpha("orange", 0.2), lty=0)
rect(anytime("2013-10-17"), min(x0_r[,2])-2, anytime("2013-10-18"), max(x0_r[,2])+2, col = alpha("orange", 0.2), lty=0)

#Fall transoceanic flight 

start <- "2013-10-05"
end <- "2013-10-25"

par(cex.lab=1.4)
par(cex.axis=1.4)
par(mfrow=c(3,1), mar = c(5,5,0.1,5))
plot(lig$Date[lig$Date > start & lig$Date < end], lig$Light[lig$Date > start & lig$Date < end], type = "o",
     ylab = "Light level", xlab = "Time")
rect(anytime("2013-10-13"), min(lig$Light)-2, anytime("2013-10-15"), max(lig$Light)+2, col = alpha("yellow", 0.2), lty=0)
rect(anytime("2013-10-17"), min(lig$Light)-2, anytime("2013-10-18"), max(lig$Light)+2, col = alpha("yellow", 0.2), lty=0)

plot(twl$Twilight[twl$Twilight> start & twl$Twilight < end], x0_r[,1][twl$Twilight > start & twl$Twilight < end],
     ylab = "Longitude", xlab = "Time")
rect(anytime("2013-10-13"), min(x0_r[,1])-2, anytime("2013-10-15"), max(x0_r[,1])+2, col = alpha("yellow", 0.2), lty=0)
rect(anytime("2013-10-17"), min(x0_r[,1])-2, anytime("2013-10-18"), max(x0_r[,1])+2, col = alpha("yellow", 0.2), lty=0)

plot(twl$Twilight[twl$Twilight > start & twl$Twilight < end], x0_r[,2][twl$Twilight > start & twl$Twilight < end],
     ylab = "Latitude", xlab = "Time")
rect(anytime("2013-10-13"), min(x0_r[,2])-2, anytime("2013-10-15"), max(x0_r[,2])+2, col = alpha("yellow", 0.2), lty=0)
rect(anytime("2013-10-17"), min(x0_r[,2])-2, anytime("2013-10-15"), max(x0_r[,2])+2, col = alpha("yellow", 0.2), lty=0)

par(cex.lab= 1)
par(cex.axis= 1)

