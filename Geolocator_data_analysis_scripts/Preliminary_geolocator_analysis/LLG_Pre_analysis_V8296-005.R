# source: unpublished data 
# tag number: V8296 005
# site: Quebec

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

geo.id <- "V8296-005"

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# geo deployment location 
lat.calib <- 47.38323	
lon.calib <- -71.09479

# time of deployment (from reference file)
deploy.start <- anytime("2019-06-18	06:45:00", asUTC = T, tz = "GMT")

# time of recovery (estimate from light data)
deploy.end <- anytime("2020-07-12 17:55:00", asUTC = T, tz = "GMT")

#Equinox times
fall.equi <- anytime("2019-09-23", asUTC = T, tz = "GMT")
spring.equi <- anytime("2020-03-19", asUTC = T, tz = "GMT")

#DATA EXTRACTION ##############################################################

# import lig data 
lig <- readLig(paste0(dir,"/ML6740 V8296 005 reconstructed_000.lig"), skip = 1)

#remove rows before and after deployment time 
lig <- lig[(lig$Date > deploy.start),]

#parameter to visualize the data 
offset <- 16 # adjusts the y-axis to put night (dark shades) in the middle

#Threshold light level 
threshold <- 1.5 

# visualize threshold over light levels  
thresholdOverLight(lig, threshold, span =c(0, 25000))

# FIND TIME SHIFT ##############################################################

#This geolocator has a a time shift visible on this plot
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

#write.csv(twl_in, paste0(dir,"/Pre_analysis_V8296_005_twl_times_initial.csv"))
twl_in <- read.csv(paste0(dir,"/Pre_analysis_V8296_005_twl_times_initial.csv"))
twl_in$Twilight <- as.POSIXct(twl_in$Twilight, tz = "UTC")

# Period over which to calculate the time shift. It should be while the bird is
# still in the breeding grounds
period <- as.POSIXct(c("2019-06-19", "2019-08-10"), tz = "UTC")

#plot the period over the light image
lightImage( tagdata = lig,
            offset = offset,
            zlim = c(0, 64))

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
#lig$Date <- lig$Date - (shift$shift)
lig$Date <- lig$Date - 1*60*60

#TWILIGHT ANNOTATION ##########################################################

#open jpeg
jpeg(paste0(dir, "/V8296-005_light_plot.png"), width = 1024, height = 990)

lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))
dev.off()

#Detect twilight times, for now do not edit twilight times  
twl <- preprocessLight(lig, 
                       threshold = threshold,
                       offset = offset, 
                       lmax = 64,         # max. light value
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
# write.csv(twl, paste0(dir,"/Pre_analysis_V8296_005_twl_times.csv"))

###############################################################################
# SGAT ANALYSIS ###############################################################
###############################################################################

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/Pre_analysis_V8296_005_twl_times.csv"))
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

#calibration period before the migration 
tm.calib <- as.POSIXct(c("2019-06-19", "2019-09-10"), tz = "UTC")

#calibration period after the migration 
tm.calib2 <- as.POSIXct(c("2020-05-30", "2020-06-15"), tz = "UTC")

abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")
abline(v = tm.calib2, lwd = 2, lty = 2, col = "orange")

# subset of twilight times that are within the calibration period
d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2] | Twilight>=tm.calib2[1] & Twilight<=tm.calib2[2])

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

# Vector with different zenith angle for the breeding and nonbreeding periods 
zenith_twl <- data.frame(Date = twl$Twilight) %>%
  mutate(zenith = case_when(Date < fall.equi ~ zenith0,
                            Date > fall.equi ~ zenith_sd))
zeniths <- zenith_twl$zenith

# adjust the zenith angles calculated from the breeding sites 
zenith0_ad <- zenith0 + abs(zenith - zenith_sd)
zenith_ad  <- zenith_sd

# use a different zenith angle for the breeding and nonbreeding periods 
zenith_twl <- data.frame(Date = twl$Twilight) %>%
  mutate(zenith = case_when(Date < fall.equi ~ zenith0,
                            Date > fall.equi ~ zenith0_ad))

zeniths <- zenith_twl$zenith

# Movement model ###############################################################

#this movement model should be based on the estimated migration speed of the blackpoll warbler 
beta  <- c(0.7, 0.05)
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

# Initial Path #################################################################
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol=0.01)

x0 <- path$x
z0 <- trackMidpts(x0)

# open jpeg
jpeg(paste0(dir, "/ML6440_V8296_005_Threshold_path.png"), width = 1024, height = 990)

data(wrld_simpl)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(path$x, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()

dev.off()

# Define known locations #######################################################

#we set the location of geolocator deployment and recovery as fixed locations for the MCMC sampler 

fixedx <- rep(F, nrow(x0))
fixedx[1:2] <- T # first two location estimates

fixedx[(nrow(x0) - 1):nrow(x0)] <- T # last two location estimates

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
                        #logp.x = log.prior, logp.z = log.prior, 
                        x0 = x0,
                        z0 = z0,
                        zenith = zeniths,
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
                        twilight.model = "Gamma",
                        alpha = alpha,
                        beta = beta,
                        #logp.x = log.prior, logp.z = log.prior, 
                        x0 = x0,
                        z0 = z0,
                        zenith = zeniths,
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

# Plot Results #################################################################

# open jpeg
jpeg(paste0(dir, "/ML6440_V8296_005_Estelle_path.png"), width = 1024 , height = 990)

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
jpeg(paste0(dir, "/ML6440_V8296_005_mean_lon_lat.png"), width = 1024 , height = 990)

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
sm <- sm %>% mutate(stationary = ifelse(abs(lead(Lon.mean) - Lon.mean) < 1 & abs(lead(Lat.mean) - Lat.mean) < 1, 1, 0)) 

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

#SGAT Groupe model analysis ####################################################

# group twilight times were birds were stationary 
geo_twl <- export2GeoLight(twl)

# Often it is necessary to play around with quantile and days
# quantile defines how many stopovers there are. the higher, the fewer there are
# days indicates the duration of the stopovers 
cL <- changeLight(twl=geo_twl, quantile=0.86, summary = F, days = 2, plot = T)

# merge site helps to put sites together that are separated by single outliers.
mS <- mergeSites(twl = geo_twl, site = cL$site, degElevation = 90-zenith0, distThreshold = 500)

#back transfer the twilight table and create a group vector with TRUE or FALSE according to which twilights to merge 
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
fixedx[1] <- TRUE
fixedx[c(1, length(fixedx))] <- TRUE

x0[fixedx,1] <- lon.calib
x0[fixedx,2] <- lat.calib

z0 <- trackMidpts(x0)

# plot stationary locations ####################################################
dtx0 <- as.data.frame(x0)
names(dtx0) <- c("x", "y")

data(wrld_simpl)
plot(dtx0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(dtx0, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
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
                               x0 = x0, # median point for each greoup (defined by twl$group)
                               z0 = z0, # middle points between the x0 points
                               zenith = zeniths,
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
                               zenith = zeniths,
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

#Save the output of the model 
#save(sm, file = paste0(dir,"/Pre_analysis_8296_005_SGAT_GroupedThreshold_summary.csv"))
#save(fit, file = paste0(dir,"/Pre_analysis_8296_005_SGAT_GroupedThreshold_fit.R"))

#create a plot of the stationary locations #####################################
par(mfrow=c(1,1))
colours <- c("black",colorRampPalette(c("blue","yellow","red"))(max(twl.rev$Site)))
data(wrld_simpl)

# empty raster of the extent
r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = proj4string(wrld_simpl))

s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- slice(s, sliceIndices(s))

plot(sk, useRaster = F,col = c("transparent", rev(viridis::viridis(50))))
plot(wrld_simpl, xlim=xlim, ylim=ylim,add = T, bg = adjustcolor("black",alpha=0.1))

with(sm[sitenum>0,], arrows(`Lon.50.`, `Lat.2.5.`, `Lon.50.`, `Lat.97.5.`, length = 0, lwd = 2.5, col = "firebrick"))
with(sm[sitenum>0,], arrows(`Lon.2.5.`, `Lat.50.`, `Lon.97.5.`, `Lat.50.`, length = 0, lwd = 2.5, col = "firebrick"))
lines(sm[,"Lon.50."], sm[,"Lat.50."], col = adjustcolor("black", alpha = 0.6), lwd = 2)
points(sm[,"Lon.50."], sm[,"Lat.50."], col = ifelse(sm$StartTime > fall.equi - days(10) & sm$StartTime < fall.equi + days(10), "blue", "darkorchid4"), lwd = 2)

points(sm[,"Lon.50."], sm[,"Lat.50."], pch=21, bg=colours[sitenum+1], 
       cex = ifelse(sitenum>0, 3, 0), col = "firebrick", lwd = 2.5)

# Use this to number the stationary locations in the order they were use by the bird 
# points(sm[,"Lon.50."], sm[,"Lat.50."], pch=as.character(sitenum),
#        cex = ifelse(sitenum>0, 1, 0))

#The text in the symbols indicates the estimated number of days spent at each stopover location 
text(sm[,"Lon.50."], sm[,"Lat.50."], ifelse(sitenum>0, as.integer(((sm$EndTime - sm$StartTime)/86400)), ""), col="black") 

#Show dates
#text(sm[,"Lon.50."], sm[,"Lat.50."], ifelse(sitenum>0, as.character(sm$StartTime), ""), col="red", pos = 1) 

#plot of longitude and latitude
par(mfrow=c(2,1))

plot(sm$StartTime, sm$"Lon.50.", ylab = "Longitude", xlab = "", yaxt = "n", type = "n", ylim = c(min(sm$Lon.50.) - 10, max(sm$Lon.50.) + 10))
axis(2, las = 2)
polygon(x=c(sm$StartTime,rev(sm$StartTime)), y=c(sm$`Lon.2.5.`,rev(sm$`Lon.97.5.`)), border="gray", col="gray")
lines(sm$StartTim,sm$"Lon.50.", lwd = 2)
abline(v = fall.equi, lwd = 2, lty = 2, col = "orange")
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
abline(v = fall.equi, lwd = 2, lty = 2, col = "orange")
abline(v = spring.equi, lwd = 2, lty = 2, col = "orange")

#Add points for stopovers 
points(sm$StartTime, sm$"Lat.50.", pch=21, bg=colours[sitenum+1], 
       cex = ifelse(sitenum>0, 3, 0), col = "firebrick", lwd = 2.5)

#The text in the symbols indicates the estimated number of days spent at each stopover location 
text(sm$StartTime, sm$"Lat.50.", ifelse(sitenum>0, as.integer(((sm$EndTime - sm$StartTime)/86400)), ""), col="black") 


################################################################################
# FLIGHTR ANALYSIS #############################################################
################################################################################

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/Pre_analysis_V8296_005_twl_times.csv"))
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

twlexp <- twGeos2TAGS(raw = lig[, c("Date", "Light")], twl = twl,
                      threshold = 1.5,
                      filename = paste0(dir, "/V8296_005_TAGS_data.csv"))

tags <- get.tags.data(paste0(dir, "/V8296_005_TAGS_data.csv"))

#Calibration ###################################################################

# plot slopes of light transition over calibration period 
plot_slopes_by_location(Proc.data=tags, location=c(lon.calib, lat.calib), ylim=c(-2, 2))

abline(v=as.POSIXct("2019-09-02"), col = "green") # end of first calibration period
abline(v=as.POSIXct("2020-05-20"), col = "green") # start of the second calibration period

Calibration.periods<-data.frame(
  calibration.start=as.POSIXct(c(NA)),
  calibration.stop=as.POSIXct(c("2019-09-02")),
  lon=lon.calib, lat=lat.calib) 
# use c() also for the geographic coordinates, 
# if you have more than one calibration location
# (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))
print(Calibration.periods)

Calibration<-make.calibration(tags, Calibration.periods, model.ageing=TRUE, plot.final = T)

# Create grid for spatial extent ###############################################

Grid <- make.grid(left=lon.calib -50, bottom=lat.calib-60, right=lon.calib+50, top= lat.calib + 10,
                  distance.from.land.allowed.to.use=c(-Inf, 1100),
                  distance.from.land.allowed.to.stay=c(-Inf, 50))

# create the model prerun object ###############################################
#all.in <- make.prerun.object(tags, Grid, start=c(lon.calib, lat.calib),
#                             Calibration=Calibration, M.mean=750)

#save(all.in, file = paste0(dir, "/V8296-005_FlightRCalib.RData"))

#run twilight filter ###########################################################

#Load results 
load(paste0(dir, "/V8296-005_FlightRCalib.RData"))

#run filter 
#nParticles=1e4 #increase to 1e6 for test 
#Result<-run.particle.filter(all.in, threads=-1,
#                            nParticles=nParticles, known.last=TRUE,
#                            precision.sd=25, check.outliers= T, 
#                            b=1700)

#save(Result, file = paste0(dir, "/V8296-005_FLightRResult.RData"))

# Plot results ################################################################# 
load(paste0(dir, "/V8296-005_FLightRResult.RData"))

#longitude and latitude
plot_lon_lat(Result)

#simple map 
library(ggmap)
ggmap::register_google("AIzaSyABANOgjTyVFpOuDOiyPlBL4geijIy6vPo")
map.FLightR.ggmap(Result, zoom=3, save = FALSE)

# find Stationary sites #############################################################
Summary <-stationary.migration.summary(Result, prob.cutoff = 0.5)
view(Summary$Stationary.periods)

# find longest stationary period 
Summary$Potential_stat_periods

# this is the period between 357-485
Result$Results$Quantiles[c(357,485),]$time

# The new calibration period ranges between "2019-12-13 22:42:48 GMT" and "2020-02-15 23:07:30 GMT"
# mean lat =  14.3838540
# mean lon =  -67.71948

# Calibration with new stationary period #######################################
# plot slopes of light transition over calibration period 

Calibration.periods<-data.frame(
  calibration.start=as.POSIXct(c("2019-12-13")),
  calibration.stop=as.POSIXct(c("2020-02-15")),
  lon=lon.calib, lat=lat.calib) 
# use c() also for the geographic coordinates, 
# if you have more than one calibration location
# (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))
print(Calibration.periods)

NewCalibration <- make.calibration(tags, Calibration.periods, model.ageing=TRUE, plot.final = T)

# create new model prerun object ###############################################
#all.in.adjusted <- make.prerun.object(tags, Grid, start=c(lon.calib, lat.calib),
#                             Calibration=NewCalibration, M.mean=750)

#save(all.in.adjusted, file = paste0(dir, "/V8296-005_FlightRCalib_nonbreed_calib.RData"))

#run twilight filter again #####################################################
#nParticles=1e6 #increase to 1e6 for test 
Result_adjusted <-run.particle.filter(all.in, threads=-1,
                                      nParticles=nParticles, known.last=TRUE,
                                      precision.sd=25, check.outliers= T, 
                                      b=1700)

#save(Result_adjusted, file = paste0(dir, "/V8296-005_FLightRResult_nonbreed_calib.RData"))

# Plot new results ################################################################# 
load(paste0(dir, "/V8296-005_FLightRResult_nonbreed_calib.RData"))

#longitude and latitude
plot_lon_lat(Result_adjusted)

#simple map 
library(ggmap)
ggmap::register_google("AIzaSyABANOgjTyVFpOuDOiyPlBL4geijIy6vPo")
map.FLightR.ggmap(Result_adjusted, zoom=3, save = FALSE)

# find Stationary sites #############################################################
Summary <-stationary.migration.summary(Result_adjusted, prob.cutoff = 0.3)
view(Summary$Stationary.periods)


# Now we want to plot the detected stationary periods on a map
Summary$Stationary.periods$stopover_duration<-as.numeric(difftime(Summary$Stationary.periods$Departure.Q.50,Summary$Stationary.periods$Arrival.Q.50, units='days'))
# Now I want to select the periods which were >=2 days and 
Main_stopovers<-Summary$Stationary.periods[is.na(Summary$Stationary.periods$stopover_duration) | Summary$Stationary.periods$stopover_duration>=2,]
# delete breeding season
Main_stopovers<-Main_stopovers[-which(is.na(Main_stopovers$stopover_duration)),]

Coords2plot<-cbind(Result$Results$Quantiles$Medianlat, Result$Results$Quantiles$Medianlon)

for (i in 1:nrow(Summary$Potential_stat_periods)) {
  Coords2plot[Summary$Potential_stat_periods[i,1]:
                Summary$Potential_stat_periods[i,2],1] =  
    Summary$Stationary.periods$Medianlat[i]
  
  Coords2plot[Summary$Potential_stat_periods[i,1]:
                Summary$Potential_stat_periods[i,2],2] =  
    Summary$Stationary.periods$Medianlon[i]
}
Coords2plot<-Coords2plot[!duplicated(Coords2plot),]

#pdf('FLightR_shrike_migration_with_stopovers.pdf', width=6, height=9)
par(mar=c(0,0,0,0))
map('worldHires', ylim=c(-20, 60), xlim=c(-120, -50), col=grey(0.7),
    fill=TRUE, border=grey(0.9), mar=rep(0.5, 4), myborder=0)

lines(Coords2plot[,1]~Coords2plot[,2], col='red', lwd=2)
points(Coords2plot[,1]~Coords2plot[,2], ,lwd=2, col='red', pch=19)

# Here we assign the colours to represent time of the year
Seasonal_palette<-grDevices::colorRampPalette(grDevices::hsv(1-((1:365)+(365/4))%%365/365,
                                                             s=0.8, v=0.8), space="Lab")
Seasonal_colors<-Seasonal_palette(12)

Main_stopovers$Main_month<-as.numeric(format(Main_stopovers$Arrival.Q.50+
                                               Main_stopovers$stopover_duration/2, format='%m'))

points(Main_stopovers$Medianlat~Main_stopovers$Medianlon, pch=21, 
       cex=log(as.numeric(Main_stopovers$stopover_duration)),
       bg=Seasonal_colors[Main_stopovers$Main_month])

# Now, for each of these points we plot the uncertainties
# Horizontal
segments(y0=Main_stopovers$Medianlat, x0=Main_stopovers$FstQu.lon,
         x1=Main_stopovers$TrdQu.lon, lwd=2)
# Vertical
segments(x0=Main_stopovers$Medianlon, y0=Main_stopovers$FstQu.lat,
         y1=Main_stopovers$TrdQu.lat, lwd=2)

# dev.off()
