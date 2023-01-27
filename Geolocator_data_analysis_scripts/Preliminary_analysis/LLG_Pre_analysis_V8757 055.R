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

geo.id <- "V8757-055"

# geo deployment location 
lat.calib <- 62.475799
lon.calib <- -114.699932

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# time of deployment
deploy.start <- anytime("2019-06-27", tz = "GMT")

#Equinox times
fall.equi <- anytime("2019-09-23", tz = "GMT")
spring.equi <- anytime("2020-03-19", tz = "GMT")

#Find number of cores available for analysis
Threads= detectCores()-1

###############################################################################
#DATA EXTRACTION ##############################################################
###############################################################################

# import lig data 
lig <- readLig(paste0(dir,"/ML6740 V8757 055 reconstructed_000.lig"), skip = 1)

#adjust year
lig$Date <- lig$Date %m+% years(7)

#remove rows before and after deployment time 
lig <- lig[(lig$Date > deploy.start),]

#adjust time 
lig$Date <- lig$Date + 5.5 *60*60

###############################################################################
#TWILIGHT ANNOTATION ##########################################################
###############################################################################

threshold <- 1.5 

# visualize threshold over light levels  
thresholdOverLight(lig, threshold, span =c(70000, 75000))

# plot light levels 
offset <- 20 # adjusts the y-axis to put night (dark shades) in the middle

# open jpeg
jpeg(paste0(dir, "/V8757_055_light_plot.png"), width = 1024, height = 990)

lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 64))

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

# This twilight file has some duplicates between rows 456 and 707
# No data appears to have been lost, so the duplicates can simply be removed 
twl[duplicated(twl$Twilight),]
twl[456:707,]
twl <- filter(twl, !duplicated(twl$Twilight))

# Save the twilight times 
# write.csv(twl, paste0(dir,"/Pre_analysis_V8757_055_twl_times.csv"))

###############################################################################
# SGAT ANALYSIS ###############################################################
# Threshold analysis and Estelle model ########################################
###############################################################################

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/Pre_analysis_V8757_055_twl_times.csv"))
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

#calibration period before the migration 
tm.calib <- as.POSIXct(c("2019-08-03", "2019-09-03"), tz = "UTC")

abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

# subset of twilight times that are within the calibration period
d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])

# perform the calibration and verify the fit with the gamma or log distribution
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

#parameters of the error distribution 
#parameters of the error distribution 
zenith  <- calib[1] 
zenith0 <- calib[2]

alpha <- calib[3:4]

# Calibration in the breeding grounds is difficult here because of very brief nights and long periods of incomplete darkness

#alternative calibration #######################################################
startDate <- "2019-11-15"
endDate   <- "2020-04-15"

start = min(which(as.Date(twl$Twilight) == startDate))
end = max(which(as.Date(twl$Twilight) == endDate))

(zenith_sd <- findHEZenith(twl, tol=0.01, range=c(start,end)))

# Movement model ###############################################################

#this movement model should be based on the estimated migration speed of the blackpoll warbler 
beta  <- c(0.45, 0.05)
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

# Initial Path #################################################################
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol=0.18)

x0 <- path$x
z0 <- trackMidpts(x0)

# open jpeg
jpeg(paste0(dir, "/V8757_055_Threshold_path.png"), width = 1024, height = 990)

# plot of threshold path  
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
  
  abundance <- raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geo_spatial_data/bkpwar_abundance_seasonal_full-year_mean_2021.tif")
  abundance_resamp <- projectRaster(abundance, mask, method = "ngb")
  abundance_resamp[is.nan(abundance_resamp) | abundance_resamp == 0] <- NA
  abundance_resamp[abundance_resamp > 0 ] <- 1
  
  mask <- mask * abundance_resamp
  
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
                        zenith = zenith,
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
                        logp.x = log.prior, logp.z = log.prior, 
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith,
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

# open jpeg
jpeg(paste0(dir, "/V8757_055__Estelle_path.png"), width = 1024 , height = 990)

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
jpeg(paste0(dir, "V8757_055__mean_lon_lat.png"), width = 1024 , height = 990)

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

################################################################################
# FLIGHTR ANALYSIS #############################################################
################################################################################

# Import file with twilight times  
twl <- read.csv(paste0(dir,"/Pre_analysis_V8757_055_twl_times.csv"))
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

twlexp <- twGeos2TAGS(raw = lig[, c("Date", "Light")], twl = twl,
                      threshold = 1.5,
                      filename = paste0(dir, "/V8757_055_TAGS_data.csv"))

tags <- get.tags.data(paste0(dir, "/V8757_055_TAGS_data.csv"))

#Calibration ###################################################################

# plot slopes of light transition over calibration period 
plot_slopes_by_location(Proc.data=tags, location=c(lon.calib, lat.calib), ylim=c(-3, 3))

abline(v=as.POSIXct(deploy.start + days(10)), col = "green") # start of first calibration period
abline(v=as.POSIXct(deploy.start + days(60)), col = "green") # end of first calibration period
abline(v=as.POSIXct("2020-06-25"), col = "green") # start of the second calibration period
abline(v=as.POSIXct("2020-07-12"), col = "green") # end of the second calibration period

Calibration.periods<-data.frame(
  calibration.start=as.POSIXct(c(deploy.start + days(10), "2020-06-25")),
  calibration.stop=as.POSIXct(c(deploy.start + days(60), "2020-07-02" )),
  lon=lon.calib, lat=lat.calib) 
print(Calibration.periods)

Calibration<-make.calibration(tags, Calibration.periods, model.ageing=TRUE, plot.final = T)

# Create grid for spatial extent ###############################################

Grid <- make.grid(left=lon.calib -30, bottom=lat.calib-80, right=lon.calib+80, top= lat.calib + 5,
                  distance.from.land.allowed.to.use=c(-Inf, 1100),
                  distance.from.land.allowed.to.stay=c(-Inf, 50))

# create the model prerun object ###############################################
all.in <- make.prerun.object(tags, Grid, start=c(lon.calib, lat.calib),
                             Calibration=Calibration, threads = min(Threads, 6))

save(all.in, file = paste0(dir, "/V8757_055_FlightRCalib.RData"))

#run twilight filter ###########################################################

#Load  model prerun object
load(paste0(dir, "/V8757_055_FlightRCalib.RData"))

#run filter 
nParticles=1e4 # start at 1e4 for initial run 
Result<-run.particle.filter(all.in, threads= min(Threads, 6),
                            nParticles=nParticles, known.last=TRUE,
                            precision.sd=25, check.outliers= T, 
                            b=2700)

save(Result, file = paste0(dir, "/V8757_055_FLightRResult.RData"))

# Plot results ################################################################# 
load(paste0(dir, "/V8757_055_FLightRResult.RData"))

#longitude and latitude
plot_lon_lat(Result)

#simple map 
library(ggmap)
ggmap::register_google("AIzaSyABANOgjTyVFpOuDOiyPlBL4geijIy6vPo")
map.FLightR.ggmap(Result, zoom=3, save = FALSE)

# find Stationary sites #############################################################
Summary <-stationary.migration.summary(Result, prob.cutoff = 0.8)
view(Summary$Stationary.periods)

# find longest stationary period 
Summary$Potential_stat_periods

# this is the period between 357-485
Result$Results$Quantiles[c(236,302),]$time

# The new calibration period ranges between "2020-02-18" and "2020-04-22"
new.calib.start <- "2019-11-10"
new.calib.end  <- "2019-12-03"
mean.lat <-  9.900000
mean.lon <- -63.05692

# Calibration with new stationary period #######################################
# plot slopes of light transition over calibration period 

Calibration.periods<-data.frame(
  calibration.start=as.POSIXct(c(new.calib.start)),
  calibration.stop=as.POSIXct(c(new.calib.end)),
  lon= mean.lon, lat= mean.lat) 
print(Calibration.periods)

NewCalibration <- make.calibration(tags, Calibration.periods, model.ageing=TRUE, plot.final = T)

# create new model prerun object ###############################################
all.in.adjusted <- make.prerun.object(tags, Grid, start=c(lon.calib, lat.calib),
                                      Calibration=NewCalibration)

save(all.in.adjusted, file = paste0(dir, "/V8757_055_FlightRCalib_nonbreed_calib.RData"))

#run twilight filter again #####################################################
nParticles=1e4 #increase to 1e6 for final run  

a <- Sys.time()
Result_adjusted <-run.particle.filter(all.in, threads= -1,
                                      nParticles=nParticles, known.last=TRUE,
                                      precision.sd=25, check.outliers= T, 
                                      b=2700)

b <- Sys.time()
d <- a-b

save(Result_adjusted, file = paste0(dir, "/V8757_055_FLightRResult_nonbreed_calib.RData"))

# Plot new results ################################################################# 
load(paste0(dir, "/V8757_055_FLightRResult_nonbreed_calib.RData"))

#longitude and latitude
plot_lon_lat(Result_adjusted)

#simple map 
library(ggmap)
ggmap::register_google("AIzaSyABANOgjTyVFpOuDOiyPlBL4geijIy6vPo")
map.FLightR.ggmap(Result_adjusted, zoom=3, save = FALSE)

# find Stationary sites #############################################################
Summary <-stationary.migration.summary(Result_adjusted, prob.cutoff = 0.2)
view(Summary$Stationary.periods)



