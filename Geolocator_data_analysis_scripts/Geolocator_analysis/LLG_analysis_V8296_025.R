# source: unpublished data 
# tag number: V8296-025
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
# rm(list=ls())

# Load helper functions 
# Load helper functions 
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")
geo.id <- "V8296_025"

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# geo deployment location 
lat.calib <- 47.37864
lon.calib <- -71.09769

# time of deployment (from reference file)
#deploy.start <- anytime("", asUTC = T, tz = "GMT")

# time of recovery (estimate from light data)
deploy.end <- anytime("", asUTC = T, tz = "GMT")

#Equinox times
fall.equi <- anytime("2019-09-23", asUTC = T, tz = "GMT")
spring.equi <- anytime("2020-03-19", asUTC = T, tz = "GMT")

#DATA EXTRACTION ##############################################################

# import lig data 
lig <- readLig(paste0(dir,"/Raw_light_data_", geo.id, ".lig"), skip = 1)

#remove rows before and after deployment time 
#lig <- lig[(lig$Date > deploy.start),]

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
            zlim = c(0, 64))

tsimageDeploymentLines(lig$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# we will do an initial twilight annotation to find identify the time interval
# by which we need to shift time
# There should be not need to edit, delete or insert twilights for this
# twl_in <- preprocessLight(lig,
#                           threshold = threshold,
#                           offset = offset,
#                           lmax = 64,         # max. light value
#                           gr.Device = "x11", # MacOS version (and windows)
#                           dark.min = 60)

#write.csv(twl_in, paste0(dir,"/", geo.id,"_twl_times_initial.csv"))
twl_in <- read.csv(paste0(dir,"/", geo.id,"_twl_times_initial.csv"))
twl_in$Twilight <- as.POSIXct(twl_in$Twilight, tz = "UTC")

# Period over which to calculate the time shift. It should be while the bird is
# still in the breeding grounds
period <- as.POSIXct(c("2019-06-25", "2019-08-14"), tz = "UTC")

# #plot the period over the light image
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


# verify the that the time shift measured makes sense 
shift

#adjust  the based on the time shift detected 
lig$Date <- lig$Date - (shift$shift)
#lig$Date <- lig$Date + 1*60*60

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

# read parameters that will be used during the analysis
niter <- ref_data$MCMC.iter[which(ref_data$geo.id == geo.id)]
nthin <- ref_data$MCMC.thin[which(ref_data$geo.id == geo.id)]
chains <- ref_data$MCMC.chains[which(ref_data$geo.id == geo.id)]
nday <- ref_data$changeLight.days[which(ref_data$geo.id == geo.id)]

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

#calibration period before the migration 
tm.calib1 <- as.POSIXct(c("2019-06-25", "2019-08-25"), tz = "UTC")

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
cL <- changeLight(twl=geo_twl, quantile=0.6, summary = F, days = 10, plot = T)
# merge site helps to put sites together that are separated by single outliers.
mS <- mergeSites(twl = geo_twl, site = cL$site, degElevation = 90-zenith, distThreshold = 250)

#specify which site is the stationary one
site           <- mS$site[mS$site>0] # get rid of movement periods
stationarySite <- which(table(site) == max(table(site))) # find the site where bird is the longest

#find the dates that the bird arrives and leaves this stationary site
start <- min(which(mS$site == stationarySite))
end   <- max(which(mS$site == stationarySite))

#calculate the zenith angle that minimizes variation in latitude during this period 
(zenith_sd <- findHEZenith(twl, tol=0.01, range=c(start,end)))

# The angles obtained with in-habitat and Hill-Ekstrom Calibration differ by less than 0.5
# we can use the same zenith angle throughout the annual cycle

# adjust the zenith angles calculated from the breeding sites for the non-breeding sites
zenith0_ad <- zenith0 + abs(zenith - zenith_sd)
zenith_ad  <- zenith_sd

# Find approximate  timing of arrival and departure from the nonbreeding grounds 
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol= 0)

x0_r<- path$x
z0 <- trackMidpts(x0_r)

#Save raw path (no linear interpolation around the equinox)
save(x0_r, file = paste0(dir,"/", geo.id, "_initial_path_raw.csv"))

# Check the following times of arrival and departure using a plot 
arr.nbr <- "2019-10-22" 

# open jpeg
jpeg(paste0(dir, "/", geo.id, "_LatLon_scatterplot.png"), width = 1024, height = 990)

par(mfrow = c(2,1))
plot(twl$Twilight, x0_r[,1], ylab = "longitude")
abline(v = anytime(arr.nbr))
abline(v = fall.equi, col = "orange")

plot(twl$Twilight, x0_r[,2], ylab = "latitude")
abline(v = anytime(arr.nbr))
abline(v = fall.equi, col = "orange")

dev.off()

# Using approximate timings of arrival and departure from the breeding grounds
zenith_twl_zero <- data.frame(Date = twl$Twilight) %>%
  mutate(zenith = case_when(Date < anytime(arr.nbr) ~ zenith0,
                            Date > anytime(arr.nbr) ~ zenith0_ad))

zeniths0 <- zenith_twl_zero$zenith

zenith_twl_med <- data.frame(Date = twl$Twilight) %>%
  mutate(zenith = case_when(Date < anytime(arr.nbr) ~ zenith,
                            Date > anytime(arr.nbr) ~ zenith_sd))

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
abline(v = fall.equi, col = "orange")

plot(twl$Twilight, x0_ad[,2], ylab = "latitude")
abline(v = anytime(arr.nbr))
abline(v = fall.equi, col = "orange")

dev.off()

# Initial Path #################################################################
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zeniths_med, tol=0.08)

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

#SGAT Group model analysis ####################################################

# group twilight times were birds were stationary 
geo_twl <- export2GeoLight(twl)

# Often it is necessary to play around with quantile and days
# quantile defines how many stopovers there are. the higher, the fewer there are
# days indicates the duration of the stopovers 
cL <- changeLight(twl=geo_twl, quantile=0.90, summary = F, days = nday, plot = T)

# merge site helps to put sites together that are separated by single outliers.
mS <- mergeSites2(geo_twl, site = cL$site, distThreshold = 250, degElevation = 90-zenith0, alpha = calib[3:4] , method = "gamma", map = wrld_simpl)

##back transfer the twilight table and create a group vector with TRUE or FALSE according to which twilights to merge 
twl.rev <- data.frame(Twilight = as.POSIXct(geo_twl[,1], geo_twl[,2], tz = "UTC"), 
                      Rise     = c(ifelse(geo_twl[,3]==1, TRUE, FALSE), ifelse(geo_twl[,3]==1, FALSE, TRUE)),
                      Site     = append(rep(mS$site,2), c(0,0)))
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
fixedx[1] <- TRUE # we only fix the breeding location because this is not a complete track 

x0[fixedx,1] <- lon.calib
x0[fixedx,2] <- lat.calib

z0 <- trackMidpts(x0)

# plot stationary locations ####################################################

par(mfrow=c(1,1), mar = c(3,3,3,3))
dtx0 <- as.data.frame(x0)
names(dtx0) <- c("x", "y")

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

#Set limits of the mask
xlim <- range(x0[,1])+c(-5,5)
ylim <- range(x0[,2])+c(-5,5)

index <- ifelse(stationary, 1, 2)
#mask <- earthseaMask(xlim, ylim, n = 10, index=index)
mask <- earthseaMask3(xlim, ylim, res = "lr", index=index, span = 4, twl = twl)

# We will give locations on land a higher prior 
## Define the log prior for x and z
logp <- function(p) {
  f <- mask(p)
  #ifelse(is.na(f), -1000, log(2))
  ifelse(is.na(f), -1000, 100*f)
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
                               zenith = zeniths0,
                               logp.x = logp,# land sea mask
                               fixedx = fixedx)


# define the error shape
x.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(x0))
z.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(z0))

# Fit the model
fit <- estelleMetropolis(model, x.proposal, z.proposal, iters = niter, thin = nthin, chain = chains)

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
                               zenith = zeniths0,
                               fixedx = fixedx)

for (k in 1:3) {
  x.proposal <- mvnorm(chainCov(fit$x), s = 0.3)
  z.proposal <- mvnorm(chainCov(fit$z), s = 0.3)
  fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x),
                           z0 = chainLast(fit$z), iters = niter, thin = nthin, chain = chains)
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
                         z0 = chainLast(fit$z), iters = niter, thin = nthin, chain = 3)

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
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = proj4string(wrld_simpl))

s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- slice(s, sliceIndices(s))

plot(sk, useRaster = F,col = c("transparent", rev(viridis::viridis(50))))
plot(wrld_simpl, xlim=xlim, ylim=ylim,add = T, bg = adjustcolor("black",alpha=0.1))

with(sm[sitenum>0,], arrows(`Lon.50.`, `Lat.2.5.`, `Lon.50.`, `Lat.97.5.`, length = 0, lwd = 2.5, col = "firebrick"))
with(sm[sitenum>0,], arrows(`Lon.2.5.`, `Lat.50.`, `Lon.97.5.`, `Lat.50.`, length = 0, lwd = 2.5, col = "firebrick"))
lines(sm[,"Lon.50."], sm[,"Lat.50."], col = adjustcolor("black", alpha = 0.6), lwd = 2)
points(sm[,"Lon.50."], sm[,"Lat.50."], col = ifelse(sm$StartTime > fall.equi - days(20) & sm$StartTime < fall.equi + days(20), "blue", "darkorchid4"), lwd = 2)

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

# add column with geolocator ID and group 
sm$geo_id <- geo.id

#add a column that categorizes the locations 
sm <- sm %>% mutate(period= case_when(StartTime < anytime("2019-10-24 10:24:07", asUTC = T, tz = "GMT")  ~ "Post-breeding migration",
                                      StartTime >= anytime("2019-10-24 10:24:07", asUTC = T, tz = "GMT") ~ "Non-breeding period"))

#Save the output of the model 
save(sm, file = paste0(dir,"/", geo.id,"_SGAT_GroupedThreshold_summary.csv"))
save(fit, file = paste0(dir,"/", geo.id,"_SGAT_GroupedThreshold_fit.R"))

#load the output of the model 
#load(file = paste0(dir,"/", geo.id,"_SGAT_GroupedThreshold_summary.csv"))
#load(file = paste0(dir,"/", geo.id,"_SGAT_GroupedThreshold_fit.R"))

# Examine twilights ############################################################

#load the adjusted threshold path path x0_ad
load(file = paste0(dir,"/", geo.id, "adjusted_initial_path_raw.csv"))

#Fall transoceanic flight
start <- "2019-09-01"
end <- "2019-11-01"

# Plot lat, lon and light transitions  
jpeg(paste0(dir, "/", geo.id,"_fall_ocean_light_transition.png"), width = 1024 , height = 990, quality = 100, res = 200)

par(cex.lab=1.4)
par(cex.axis=1.4)
par(mfrow=c(3,1), mar = c(5,5,0.1,5))
plot(lig$Date[lig$Date > start & lig$Date < end], lig$Light[lig$Date > start & lig$Date < end], type = "o",
     ylab = "Light level", xlab = "Time")
rect(anytime("2019-10-03"), min(lig$Light)-2, anytime("2019-10-06"), max(lig$Light)+2, col = alpha("yellow", 0.2), lty=0)
rect(anytime("2019-10-16"), min(lig$Light)-2, anytime("2019-10-17"), max(lig$Light)+2, col = alpha("yellow", 0.2), lty=0)

plot(twl$Twilight[twl$Twilight> start & twl$Twilight < end], x0_ad[,1][twl$Twilight > start & twl$Twilight < end],
     ylab = "Longitude", xlab = "Time")
rect(anytime("2019-10-03"), min(x0_ad[,1])-2, anytime("2019-10-06"), max(x0_ad[,1])+2, col = alpha("yellow", 0.2), lty=0)
rect(anytime("2019-10-16"), min(x0_ad[,1])-2, anytime("2019-10-17"), max(x0_ad[,1])+2, col = alpha("yellow", 0.2), lty=0)

plot(twl$Twilight[twl$Twilight > start & twl$Twilight < end], x0_ad[,2][twl$Twilight > start & twl$Twilight < end],
     ylab = "Latitude", xlab = "Time")
rect(anytime("2019-10-03"), min(x0_ad[,2])-2, anytime("2019-10-06"), max(x0_ad[,2])+2, col = alpha("yellow", 0.2), lty=0)
rect(anytime("2019-10-16"), min(x0_ad[,2])-2, anytime("2019-10-17"), max(x0_ad[,2])+2, col = alpha("yellow", 0.2), lty=0)
par(cex.lab= 1)
par(cex.axis= 1)

dev.off()

# A stopover over the Caribbean was detected between 2019-10-06 and 2019-10-16
# Add the new stopover to the location summary obtained at the end of the geolocator analysis
sm.fall.edit <- insertLoc(data = sm,
                          lat.at.loc = 18.747636,
                          start.date = "2019-10-06" ,
                          end.date = "2019-10-16" , 
                          period = "Post-breeding migration",
                          thresh.locs = x0_ad,
                          twl = twl,
                          geo_id = geo.id,
                          sep1 = days(3),
                          sep2 = days (1))

#plot the final stationary locations 
sm.fall.stat <- sm.fall.edit[(sm.fall.edit$sitenum > 0), ]

par(mfrow=c(1,1))

data(wrld_simpl)
plot(wrld_simpl, xlim=xlim, ylim=ylim, col = "grey95")
points(sm.fall.edit$Lon.50., sm.fall.edit$Lat.50., pch = 16, cex = 0, col = "firebrick", type = "o")
points(sm.fall.stat$Lon.50., sm.fall.stat$Lat.50., pch = 16, cex = 1.5, col = "firebrick")

#Save the final location summary
save(sm.fall.edit , file = paste0(dir,"/", geo.id,"_SGAT_GroupedThreshold_summary_fall_edit.csv"))

# Record details for the geolocator analysis ###################################
geo.ref <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv") 
geo.ref[(geo.ref$geo.id == geo.id),]$In_habitat_median_zenith_angle <- zenith
geo.ref[(geo.ref$geo.id == geo.id),]$Hill_Ekstrom_median_angle <- zenith_sd 
geo.ref[(geo.ref$geo.id == geo.id),]$Fall_carrib_edits <- TRUE
geo.ref[(geo.ref$geo.id == geo.id),]$Time_shift_hours <- shift$shift
geo.ref[(geo.ref$geo.id == geo.id),]$nbr.arrival <- arr.nbr
write.csv(geo.ref, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv", row.names=FALSE) 

