# source: unpublished data 
# tag number: 4068-014
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
# install_github("SLisovski/GeoLocTools")
library(GeoLocTools)
#setupGeolocation()

# clear object from workspace
rm(list=ls())

# Load helper functions 
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

geo.id <- "4068_014"

# data directory
dir <- paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data/", geo.id)

# read file with consolidated geolocator data
ref_data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/Geolocator_reference_data_consolidated.csv")

#Assign geolocator deployment site 
lat.calib <- ref_data$deploy.latitude[which(ref_data$geo.id == geo.id)]
lon.calib <- ref_data$deploy.longitude[which(ref_data$geo.id == geo.id)]

# time of deployment (from reference file)
#deploy.start <- anytime("", asUTC = T, tz = "GMT")

# time of recovery (estimate from light data)
#deploy.end <- anytime("", asUTC = T, tz = "GMT")

#Equinox times
fall.equi <- anytime("2016-09-22", asUTC = T, tz = "GMT")
spring.equi <- anytime("2017-03-20", asUTC = T, tz = "GMT")

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

# This geolocator has no time shift 

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

# load some parameters for the analysis 
geo.ref <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/Geolocator_reference_data_consolidated.csv") 
days <- geo.ref[(geo.ref$geo.id == geo.id),]$changeLight.days
dist <- geo.ref[(geo.ref$geo.id == geo.id),]$mergesites.distance
stat.nbr.lim <- geo.ref[(geo.ref$geo.id == geo.id),]$stat.nbr.limit 

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

# Save edited twilights 
write.csv(twl, paste0(dir,"/",geo.id , "_twl_times_edited.csv"))

# Calibration ##################################################################

# We start with calibration based on the stationary periods before and after the migration
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

#calibration period before the migration 
tm.calib <- as.POSIXct(c("2016-08-01", "2016-08-20"), tz = "UTC")

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

#convert to geolight format
geo_twl <- export2GeoLight(twl)

# this is just to find places where birds have been for a long time, would not use these parameters for stopover identification, detailed can be found in grouped model section
cL <- changeLight(twl=geo_twl, quantile=0.8, summary = F, days = 10, plot = T)
# merge site helps to put sites together that are separated by single outliers.
mS <- mergeSites(twl = geo_twl, site = cL$site, degElevation = 90-zenith0, distThreshold = 250)

#specify which site is the stationary one
site           <- mS$site[mS$site>0] # get rid of movement periods
stationarySite <- which(table(site) == max(table(site))) # find the site where bird is the longest

#find the dates that the bird arrives and leaves this stationary site
start <- min(which(mS$site == stationarySite))
end   <- max(which(mS$site == stationarySite))

(zenith_sd <- findHEZenith(twl, tol=0.01, range=c(start,end)))

# The Hill-ekstrom zenith and in-habitat zenith angles differ by less than 0.5 degrees
# Alternative calibration is not necessary 

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
arr.nbr <- "2016-10-12" 
dep.nbr <- "2017-04-18" 

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
                            Date > anytime(arr.nbr) & Date < anytime(dep.nbr) ~ zenith0_ad+2.5,
                            Date > anytime(dep.nbr) ~ zenith0_ad + 1.5))

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
tol_ini <- 0.1
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zeniths_med, tol = tol_ini)

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
q = 0.86
cL <- changeLight(twl=geo_twl, quantile= q, summary = F, days = days, plot = T)

# merge site helps to put sites together that are separated by single outliers.
mS <- mergeSites(twl = geo_twl, site = cL$site, degElevation = 90-zenith0, distThreshold = dist)

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
fixedx[1] <- TRUE
fixedx[c(1, length(fixedx))] <- TRUE

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
#Set limits of the mask
xlim <- range(x0[,1])+c(-5,5)
ylim <- range(x0[,2])+c(-5,5)

index <- ifelse(stationary, 1, 2)
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
                               zenith = zeniths0,
                               logp.x = logp,# land sea mask
                               fixedx = fixedx)


# define the error shape
x.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(x0))
z.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(z0))

#Tuning ########################################################################

# Fit a first chain for tuning
fit <- estelleMetropolis(model, x.proposal, z.proposal, iters = 1000, thin = 20)

# Fit additional chains for tuning
for (k in 1:2) {
  x.proposal <- mvnorm(chainCov(fit$x), s = 0.3)
  z.proposal <- mvnorm(chainCov(fit$z), s = 0.3)
  
  # get Median of chains
  chain.sm.x <- SGAT2Movebank(fit$x)
  x.med <- list(as.matrix(chain.sm.x[,c("Lon.50%","Lat.50%")]))
  chain.sm.z <- SGAT2Movebank(fit$z)
  z.med <- list(as.matrix(chain.sm.z[,c("Lon.50%","Lat.50%")]))
  
  fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = x.med,
                           z0 = z.med, iters = 1000, thin = 20)
}

## Check if chains mix
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!fixedx, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!fixedx, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)

#Final model run ###############################################################

# get Median of chains
chain.sm.x <- SGAT2Movebank(fit$x)
x.med <- list(as.matrix(chain.sm.x[,c("Lon.50%","Lat.50%")]))
chain.sm.z <- SGAT2Movebank(fit$z)
z.med <- list(as.matrix(chain.sm.z[,c("Lon.50%","Lat.50%")]))

# get proposals 
x.proposal <- mvnorm(chainCov(fit$x), s = 0.3)
z.proposal <- mvnorm(chainCov(fit$z), s = 0.3)

fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = x.med,
                         z0 = z.med , iters = 3000, thin = 20, chain = 1)

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

# Stationary locations from the threshold estimates 
points(dtx0[sitenum > 0,], pch = 16, cex = 1, col = "green")

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

# add column with geolocator ID
sm$geo_id <- geo.id

# find time of establishment and departure from the nonbreeding grounds 
arr.nbr.sgat <- sm %>% filter(Lat.50. < 12 & sitenum > 0 & duration > stat.nbr.lim) %>% 
  first(.$StartTime) %>% .$StartTime
dep.nbr.sgat <- sm %>% filter(Lat.50. < 12 & sitenum > 0 & duration > stat.nbr.lim) %>%
  last(.$EndTime)%>% .$EndTime

#add a column that categorizes the locations (based on the groupthreshold model output)
sm <- sm %>% rowwise() %>% mutate(period= case_when(StartTime < anytime(arr.nbr.sgat, asUTC = T, tz = "GMT")  ~ "Post-breeding migration",
                                                    StartTime >= anytime(arr.nbr.sgat, asUTC = T, tz = "GMT") & StartTime < anytime(dep.nbr.sgat, asUTC = T, tz = "GMT") ~ "Non-breeding period",
                                                    StartTime > anytime(dep.nbr.sgat , asUTC = T, tz = "GMT") ~ "Pre-breeding migration"))
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
start <- "2016-09-20"
end <- "2016-10-20"

# flight
f1.start <- "2016-10-10"
f1.end <- "2016-10-12"

# Plot lat, lon and light transitions  
jpeg(paste0(dir, "/", geo.id,"_fall_ocean_light_transition.png"), width = 1024 , height = 990, quality = 100, res = 200)

par(cex.lab=1.4)
par(cex.axis=1.4)
par(mfrow=c(3,1), mar = c(5,5,0.1,5))
plot(lig$Date[lig$Date > start & lig$Date < end], lig$Light[lig$Date > start & lig$Date < end], type = "o",
     ylab = "Light level", xlab = "Time")
rect(anytime(f1.start), min(lig$Light)-2, anytime(f1.end), max(lig$Light)+2, col = alpha("yellow", 0.2), lty=0)

plot(twl$Twilight[twl$Twilight> start & twl$Twilight < end], x0_ad[,1][twl$Twilight > start & twl$Twilight < end],
     ylab = "Longitude", xlab = "Time")
rect(anytime(f1.start), min(x0_ad[,1])-2, anytime(f1.end), max(x0_ad[,1])+2, col = alpha("yellow", 0.2), lty=0)

plot(twl$Twilight[twl$Twilight > start & twl$Twilight < end], x0_ad[,2][twl$Twilight > start & twl$Twilight < end],
     ylab = "Latitude", xlab = "Time")
rect(anytime(f1.start), min(x0_ad[,2])-2, anytime(f1.end), max(x0_ad[,2])+2, col = alpha("yellow", 0.2), lty=0)
par(cex.lab= 1)
par(cex.axis= 1)

dev.off()

# A transoceanic flight occured between 2016-10-10 and 2019-10-12 

# Estimate timing of departure and arrival from the breeding and nonbreeding grounds ############################################################
dep.br <- "2016-08-18 13:22"
arr.br <- "2017-05-29 8:59"

par(mfrow=c(2,1))
plot(twl$Twilight, type  = "l", x0_ad[,1])
abline(v = anytime(dep.br))
abline(v = anytime(arr.br))
plot(twl$Twilight, type  = "l", x0_ad[,2])
abline(v = anytime(dep.br))
abline(v = anytime(arr.br))
par(mfrow=c(1,1))

# Record details for the geolocator analysis ###################################
geo.ref <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/Geolocator_reference_data_consolidated.csv") 
geo.ref[(geo.ref$geo.id == geo.id),]$In_habitat_median_zenith_angle <- zenith
geo.ref[(geo.ref$geo.id == geo.id),]$Hill_Ekstrom_median_angle <- zenith_sd 
geo.ref[(geo.ref$geo.id == geo.id),]$nbr.arrival <- arr.nbr
geo.ref[(geo.ref$geo.id == geo.id),]$nbr.departure <- dep.nbr
geo.ref[(geo.ref$geo.id == geo.id),]$IH.calib.start <- as.character(tm.calib[1])
geo.ref[(geo.ref$geo.id == geo.id),]$IH.calib.end <- as.character(tm.calib[2])
geo.ref[(geo.ref$geo.id == geo.id),]$tol <-tol_ini
geo.ref[(geo.ref$geo.id == geo.id),]$nbr.arrival <- as.Date(arr.nbr.sgat)
geo.ref[(geo.ref$geo.id == geo.id),]$nbr.departure <- as.Date(dep.nbr.sgat)
geo.ref[(geo.ref$geo.id == geo.id),]$br.departure <- as.Date(dep.br)
geo.ref[(geo.ref$geo.id == geo.id),]$br.arrival <- as.Date(arr.br)
geo.ref[(geo.ref$geo.id == geo.id),]$changelight.quantile <- q
write.csv(geo.ref, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/Geolocator_reference_data_consolidated.csv", row.names=FALSE) 

# Examine twilights ############################################################

#################################################################################
# Simple thresold estimate
#################################################################################

# movement model with slower speed ##############################################
beta  <- c(2.2, 0.08)
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

# Perform threshold model without stationary location ###########################
twl$index <- seq(1, length(twl$group))
x0 <- path$x
z0 <- trackMidpts(x0)

model <- groupedThresholdModel(twl$Twilight,
                               twl$Rise,
                               group = twl$index, 
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

#Tuning ########################################################################

# fit.s a first chain for tuning
fit.s <- estelleMetropolis(model, x.proposal, z.proposal, iters = 300, thin = 20)

# fit.s additional chains for tuning
for (k in 1:2) {
  x.proposal <- mvnorm(chainCov(fit.s$x), s = 0.3)
  z.proposal <- mvnorm(chainCov(fit.s$z), s = 0.3)
  
  # get Median of chains
  chain.sm.s.x <- SGAT2Movebank(fit.s$x)
  x.med <- list(as.matrix(chain.sm.s.x[,c("Lon.50%","Lat.50%")]))
  chain.sm.s.z <- SGAT2Movebank(fit.s$z)
  z.med <- list(as.matrix(chain.sm.s.z[,c("Lon.50%","Lat.50%")]))
  
  fit.s <- estelleMetropolis(model, x.proposal, z.proposal, x0 = x.med,
                             z0 = z.med, iters = 300, thin = 20)
}

## Check if chains mix
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit.s$x[[1]][!fixedx, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit.s$x[[1]][!fixedx, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)

#Final model run ###############################################################

# get Median of chains
chain.sm.s.x <- SGAT2Movebank(fit.s$x)
x.med <- list(as.matrix(chain.sm.s.x[,c("Lon.50%","Lat.50%")]))
chain.sm.s.z <- SGAT2Movebank(fit.s$z)
z.med <- list(as.matrix(chain.sm.s.z[,c("Lon.50%","Lat.50%")]))

# get proposals 
x.proposal <- mvnorm(chainCov(fit.s$x), s = 0.3)
z.proposal <- mvnorm(chainCov(fit.s$z), s = 0.3)

fit.s <- estelleMetropolis(model, x.proposal, z.proposal, x0 = x.med,
                           z0 = z.med , iters = 1000, thin = 20, chain = 1)

# Summarize and plot the results ########################################################
sm.s <- locationSummary(fit.s$z, time=fit.s$model$time)

# open jpeg
jpeg(paste0(dir, "/", geo.id,"_simple_threshold_model_map.png"), width = 1024 , height = 990)

# empty raster of the extent
r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = proj4string(wrld_simpl))

s <- slices(type = "intermediate", breaks = "week", mcmc = fit.s, grid = r)
sk <- slice(s, sliceIndices(s))

plot(sk, useRaster = F,col = rev(viridis::viridis(50)))
plot(wrld_simpl, xlim=xlim, ylim=ylim,add = T, bg = adjustcolor("black",alpha=0.1))

lines(sm.s[,"Lon.50%"], sm.s[,"Lat.50%"], col = adjustcolor("firebrick", alpha.f = 0.6), type = "o", pch = 16)

dev.off()

# open jpeg
jpeg(paste0(dir, "/", geo.id,"_threshold_model_lon_lat.png"), width = 1024 , height = 990)

par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(sm.s$Time1, sm.s$"Lon.50%", ylab = "Longitude", xlab = "", yaxt = "n", type = "n", ylim = xlim)
axis(2, las = 2)
polygon(x=c(sm.s$Time1,rev(sm.s$Time1)), y=c(sm.s$`Lon.2.5%`,rev(sm.s$`Lon.97.5%`)), border="gray", col="gray")
lines(sm.s$Time1,sm.s$"Lon.50%", lwd = 2)

plot(sm.s$Time1,sm.s$"Lat.50%", type="n", ylab = "Latitude", xlab = "", yaxt = "n", ylim = ylim)
axis(2, las = 2)
polygon(x=c(sm.s$Time1,rev(sm.s$Time1)), y=c(sm.s$`Lat.2.5%`,rev(sm.s$`Lat.97.5%`)), border="gray", col="gray")
lines(sm.s$Time1,sm.s$"Lat.50%", lwd = 2)

par(mfrow=c(1,1),mar=c(1,1,1,1))

dev.off()


# Save the threshold model estimate ############################################
save(sm.s, file = paste0(dir,"/", geo.id,"_SGAT_Threshold_summary.csv"))
save(fit.s, file = paste0(dir,"/", geo.id,"_SGAT_Threshold_fit.s.R"))
