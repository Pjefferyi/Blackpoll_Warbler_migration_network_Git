require(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(devtools)
library(remotes)
library(anytime)
library(ggmap)
library(FLightR)

#Sys.getenv("GITHUB_PAT")
#Sys.unsetenv("GITHUB_PAT")

#instal SGAT (dependency for TWGeos)
#devtools::install_github("SWotherspoon/SGAT")
library(SGAT)

#install the R package TwGeos from github to define twilight events 
#devtools::install_github("SLisovski/TwGeos")
library(TwGeos)

# TwGeos manual: https://rdrr.io/github/slisovski/TwGeos/

#other important packages 
#install_github("SLisovski/GeoLocTools")
library(GeoLocTools)
setupGeolocation()

#generate migration schedule
#install.packages("remotes")
#remotes::install_github("MTHallworth/LLmig")
library(LLmig)

#get data for ML6740 V8296 017 
ligdata <- readLig("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Data_raw/BLPW_geo_2020/breed/ML6740 V8296 017/ML6740 V8296 017 reconstructed_000.lig", skip = 0)

ligdata <- ligdata[2:nrow(ligdata),]

#latitude and longitude at the location where bird A was tagged: 
lat.calib <- 47.38355
lon.calib <- -71.08777

#rename columns for dates and light level 
#Aldata <- rename(ligdata, Date = timestamp, Light = gls.light.level) %>%
# convert Dates from character to POSIXct
# set asUTC to true to avoid dayligth saving times 
#    mutate(Date = anytime(Date, asUTC = TRUE, tz = "UTC") + 4*60*60)

# convert  dates to UTC
ligdata$Date <- anytime(ligdata$Date, asUTC = TRUE, tz = "UTC") - 1*60*60

#plot light levels for a few days
threshold <- 3 # slightly above night recordings
with(ligdata[45000:50000,], plot(Date, Light, type = "o"))
abline(h=threshold, col="orange", lty = 2, lwd = 2)

#another plot 
offset = -4# adjusts the y-axis to put night (dark shades) in the middle

lightImage( tagdata = ligdata,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(ligdata$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

#Twilight annotation
#twl <- preprocessLight(tagdata = ligdata, threshold = threshold, offset = 15)

# twledit <- twilightEdit(twilights = twl,
#                         offset = offset,
#                         window = 4,           # two days before and two days after
#                         outlier.mins = 45,    # difference in mins
#                         stationary.mins = 15, # are the other surrounding twilights within 25 mins of one another
#                         plot = TRUE)

#save output 
#write.csv(twl,"C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Blackpoll_twl_data/Twilight_times/Twilight_times_ML6740_V8296_017.csv")
twlV8296_017 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_analysis_intermediate_data/Twilight_times/Twilight_times_ML6740_V8296_017.csv")

# SGAT analysis ================================================================
twlV8296_017$Twilight <- as.POSIXct(twlV8296_017$Twilight, tz = "UTC")

#filtering twilight times
twlV8296_017 <- twilightEdit(twilights = twlV8296_017,
                      offset = offset,
                      window = 4,           # two days before and two days after
                      outlier.mins = 45,    # difference in mins
                      stationary.mins = 15, # are the other surrounding twilights within 25 mins of one another
                      plot = TRUE)

#inspect the data 
lightImage( tagdata = ligdata,
            offset = offset,     
            zlim = c(0, 20))

tsimagePoints(twlV8296_017$Twilight, offset = offset, pch = 16, cex = 1.2,
              col = ifelse(twlV8296_017$Rise, "firebrick", "cornflowerblue"))

#Select a calibration period====================================================
lightImage( tagdata = ligdata,
            offset = offset,     
            zlim = c(0, 20))

#twilight times in the deployment area
tsimageDeploymentLines(twlV8296_017$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

tm.calib <- as.POSIXct(c("2019-07-14", "2019-09-30"), tz = "UTC")
abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

d_calib <- subset(twlV8296_017, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])

#perform calibration 
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

#parameters for the analysis====================================================
zenith  <- calib[1]
zenith0 <- calib[2]
alpha <- calib[3:4]

#movement model ================================================================
beta  <- c(2.2, 0.08) #parameters for a gamma distribution of flight speeds
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

#initial path ===================================================================
path <- thresholdPath(twlV8296_017$Twilight, twlV8296_017$Rise, zenith = zenith, tol=0.01)

x0 <- path$x
z0 <- trackMidpts(x0)

data(wrld_simpl)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(path$x, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()

#define known locations =========================================================
fixedx <- rep(F, nrow(x0))
fixedx[1:2] <- T # first two location estimates

fixedx[(nrow(x0) - 1):nrow(x0)] <- T # last two location estimates

x0[fixedx, 1] <- lon.calib
x0[fixedx, 2] <- lat.calib

z0 <- trackMidpts(x0) # we need to update the z0 locations

#Create a Land mask ============================================================

earthseaMask <- function(xlim, ylim, n = 2, pacific=FALSE) {
  
  if (pacific) { wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1 and sea NA
  mask = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
               rasterize(wrld_simpl, r, 1, silent = TRUE), 
               rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  
  xbin = seq(xmin(mask),xmax(mask),length=ncol(mask)+1)
  ybin = seq(ymin(mask),ymax(mask),length=nrow(mask)+1)
  
  function(p) mask[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin))]
}

xlim <- range(x0[,1]+c(-5,5))
ylim <- range(x0[,2]+c(-5,5))

mask <- earthseaMask(xlim, ylim, n = 1)


## Define the log prior for x and z (give values on land a higher prior). The prior is defined on the logscale
log.prior <- function(p) {
  f <- mask(p)
  ifelse(f | is.na(f), log(2), log(1))
}

#The Estelle model (Uses teh threshold method) =================================
model <- thresholdModel(twilight = twlV8296_017$Twilight,
                        rise = twlV8296_017$Rise,
                        twilight.model = "ModifiedGamma",
                        alpha = alpha,
                        beta = beta,
                        logp.x = log.prior, logp.z = log.prior,
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = fixedx)

proposal.x <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0))

fit <- estelleMetropolis(model, proposal.x, proposal.z, iters = 1000, thin = 20)

#Tuning the proposals ==========================================================
x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

model <- thresholdModel(twilight = twlV8296_017$Twilight,
                        rise = twlV8296_017$Rise,
                        twilight.model = "Gamma",
                        alpha = alpha,
                        beta = beta,
                        logp.x = log.prior, logp.z = log.prior, 
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = fixedx)

x.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twlV8296_017))
z.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twlV8296_017) - 1)

#short runs
for (k in 1:3) {
  fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x), 
                           z0 = chainLast(fit$z), iters = 300, thin = 20)
  
  x.proposal <- mvnorm(chainCov(fit$x), s = 0.2)
  z.proposal <- mvnorm(chainCov(fit$z), s = 0.2)
}

#examine chain mixing 
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!fixedx, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!fixedx, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)

#Final run =====================================================================
x.proposal <- mvnorm(chainCov(fit$x), s = 0.25)
z.proposal <- mvnorm(chainCov(fit$z), s = 0.25)

fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x), 
                         z0 = chainLast(fit$z), iters = 1000, thin = 20)

#Summarize the results =========================================================
sm <- locationSummary(fit$z, time=fit$model$time)
head(sm)

# Plot the results =============================================================

# empty raster of the extent
r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = proj4string(wrld_simpl))

s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- slice(s, sliceIndices(s))

plot(sk, useRaster = F,col = rev(viridis::viridis(50)))
plot(wrld_simpl, xlim=xlim, ylim=ylim,add = T, bg = adjustcolor("black",alpha=0.1))

lines(sm[,"Lon.50%"], sm[,"Lat.50%"], col = adjustcolor("firebrick", alpha.f = 0.6), type = "o", pch = 16)

# MigSchedule ==================================================================
library(LLmig)
s2 <- slices(type = "intermediate", breaks = "day", mcmc = fit, grid = r)
sites <- MigSchedule(s2, plot = FALSE)

# Plot longitude and latitude separately =======================================

par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(sm$Time1, sm$"Lon.50%", ylab = "Longitude", xlab = "", yaxt = "n", type = "n", ylim = c(-85, -55))
axis(2, las = 2)
polygon(x=c(sm$Time1,rev(sm$Time1)), y=c(sm$`Lon.2.5%`,rev(sm$`Lon.97.5%`)), border="gray", col="gray")
lines(sm$Time1,sm$"Lon.50%", lwd = 2)

plot(sm$Time1,sm$"Lat.50%", type="n", ylab = "Latitude", xlab = "", yaxt = "n", ylim = c(-20,90))
axis(2, las = 2)
polygon(x=c(sm$Time1,rev(sm$Time1)), y=c(sm$`Lat.2.5%`,rev(sm$`Lat.97.5%`)), border="gray", col="gray")
lines(sm$Time1,sm$"Lat.50%", lwd = 2)
par(mfrow=c(1,1),mar=c(4,4,1,1))


# Groupe Model #================================================================

geo_twl <- export2GeoLight(twl)

# Often it is necessary to play around with quantile and days
# quantile defines how many stopovers there are. the higher, the fewer there are
# days indicates the duration of the stopovers 
cL <- changeLight(twl=geo_twl, quantile=0.86, summary = F, days = 2, plot = T)

# merge site helps to put sites together that are separated by single outliers.
mS <- mergeSites(twl = geo_twl, site = cL$site, degElevation = 90-zenith0, distThreshold = 500)



