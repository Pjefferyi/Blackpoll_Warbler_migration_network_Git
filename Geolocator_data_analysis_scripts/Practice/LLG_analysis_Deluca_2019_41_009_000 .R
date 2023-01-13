require(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(devtools)
library(remotes)
library(anytime)
library(ggmap)
library(FLightR)
library(LLmig)

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

# Twilight annotation ==========================================================

#import lig data and remove first row 
ligdata <- readLig("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Data_raw/Deluca_et_al_2019/churchill/Geo_41_009_000.lig")
ligdata <- ligdata[2:nrow(ligdata),]

#latitude and longitude at the location where bird A was tagged: 
lat.calib <- 	58.731288
lon.calib <- -93.82587
  
# If needed, convert dates to POSIXct (UTC/GMT) and rename columns 
ligdata <- ligdata %>% mutate(Date = anytime(Date))
  
# plot light levels for a few days
threshold <- 2 #Note, this threshold is also used for to annotate twilights woth preprocessLight 
with(ligdata[45000:50000,], plot(Date, Light, type = "o"))
abline(h=threshold, col="orange", lty = 2, lwd = 2)

# plot light levels over the entire deployment period
offset = 8# adjusts the y-axis to put night (dark shades) in the middle
  
  lightImage( tagdata = ligdata,
              offset = offset,     
              zlim = c(0, 20),
              xlab = "Date")

tsimageDeploymentLines(ligdata$Date, lon = lon.calib, lat = lat.calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

#Twilight annotation
protwl <- preprocessLight(tagdata = ligdata, threshold = threshold, offset = 15)

twledit <- twilightEdit(twilights = protwl,
                        offset = offset,
                        window = 4,           # two days before and two days after
                        outlier.mins = 35,    # difference in mins
                        stationary.mins = 15, # are the other surrounding twilights within 25 mins of one another
                        plot = TRUE)

#save output (fill in the end of the filepaths)
write.csv(twledit, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_analysis_intermediate_data/Twilight_times/Twilight_times_Deluca_2019_Geo_41_009_000.csv")

# SGAT analysis Starts =========================================================
twl<- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_analysis_intermediate_data/Twilight_times/Twilight_times_Deluca_2019_Geo_41_009_000.csv")
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")

#inspect the data 
offset <- 16 # adjusts the y-axis to put night (dark shades) in the middle
  
  lightImage( tagdata = ligdata,
              offset = offset,     
              zlim = c(0, 4))

tsimagePoints(twl$Twilight, offset = offset, pch = 16, cex = 1.2,
              col = ifelse(twl$Deleted, "grey20", ifelse(twl$Rise, "firebrick", "cornflowerblue")))

#remove rows that have been marked for deletion
twl <- subset(twl, !Deleted)

#Calibration====================================================================
lightImage( tagdata = ligdata,
            offset = offset,     
            zlim = c(0, 20))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

tm.calib <- as.POSIXct(c("2016-06-20", "2016-09-01"), tz = "UTC")
abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2]) #subset of data used for the calibration 

calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

#Zenith angle at the zero deviation, the median deviation and parameters of the error distribution 
zenith  <- calib[1]
zenith0 <- calib[2]

alpha <- calib[3:4]

#movement model
beta  <- c(2.2, 0.08)
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

#initial path===================================================================
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol=0.01)

x0 <- path$x
z0 <- trackMidpts(x0)

data(wrld_simpl)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey95", add = T)

points(path$x, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()

#Define known locations=========================================================
fixedx <- rep(F, nrow(x0))
fixedx[1:2] <- T # first two location estimates

fixedx[(nrow(x0) - 1):nrow(x0)] <- T # last two location estimates

x0[fixedx, 1] <- lon.calib
x0[fixedx, 2] <- lat.calib

z0 <- trackMidpts(x0) # we need to update the z0 locations

#Land mask =====================================================================
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

## Define the log prior for x and z
log.prior <- function(p) {
  f <- mask(p)
  ifelse(f | is.na(f), log(2), log(1))
}

#Estelle model =================================================================
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

proposal.x <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0))

fit <- estelleMetropolis(model, proposal.x, proposal.z, iters = 1000, thin = 20)

#Tuning the proposals===========================================================
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
                        zenith = zenith0,
                        fixedx = fixedx)

x.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl))
z.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl) - 1)

for (k in 1:3) {
  fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x), 
                           z0 = chainLast(fit$z), iters = 300, thin = 20)
  
  x.proposal <- mvnorm(chainCov(fit$x), s = 0.2)
  z.proposal <- mvnorm(chainCov(fit$z), s = 0.2)
}

#examine the samples
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!fixedx, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!fixedx, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)

#Final run
x.proposal <- mvnorm(chainCov(fit$x), s = 0.25)
z.proposal <- mvnorm(chainCov(fit$z), s = 0.25)

fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x), 
                         z0 = chainLast(fit$z), iters = 1000, thin = 20)

#Summarize the results==========================================================
#all locations 
sm <- locationSummary(fit$z, time=fit$model$time)
head(sm)

#Plot the results ==============================================================
# empty raster of the extent
r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = proj4string(wrld_simpl))

s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- slice(s, sliceIndices(s))

plot(sk, useRaster = F,col = rev(viridis::viridis(50)))
plot(wrld_simpl, xlim=xlim, ylim=ylim,add = T, bg = adjustcolor("black",alpha=0.1))

lines(sm[,"Lon.50%"], sm[,"Lat.50%"], col = adjustcolor("firebrick", alpha.f = 0.6), type = "o", pch = 16)

# MigSchedule ==================================================================
s <- slices(type = "intermediate", breaks = "day", mcmc = fit, grid = r)
sites <- MigSchedule(s, rm.lat.equinox = TRUE, collapseSites = TRUE)
