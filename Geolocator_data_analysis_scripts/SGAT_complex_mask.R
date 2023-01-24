


n = 4

# create empty raster with desired resolution
r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
           xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))

# create a raster for the stationary period, in this case by giving land a value of 1 and sea NA
 maski = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
              rasterize(wrld_simpl, r, 1, silent = TRUE), 
              rasterize(elide(wrld_simpl, shift = c(360, 0)), r, 1, silent = TRUE))

#abundance <- raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geo_spatial_data/bkpwar_abundance_seasonal_full-year_mean_2021.tif")
#abundance_resamp <- projectRaster(abundance, maski, method = "ngb")
#abundance_resamp[abundance_resamp == 0 | is.nan(abundance_resamp)] <- NA
#abundance_resamp[is.nan(abundance_resamp)] <- NA
#abundance_resamp[abundance_resamp > 0 ] <- 1

#maski <- abundance_resamp * maski
 
xbin = seq(xmin(maski),xmax(maski),length=ncol(maski)+1)
ybin = seq(ymin(maski),ymax(maski),length=nrow(maski)+1)

maskx <- as.matrix(maski)

## Define the log prior for x and z
log.prior <- function(p) {
  f <- mask(p)
  ifelse(is.na(f), log(1), f) 
}

#modified matrix subset approach 
maski[cbind(.bincode(lat.calib,ybin), .bincode(lon.calib,xbin))]

#original matrix subset approach
maski[cbind(length(ybin) -.bincode(lat.calib,ybin), .bincode(lon.calib,xbin))]

#breeding site location 
maski[cbind(length(ybin) -.bincode(x0[,2],ybin), .bincode(x0[,1],xbin))] <-2
maski[cbind(.bincode(x0[,2],ybin), .bincode(x0[,1],xbin))] <- 1

plot(maski)

getValues(mask)
