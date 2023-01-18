n = 4

# create empty raster with desired resolution
r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
           xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))

# create a raster for the stationary period, in this case by giving land a value of 1 and sea NA
 mask = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
              rasterize(wrld_simpl, r, 1, silent = TRUE), 
              rasterize(elide(wrld_simpl, shift = c(360, 0)), r, 1, silent = TRUE))

abundance <- raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geo_spatial_data/bkpwar_abundance_seasonal_full-year_mean_2021.tif")
abundance_resamp <- projectRaster(abundance, mask, method = "ngb")
abundance_resamp[abundance_resamp == 0 | is.nan(abundance_resamp)] <- NA

mask <- abundance_resamp * mask
 
xbin = seq(xmin(mask),xmax(mask),length=ncol(mask)+1)
ybin = seq(ymin(mask),ymax(mask),length=nrow(mask)+1)


## Define the log prior for x and z
log.prior <- function(p) {
  f <- mask(p)
  ifelse(is.na(f), log(1), f) 
}

function(p) mask[cbind(.bincode(58.4176050,ybin),.bincode(-93.74039,xbin))]

getValues(mask)
