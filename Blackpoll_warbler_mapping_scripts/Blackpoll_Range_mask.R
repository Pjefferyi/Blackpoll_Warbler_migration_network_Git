xlim <- range(x0[,1]+c(-5,5))
ylim <- range(x0[,2]+c(-5,5))
n = 2
pacific=FALSE

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

