# Land mask for the SGAT threshold model 

# Land mask ####################################################################
#earthseaMask <- function(xlim, ylim, n = 2, pacific=FALSE) {

  n <-  2 
  pacific <- FALSE 
  xlim <- range(x0[,1])+c(-5,5)
  ylim <- range(x0[,2])+c(-5,5)
  
  if (pacific) { wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  #load weekly rasters of blackpoll warbler abundance
  ab.ras <- load_raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports/2021/bkpwar",
                        product = "abundance",
                        period = "weekly",
                        resolution = "lr")
  
  names(ab.ras) <- as.numeric(strftime(names(ab.ras), format = "%j"))
  names(ab.ras)[1] <- 0
  names(ab.ras)[length(names(ab.ras))] <- 366
  
  #project the abundance rasters
  ab.ras.pr <- project(ab.ras, as.character(crs(r)), method = "near") 
  
  values(ab.ras.pr)[is.nan(values(ab.ras.pr))] <- 0
  
  xbin = seq(xmin(ab.ras.pr),xmax(ab.ras.pr),length=ncol(ab.ras.pr)+1)
  ybin = seq(ymin(ab.ras.pr),ymax(ab.ras.pr),length=nrow(ab.ras.pr)+1)
  ab.arr <- as.array(ab.ras.pr)
  
  function(p) {
    
    # get bincodes linking geolocator twilight measurement times to the weeks of  
    # each abundance layer
    # This has to be done here because x0 and z0 (midpoints) can have different lengths 
    times <- seq(from = min(twl$Twilight), to = max(twl$Twilight), length.out = nrow(p))
    doy <- as.numeric(strftime(times, format = "%j"))
    t.code <- .bincode(doy, as.numeric(names(ab.ras)))
    
    ab.arr[cbind(length(ybin) -.bincode(p[,2],ybin),.bincode(p[,1],xbin), t.code)]
    
  } 
  
#}

xlim <- range(x0[,1]+c(-5,5))
ylim <- range(x0[,2]+c(-5,5))

mask <- earthseaMask(xlim, ylim, n = 4)

## Define the log prior for x and z
log.prior <- function(p) {
  f <- mask(p)
}


# Land mask for the SGAT group threshold model #################################
#earthseaMask <- function(
    n <-  2 
    pacific <- FALSE 
    xlim <- range(x0[,1])+c(-5,5)
    ylim <- range(x0[,2])+c(-5,5)
    
    index <- ifelse(stationary, T, F)
  
  if (pacific) { wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1
  rs = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
             rasterize(wrld_simpl, r, 1, silent = TRUE), 
             rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  
  #load weekly rasters of blackpoll warbler abundance
  ab.ras <- load_raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports/2021/bkpwar",
                        product = "abundance",
                        period = "weekly",
                        resolution = "lr")
  
  names(ab.ras) <- as.numeric(strftime(names(ab.ras), format = "%W"))
  
  #project the abundance rasters
  ab.ras.pr <- project(ab.ras, as.character(crs(rs)), method = "near") 
  
  values(ab.ras.pr)[is.nan(values(ab.ras.pr))] <- NA
  
  # get bincodes linking geolocator twilight measurement times to the weeks of  
  # each abundance layer 
  t.datetime <- (twl %>% filter(group != lag(group, default = -1)))$Twilight
  t.weeks <- week(t.datetime)
  t.code <- .bincode(t.weeks, as.numeric(names(ab.ras)))
  
  xbin = seq(xmin(ab.ras.pr),xmax(ab.ras.pr),length=ncol(ab.ras.pr)+1)
  ybin = seq(ymin(ab.ras.pr),ymax(ab.ras.pr),length=nrow(ab.ras.pr)+1)
  ab.arr <- as.array(ab.ras.pr)
  
  # If the bird becomes stationary, the prior for a location is selected from an eBird 
  # relative abundance raster for the week during which the bird arrived at the location 
  # While the bird is moving, the prior always has a value of 0
  function(p){
    
    ifelse(stationary, 
           ab.arr[cbind(length(ybin)-.bincode(p[,2],ybin), .bincode(p[,1],xbin), t.code)],
           0)
  }
  #}
  
  mask <- earthseaMask(xlim, ylim, n = 10, index=index)
  
  # We will give locations on land a higher prior 
  ## Define the log prior for x and z
  logp <- function(p) {
    f <- mask(p)
    ifelse(is.na(f), -1000, f)
  }

  
  
# code trials ##################################################################
.bincode(week(anytime("2019-02-04")), as.numeric(names(ab.ras)))

ab.arr[cbind(length(ybin)-.bincode(-72.959498, ybin), .bincode(-72.959498,xbin), 4)]

x <- raster(ab.arr[,,4], crs = rs)
values(x)[values(x) > 0 & !is.na(values(x))] <- 1
plot(x)


p <- x0
coord <- cbind(length(ybin)-.bincode(p[,2],ybin), .bincode(p[,1],xbin))
               
ab.ras.pr2 <- ab.ras.pr$"41"         
values(ab.ras.pr2)[values(ab.ras.pr2) > 0 & !is.na(values(ab.ras.pr2))] <- 1
plot(ab.ras.pr2)

ab.ras.pr2[cbind(length(ybin)-.bincode(p[,2],ybin), .bincode(p[,1],xbin))] <- 10 
plot(ab.ras.pr2)