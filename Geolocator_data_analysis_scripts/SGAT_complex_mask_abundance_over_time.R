# Create a Land mask for the group model #######################################
#earthseaMask <- function(
    n <-  2 
    pacific <- FALSE 
    xlim <- range(x0[,1])+c(-5,5)
    ylim <- range(x0[,2])+c(-5,5)
    
    index <- ifelse(stationary, 1, 2)
    
  
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
  
  names(ab.ras) <- as.numeric(strftime(names(ab.ras), format = "%j"))
  names(ab.ras)[1] <- 0
  names(ab.ras)[length(names(ab.ras))] <- 366
  
  #project the abundance rasters
  ab.ras.pr <- project(ab.ras, crs(rs), method = "near") 
  
  # get bincodes linking geolocator twilight measurment times to the weeks of  
  # each abundance layer 
  doy <- as.numeric(strftime(twl$Twilight, format = "%j"))
  t.code <- .bincode(doy, as.numeric(names(ab.ras)))
  
  xbin = seq(xmin(ab.ras.pr),xmax(ab.ras.pr),length=ncol(ab.ras.pr)+1)
  ybin = seq(ymin(ab.ras.pr),ymax(ab.ras.pr),length=nrow(ab.ras.pr)+1)
  ab.arr <- as.array(ab.ras.pr)
  
  ab.arr[cbind(length(ybin)-.bincode(dtx0[,2],ybin), .bincode(dtx0[,1],xbin), t.code)]

  
  
  twl.rev2 <- twl.rev %>% filter(Site == 0 | Site != lag(Site))
  
  
  
  # make the movement raster the same resolution as the stationary raster, but allow the bird to go anywhere by giving all cells a value of 1
  rm = rs; rm[] = 1
  
  # stack the movement and stationary rasters on top of each other
  mask = stack(rs, rm)
  
  xbin = seq(xmin(ab.ras.pr ),xmax(ab.ras.pr ),length=ncol(ab.ras.pr)+1)
  ybin = seq(ymin(ab.ras.pr ),ymax(ab.ras.pr ),length=nrow(ab.ras.pr)+1)
  mask = as.array(mask)[,,sort(unique(index)),drop=FALSE]
  
  p = dtx0
  
  function(p) mask[cbind(length(ybin)-.bincode(p[,2],ybin), .bincode(p[,1],xbin), index)]

  
  
  #}

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






a <- anytime(names(ab.ras), asUTC = T, tz = 'UTC')
a[1] <- "2021-01-01"
a[length(a)] <- "2021-12-31"

year(twl$Twilight) <- 2021

b <- .bincode(twl$Twilight, a)
