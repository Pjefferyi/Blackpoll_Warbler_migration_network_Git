# mapping Blackpoll warbler Distributions 

#Load the helper functions (and associated packages)
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis/Geolocator_analysis_helper_functions.R")

# Plot of nonbreeding probability densities ###################################

#limits for plotting 
xl <- c(-100, -40)
yl <- c(-40, 75)

# extract density rasters from geolocator analysis output folders
periods = c("Non-breeding period")
nbr.dens <- findSlicesData(periods, xlim = xl, ylim = yl) #helper function

for (i in names(nbr.dens)){

  #values(nbr.dens[[i]])[!is.na(values(nbr.dens[[i]]))] <- 1
  
  values(nbr.dens[[i]])[is.na(values(nbr.dens[[i]]))] <- 0

  #sum of non-NA values in raster i
  ras.sum <- sum(values(nbr.dens[[i]])[!is.na(values(nbr.dens[[i]]))])

  #we normalize the values in the raster
  values(nbr.dens[[i]])[!is.na(values(nbr.dens[[i]]))] <- values(nbr.dens[[i]])[!is.na(values(nbr.dens[[i]]))]/ras.sum
}

#Verify that each raster extracted makes sense.
# for (i in names(nbr.dens)){
#   plot(nbr.dens[[i]], useRaster = F,col = c("transparent", rev(viridis::viridis(50))),
#        main = i)
#   plot(wrld_simpl, xlim=xl, ylim=yl,add = T, bg = adjustcolor("black",alpha=0.1))
# }

#Extract the geolocator reference data to create several groups
ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")

# Create a raster for eastern populations ###################################### 
east.geo <- ref.data[(ref.data$Range_region == "Eastern"),]$geo.id
east.dens <- nbr.dens[names(nbr.dens) %in% east.geo] # remove any geolocators taht have no density rasters 

#merge the rasters
names(east.dens)[1:2] <- c('x', 'y')
east.dens$fun <- mean
east.dens$na.rm <- TRUE

east.ras <- do.call(mosaic, east.dens)

#plot probability density for eastern populations
east.col <- colorRampPalette(c("lightgray", "#D55E00"))

plot(east.ras , col = east.col(10))
plot(wrld_simpl, xlim=xl, ylim=yl,add = T, bg = adjustcolor("black",alpha=0.1))

# Create a raster for central populations ###################################### 
cent.geo <- ref.data[(ref.data$Range_region == "Central"),]$geo.id
cent.dens <- nbr.dens[names(nbr.dens) %in% cent.geo] # remove any geolocators taht have no density rasters 

#merge the rasters
names(cent.dens)[1:2] <- c('x', 'y')
cent.dens$fun <- mean
cent.dens$na.rm <- TRUE

cent.ras <- do.call(mosaic,cent.dens)

#plot probability density for eastern populations
cent.col <- colorRampPalette(c("lightgray", "#009E73"))

plot(cent.ras , col = cent.col(10))
plot(wrld_simpl, xlim=xl, ylim=yl,add = T, bg = adjustcolor("black",alpha=0.1))

# Create a raster for Western populations ###################################### 
west.geo <- ref.data[(ref.data$Range_region == "West"),]$geo.id
west.dens <- nbr.dens[names(nbr.dens) %in% west.geo] # remove any geolocators taht have no density rasters 

#merge the rasters
names(west.dens)[1:2] <- c('x', 'y')
west.dens$fun <-mean
west.dens$na.rm <- TRUE

west.ras <- do.call(mosaic, west.dens)

#plot probability density for eastern populations
west.col <- colorRampPalette(c("lightgray", "#0072B2"))

plot(west.ras , col = west.col(10))
plot(wrld_simpl, xlim=xl, ylim=yl,add = T, bg = adjustcolor("black",alpha=0.1))

# Create masking raster for areas over water ###################################

# empty raster 
r <- raster(nrows = 2 * diff(yl), ncols = 2 * diff(xl), xmn = xl[1],
            xmx = xl[2], ymn = yl[1], ymx = yl[2], crs = proj4string(wrld_simpl))

# create a raster
rs = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
           rasterize(wrld_simpl, r, 1, silent = TRUE), 
           rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))

rs.repro <-projectRaster(rs, east.ras, method = 'ngb')

# Apply the mask to the other rasters
east.ras.m <- east.ras * rs.repro
cent.ras.m <- cent.ras * rs.repro
west.ras.m <- west.ras * rs.repro

# plot each raster after masking
plot(east.ras.m , col = east.col(4))
plot(wrld_simpl, xlim=xl, ylim=yl,add = T, bg = adjustcolor("black",alpha=0.1))

plot(cent.ras.m , col = cent.col(4))
plot(wrld_simpl, xlim=xl, ylim=yl,add = T, bg = adjustcolor("black",alpha=0.1))

plot(west.ras.m , col = west.col(4))
plot(wrld_simpl, xlim=xl, ylim=yl,add = T, bg = adjustcolor("black",alpha=0.1))


# Create a plot with all three raster distributions (with generic plot function) ################

#Create a layout for the plot
layout.matrix <- matrix(c(1, 2, 1, 3, 1, 4), nrow = 2, ncol = 3)

layout(mat = layout.matrix,
       heights = c(1, 1), # Heights of the two rows
       widths = c(2, 2, 2)) # Widths of the two columns

# Plot 1: Nonbreeding sites

par (mar = c(1,1,1,1), mai = c(0, 0, 0, 0))

world <- spData::world

plot(wrld_simpl, xlim= c(-170, -40), ylim= c(42, 70), bg = "white", col = "lightgray", border ="darkgray", lwd=0.5)
points(ref.data[(ref.data$Range_region == "Eastern"),]$deploy.longitude,
       ref.data[(ref.data$Range_region == "Eastern"),]$deploy.latitude,
     col = "black", bg = "#D55E00", pch = 21, cex = 2)
points(ref.data[(ref.data$Range_region == "Central"),]$deploy.longitude,
       ref.data[(ref.data$Range_region == "Central"),]$deploy.latitude,
       col = "black", bg = "#009E73", pch = 21, cex = 2)
points(ref.data[(ref.data$Range_region == "West"),]$deploy.longitude,
       ref.data[(ref.data$Range_region == "West"),]$deploy.latitude,
       col = "black", bg = "#0072B2", pch = 21, cex = 2)
box(col = "darkgray", lwd = 0.5)

# Plot 2: Breeding areas for the Western populations

# limits for the plot
xl.p <- c(-90, -30)
yl.p  <- c(-1, 1)

par (mar = c(1,1,1,1), mai = c(0, 0, 0, 0))

crs(east.ras.m) <- crs(world)
plot(wrld_simpl, xlim= xl.p, ylim= yl.p, bg = "white", col = "lightgray", border = "darkgray", lwd=0.5)
plot(west.ras.m , col = west.col(100), alpha = 0.9, legend = F, axes=FALSE, add = T)
plot(wrld_simpl, xlim= xl.p, ylim= yl.p, bg = "white", add = T, border = "darkgray", lwd=0.5)

west.samp <- paste( "N = ", length(west.geo) -2) # number of geolocators (sample size)
#Added -2 because some geolocators were excluded from the dataset 

text(-45, 10, cex = 1.6, (bquote(paste(bold(.(west.samp))))))



box(col = "darkgray", lwd = 0.5)

# Plot 3: Breeding areas for the central populations

par (mar = c(1,1,1,1), mai = c(0, 0, 0, 0))

crs(east.ras.m) <- crs(world)
plot(wrld_simpl, xlim= xl.p, ylim= yl.p, bg = "white", col = "lightgray", border = "darkgray", lwd=0.5)
plot(cent.ras.m , col = cent.col(100), alpha = 0.9, legend = F, axes=FALSE, add = T)
plot(wrld_simpl, xlim= xl.p, ylim= yl.p, bg = "white", add = T, border = "darkgray", lwd=0.5)

cent.samp <- paste( "N = ", length(cent.geo)) # number of geolocators (sample size)

text(-45, 10, cex = 1.6, (bquote(paste(bold(.(cent.samp))))))

box(col = "darkgray", lwd = 0.5)

# Plot 4: Breeding areas for the Eastern populations

par (mar = c(1,1,1,1), mai = c(0, 0, 0, 0))

crs(east.ras.m) <- crs(world)
plot(wrld_simpl, xlim= xl.p, ylim= yl.p, bg = "white", col = "lightgray", border = "darkgray", lwd=0.5)
plot(east.ras.m , col = east.col(100), alpha = 0.9, legend = F, axes=FALSE, add = T)
plot(wrld_simpl, xlim= xl.p, ylim= yl.p, bg = "white", add = T, border = "darkgray", lwd=0.5)

east.samp <- paste( "N = ", length(east.geo)-1) # number of geolocators (sample size)
#Substract 1 due to the geolocator located in Columbia 

text(-45, 10, cex = 1.6, (bquote(paste(bold(.(east.samp))))))

box(col = "darkgray", lwd = 0.5)

# Create a plot with all three raster distributions (with ggplot2) #############

#Convert all rasters to dataframes 
east.ras.p <- rasterToPoints(east.ras.m)
east.ras.df <- data.frame(east.ras.p)

cent.ras.p <- rasterToPoints(cent.ras.m)
cent.ras.df <- data.frame(cent.ras.p)

west.ras.p <- rasterToPoints(west.ras.m)
west.ras.df <- data.frame(west.ras.p)

# ggplot for of the onbreeding distribution for each region 

east.nbr.dist <- ggplot(st_as_sf(wrld_simpl))+
                      geom_sf(fill = "white", colour = NA) +
                      coord_sf(xlim = c(-88, -30),ylim = c(30, -15)) +
                      geom_tile(east.ras.df, mapping = aes(x, y, fill = layer)) +
                      scale_fill_gradientn(colours = east.col(3), guide = "none")+
                      geom_sf(fill = NA, colour = "darkgray") +
                      coord_sf(xlim = c(-88, -30),ylim = c(30, -15)) +
                      theme_bw() +
                      theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank()) +
                      theme(text = element_text(size = 14))

cent.nbr.dist <- ggplot(st_as_sf(wrld_simpl))+
                      geom_sf(fill = "white", colour = NA) +
                      coord_sf(xlim = c(-88, -30),ylim = c(30, -15)) +
                      geom_tile(cent.ras.df, mapping = aes(x, y, fill = layer)) +
                      scale_fill_gradientn(colours = cent.col(3), guide = "none")+
                      geom_sf(fill = NA, colour = "darkgray") +
                      coord_sf(xlim = c(-88, -30),ylim = c(30, -15)) +
                      theme_bw() +
                      theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank()) +
                      theme(text = element_text(size = 14))


west.nbr.dist <- ggplot(st_as_sf(wrld_simpl))+
                      geom_sf(fill = "white", colour = NA) +
                      coord_sf(xlim = c(-88, -30),ylim = c(30, -15)) +
                      geom_tile(west.ras.df, mapping = aes(x, y, fill = layer)) +
                      scale_fill_gradientn(colours = west.col(3), guide = "none") +
                      geom_sf(fill = NA, colour = "darkgray") +
                      coord_sf(xlim = c(-88, -30),ylim = c(30, -15)) +
                      theme_bw() +
                      theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank()) +
                      theme(text = element_text(size = 14))

# Arrange the plots together 

#Sys.getenv("GITHUB_PAT")
#Sys.unsetenv("GITHUB_PAT")

# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")

#library (ggpubr)

ggarrange(east.nbr.dist, cent.nbr.dist, west.nbr.dist, ncol = 3, nrow = 3)
