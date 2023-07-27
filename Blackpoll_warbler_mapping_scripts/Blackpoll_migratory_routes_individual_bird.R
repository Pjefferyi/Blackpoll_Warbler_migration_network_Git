#Load libraries
library(tidyverse)
library(igraph)
library(lubridate)
library(anytime)
library(ggplot2)
library(geosphere)
library(terra)
library(sf)
library(maptools)
library(ebirdst)


# Create a list of path to all files with fit data
raster.paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data",
                    pattern = "SGAT_GroupedThreshold_fit.R", recursive = T, full.names = T )

# Create a list of geolocator names for the slices
spl.paths <- strsplit(raster.paths , "/")
names(raster.paths) <- sapply(spl.paths , "[[", 11)

# location data 
geo.all <- findLocData(geo.ids = c("V8757_010",
                                   "V8296_004",
                                   "V8296_005",
                                   "V8296_006",
                                   "V8757_055",
                                   "V8757_018",
                                   "V8757_021",
                                   "V8296_015",
                                   "V8296_017",
                                   "V8296_026",
                                   "V8296_025",
                                   "V8296_007",
                                   "V8296_008",
                                   "V8757_019",
                                   "V8757_096",
                                   "V8757_134",
                                   "V8757_029",
                                   "V8757_078",
                                   "blpw09",
                                   "blpw12",
                                   "3254_001",
                                   "4068_014",
                                   "blpw14",
                                   "3254_003",
                                   "3254_008",
                                   "3254_011",
                                   "3254_057",
                                   "blpw15",
                                   "blpw25",
                                   "4105_008",
                                   "4105_009",
                                   "4105_016",
                                   "4105_017",
                                   "4210_002",
                                   "4210_004",
                                   "4210_006",
                                   "4210_010",
                                   "WRMA04173",
                                   "A",
                                   "B",
                                   "C",
                                   #"E",
                                   "D"), check_col_length = F)

data(wrld_simpl)


#Create a grid to project density rasters
xlim <- c(-170, -30)
ylim <- c(-15, 70)
r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1],
            xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))


for (i in unique(geo.all$geo_id)){
  
  # open jpeg
  jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/Appendix_figures/",i,"_locations.jpg"),units = "cm", width = 15, height = 15,  quality = 100, res = 600)
  
  geo.data <- geo.all %>% filter(geo_id == i) 
  path.palette <- colorRampPalette(c("yellow", "red"))(max(geo.data$sitenum) * 2)
  stat.position = geo.data %>% filter(sitenum > 0)
  
  equinox.time <- anytime(paste0(year(as.character(geo.data$StartTime[1])),
                                 "-09-22"))
  
  # extract raster
  load(raster.paths[i])
  s <- slices(type = "intermediate", breaks = "weeks", mcmc = fit, grid = r)
  sk <- slice(s, sliceIndices(s))
  
  #add column to indicate whether the locations are within 15 days of the equinox 
  geo.data <- geo.data %>% 
    mutate(equinox = case_when(StartTime > equinox.time - days(15) & EndTime < equinox.time + days(15) ~ "equinox",
           EndTime < equinox.time - days(15) ~ "pre.equinox", 
           EndTime  > equinox.time + days(15) ~ "post.equinox")) 
  
  #plot raster 
  plot(sk, useRaster = F, col = c("transparent", rev(viridis::viridis(50))))
  
  #plot map
  plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
       xlim = c(-170, -30), ylim = c(-15, 70), col = NA, lwd = 0.5, add = T)
  
  #line for the path used by the bird
  lines(geo.data$Lon.50.[geo.data$equinox == "pre.equinox"], geo.data$Lat.50.[geo.data$equinox== "pre.equinox"], col = "blue", lty = "solid")
  lines(geo.data$Lon.50., geo.data$Lat.50., col = "blue", lty = "dashed")
  lines(geo.data$Lon.50.[geo.data$equinox== "post.equinox"], geo.data$Lat.50.[geo.data$equinox == "post.equinox"], col = "blue", lty = "solid")
  
  # information for stationary points 
  with(geo.data[geo.data$sitenum>0,], arrows(`Lon.50.`, `Lat.2.5.`, `Lon.50.`, `Lat.97.5.`, length = 0, lwd = 1.5, col = "firebrick"))
  with(geo.data[geo.data$sitenum>0,], arrows(`Lon.2.5.`, `Lat.50.`, `Lon.97.5.`, `Lat.50.`, length = 0, lwd = 1.5, col = "firebrick"))
  points(stat.position$Lon.50., stat.position$Lat.50., cex = 2, col =  path.palette, pch = 19)
  text(stat.position$Lon.50., stat.position$Lat.50., round(stat.position$duration))
  
  #Add plot tile
  title(main = paste0(i," - ", geo.data$study.site[1]))
  
  dev.off()
  
}

     