#load packages
require(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(remotes)
library(anytime)
library(lubridate)
library(ebirdst)

#load spatial packages 
library(ggmap)
library(terra)
library(spData)
library(sf)
library(cluster)

# Function that takes a point of longitude and returns the intersection between 
# this point and the coast of North America (lon, lat)

  lon.pt = -73.43 
  
  # Transfrom the point into a line object in terra. 
  d <- data.frame(id = 1:2, name = c("a", "b"))
  d$lon <- c(lon.pt, lon.pt)
  d$lat <- c(60, -60)
  
  lon.line <- vect(cbind(d$lon, d$lat),  type="lines",  crs= crs(st_as_sf(wrld_simpl)))
  
  # Coast of Northern south america 
  nsa.c <- (wrld_simpl[(wrld_simpl$NAME %in% c("Brazil",
                                               "Colombia",
                                               "Venezuela",
                                               "Guyana",
                                               "Suriname",
                                               "French Guiana",
                                               "Panama",
                                               "Ecuador")),]
  
  plot(nsa.c, xlim = c(-165, -35), ylim = c(-10, 65))
  lines(lon.line, col = "red")
  
  intersect()

  
  