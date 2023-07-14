# Calculate the strength of migratory connectivity for the blackpoll warbler

# migconnectivity package 
#install.packages("devtools")
# usethis::create_github_token()
# usethis::edit_r_environ()
#devtools::install_github("SMBC-NZP/MigConnectivity")
library(MigConnectivity)

# load required packages 
library(tidyverse)
library(lubridate)
library(anytime)
library(ggplot2)
library(geosphere)
library(terra)
library(sf)
library(maptools)
library(ebirdst)


# data preparation############################################################## 

# run the network preparation script
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Network_construction.R")

# fall network breeding sites
fall.br <- fall.breed 
fall.br.sf <- st_as_sf(fall.br, coords = c("Lon.50.", "Lat.50."))
st_crs(fall.nbr.sf ) <- st_crs(wrld_simpl) 

# fall network nonbreeding sites 
fall.nbr <- fall.stat %>% group_by(geo_id) %>% filter(sitenum == max(sitenum))
fall.nbr.sf <- st_as_sf(fall.nbr, coords = c("Lon.50.", "Lat.50."))
st_crs(fall.nbr.sf ) <- st_crs(wrld_simpl) 
#st_write(fall.nbr.sf, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data_exports/bpw_nbr_sites.shp")

# Origin sites 
fall.br.regions <- abundance.regions 

#Target sites 
fall.nbr.regions <-  read_sf(


## migratory connectivity between the blackpoll warbler's breeding and first fall nonbreeding sites ####









