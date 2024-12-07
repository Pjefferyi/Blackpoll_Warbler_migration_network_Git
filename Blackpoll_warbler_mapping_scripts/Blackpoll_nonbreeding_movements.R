# Processing of nonbreeding movement data for plotting 

#Load libraries
library(tidyverse)
library(igraph)
library(lubridate)
library(anytime)
library(ggplot2)
library(geosphere)
library(terra)
library(sf)
library(ebirdst)
library(fmdates)

source("Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Import data ##################################################################

# path to reference data file 
ref_path <- "Analysis_Output_Data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

# location data 
geo.all <- read.csv("Network_construction/All.locations.csv")

# data preparation ####################################################

# extract Nonbreeding stationary sites
NB.stat <- geo.all %>% dplyr::filter(sitenum > 0, site_type %in% c("Nonbreeding"),
                   period %in% c("Non-breeding period"),
                   Recorded_North_South_mig %in% c("Both", "South and partial North") | geo_id == "WRMA04173") %>%
  ungroup()

# Modify the data so we have two points for each movement between two sites (which means that we will have duplicate points, for ease of plotting)
NB.stat <- rbind(NB.stat, NB.stat)
NB.stat <- NB.stat %>% arrange(geo_id, StartTime)

# Remove movements where the distance is less than 250 km
NB.stat <- NB.stat %>% group_by(geo_id) %>% mutate(move.start = lag(EndTime),
                                                   move.end = StartTime,
                                                   start.lon = lag(Lon.50.),
                                                   start.lat = lag(Lat.50.),
                                                   end.lon = Lon.50.,
                                                   end.lat = Lat.50.,
                                                   start.lon.2.5 = lag(Lon.2.5.),
                                                   start.lat.2.5 = lag(Lat.2.5.),
                                                   start.lon.97.5 = lag(Lon.97.5.),
                                                   start.lat.97.5 = lag(Lat.97.5.))%>%
  rowwise() %>%
  mutate(dist = distHaversine(c(Lon.50.,Lat.50.), c(start.lon, start.lat)),.after = geo_id) %>%
  filter(dist >= 250000)

# Add a new column specifying whether individuals started making their nonbreeding movements during the equinox, and during the late or early winter 
NB.stat <- NB.stat %>% mutate(spring.equinox.date = anytime(spring.equinox.date),
                               equinox.nbr.move = ifelse(move.start > anytime(spring.equinox.date) - days(14) & move.start< anytime(spring.equinox.date) + days(14),
                                      "equinox affected", "equinox free"),
                               timing.nbr.move =ifelse(month(anytime(move.start)) >= 9,
                                      "fall.nbr.movement", "winter.nbr.movements"))


# Create a dataset wit the mean location of individuals that did not move 
NB.stat.mean <-  geo.all %>% dplyr::filter(sitenum > 0, site_type %in% c("Nonbreeding"),
                                             period %in% c("Non-breeding period"),
                                             Recorded_North_South_mig %in% c("Both", "South and partial North") | geo_id == "WRMA04173") %>%
  filter(!(geo_id %in% NB.stat$geo_id)) %>% group_by(geo_id) %>%
  summarize(mean.lon = mean(Lon.50.), mean.lat = mean(Lat.50.),
            mean.lon2.5 = mean(Lon.2.5.), mean.lat2.5  = mean(Lat.2.5.),
            mean.lon97.5 = mean(Lon.97.5.), mean.lat97.5 = mean(Lat.97.5.))

# Save nonbreeding movement data ###############################################
write.csv(NB.stat, "Analysis_input_Data/nonbreeding.movements.csv")
write.csv(NB.stat.mean, "Analysis_input_Data/nonbreeding.mean.csv")

  







