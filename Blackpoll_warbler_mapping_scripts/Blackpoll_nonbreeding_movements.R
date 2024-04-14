# Plot of of nonbreeding movements 

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
library(fmdates)

source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Import data ##################################################################

# path to reference data file 
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

# location data 
geo.all <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/All.locations.csv")

# data preparation ####################################################

# extract Nonbreeding stationary sites
NB.stat <- geo.all %>% dplyr::filter(sitenum > 0, site_type %in% c("Nonbreeding"),
                   period %in% c("Non-breeding period"),
                   Recorded_North_South_mig %in% c("Both", "South and partial North") | geo_id == "WRMA04173") %>%
  ungroup()

# # remove points for individuals that moved less than 250 km over the course of the wintering season 
# NB.stat <- NB.stat %>% group_by(geo_id) %>% mutate(Lon.50.start = first(Lon.50.),
#                                                    Lat.50.start = first(Lat.50.),
#                                                    Lon.50.end = last(Lon.50.),
#                                                    Lat.50.end = last(Lat.50.))%>%
#   rowwise() %>%
#   mutate(dist = distHaversine(c(Lon.50.start ,Lat.50.start), c(Lon.50.end, Lat.50.end)),.after = geo_id) %>%
#   filter(dist > 250000 | dist == 0)

# Modify the data so we have two points for each movement between two sites (which means that we will have duplicate points)
NB.stat <- rbind(NB.stat, NB.stat)
NB.stat <- NB.stat %>% arrange(geo_id, StartTime)

# Remove movements where the distance is less than 250 km
NB.stat <- NB.stat %>% group_by(geo_id) %>% mutate(move.start = lag(EndTime),
                                                   move.end = StartTime,
                                                   start.lon = lag(Lon.50.),
                                                   start.lat = lag(Lat.50.),
                                                   end.lon = Lon.50.,
                                                   end.lat = Lat.50.)%>%
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
            mean.lon2.5 = mean(Lon.2.5.), mean.lat2.5  = mean(Lat.97.5.),
            mean.lon97.5 = mean(Lon.2.5.), mean.lat97.5 = mean(Lat.97.5.),)

# Save nonbreeding movement data ###############################################
write.csv(NB.stat, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/nonbreeding.movements.csv")
write.csv(NB.stat.mean, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/nonbreeding.mean.csv")

  







