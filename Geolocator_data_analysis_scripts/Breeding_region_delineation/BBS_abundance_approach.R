# Using BBS abundance trends to delineate populations following Rushing et al. 2015
# https://besjournals-onlinelibrary-wiley-com.subzero.lib.uoguelph.ca/doi/full/10.1111/1365-2664.12579

# clear objects from workspace
rm(list=ls())

#Load libraries
library(tidyverse)
library(igraph)
library(lubridate)
library(anytime)
library(ggplot2)
library(geosphere)
library(terra)
library(sf)
library(spData)
library(maptools)
library(ebirdst)
library(maptree)

# visualize BBS routes where blackpoll warblers were encountered between 2000 - 2021
BBS_routes_bpw <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/BBS_data/bpw_route_profiles.csv")
bpw_fall_ab <- load_raster("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports/2021/bkpwar",
                             product = "abundance",
                             period = "seasonal",
                             resolution = "mr")

bpw_fall_ab  <- terra::project(bpw_fall_ab, crs(wrld_simpl))

data(wrld_simpl)
plot(bpw_fall_ab$breeding, xlim = c(-170, -50), ylim = c(40, 75))
plot(wrld_simpl,add = T)
points(BBS_routes_bpw$Longitude, BBS_routes_bpw$Latitude)

# Load BBS yearly count data 
BBS_counts <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/BBS_data/bpw_total_counts.csv") %>%
  arrange(Route)
BBS_counts$Route <- as.factor(BBS_counts$Route)

# fit poisson model for abundance by year
ab.model <- glm(SpeciesTotal ~ Year + Route, family = poisson(link = "log"), data = BBS_counts)
summary(ab.model)

# create a dataset of counts per route
route_trends <- data.frame(Route = as.numeric(unique(BBS_counts$Route)),
                             Trend = c(1, ab.model$coefficients[3:length(ab.model$coefficients)]))
route_trends <- merge(route_trends, BBS_routes_bpw, by = "Route")

# add an abundance estimate
ab.est <- BBS_counts %>% group_by(Route) %>% summarise(Abundance = mean(SpeciesTotal))
route_trends <- merge(route_trends, ab.est, by = "Route")

# plot trends
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -45),ylim = c(40, 70))+
  geom_point(data = route_trends, mapping = aes(x = Longitude, y = Latitude, colour = Trend))


# cluster points by trend and position 
dist.mat <- dist(route_trends[, c("Trend", "Abundance", "Longitude", "Latitude")])
clust.tree <- hclust(dist.mat)
#min(kgs(clust.tree, dist.mat))
route_trends$Cluster <- as.factor(cutree(clust.tree, k = 13))


# plot clusters
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -45), ylim = c(40, 70))+
  geom_point(data = route_trends, mapping = aes(x = Longitude, y = Latitude, colour = Cluster))


