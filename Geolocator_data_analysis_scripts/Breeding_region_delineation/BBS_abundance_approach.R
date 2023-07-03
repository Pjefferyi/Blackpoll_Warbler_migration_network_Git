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

# visualize BBS routes where blackpoll warblers were encountered between 2000 - 2021
BBS_routes_bpw <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/BBS_data/bpw_route_profiles.csv")

data(wrld_simpl)
plot(wrld_simpl, xlim = c(-170, -55), ylim = c(40, 50))
points(BBS_routes_bpw$Longitude, BBS_routes_bpw$Latitude)
