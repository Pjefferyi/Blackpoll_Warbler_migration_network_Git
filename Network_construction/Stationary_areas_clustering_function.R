# clear objects from workspace
#rm(list=ls())

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
library(cluster)

# Load helper functions 
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Import data ##################################################################

# path to reference data file 
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

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


# function to iteratively cluster geolocator stationary location estimates
# This script is based on the methods described by Lagasse et al. 2022: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0270957#sec015

# First extract stationary locations for the fall 
geo.all <- geo.all %>% group_by(geo_id) %>% mutate(site_type = case_when(
  (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ "Breeding",
  sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ "Breeding",
  (sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ "Breeding",
  (sitenum == max(sitenum)) & Recorded_North_South_mig == "South and partial North" & path.elongated == T ~ "Breeding",
  period == "Non-breeding period" & (duration >= 14 | sitenum == 1 | sitenum == max(sitenum)) ~ "Nonbreeding",
  .default = "Stopover"))

fall.stat <- geo.all %>% filter(sitenum > 0, site_type %in% c("Stopover","Nonbreeding"),
                                period %in% c("Post-breeding migration","Non-breeding period"),
                                Recorded_North_South_mig %in% c("Both", "South and partial North", "South"))


# Generate an initial distance matrix

# function to Calculate the geodesic distance between points and creates a distance matrix
geo.dist = function(df) {
  require(geosphere)
  d <- function(i,z){         # z[1:2] contain long, lat
    dist <- rep(0,nrow(z))
    dist[i:nrow(z)] <- distHaversine(z[i:nrow(z),1:2],z[i,1:2])
    return(dist/1000)
  }
  dm <- do.call(cbind,lapply(1:nrow(df),d,df))
  return(as.dist(dm))
}

#calculate distance matrix  
dist.matrix <- geo.dist(geo.stat[,c("Lon.50.","Lat.50.")])

#cluster stationary locations using the partition PAM algorithm 
statCluster <- function(geo.stat, maxdiam =  600, k = 1, dist.matrix){
  
   geo.stat = fall.stat
   maxdiam =  700
   k = 1
  
  # The cluster diameter is in km 
  diam = Inf
  runs = 1
  
  while (diam > maxdiam){
    
    # cluster stationary location and add locations to dataframe 
    clust <- pam(dist.matrix, k, diss  = T, cluster.only =  T)
    geo.stat$cluster <- clust
    
    # for each cluster, calculate the diameter
    diam.list <- c()
    
    for (i in unique(clust)){
    # Get the stationary locations in each subset 
    clust.subset <- geo.stat %>% dplyr::filter(cluster == i)
    
    if (nrow(clust.subset) == 1){
        subset.dist.matrix <- 0
      } else {
    
      # Generate distance matrix
      subset.dist.matrix <- geo.dist(clust.subset [,c("Lon.50.","Lat.50.")])
  
      # Calculate the cluster diameter (here defined as the distance between the two furtherst locations in the cluster)
      diam.list <- append(diam.list, max(subset.dist.matrix))
      }
    }
    
    diam <- mean(diam.list)
    print(paste("mean diameter =", as.character(diam), "km"))
    print(paste("Run #", as.character(runs)))
    print(paste("k =", as.character(k)))
    runs <-  runs + 1 
    k = k + 1 
  }
  
  return(k-1)
}
  

# A look at the output cluster parameter ########################################

# generate clusters with the value of k output form the function 
clust <- pam(dist.matrix, k, diss  = T, cluster.only =  T)
fall.stat$cluster <- clust

ggplot(spData::world[(spData::world$continent %in% c("North America", "South America")),]) +
  geom_sf() +
  coord_sf() +
  geom_point(data = fall.stat[(geo.stat$sitenum) > 0,], aes(x = Lon.50., y = Lat.50., colour = as.factor(cluster)), size = 0.3)+
  theme(legend.position = "None")
